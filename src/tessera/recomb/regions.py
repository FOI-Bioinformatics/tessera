"""Recombination-region calling.

The query is most similar to one dataset across the genome -- the **major
parent** (the backbone donor). A recombination is indicated where a stretch of
the query is instead closest to another dataset (a **minor parent**).

The default caller (``method="hmm"``) segments the query against the reference
panel with a hidden Markov model (:mod:`tessera.recomb.hmm`): emissions are the
per-window binomial copying likelihoods and a single jump rate penalises
switching reference, so near-identical references do not flip and thin windows
cannot drive a call. It yields breakpoint uncertainty intervals and a posterior
support per segment.

The legacy heuristic (``method="heuristic"``) is kept for comparison: major =
most windows won; a window is recombinant when a minor beats the major by
``margin``; consecutive recombinant windows are merged within ``merge_gap`` and
runs shorter than ``min_region`` dropped. Both report coordinates in MSA columns
and query bases.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isnan
from statistics import mean

import numpy as np

from ..core.errors import UserInputError
from .analyze import AnalysisResult, rank_datasets
from .clusters import _window_bounds, clustered_view, lineage_clusters
from .hmm import DEFAULT_JUMP_RATE, Segment, segment_query
from .siblings import SiblingEvidence, sibling_aware_states
from .similarity import GAP, WindowSimilarity, discordant_counts, region_identity
from .stats import benjamini_hochberg, sign_test_pvalue
from .threeseq import triplet_steps

# The region callers, in canonical (display) order; the single source of truth for the
# valid ``--method`` set. The default ensemble runs the four statistical callers and
# merges them; heuristic is legacy and opt-in (or via ``all``).
CALLERS = ("hmm", "3seq", "maxchi", "bootscan", "geneconv", "barcode", "heuristic")
DEFAULT_METHODS = ("hmm", "3seq", "maxchi", "bootscan")


def parse_methods(spec: str) -> tuple[str, ...]:
    """Parse a ``--method`` value into an ordered, de-duplicated tuple of callers.

    Accepts a comma-separated list of :data:`CALLERS` or the keyword ``all``. Order
    follows :data:`CALLERS` regardless of the input order, so the report is stable.
    Raises :class:`UserInputError` on an unknown token.
    """
    tokens = [t.strip() for t in spec.split(",") if t.strip()]
    if not tokens:
        raise UserInputError("No --method given. Choose from: all, " + ", ".join(CALLERS))
    if "all" in tokens:
        return CALLERS
    unknown = [t for t in tokens if t not in CALLERS]
    if unknown:
        raise UserInputError(
            f"Invalid --method {', '.join(unknown)}. Choose from: all, " + ", ".join(CALLERS)
        )
    return tuple(c for c in CALLERS if c in tokens)


@dataclass
class RegionParams:
    """Region-calling parameters. Heuristic thresholds default from the window size."""

    min_region: int
    margin: float
    merge_gap: int
    method: str = "hmm"  # "hmm" (default) or "heuristic"
    jump_rate: float = DEFAULT_JUMP_RATE  # HMM: prior switch probability per window
    identity: float | None = None  # HMM emission identity; None -> estimate from data
    alpha: float = 0.05  # HMM: significance level for the donor-vs-major test
    # Exclude the query's whole-genome siblings (its own lineage) from the parental
    # competition before segmenting, so a sibling cannot win every window and mask
    # the event. See :mod:`tessera.recomb.siblings`. Defaults from first principles.
    exclude_siblings: bool = True
    sibling_coverage: float = 0.90  # a sibling is comparable over >= this of the query
    sibling_beaten: float = 0.15  # ... significantly beaten in <= this fraction of windows
    sibling_alpha: float = 0.05  # significance level for the per-window beaten test
    # ... and never beaten over a region-sized contiguous run (the primary signal; a
    # parent is beaten over the co-parent's tract). Derived from min_region / step.
    # Cluster near-duplicate references into lineages and compete pooled lineages, so
    # duplicates do not tie every window and fragment the call. See recomb/clusters.py.
    cluster_lineages: bool = True
    # Two references are one lineage when they agree at least this much in every window
    # with no region-sized run below it. Absolute and conservative: it merges only
    # clearly near-identical duplicates (the genomes that tie and fragment the call),
    # leaving distinct parents -- even diffusely-divergent ones a few % apart -- separate.
    cluster_merge_identity: float = 0.985
    cluster_min_refs: int = 4  # skip clustering below this many references (curated panels)
    # ... and an upper bound. Clustering compares every pair, allocating two full-width
    # cumulative-sum arrays per pair, so the cost grows quadratically with the panel: a
    # 1000-genome panel is ~500k pairs. Above this, competing genomes individually is
    # the lesser evil -- clustering exists to stop near-duplicates fragmenting a call,
    # and a panel this large has been dereplicated upstream.
    cluster_max_refs: int = 200
    # Reference genotypes (label -> clade) for the barcode caller's per-clade markers. Used
    # ONLY by that caller, never the genome competition, so it cannot reproduce the reverted
    # competition-layer pooling regression.
    lineage_map: dict[str, str] | None = None
    # Lineage clusters computed once for the run (clusters.lineage_clusters): multi-member
    # clusters only, ``[]`` for none. ``None`` means "not computed": the HMM caller then
    # computes them itself, as it did before they were shared.
    reference_clusters: list[list[str]] | None = None

    @classmethod
    def with_defaults(
        cls,
        window_size: int,
        min_region: int | None = None,
        margin: float = 0.0,
        merge_gap: int | None = None,
        method: str = "hmm",
        jump_rate: float = DEFAULT_JUMP_RATE,
        identity: float | None = None,
        alpha: float = 0.05,
        exclude_siblings: bool = True,
        cluster_lineages: bool = True,
        lineage_map: dict[str, str] | None = None,
        reference_clusters: list[list[str]] | None = None,
    ) -> RegionParams:
        return cls(
            min_region=window_size if min_region is None else min_region,
            margin=margin,
            merge_gap=window_size if merge_gap is None else merge_gap,
            method=method,
            jump_rate=jump_rate,
            identity=identity,
            alpha=alpha,
            exclude_siblings=exclude_siblings,
            cluster_lineages=cluster_lineages,
            lineage_map=lineage_map,
            reference_clusters=reference_clusters,
        )


@dataclass
class Region:
    minor_parent: str
    major_parent: str
    msa_start: int
    msa_end: int
    query_start: int
    query_end: int
    n_windows: int
    mean_sim_minor: float
    mean_sim_major: float
    margin: float
    # HMM caller: posterior support of the donor over the segment, and the query
    # coordinate uncertainty of the segment's breakpoint.
    posterior_support: float = 1.0
    breakpoint_lo: int | None = None
    breakpoint_hi: int | None = None
    # Share of distinguishing (discordant) sites that favour the donor, and the
    # sign-test p-value with its Benjamini-Hochberg q-value (FDR across segments).
    # What `support` counts differs by caller, and `statistic` names the test rather
    # than the unit, so one case needs saying here: Bootscan's value is the support of
    # the donor's lineage cluster, pooled over its members, where the donor is one of a
    # cluster of near-duplicates, and the support of the single genome otherwise. Its
    # `statistic` reads "mean bootstrap support for the donor" either way, so the row
    # does not say which of the two it is.
    support: float | None = None
    pvalue: float | None = None
    qvalue: float | None = None
    # Set by coverage analysis: the donor's own similarity is poor, so it may be a
    # stand-in for an absent true donor. Default False (not yet evaluated).
    donor_undercovered: bool = False
    # The backbone is locally gap-degraded over this region: gappier here than it is
    # across the rest of its own genome, so the comparison that produced the region
    # was not between comparably aligned sequences. Deliberately NOT a low
    # `mean_sim_major`: over a true donor region the backbone is the poorer match by
    # definition, so that would fire on every genuine recombinant. The gap-rate test
    # is issue #81's measured discriminator, applied here below the bar at which #81
    # withholds the region outright (issue #82).
    major_degraded: bool = False
    # A putative recombination where the query is locally far from EVERY reference:
    # the true donor is likely missing from the collection (bridged from a gap).
    donor_absent: bool = False
    # Number of genomes in the donor's / backbone's lineage cluster (1 = a single
    # genome competed; >1 = a pooled lineage), for the report.
    minor_cluster_size: int = 1
    major_cluster_size: int = 1
    # The caller(s) that produced this region. A single-method run sets one entry; the
    # ensemble merge (recomb/ensemble.py) sets the set of agreeing callers.
    methods: tuple[str, ...] = ()
    # Set by the ensemble merge when the region overlaps a Hudson-Kaplan Rmin interval
    # (parent-free corroboration from the PHI/Rmin signal).
    parent_free_support: bool = False
    # Which test produced `pvalue`, and what `support` counts. The callers do not share a
    # statistic -- a sign test on discordant sites, a max-descent test, a scan-aware
    # chi-square and a block-permutation run length all land in the same column -- so a
    # row is not interpretable without naming its own test. The ensemble merge keeps the
    # most significant member's p-value, and carries that member's names with it.
    test: str = ""
    statistic: str = ""

    @property
    def length_bp(self) -> int:
        """Region length in query bases.

        Derived rather than stored: it was set from the MSA span by every caller and
        from the query span by the ensemble merge, so one column meant two different
        things depending on which caller produced the row.
        """
        return self.query_end - self.query_start

    @property
    def length_msa(self) -> int:
        """Region width in alignment columns (includes gaps, so >= ``length_bp``)."""
        return self.msa_end - self.msa_start


@dataclass
class _Run:
    minor_counts: Counter
    indices: list[int]
    msa_start: int
    msa_end: int


def _window_minor(
    result: WindowSimilarity, major: str, i: int, margin: float
) -> str | None:
    """Minor parent that beats ``major`` by >= ``margin`` at window ``i``, else None."""
    sim_major = result.similarities[major][i]
    if isnan(sim_major):  # uninformative window for the major parent
        return None
    best_label = None
    best_sim = sim_major
    for label, sims in result.similarities.items():
        if label == major or isnan(sims[i]):
            continue
        if sims[i] > best_sim:
            best_sim = sims[i]
            best_label = label
    if best_label is None or (best_sim - sim_major) < margin:
        return None
    return best_label


# A region is withheld only when the backbone's gap rate over it is both this many
# times its own genome-wide rate and this high in absolute terms. Measured over the
# five seeds of the gap-degraded fixture, the mpox panel of issue #81 and every
# region of the eight curated datasets in validation/: a backbone locally degraded
# enough to lose the comparison sits at 2.27x to 19.06x its genome-wide rate with a
# local rate of 5.02 % to 31.21 %, while the backbone under a real recombinant region
# stays at or below 0.74x and 3.09 %. Both conditions are required, so a backbone
# with no gaps at all elsewhere is not withheld for a single gap column, and a draft
# assembly gappy everywhere is not withheld for being uniformly gappy.
MIN_LOCAL_GAP_RATE = 0.05


def withhold_misaligned_backbone(
    regions: list[Region],
    result: WindowSimilarity,
    max_ratio: float,
    min_local_gap: float = MIN_LOCAL_GAP_RATE,
) -> tuple[list[Region], list[Region]]:
    """Split ``regions`` into those to report and those resting on a misaligned
    backbone; returns ``(kept, withheld)``.

    Region calling is a comparison between the backbone and a candidate donor over
    the same columns. Where the backbone is locally gap-degraded that comparison is
    not between comparably aligned genomes: the backbone loses the columns it was
    misaligned over, any reference that stayed intact matches the query better
    there, and the callers read the switch as a minor parent (issue #81).

    The test is local and relative to the backbone itself -- its gap rate over the
    region against its genome-wide rate -- so a uniformly gappy backbone is not
    penalised for being gappy, only for being worse here than it is elsewhere.
    ``max_ratio`` of 0 disables it.

    A donor-absent row is exempt. It says the query is far from every reference, so
    the true donor is likely missing, which is not a comparison the backbone can
    lose; one curated dataset carries such a row at twice its backbone's
    genome-wide rate.
    """
    if max_ratio <= 0:
        return list(regions), []
    kept: list[Region] = []
    withheld: list[Region] = []
    for region in regions:
        row = result.rows.get(region.major_parent)
        if region.donor_absent or row is None or row.size == 0:
            kept.append(region)
            continue
        genome = float(np.count_nonzero(row == GAP)) / row.size
        span = row[region.msa_start:region.msa_end]
        local = float(np.count_nonzero(span == GAP)) / span.size if span.size else 0.0
        degraded = local >= min_local_gap and (
            genome <= 0.0 or local / genome >= max_ratio
        )
        (withheld if degraded else kept).append(region)
    return kept, withheld


def backbone_gap_rates(
    region: Region, result: WindowSimilarity
) -> tuple[float, float] | None:
    """``(local, genome_wide)`` gap rate for ``region``'s backbone, or None.

    None when the backbone has no row in the alignment, so nothing was measured.
    Shared by the withholding above and the caveat below, which differ only in where
    they put the bar.
    """
    row = result.rows.get(region.major_parent)
    if row is None or row.size == 0:
        return None
    genome = float(np.count_nonzero(row == GAP)) / row.size
    span = row[region.msa_start:region.msa_end]
    local = float(np.count_nonzero(span == GAP)) / span.size if span.size else 0.0
    return local, genome


# The caveat bar, below the withholding bar above. Issue #81 measured a backbone
# degraded enough to lose the comparison at 2.27x to 19.06x its genome-wide rate with a
# local rate of 5.02 % to 31.21 %, against at most 0.74x and 3.09 % for the backbone
# under a real recombinant region. These sit in the empty band between the two, so the
# flag covers a backbone that is degraded but not far enough for #81 to remove the
# region, without firing on anything a real region was measured to show.
CAVEAT_LOCAL_GAP_RATE = 0.035
CAVEAT_GAP_RATIO = 1.5


def flag_degraded_backbone(
    regions: list[Region],
    result: WindowSimilarity,
    min_local_gap: float = CAVEAT_LOCAL_GAP_RATE,
    min_ratio: float = CAVEAT_GAP_RATIO,
) -> None:
    """Mark regions whose backbone is locally gap-degraded (``major_degraded``).

    Applied to the regions that survive ``withhold_misaligned_backbone``, so it
    describes the band below that bar: degraded enough to be worth saying, not enough
    to remove the region. A donor-absent row is exempt for the same reason it is
    exempt there -- it reports that every reference is far from the query, which is not
    a comparison the backbone can lose.
    """
    for region in regions:
        rates = backbone_gap_rates(region, result)
        if region.donor_absent or rates is None:
            continue
        local, genome = rates
        region.major_degraded = local >= min_local_gap and (
            genome <= 0.0 or local / genome >= min_ratio
        )


def call_regions(
    result: WindowSimilarity,
    analysis: AnalysisResult,
    window_size: int,
    params: RegionParams,
    *,
    not_run: dict[str, str] | None = None,
) -> tuple[list[Region], str | None, list[SiblingEvidence]]:
    """Call recombinant regions.

    Returns ``(regions, major_parent, excluded_siblings)`` -- the third element lists
    the query's whole-genome siblings set aside before segmenting (HMM caller only).

    A caller that ran but could not test (the HMM left with one state) records why in
    ``not_run`` under its own name, when a mapping is given; its result is unchanged.
    """
    if params.method == "hmm":
        regions, major, siblings = _call_regions_hmm(result, window_size, params, not_run)
    elif params.method == "3seq":
        from .threeseq import call_regions_3seq
        regions, major, siblings = call_regions_3seq(result, analysis, params)
    elif params.method == "maxchi":
        from .maxchi import call_regions_maxchi
        regions, major, siblings = call_regions_maxchi(result, analysis, params)
    elif params.method == "bootscan":
        from .bootscan import call_regions_bootscan
        regions, major, siblings = call_regions_bootscan(result, analysis, window_size, params)
    elif params.method == "geneconv":
        from .geneconv import call_regions_geneconv
        regions, major, siblings = call_regions_geneconv(result, analysis, params)
    elif params.method == "barcode":
        from .barcode import call_regions_barcode
        regions, major, siblings = call_regions_barcode(result, analysis, window_size, params)
    else:
        regions, major = _call_regions_heuristic(result, analysis, window_size, params)
        siblings = []
    for r in regions:  # tag provenance so the ensemble merge knows the source caller
        r.methods = (params.method,)
    return regions, major, siblings


# Why the HMM could not test when lineage clustering left it a single state.
HMM_ONE_STATE = (
    "lineage clustering pooled every reference into one lineage, so the HMM had no "
    "second state to switch to; rerun with --no-cluster-lineages"
)


def _call_regions_hmm(
    result: WindowSimilarity, window_size: int, params: RegionParams,
    not_run: dict[str, str] | None = None,
) -> tuple[list[Region], str | None, list[SiblingEvidence]]:
    """HMM segmentation; a non-major segment is reported only when its donor is a
    significantly better match than the major on the sites distinguishing them.

    Whole-genome siblings of the query (its own lineage) are excluded from the
    competition first, so a sibling cannot win every window and mask the event.

    Near-duplicate references are first pooled into lineage clusters, so duplicates
    do not tie every window and fragment the call; the rest runs on the pooled view.
    """
    if len(result.similarities) < 2:
        only = next(iter(result.similarities), None)
        return [], only, []

    # Front-end: pool near-duplicate references into lineage clusters and compete the
    # pooled lineages. A no-op for a curated panel (few refs, or all singletons).
    work = result
    cluster_size: dict[str, int] = {}
    pooled = params.reference_clusters
    if pooled is None:
        pooled = lineage_clusters(result, window_size, params)
    if pooled:
        work, cluster_size = clustered_view(result, _with_singletons(result, pooled))
    if len(work.similarities) < 2 and not_run is not None:
        # One state cannot switch, so the decode below cannot call a region whatever
        # the data. The run continues unchanged; the caller is reported as not run.
        not_run["hmm"] = HMM_ONE_STATE

    states: list[str] | None = None
    dropped: list[SiblingEvidence] = []
    if params.exclude_siblings:
        states, dropped, _ = sibling_aware_states(work, params)
    segments, major = segment_query(
        work, identity=params.identity, jump_rate=params.jump_rate, states=states
    )
    if major is None:
        return [], None, dropped

    # Pass 1: every non-major segment that clears the length filter and has a donor
    # genuinely closer than the major (favor_minor > favor_major on discordant
    # sites) is a candidate; collect its one-sided sign-test p-value.
    candidates = []
    for seg in segments:
        if seg.state == major or seg.msa_end - seg.msa_start < params.min_region:
            continue
        favor_minor, favor_major = discordant_counts(
            work.rows, work.query, major, seg.state, seg.msa_start, seg.msa_end
        )
        if favor_minor <= favor_major:
            continue
        candidates.append((seg, favor_minor, favor_major,
                           sign_test_pvalue(favor_minor, favor_major)))

    # Pass 2: a segment is reported when its Benjamini-Hochberg q-value clears alpha --
    # the multiplicity-adjusted gate, so the number of candidate segments scanned does not
    # inflate the false-call rate. (A low-divergence event between near-identical parents
    # has few distinguishing sites but, when real, still clears the FDR gate.)
    qvalues = benjamini_hochberg([c[3] for c in candidates])
    regions: list[Region] = []
    for (seg, favor_minor, favor_major, pvalue), q in zip(candidates, qvalues, strict=True):
        if q > params.alpha:
            continue
        support = favor_minor / (favor_minor + favor_major)
        # The segment is where the HMM path sat in the donor state, which is not the
        # same as where the donor is the better match: once in the donor state the path
        # has no reason to leave while donor and major are identical, and its ends are
        # only as fine as the window grid. The reported span is the part of the segment
        # that the distinguishing sites support (see _tract_span). The sign test above
        # stays on the segment, so which regions are called does not depend on this
        # step -- only their coordinates do.
        lo, hi = _tract_span(work, seg, major)
        # Identity over the reported columns, the same quantity the site-based callers
        # report, on the same scale under either windowing (the per-window values of
        # informative-site windowing are identity at polymorphic columns only).
        mean_minor = region_identity(work.rows, work.query, seg.state, lo, hi)
        mean_major = region_identity(work.rows, work.query, major, lo, hi)
        # The window count and the posterior describe the reported span too, not the
        # segment it was trimmed from: the segment's windows that overlap the span.
        n_windows, posterior = _span_windows(work, window_size, seg, lo, hi)
        regions.append(
            Region(
                minor_parent=seg.state, major_parent=major,
                msa_start=lo, msa_end=hi,
                query_start=work.column_to_query(lo), query_end=work.column_to_query(hi),
                n_windows=n_windows,
                mean_sim_minor=round(mean_minor, 4), mean_sim_major=round(mean_major, 4),
                margin=round(mean_minor - mean_major, 4),
                posterior_support=posterior,
                breakpoint_lo=seg.breakpoint_lo, breakpoint_hi=seg.breakpoint_hi,
                support=round(support, 3),
                pvalue=pvalue, qvalue=q,
                test="sign test on discordant sites (one-sided)",
                statistic="fraction of discordant sites favouring the donor",
                minor_cluster_size=cluster_size.get(seg.state, 1),
                major_cluster_size=cluster_size.get(major, 1),
            )
        )
    return regions, major, dropped


def _with_singletons(result: WindowSimilarity, pooled: list[list[str]]) -> list[list[str]]:
    """Every reference in one cluster: the multi-member ``pooled`` clusters plus each
    remaining label alone, in the order the labels first appear in ``result`` (the
    order ``cluster_references`` returns, so the HMM's state order is unchanged)."""
    index = {label: i for i, members in enumerate(pooled) for label in members}
    out: list[list[str]] = []
    seen: set[int] = set()
    for label in result.similarities:
        i = index.get(label)
        if i is None:
            out.append([label])
        elif i not in seen:
            seen.add(i)
            out.append(pooled[i])
    return out


def _span_windows(
    work: WindowSimilarity, window_size: int, seg: Segment, lo: int, hi: int
) -> tuple[int, float]:
    """``(n_windows, mean posterior)`` over the segment's windows overlapping ``[lo, hi)``.

    Falls back to the whole segment when no window of it overlaps the span (a span
    can only lie outside its segment's windows at an alignment end).
    """
    bounds = _window_bounds(work, window_size)
    idx = [
        i for i in range(seg.start_window, seg.end_window + 1)
        if bounds[i][0] < hi and bounds[i][1] > lo
    ]
    if not idx or not seg.window_posteriors:
        return seg.n_windows, seg.mean_posterior
    post = [seg.window_posteriors[i - seg.start_window] for i in idx]
    return len(idx), round(float(sum(post) / len(post)), 4)


def _tract_span(work: WindowSimilarity, seg: Segment, major: str) -> tuple[int, int]:
    """MSA columns ``[start, end)`` of an HMM segment that its donor actually explains.

    The span from the first to the last site where the query matches the donor and not
    the major. Columns where the query matches both parents or neither carry no
    information about which one it was copied from, so they cannot hold a region open:
    a run of windows in which donor and major are identical is left out, and so is the
    padding the window grid adds at either end.

    A segment that includes the first (last) window is searched from the alignment
    start (to its end): window centres stop half a window short of the ends, so a tract
    reaching an end would otherwise be reported from the first centre.

    Falls back to the segment's own span when no site in it favours the donor.
    """
    lo = 0 if seg.start_window == 0 else seg.msa_start
    hi = work.width if seg.end_window == len(work.positions) - 1 else seg.msa_end
    steps, cols = triplet_steps(work.rows, work.query, major, seg.state)
    donor_cols = cols[(steps == -1) & (cols >= lo) & (cols < hi)]
    if donor_cols.size == 0:
        return seg.msa_start, seg.msa_end
    return int(donor_cols[0]), int(donor_cols[-1]) + 1


def _call_regions_heuristic(
    result: WindowSimilarity,
    analysis: AnalysisResult,
    window_size: int,
    params: RegionParams,
) -> tuple[list[Region], str | None]:
    """Legacy greedy run-merge caller (kept for comparison)."""
    top = rank_datasets(analysis, 1)
    if not top or len(result.similarities) < 2:
        return [], (top[0] if top else None)
    major = top[0]
    half = window_size // 2

    def span(index: int) -> tuple[int, int]:
        center = result.positions[index]
        return max(0, center - half), min(result.width, center + half)

    # Build runs of consecutive recombinant windows.
    runs: list[_Run] = []
    for i in range(len(result.positions)):
        minor = _window_minor(result, major, i, params.margin)
        if minor is None:
            continue
        lo, hi = span(i)
        if runs and i == runs[-1].indices[-1] + 1:
            runs[-1].indices.append(i)
            runs[-1].minor_counts[minor] += 1
            runs[-1].msa_end = hi
        else:
            runs.append(_Run(Counter({minor: 1}), [i], lo, hi))

    # Merge adjacent runs sharing the dominant minor parent within merge_gap.
    merged: list[_Run] = []
    for run in runs:
        run_minor = run.minor_counts.most_common(1)[0][0]
        if merged:
            prev = merged[-1]
            prev_minor = prev.minor_counts.most_common(1)[0][0]
            if run_minor == prev_minor and run.msa_start - prev.msa_end <= params.merge_gap:
                prev.indices.extend(run.indices)
                prev.minor_counts.update(run.minor_counts)
                prev.msa_end = run.msa_end
                continue
        merged.append(run)

    regions: list[Region] = []
    for run in merged:
        length = run.msa_end - run.msa_start
        if length < params.min_region:
            continue
        minor = run.minor_counts.most_common(1)[0][0]
        minor_sims = [result.similarities[minor][i] for i in run.indices
                      if not isnan(result.similarities[minor][i])]
        major_sims = [result.similarities[major][i] for i in run.indices
                      if not isnan(result.similarities[major][i])]
        mean_minor = mean(minor_sims) if minor_sims else float("nan")
        mean_major = mean(major_sims) if major_sims else float("nan")
        regions.append(
            Region(
                minor_parent=minor,
                major_parent=major,
                msa_start=run.msa_start,
                msa_end=run.msa_end,
                query_start=result.column_to_query(run.msa_start),
                query_end=result.column_to_query(run.msa_end),
                n_windows=len(run.indices),
                mean_sim_minor=round(mean_minor, 4),
                mean_sim_major=round(mean_major, 4),
                margin=round(mean_minor - mean_major, 4),
            )
        )
    return regions, major
