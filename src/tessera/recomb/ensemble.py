"""Combine several callers' regions into one transparent consensus.

The HMM and 3SEQ callers are complementary -- 3SEQ recovers short low-divergence tracts
the HMM segmentation dilutes, the HMM localizes breakpoints better on divergent parents.
Running them together raises recall (the union of what each finds) and, where they agree,
confidence. This module merges the per-method region sets into one annotated list, the
ensemble's output.

The merge is deliberately transparent: there is no fabricated combined score. Two regions
are *the same event* when they overlap in query coordinates and name the same minor
parent; a merged region records exactly **which** callers found it (``Region.methods``)
and whether the parent-free Hudson-Kaplan Rmin signal corroborates it
(``Region.parent_free_support``). Overlapping regions that name different minor parents
are a genuine disagreement and stay separate.
"""

from __future__ import annotations

from .clusters import _window_bounds
from .regions import CALLERS, Region
from .similarity import WindowSimilarity, region_identity
from .typing import LineageMap, lineage_of


def _rank(method: str) -> int:
    return CALLERS.index(method) if method in CALLERS else len(CALLERS)


def _overlap(a: Region, b: Region) -> bool:
    """Do two regions overlap in query coordinates?"""
    return a.query_start < b.query_end and b.query_start < a.query_end


def _same_donor(
    a: Region, b: Region, lineage_map: LineageMap | None,
    cluster_of: dict[str, int] | None = None,
) -> bool:
    """Do two regions name the same donor? Equal genome labels; or members of one
    lineage cluster (near-duplicates the HMM pools; see ``clusters.lineage_clusters``);
    or -- when references are typed -- the same lineage. Two callers that pick different
    genomes of one lineage then still count as the same event (and as agreement)."""
    if a.minor_parent == b.minor_parent:
        return True
    if cluster_of:
        ca = cluster_of.get(a.minor_parent)
        if ca is not None and ca == cluster_of.get(b.minor_parent):
            return True
    if lineage_map:
        la = lineage_of(a.minor_parent, lineage_map)
        return la is not None and la == lineage_of(b.minor_parent, lineage_map)
    return False


def reconcile_major(
    majors: dict[str, str | None],
    *,
    window_wins: dict[str, int] | None = None,
) -> tuple[str | None, dict[str, str]]:
    """The canonical backbone (major parent) plus each method's own major.

    ``majors`` maps each method to the major parent its caller reported -- which exists
    even when the caller found no recombinant region (the no-recombination case), so the
    backbone is always available for the verdict. Callers usually agree on the backbone
    (the genome-wide closest reference); when they differ the HMM's major wins if HMM ran,
    else the most frequently reported one. The per-method mapping is returned so the
    report can surface a disagreement.

    ``window_wins`` (label -> windows won) guards one failure mode: the HMM excludes the
    query's whole-genome siblings before segmenting, and for a very short recombinant
    tract the true backbone is ~identical to the query and can be dropped as a sibling,
    leaving the HMM with a major that wins no window. When that happens and the windowed
    vote has a clear winner the HMM's major is not, the winner is trusted instead -- the
    genome-wide closest reference is the backbone, and a zero-win HMM major is an artefact
    of the exclusion. A HMM major that wins any window (the usual case, and the HIV/RSV
    masking-twin case where the true parent wins its own segment) is left untouched.
    """
    per = {m: major for m, major in majors.items() if major is not None}
    if not per:
        return None, per
    if "hmm" in per:
        hmm_major = per["hmm"]
        if window_wins and window_wins.get(hmm_major, 0) == 0:
            top = max(window_wins, key=lambda label: window_wins[label])
            if window_wins.get(top, 0) > 0 and top != hmm_major:
                return top, per
        return hmm_major, per
    values = list(per.values())
    canonical = max(values, key=lambda m: (values.count(m), -values.index(m)))
    return canonical, per


def is_inverted_backbone(
    region: Region, *, backbone: str | None, window_wins: dict[str, int] | None
) -> bool:
    """Does ``region`` rest on a backbone the run has rejected? (issue #86)

    ``reconcile_major`` above already recognises one artefact: a caller segments on a
    view that can disagree with the genome-wide evidence, and a major that is the closest
    reference in not a single base-pair window is an artefact of that view rather than
    the query's backbone. It corrects the run's backbone and stops there, which left the
    regions the caller had already built on the rejected major in the table. The report
    could then name one backbone in its headline and the opposite in its rows, casting
    the genome-wide closest reference as a donor into a reference closest nowhere -- the
    mosaic read the wrong way round.

    Both conditions are required:

    * the reconciled ``backbone`` is the region's *minor* parent, so the row asserts the
      inverse of what the run concluded, and
    * the region's major parent won no window in the base-pair vote, the same condition
      ``reconcile_major`` used to reject it.

    The second ties the guard to that artefact. Callers may disagree about the backbone
    for real reasons, and a disagreement backed by window evidence of its own is left
    alone rather than overruled by the headline. An empty vote matches nothing, since it
    cannot establish that a major won no window.
    """
    if backbone is None or not window_wins:
        return False
    return (
        region.minor_parent == backbone
        and window_wins.get(region.major_parent, 0) == 0
    )


def withhold_inverted_backbone(
    regions: list[Region],
    breakdown: list[dict],
    *,
    backbone: str | None,
    window_wins: dict[str, int] | None,
) -> tuple[list[Region], list[dict], list[Region]]:
    """Drop the regions ``is_inverted_backbone`` matches, keeping ``breakdown`` in step.

    ``regions`` and ``breakdown`` are parallel, as in ``filter_by_length``. Returns the
    kept regions, their breakdown rows, and the withheld regions so the run can report
    them rather than let the table come up short for no stated reason.
    """
    kept: list[tuple[Region, dict]] = []
    withheld: list[Region] = []
    for region, row in zip(regions, breakdown, strict=True):
        if is_inverted_backbone(region, backbone=backbone, window_wins=window_wins):
            withheld.append(region)
        else:
            kept.append((region, row))
    return [r for r, _ in kept], [b for _, b in kept], withheld


def _best(group: list[Region]) -> Region:
    """The most significant member of a group (smallest q-value, then highest support)."""
    return min(group, key=lambda r: (r.qvalue if r.qvalue is not None else 1.0,
                                      -(r.support if r.support is not None else 0.0)))


def _merge_group(
    group: list[Region], major: str | None,
    scan: tuple[WindowSimilarity, int] | None = None,
) -> Region:
    """One consensus region from a set of overlapping same-minor-parent regions.

    With ``scan`` (the base-pair similarity result and its window size) the merged row's
    identities and window count are restated over the merged span; without it they are
    the most significant member's, as before.
    """
    best = _best(group)
    methods = tuple(sorted({m for r in group for m in r.methods}, key=_rank))
    q_start = min(r.query_start for r in group)
    q_end = max(r.query_end for r in group)
    # Prefer the HMM member's breakpoint (best localization) when present.
    hmm_member = next((r for r in group if "hmm" in r.methods), None)
    bp_lo = (hmm_member or best).breakpoint_lo
    bp_hi = (hmm_member or best).breakpoint_hi
    # The canonical backbone replaces each caller's own, so the table reads against one
    # backbone -- except when it IS this region's donor. That happens when the callers
    # disagree on the backbone (the HMM set a sibling aside, a site caller did not): the
    # region was called against the other backbone, and relabelling it would report a
    # genome as a donor into itself. Keep the pair the caller actually tested.
    relabel = bool(major) and major != best.minor_parent
    major_label = major if relabel and major else best.major_parent
    msa_start = min(r.msa_start for r in group)
    msa_end = max(r.msa_end for r in group)
    n_windows = max(r.n_windows for r in group)
    sim_minor, sim_major, margin = best.mean_sim_minor, best.mean_sim_major, best.margin
    if scan is not None:
        # The span is the union of the members, so the identities must be measured over
        # that union: the best member's own values describe its own (often much shorter)
        # span. The window count is restated in one unit -- base-pair scan windows
        # overlapping the span -- because the callers count different things (windows,
        # discriminating sites, run length) and a maximum across them has no unit.
        result, window_size = scan
        q = result.query
        raw_minor = region_identity(result.rows, q, best.minor_parent, msa_start, msa_end)
        raw_major = region_identity(result.rows, q, major_label, msa_start, msa_end)
        sim_minor, sim_major = round(raw_minor, 4), round(raw_major, 4)
        margin = round(raw_minor - raw_major, 4)
        n_windows = sum(
            1 for lo, hi in _window_bounds(result, window_size) if lo < msa_end and hi > msa_start
        )
    return Region(
        minor_parent=best.minor_parent,
        major_parent=major_label,
        msa_start=msa_start,
        msa_end=msa_end,
        query_start=q_start,
        query_end=q_end,
        n_windows=n_windows,
        mean_sim_minor=sim_minor,
        mean_sim_major=sim_major,
        margin=margin,
        # The posterior is an HMM quantity: the HMM member's where there is one (as for
        # the breakpoint), absent otherwise. `best` is often a site caller, whose
        # value here would be the dataclass default rather than a measurement.
        posterior_support=(hmm_member or best).posterior_support,
        breakpoint_lo=bp_lo,
        breakpoint_hi=bp_hi,
        support=best.support,
        pvalue=best.pvalue,
        qvalue=best.qvalue,
        donor_undercovered=any(r.donor_undercovered for r in group),
        major_degraded=any(r.major_degraded for r in group),
        donor_absent=any(r.donor_absent for r in group),
        minor_cluster_size=best.minor_cluster_size,
        major_cluster_size=best.major_cluster_size,
        methods=methods,
        # The p/q kept above are `best`'s, so its test names travel with them. Without
        # this the row would list every agreeing caller beside a p-value from just one
        # of them, and a reader could not tell which.
        test=best.test,
        statistic=best.statistic,
    )


def _group(
    regions: list[Region], lineage_map: LineageMap | None,
    cluster_of: dict[str, int] | None = None,
) -> list[list[Region]]:
    """Group regions that overlap (transitively) and name the same donor (genome,
    lineage cluster or, when typed, lineage)."""
    groups: list[list[Region]] = []
    for region in sorted(regions, key=lambda r: (r.minor_parent, r.query_start)):
        hits = [
            group for group in groups
            if _same_donor(group[0], region, lineage_map, cluster_of)
            and any(_overlap(region, member) for member in group)
        ]
        if not hits:
            groups.append([region])
            continue
        # A region can bridge groups that did not overlap each other -- on a typed panel
        # the sort is by genome label, so two labels of one lineage arrive out of
        # coordinate order. Fold every group it touches into the first, or the same
        # donor would be reported as two overlapping regions.
        target = hits[0]
        target.append(region)
        for other in hits[1:]:
            target.extend(other)
        groups = [g for g in groups if all(g is not other for other in hits[1:])]
    return groups


def filter_by_agreement(
    regions: list[Region], breakdown: list[dict], min_methods: int
) -> tuple[list[Region], list[dict], int]:
    """Keep only regions at least ``min_methods`` callers found; returns the kept
    regions, their breakdown rows and how many were suppressed.

    The union of several callers raises recall but inherits every caller's false
    positives, so a region resting on one caller alone is the weakest thing the
    ensemble can report. RDP5 ships the same control ("list events detected by more
    than N methods"). ``regions`` and ``breakdown`` are built in one pass and sorted
    alike, so they are filtered together to stay in step.
    """
    if min_methods <= 1:
        return regions, breakdown, 0
    kept = [
        (region, row)
        for region, row in zip(regions, breakdown, strict=True)
        if len(region.methods) >= min_methods
    ]
    return [r for r, _ in kept], [b for _, b in kept], len(regions) - len(kept)


def filter_by_length(
    regions: list[Region], breakdown: list[dict], min_length: int
) -> tuple[list[Region], list[dict], int]:
    """Keep only regions at least ``min_length`` query bases long (``--min-region``).

    Applied to the merged regions, so it bounds what is reported whichever caller found
    a region; ``regions`` and ``breakdown`` are parallel and filtered together. Returns
    the kept regions, their breakdown rows and how many were removed.
    """
    kept = [
        (region, row) for region, row in zip(regions, breakdown, strict=True)
        if region.length_bp >= min_length
    ]
    return [r for r, _ in kept], [b for _, b in kept], len(regions) - len(kept)


def consensus_regions(
    per_method: dict[str, list[Region]],
    *,
    major: str | None,
    rmin_intervals: list[tuple[int, int]] | None = None,
    lineage_map: LineageMap | None = None,
    scan: tuple[WindowSimilarity, int] | None = None,
    clusters: list[list[str]] | None = None,
) -> tuple[list[Region], list[dict]]:
    """Merge per-method regions into one consensus list plus a per-region breakdown.

    ``rmin_intervals`` are the Hudson-Kaplan recombination intervals in query
    coordinates; a consensus region overlapping one is flagged ``parent_free_support``.
    When ``lineage_map`` is given, two callers' regions over the same span merge if they
    name the same *lineage* (not only the same genome), so picking different
    representative genomes of one lineage still counts as agreement. The breakdown lists,
    per consensus region, which callers found it and their support, for the report's
    "method comparison" section.

    ``scan`` is ``(base-pair result, window size)``. When given, each merged row's
    ``mean_sim_minor`` / ``mean_sim_major`` / ``margin`` are identity over the merged span
    and ``n_windows`` counts the base-pair windows overlapping it.

    ``clusters`` are the run's lineage clusters (``clusters.lineage_clusters``). Regions
    naming two members of one cluster are the same event: the HMM reports a cluster under
    its representative, while 3SEQ and MaxChi rank individual genomes and may name
    another member.
    """
    all_regions = [r for regions in per_method.values() for r in regions]
    intervals = rmin_intervals or []
    consensus: list[Region] = []
    breakdown: list[dict] = []
    cluster_of = {label: i for i, members in enumerate(clusters or []) for label in members}
    for group in _group(all_regions, lineage_map, cluster_of):
        region = _merge_group(group, major, scan)
        region.parent_free_support = any(
            region.query_start < hi and lo < region.query_end for lo, hi in intervals
        )
        consensus.append(region)
        per: dict[str, float | None] = {}
        for member in group:
            for m in member.methods:
                per[m] = member.support
        breakdown.append({
            "minor_parent": region.minor_parent,
            "query_start": region.query_start,
            "query_end": region.query_end,
            "methods": region.methods,
            "per_method_support": per,
            "parent_free_support": region.parent_free_support,
        })
    consensus.sort(key=lambda r: r.query_start)
    breakdown.sort(key=lambda b: b["query_start"])
    return consensus, breakdown
