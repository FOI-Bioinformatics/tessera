"""Recombination-region calling.

The query is most similar to one dataset across the genome -- the **major
parent** (the backbone donor). A recombination is indicated where a stretch of
the query is instead closest to another dataset (a **minor parent**).

The default caller (``method="hmm"``) segments the query against the reference
panel with a hidden Markov model (:mod:`recomfi.recomb.hmm`): emissions are the
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

from .analyze import AnalysisResult, rank_datasets
from .hmm import DEFAULT_JUMP_RATE, segment_query
from .siblings import SiblingEvidence, sibling_aware_states
from .similarity import WindowSimilarity, discordant_counts
from .stats import benjamini_hochberg, sign_test_pvalue


def _signif(x: float) -> float:
    """Round a (possibly tiny) p/q-value to two significant figures for display."""
    return 0.0 if x == 0 else float(f"{x:.2g}")


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
    # the event. See :mod:`recomfi.recomb.siblings`. Defaults from first principles.
    exclude_siblings: bool = True
    sibling_coverage: float = 0.90  # a sibling is comparable over >= this of the query
    sibling_beaten: float = 0.15  # ... significantly beaten in <= this fraction of windows
    sibling_alpha: float = 0.05  # significance level for the per-window beaten test
    # ... and never beaten over a region-sized contiguous run (the primary signal; a
    # parent is beaten over the co-parent's tract). Derived from min_region / step.

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
        )


@dataclass
class Region:
    minor_parent: str
    major_parent: str
    msa_start: int
    msa_end: int
    query_start: int
    query_end: int
    length_bp: int
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
    support: float | None = None
    pvalue: float | None = None
    qvalue: float | None = None
    # Set by coverage analysis: the donor's own similarity is poor, so it may be a
    # stand-in for an absent true donor. Default False (not yet evaluated).
    donor_undercovered: bool = False
    # A putative recombination where the query is locally far from EVERY reference:
    # the true donor is likely missing from the collection (bridged from a gap).
    donor_absent: bool = False


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


def call_regions(
    result: WindowSimilarity,
    analysis: AnalysisResult,
    window_size: int,
    params: RegionParams,
) -> tuple[list[Region], str | None, list[SiblingEvidence]]:
    """Call recombinant regions.

    Returns ``(regions, major_parent, excluded_siblings)`` -- the third element lists
    the query's whole-genome siblings set aside before segmenting (HMM caller only).
    """
    if params.method == "hmm":
        return _call_regions_hmm(result, window_size, params)
    regions, major = _call_regions_heuristic(result, analysis, window_size, params)
    return regions, major, []


def _call_regions_hmm(
    result: WindowSimilarity, window_size: int, params: RegionParams
) -> tuple[list[Region], str | None, list[SiblingEvidence]]:
    """HMM segmentation; a non-major segment is reported only when its donor is a
    significantly better match than the major on the sites distinguishing them.

    Whole-genome siblings of the query (its own lineage) are excluded from the
    competition first, so a sibling cannot win every window and mask the event.
    """
    if len(result.similarities) < 2:
        only = next(iter(result.similarities), None)
        return [], only, []
    states: list[str] | None = None
    dropped: list[SiblingEvidence] = []
    if params.exclude_siblings:
        states, dropped, _ = sibling_aware_states(result, params)
    segments, major = segment_query(
        result, identity=params.identity, jump_rate=params.jump_rate, states=states
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
            result.rows, result.query, major, seg.state, seg.msa_start, seg.msa_end
        )
        if favor_minor <= favor_major:
            continue
        candidates.append((seg, favor_minor, favor_major,
                           sign_test_pvalue(favor_minor, favor_major)))

    # Pass 2: a segment is reported when it is individually significant
    # (p <= alpha); the Benjamini-Hochberg q-value is computed across all candidates
    # and reported alongside, so the multiplicity-adjusted significance is visible
    # without discarding marginal-but-real events (e.g. a recombination between
    # near-identical parents, which has few distinguishing sites).
    qvalues = benjamini_hochberg([c[3] for c in candidates])
    regions: list[Region] = []
    for (seg, favor_minor, favor_major, pvalue), q in zip(candidates, qvalues, strict=True):
        if pvalue > params.alpha:
            continue
        support = favor_minor / (favor_minor + favor_major)
        idx = range(seg.start_window, seg.end_window + 1)
        minor_sims = [result.similarities[seg.state][i] for i in idx
                      if not isnan(result.similarities[seg.state][i])]
        major_sims = [result.similarities[major][i] for i in idx
                      if not isnan(result.similarities[major][i])]
        mean_minor = mean(minor_sims) if minor_sims else float("nan")
        mean_major = mean(major_sims) if major_sims else float("nan")
        regions.append(
            Region(
                minor_parent=seg.state, major_parent=major,
                msa_start=seg.msa_start, msa_end=seg.msa_end,
                query_start=seg.query_start, query_end=seg.query_end,
                length_bp=seg.msa_end - seg.msa_start, n_windows=seg.n_windows,
                mean_sim_minor=round(mean_minor, 4), mean_sim_major=round(mean_major, 4),
                margin=round(mean_minor - mean_major, 4),
                posterior_support=seg.mean_posterior,
                breakpoint_lo=seg.breakpoint_lo, breakpoint_hi=seg.breakpoint_hi,
                support=round(support, 3),
                pvalue=_signif(pvalue), qvalue=_signif(q),
            )
        )
    return regions, major, dropped


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
                length_bp=length,
                n_windows=len(run.indices),
                mean_sim_minor=round(mean_minor, 4),
                mean_sim_major=round(mean_major, 4),
                margin=round(mean_minor - mean_major, 4),
            )
        )
    return regions, major
