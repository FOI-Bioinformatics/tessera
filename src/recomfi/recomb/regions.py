"""Heuristic recombination-region calling.

The query is most similar to one dataset across the genome -- the **major
parent** (the backbone donor). A recombination is indicated where a stretch of
the query is instead closest to another dataset (a **minor parent**). This
module turns the per-window similarities into called regions.

The heuristic is deliberately simple and transparent -- it is an indicative
screen, not a statistical significance test (e.g. 3SEQ/RDP):

1. The major parent is the dataset winning the most windows overall.
2. A window is *recombinant* when its closest dataset is some minor parent whose
   similarity exceeds the major parent's by at least ``margin``.
3. Consecutive recombinant windows form a run; runs of the same minor parent
   separated by at most ``merge_gap`` bases are merged; runs shorter than
   ``min_region`` bases are dropped.

Region coordinates are reported in both MSA columns and query bases.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isnan
from statistics import mean

from .analyze import AnalysisResult, rank_datasets
from .similarity import WindowSimilarity


@dataclass
class RegionParams:
    """Region-calling thresholds. Defaults derive from the window size."""

    min_region: int
    margin: float
    merge_gap: int

    @classmethod
    def with_defaults(
        cls,
        window_size: int,
        min_region: int | None = None,
        margin: float = 0.0,
        merge_gap: int | None = None,
    ) -> RegionParams:
        return cls(
            min_region=window_size if min_region is None else min_region,
            margin=margin,
            merge_gap=window_size if merge_gap is None else merge_gap,
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
    # Set by coverage analysis: the donor's own similarity is poor, so it may be a
    # stand-in for an absent true donor. Default False (not yet evaluated).
    donor_undercovered: bool = False


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
) -> tuple[list[Region], str | None]:
    """Call recombinant regions. Returns ``(regions, major_parent)``."""
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
