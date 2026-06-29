"""Reference-coverage diagnostics: where might the correct reference be missing?

Tessera's region calling is purely relative -- it always crowns the closest
reference, even when *every* reference is far from the query. When the query's
true parent or donor is not in the collection, the giveaway is a low **best**
similarity: the closest reference can still only manage, say, 0.82 identity.

This module flags those stretches. For each window it uses the per-window
``best_sim`` (closest reference) already computed in
:mod:`tessera.recomb.similarity`, and reports contiguous **coverage gaps** where
even the best reference falls below a threshold. The threshold is data-adaptive
by default (relative to the query's own typical best match), so it does not flag
everything in a genuinely divergent collection; an absolute ``floor`` overrides it.

A gap is labelled ``divergent`` (ample comparable bases, the query really is far
from all references -- a likely missing reference) or ``low_information`` (few
comparable bases, so the gap is uncertain rather than informative).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isnan
from statistics import mean, median

import numpy as np

from .regions import Region
from .similarity import WindowSimilarity


@dataclass
class CoverageParams:
    """Thresholds for flagging under-covered windows. Defaults derive from window size."""

    floor: float | None  # absolute best-similarity floor; None -> adaptive baseline
    rel_drop: float  # adaptive: flag below (90th percentile of best_sim - rel_drop)
    min_gap: int  # drop gaps shorter than this many bp
    merge_gap: int  # merge gaps separated by at most this many bp

    @classmethod
    def with_defaults(
        cls,
        window_size: int,
        floor: float | None = None,
        rel_drop: float = 0.05,
        min_gap: int | None = None,
        merge_gap: int | None = None,
    ) -> CoverageParams:
        return cls(
            floor=floor,
            rel_drop=rel_drop,
            min_gap=window_size if min_gap is None else min_gap,
            merge_gap=window_size if merge_gap is None else merge_gap,
        )


@dataclass
class CoverageGap:
    """A stretch where even the closest reference matches the query poorly."""

    msa_start: int
    msa_end: int
    query_start: int
    query_end: int
    length_bp: int
    n_windows: int
    best_label: str  # the closest reference across the gap (still a poor match)
    mean_best: float  # its mean similarity over the gap
    kind: str  # "divergent" | "low_information"


def coverage_threshold(result: WindowSimilarity, params: CoverageParams) -> float:
    """The best-similarity threshold below which a window is under-covered."""
    if params.floor is not None:
        return params.floor
    finite = [b for b in result.best_sim if not isnan(b)]
    if not finite:
        return 0.0
    baseline = float(np.percentile(finite, 90)) - params.rel_drop
    return max(0.0, min(1.0, baseline))


def call_coverage_gaps(
    result: WindowSimilarity, window_size: int, params: CoverageParams
) -> tuple[list[CoverageGap], float]:
    """Return ``(gaps, threshold)`` for windows whose best reference is poor."""
    threshold = coverage_threshold(result, params)
    best, info = result.best_sim, result.informative
    half = window_size // 2

    def span(index: int) -> tuple[int, int]:
        center = result.positions[index]
        return max(0, center - half), min(result.width, center + half)

    # Runs of consecutive under-covered windows.
    runs: list[list[int]] = []
    for i in range(len(result.positions)):
        if isnan(best[i]) or best[i] >= threshold:
            continue
        if runs and i == runs[-1][-1] + 1:
            runs[-1].append(i)
        else:
            runs.append([i])

    # Merge runs separated by at most merge_gap bp.
    merged: list[dict] = []
    for run in runs:
        lo, _ = span(run[0])
        _, hi = span(run[-1])
        if merged and lo - merged[-1]["end"] <= params.merge_gap:
            merged[-1]["idx"].extend(run)
            merged[-1]["end"] = hi
        else:
            merged.append({"start": lo, "end": hi, "idx": list(run)})

    low_info = max(20, int(0.05 * window_size))
    gaps: list[CoverageGap] = []
    for run in merged:
        length = run["end"] - run["start"]
        if length < params.min_gap:
            continue
        idx = run["idx"]
        labels = [result.best_label[i] for i in idx if result.best_label[i]]
        best_label = Counter(labels).most_common(1)[0][0] if labels else "n/a"
        sims = [best[i] for i in idx if not isnan(best[i])]
        med_info = median(info[i] for i in idx) if idx else 0
        gaps.append(
            CoverageGap(
                msa_start=run["start"], msa_end=run["end"],
                query_start=result.column_to_query(run["start"]),
                query_end=result.column_to_query(run["end"]),
                length_bp=length, n_windows=len(idx),
                best_label=best_label,
                mean_best=round(mean(sims), 4) if sims else float("nan"),
                kind="low_information" if med_info < low_info else "divergent",
            )
        )
    return gaps, threshold


def flag_undercovered_regions(regions: list[Region], threshold: float) -> None:
    """Mark called regions whose donor is itself a poor match (possible stand-in)."""
    for r in regions:
        r.donor_undercovered = (
            not isnan(r.mean_sim_minor) and r.mean_sim_minor < threshold
        )


def gaps_as_regions(
    gaps: list[CoverageGap], result: WindowSimilarity, major: str | None
) -> list[Region]:
    """Turn ``divergent`` coverage gaps into 'donor-absent' putative recombinations.

    A gap where the query is locally far from *every* reference -- including its
    own major parent -- is a recombination signature whose donor is missing from
    the collection (e.g. the HIV CRF01_AE env, which has no pure subtype-E
    reference). ``low_information`` gaps carry too little evidence and are skipped.
    """
    out: list[Region] = []
    for g in gaps:
        if g.kind != "divergent":
            continue
        idx = [i for i, pos in enumerate(result.positions)
               if g.msa_start <= pos <= g.msa_end]
        major_sims = (
            [result.similarities[major][i] for i in idx
             if not isnan(result.similarities[major][i])]
            if major and major in result.similarities else []
        )
        mean_major = round(mean(major_sims), 4) if major_sims else float("nan")
        margin = round(g.mean_best - mean_major, 4) if not isnan(mean_major) else float("nan")
        out.append(Region(
            minor_parent=g.best_label, major_parent=major or "n/a",
            msa_start=g.msa_start, msa_end=g.msa_end,
            query_start=g.query_start, query_end=g.query_end,
            length_bp=g.length_bp, n_windows=g.n_windows,
            mean_sim_minor=g.mean_best, mean_sim_major=mean_major, margin=margin,
            donor_undercovered=True, donor_absent=True,
        ))
    return out
