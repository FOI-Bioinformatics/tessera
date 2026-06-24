"""Summarise per-window similarities into winners and per-dataset statistics.

For every window the closest reference (or references, on a tie) is the
"winner"; datasets are also ranked by median similarity and by how many windows
exceed a set of identity thresholds.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import isnan
from statistics import median

from .similarity import WindowSimilarity

# Similarity thresholds reported per dataset (fraction identical). The first (0)
# counts windows with any similarity; the rest are the high-identity buckets.
DEFAULT_THRESHOLDS: tuple[float, ...] = (0.0, 0.99, 0.95, 0.90, 0.80, 0.70)

STATS_HEADER: list[str] = [
    "Dataset",
    "Tot windows",
    "Median similarity",
    "Windows >0",
    "Windows >99%",
    "Windows >95%",
    "Windows >90%",
    "Windows >80%",
    "Windows >70%",
]


@dataclass
class AnalysisResult:
    winners_no_ties: dict[str, int]
    winners_with_ties: dict[str, int]
    stats: dict[str, list] = field(default_factory=dict)
    stats_header: list[str] = field(default_factory=lambda: list(STATS_HEADER))


def winners_per_window(result: WindowSimilarity) -> list[list[str]]:
    """Return, for each window, the list of datasets tied for highest similarity."""
    datasets = list(result.similarities)
    out: list[list[str]] = []
    for i in range(len(result.positions)):
        pos_sims = {
            ds: result.similarities[ds][i]
            for ds in datasets
            if not isnan(result.similarities[ds][i])
        }
        if not pos_sims:
            out.append([])  # uninformative window (no comparable position)
            continue
        best = max(pos_sims.values())
        out.append([ds for ds, sim in pos_sims.items() if sim == best])
    return out


def winner_label(winners: list[str]) -> str:
    """A single cell value for a window's winner(s); ties joined with '|'."""
    return "|".join(sorted(winners)) if winners else ""


def _tally_winners(
    per_window: list[list[str]],
) -> tuple[dict[str, int], dict[str, int]]:
    """Count, per dataset, the windows it is closest in (without and with ties)."""
    no_ties: dict[str, int] = {}
    with_ties: dict[str, int] = {}
    for winners in per_window:
        if not winners:
            continue
        if len(winners) == 1:
            no_ties[winners[0]] = no_ties.get(winners[0], 0) + 1
        for ds in winners:
            with_ties[ds] = with_ties.get(ds, 0) + 1
    return no_ties, with_ties


def _dataset_stats(result: WindowSimilarity, thresholds: tuple[float, ...]) -> dict[str, list]:
    """Median similarity plus per-threshold window counts for each dataset."""
    stats: dict[str, list] = {}
    for dataset, sims in result.similarities.items():
        finite = [s for s in sims if not isnan(s)]
        tot_windows = len(finite)
        sim_median = round(median(finite), 2) if finite else 0.0
        counts = [sum(1 for s in finite if s > t) for t in thresholds]
        stats[dataset] = [tot_windows, sim_median, *counts]
    return stats


def analyze(
    result: WindowSimilarity, thresholds: tuple[float, ...] = DEFAULT_THRESHOLDS
) -> AnalysisResult:
    """Tally window winners and per-dataset statistics for ``result``."""
    per_window = winners_per_window(result)
    no_ties, with_ties = _tally_winners(per_window)
    stats = _dataset_stats(result, thresholds)
    return AnalysisResult(
        winners_no_ties=no_ties,
        winners_with_ties=with_ties,
        stats=stats,
    )


def rank_by_wins(winners_with_ties: dict[str, int], top_n: int) -> list[str]:
    """Return up to ``top_n`` dataset labels ordered by window wins (with ties)."""
    ranked = sorted(winners_with_ties.items(), key=lambda kv: kv[1], reverse=True)
    return [dataset for dataset, _ in ranked[:top_n]]


def rank_datasets(analysis: AnalysisResult, top_n: int) -> list[str]:
    """Return up to ``top_n`` datasets, by window wins then by high-identity stats.

    Window winners come first (the datasets the query is closest to most often);
    when fewer than ``top_n`` datasets won any window, the list is filled from the
    per-dataset statistics ranking so the plots still show several candidates.
    """
    ordered = rank_by_wins(analysis.winners_with_ties, top_n)
    if len(ordered) >= top_n:
        return ordered
    seen = set(ordered)
    by_stats = sorted(analysis.stats.items(), key=stats_sort_key, reverse=True)
    for dataset, _ in by_stats:
        if len(ordered) >= top_n:
            break
        if dataset not in seen:
            ordered.append(dataset)
            seen.add(dataset)
    return ordered


def stats_sort_key(item: tuple[str, list]) -> tuple:
    """Sort key for the stats table: by the >99/95/90/80 buckets in order."""
    values = item[1]
    # values = [tot_windows, median, >0, >99, >95, >90, >80, >70]
    return (values[3], values[4], values[5], values[6])
