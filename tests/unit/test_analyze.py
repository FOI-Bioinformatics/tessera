"""Winner tallies and per-dataset statistics."""

from __future__ import annotations

import numpy as np

from tessera.recomb.analyze import analyze, rank_datasets, winner_label, winners_per_window
from tessera.recomb.similarity import WindowSimilarity


def _result(similarities: dict[str, list[float]]) -> WindowSimilarity:
    n = len(next(iter(similarities.values())))
    return WindowSimilarity(
        positions=list(range(n)),
        query_positions=list(range(n)),
        similarities=similarities,
        query="query",
        width=n,
        query_cumulative=np.arange(n + 1),
    )


def _example() -> WindowSimilarity:
    # refA wins windows 0 and 2 outright; window 1 is a refA/refB tie.
    return _result({"refA": [1.0, 0.9, 1.0], "refB": [0.8, 0.9, 0.7]})


def test_winners_with_and_without_ties() -> None:
    analysis = analyze(_example())
    assert analysis.winners_with_ties == {"refA": 3, "refB": 1}
    assert analysis.winners_no_ties == {"refA": 2}


def test_winners_per_window_and_label() -> None:
    per_window = winners_per_window(_example())
    assert per_window == [["refA"], ["refA", "refB"], ["refA"]]
    assert winner_label(per_window[1]) == "refA|refB"


def test_stats_shape_and_median() -> None:
    analysis = analyze(_example())
    stats = analysis.stats["refA"]
    # [tot_windows, median, >0, >99, >95, >90, >80, >70]
    assert stats[0] == 3
    assert stats[1] == 1.0
    assert stats[2] == 3
    assert "Median similarity" in analysis.stats_header


def test_rank_datasets_fills_from_stats() -> None:
    # Only refA ever wins a window, but top_n=2 should still surface refB.
    analysis = analyze(_result({"refA": [1.0, 1.0], "refB": [0.9, 0.9]}))
    ranked = rank_datasets(analysis, top_n=2)
    assert ranked[0] == "refA"
    assert set(ranked) == {"refA", "refB"}


def test_nan_windows_ignored_in_winners_and_stats() -> None:
    nan = float("nan")
    # window 1 is uninformative (nan for both): no winner, excluded from stats.
    analysis = analyze(_result({"refA": [1.0, nan, 0.8], "refB": [0.9, nan, 0.7]}))
    per_window = winners_per_window(_result({"refA": [1.0, nan, 0.8], "refB": [0.9, nan, 0.7]}))
    assert per_window[1] == []  # no winner in the uninformative window
    assert analysis.winners_with_ties == {"refA": 2}
    # tot_windows counts only the two finite windows
    assert analysis.stats["refA"][0] == 2
    assert analysis.stats["refB"][0] == 2
