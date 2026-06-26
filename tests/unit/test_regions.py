"""Heuristic recombination-region calling (the legacy method="heuristic" path)."""

from __future__ import annotations

import numpy as np

from recomfi.recomb.analyze import analyze
from recomfi.recomb.regions import RegionParams, call_regions
from recomfi.recomb.similarity import WindowSimilarity


def _result(similarities: dict[str, list[float]], step: int = 1) -> WindowSimilarity:
    n = len(next(iter(similarities.values())))
    # Place window centers at 0, step, 2*step, ... and map query 1:1 (no gaps).
    positions = [i * step for i in range(n)]
    width = positions[-1] + 1 if positions else 0
    return WindowSimilarity(
        positions=positions,
        query_positions=list(positions),
        similarities=similarities,
        query="query",
        width=width,
        query_cumulative=np.arange(width + 1),
    )


def test_calls_region_for_minor_parent_run() -> None:
    # major (cowpox) wins everywhere except a central run where variola wins.
    cowpox = [0.99, 0.99, 0.80, 0.80, 0.80, 0.99, 0.99]
    variola = [0.90, 0.90, 0.99, 0.99, 0.99, 0.90, 0.90]
    result = _result({"cowpox": cowpox, "variola": variola}, step=100)
    analysis = analyze(result)
    params = RegionParams.with_defaults(window_size=100, min_region=100, method="heuristic")
    regions, major = call_regions(result, analysis, window_size=100, params=params)

    assert major == "cowpox"
    assert len(regions) == 1
    region = regions[0]
    assert region.minor_parent == "variola"
    assert region.major_parent == "cowpox"
    assert region.n_windows == 3
    assert region.mean_sim_minor > region.mean_sim_major


def test_min_region_filters_short_runs() -> None:
    cowpox = [0.99, 0.80, 0.99]
    variola = [0.90, 0.99, 0.90]  # a single recombinant window
    result = _result({"cowpox": cowpox, "variola": variola}, step=100)
    analysis = analyze(result)
    # window=100 -> the lone window spans 100 bp; require 1000 bp to keep it
    params = RegionParams.with_defaults(window_size=100, min_region=1000, method="heuristic")
    regions, _ = call_regions(result, analysis, window_size=100, params=params)
    assert regions == []


def test_margin_suppresses_marginal_windows() -> None:
    cowpox = [0.99, 0.95, 0.99]
    variola = [0.90, 0.96, 0.90]  # minor beats major by only 0.01 in window 1
    result = _result({"cowpox": cowpox, "variola": variola}, step=100)
    analysis = analyze(result)
    params = RegionParams.with_defaults(
        window_size=100, min_region=1, margin=0.05, method="heuristic"
    )
    regions, _ = call_regions(result, analysis, window_size=100, params=params)
    assert regions == []


def test_nan_windows_break_runs_and_are_skipped() -> None:
    nan = float("nan")
    # cowpox is the major parent; variola wins two short runs separated by an
    # uninformative (nan) window, which must break the run rather than merge it.
    cowpox = [0.99, 0.99, 0.80, nan, 0.80, 0.99, 0.99]
    variola = [0.90, 0.90, 0.99, nan, 0.99, 0.90, 0.90]
    result = _result({"cowpox": cowpox, "variola": variola}, step=100)
    analysis = analyze(result)
    params = RegionParams.with_defaults(
        window_size=100, min_region=1, merge_gap=0, method="heuristic"
    )
    regions, major = call_regions(result, analysis, window_size=100, params=params)
    assert major == "cowpox"
    assert len(regions) == 2
    assert all(r.minor_parent == "variola" for r in regions)
    # means are computed over finite windows only (never nan)
    assert all(r.mean_sim_minor == r.mean_sim_minor for r in regions)


def test_merge_gap_joins_same_parent_runs() -> None:
    # cowpox is the major parent (wins most windows); variola wins two short runs
    # separated by one cowpox window. merge_gap should join them into one region.
    cowpox = [0.99, 0.99, 0.80, 0.99, 0.80, 0.99, 0.99]
    variola = [0.90, 0.90, 0.99, 0.90, 0.99, 0.90, 0.90]
    result = _result({"cowpox": cowpox, "variola": variola}, step=100)
    analysis = analyze(result)
    params = RegionParams.with_defaults(
        window_size=100, min_region=1, merge_gap=1000, method="heuristic"
    )
    regions, _ = call_regions(result, analysis, window_size=100, params=params)
    assert len(regions) == 1
    assert regions[0].minor_parent == "variola"
