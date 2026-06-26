"""Reference-coverage diagnostics: flagging possible missing references."""

from __future__ import annotations

import numpy as np

from recomfi.recomb.coverage import (
    CoverageParams,
    call_coverage_gaps,
    flag_undercovered_regions,
)
from recomfi.recomb.regions import Region
from recomfi.recomb.similarity import WindowSimilarity


def _result(best_sim, best_label, informative, *, window=100):
    """A minimal WindowSimilarity with non-overlapping windows over an ungapped query."""
    n = len(best_sim)
    width = n * window
    positions = [i * window + window // 2 for i in range(n)]
    return WindowSimilarity(
        positions=positions,
        query_positions=positions,
        similarities={},
        query="q",
        width=width,
        query_cumulative=np.arange(width + 1),  # ungapped: column == query base
        best_sim=best_sim,
        best_label=best_label,
        informative=informative,
    )


def test_divergent_gap_is_flagged():
    # high coverage everywhere except a contiguous block where the best ref is poor
    best = [0.98] * 10
    labels = ["refA"] * 10
    info = [100] * 10
    for i in (4, 5, 6):
        best[i], labels[i] = 0.80, "refB"
    result = _result(best, labels, info)

    gaps, threshold = call_coverage_gaps(result, 100, CoverageParams.with_defaults(100))

    assert 0.90 < threshold < 0.95  # adaptive: ~p90(best) - 0.05
    assert len(gaps) == 1
    gap = gaps[0]
    assert (gap.query_start, gap.query_end) == (400, 700)
    assert gap.best_label == "refB"
    assert gap.kind == "divergent"
    assert gap.mean_best == 0.80


def test_low_information_gap_is_labelled_distinctly():
    best = [0.98] * 10
    labels = ["refA"] * 10
    info = [100] * 10
    for i in (4, 5, 6):
        best[i], info[i] = 0.80, 3  # poor best, but almost no comparable bases
    result = _result(best, labels, info)

    gaps, _ = call_coverage_gaps(result, 100, CoverageParams.with_defaults(100))

    assert len(gaps) == 1
    assert gaps[0].kind == "low_information"


def test_full_coverage_yields_no_gaps():
    result = _result([0.99] * 8, ["refA"] * 8, [100] * 8)
    gaps, _ = call_coverage_gaps(result, 100, CoverageParams.with_defaults(100))
    assert gaps == []


def test_absolute_floor_overrides_adaptive():
    best = [0.90] * 6  # uniformly moderate: adaptive would flag nothing
    result = _result(best, ["refA"] * 6, [100] * 6)

    none_flagged, thr_adaptive = call_coverage_gaps(
        result, 100, CoverageParams.with_defaults(100)
    )
    assert none_flagged == []

    flagged, thr_floor = call_coverage_gaps(
        result, 100, CoverageParams.with_defaults(100, floor=0.95)
    )
    assert thr_floor == 0.95
    assert len(flagged) == 1  # all six windows now under the 0.95 floor -> one gap


def test_flag_undercovered_regions():
    regions = [
        Region("variola", "cowpox", 0, 100, 0, 100, 100, 5, 0.998, 0.97, 0.028),
        Region("env", "A1", 200, 300, 200, 300, 100, 5, 0.80, 0.79, 0.01),
    ]
    flag_undercovered_regions(regions, threshold=0.90)
    assert regions[0].donor_undercovered is False  # donor matches well
    assert regions[1].donor_undercovered is True  # donor itself is a poor match
