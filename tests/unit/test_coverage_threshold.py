"""The best-similarity threshold on a panel that is nearly identical to the query.

Issue #82, first half: the threshold is a fixed drop below the query's typical best
match, so on a panel whose references all sit within a percent of the query it lands
below everything the panel can show and cannot resolve anything. The mpox run that
prompted the issue recorded its threshold as 0.9500 with a minimum best similarity of
0.9950 and reported no coverage gaps at all, although two loci there are locally
misaligned.

These tests pin a threshold that falls back to the panel's own scale when the fixed
drop is vacuous, and pin that what it then finds is not reported as a missing
reference: at 99.7 % identity the closest reference is plainly present, so such a
stretch must not be bridged into a donor-absent region or searched for over the
network.
"""

from __future__ import annotations

from math import isclose

import numpy as np

from tessera.recomb.coverage import (
    LOCAL_DIP_KIND,
    CoverageParams,
    call_coverage_gaps,
    coverage_levels,
    gaps_as_regions,
    search_targets,
)
from tessera.recomb.similarity import WindowSimilarity

WINDOW = 1000


def _result(best_sim: list[float], informative: list[int]) -> WindowSimilarity:
    """Non-overlapping windows over an ungapped query, one best match per window."""
    n = len(best_sim)
    width = n * WINDOW
    positions = [i * WINDOW + WINDOW // 2 for i in range(n)]
    return WindowSimilarity(
        positions=positions,
        query_positions=positions,
        similarities={},
        query="q",
        width=width,
        query_cumulative=np.arange(width + 1),
        best_sim=best_sim,
        best_label=["refA"] * n,
        informative=informative,
    )


def _near_identical(dip: float = 0.997, n_dip: int = 2, informative: int = WINDOW):
    """A panel matching the query almost perfectly, with a short dip in the middle."""
    best = [1.0] * 25
    for i in range(10, 10 + n_dip):
        best[i] = dip
    return _result(best, [informative] * 25)


def test_a_vacuous_fixed_drop_falls_back_to_the_panels_own_scale() -> None:
    """The threshold must sit where this panel's similarities actually are.

    A fixed 0.05 drop below a typical best match of 1.0 gives 0.95, below every
    value the panel shows, so no window can ever fall under it. The fallback puts
    the threshold within reach of the observed dip.
    """
    result = _near_identical()
    threshold, coarse = coverage_levels(result, CoverageParams.with_defaults(WINDOW))

    assert isclose(coarse, 0.95)  # the fixed standard, unchanged
    assert 0.997 < threshold < 1.0
    assert threshold > min(result.best_sim)


def test_the_fallback_resolves_a_dip_the_fixed_drop_cannot_see() -> None:
    """The locus the mpox run reported nothing about becomes one coverage gap."""
    result = _near_identical()

    gaps, threshold = call_coverage_gaps(result, WINDOW, CoverageParams.with_defaults(WINDOW))

    assert len(gaps) == 1
    assert gaps[0].n_windows == 2
    assert isclose(gaps[0].mean_best, 0.997)


def test_the_fallback_is_no_finer_than_a_window_can_resolve() -> None:
    """The floor is two mismatches over the window's comparable bases.

    A window carrying 100 comparable bases cannot resolve an identity difference
    below one part in a hundred, so the threshold must not be placed there. With
    1000 comparable bases it can, and sits ten times closer to the typical match.
    """
    coarse_evidence = _near_identical(dip=0.975, informative=100)
    ample_evidence = _near_identical(dip=0.975, informative=WINDOW)

    thin, _ = coverage_levels(coarse_evidence, CoverageParams.with_defaults(WINDOW))
    ample, _ = coverage_levels(ample_evidence, CoverageParams.with_defaults(WINDOW))

    assert isclose(thin, 0.98)
    assert isclose(ample, 0.998)


def test_a_panel_with_spread_keeps_the_fixed_drop() -> None:
    """The fallback triggers only where the fixed drop is unresolvable.

    Here the worst window is already below the fixed standard, so that standard
    resolves the panel and the threshold is left exactly where it was.
    """
    best = [1.0] * 20
    best[5], best[6] = 0.94, 0.93
    result = _result(best, [WINDOW] * 20)

    threshold, coarse = coverage_levels(result, CoverageParams.with_defaults(WINDOW))

    assert isclose(threshold, 0.95)
    assert isclose(coarse, 0.95)


def test_an_absolute_floor_still_overrides_both_levels() -> None:
    """``--coverage-floor`` is an absolute statement and keeps its meaning."""
    result = _near_identical()

    threshold, coarse = coverage_levels(
        result, CoverageParams.with_defaults(WINDOW, floor=0.80)
    )

    assert isclose(threshold, 0.80)
    assert isclose(coarse, 0.80)


def test_a_dip_above_the_fixed_standard_is_not_a_missing_reference() -> None:
    """At 99.7 % identity the closest reference is present, whatever the dip.

    So the stretch is labelled as a local dip rather than ``divergent``, and the
    bridge that turns a divergent gap into a donor-absent putative recombination
    leaves it alone.
    """
    result = _near_identical()

    gaps, _ = call_coverage_gaps(result, WINDOW, CoverageParams.with_defaults(WINDOW))

    assert [g.kind for g in gaps] == [LOCAL_DIP_KIND]
    assert gaps_as_regions(gaps, result, "major") == []


def test_recruitment_does_not_search_for_a_local_dip() -> None:
    """A locally poor alignment is not a genome to go and download.

    ``find-references`` and the fill loop treat every gap as a search target, so a
    local dip has to be withheld from them or the loop chases an absent reference
    that is not absent.
    """
    result = _near_identical()
    gaps, _ = call_coverage_gaps(result, WINDOW, CoverageParams.with_defaults(WINDOW))

    assert gaps
    assert search_targets(gaps) == []


def test_recruitment_still_searches_for_a_divergent_gap() -> None:
    """The withholding is narrow: a genuinely divergent gap stays a search target."""
    best = [1.0] * 20
    best[5], best[6] = 0.80, 0.81
    result = _result(best, [WINDOW] * 20)
    gaps, _ = call_coverage_gaps(result, WINDOW, CoverageParams.with_defaults(WINDOW))

    assert len(gaps) == 1
    assert search_targets(gaps) == gaps


def test_a_panel_a_few_percent_from_the_query_keeps_the_fixed_drop() -> None:
    """Unreachable is not the same as vacuous.

    Measured on the cached rubella pool: the typical best match is 0.9592, the
    fixed standard 0.9092, and the worst window 0.9152 -- so no window reaches the
    standard, but only because the panel covers the query well everywhere. That is
    a legitimate negative, not a threshold that carries no information, and
    tightening it reports ordinary variation in a panel four percent from the
    query. The fallback is for a panel whose closest reference is within one
    percent of the query across the genome, which is the regime issues #81 and #82
    describe.
    """
    best = [0.9592] * 80 + [0.9278] * 5 + [0.9152] * 4
    result = _result(best, [WINDOW] * len(best))

    threshold, coarse = coverage_levels(result, CoverageParams.with_defaults(WINDOW))
    gaps, _ = call_coverage_gaps(result, WINDOW, CoverageParams.with_defaults(WINDOW))

    assert isclose(threshold, coarse)
    assert isclose(threshold, 0.9592 - 0.05)
    assert gaps == []
