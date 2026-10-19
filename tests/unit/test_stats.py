"""Confidence-aware scoring primitives."""

from __future__ import annotations

import math

from tessera.recomb.stats import (
    benjamini_hochberg,
    emission_loglik,
    proportion_diff_significant,
    sign_test_greater,
    sign_test_pvalue,
    wilson_ci,
)


def test_wilson_ci_widens_with_less_evidence():
    lo_big, hi_big = wilson_ci(950, 1000)
    lo_small, hi_small = wilson_ci(5, 5)
    # 950/1000 is a tight band around 0.95; 5/5 is wide despite ratio 1.0
    assert hi_big - lo_big < 0.05
    assert hi_small - lo_small > 0.4
    assert lo_small < 1.0  # 5/5 does not imply certainty of identity


def test_wilson_ci_empty_window_is_uninformative():
    assert wilson_ci(0, 0) == (0.0, 1.0)


def test_proportion_diff_significant_needs_evidence():
    # a tiny ratio edge on near-identical big windows is NOT significant
    assert proportion_diff_significant(9990, 10000, 9985, 10000) is False
    # a clear difference on well-sampled windows IS significant
    assert proportion_diff_significant(950, 1000, 850, 1000) is True
    # never significant when the first proportion is not larger
    assert proportion_diff_significant(800, 1000, 900, 1000) is False
    # thin windows: not significant
    assert proportion_diff_significant(5, 5, 4, 5) is False


def test_sign_test_greater():
    # a clear majority of discordant sites favouring the donor is significant
    assert sign_test_greater(20, 3) is True
    # an even split is not
    assert sign_test_greater(10, 10) is False
    # never significant when the donor is not favoured
    assert sign_test_greater(3, 20) is False
    # no distinguishing sites -> not significant
    assert sign_test_greater(0, 0) is False
    # large-n normal-approximation branch (no overflow)
    assert sign_test_greater(900, 700) is True


def test_sign_test_pvalue():
    assert math.isclose(sign_test_pvalue(10, 0), 0.5**10)  # exact upper tail
    assert sign_test_pvalue(20, 3) < 0.01
    assert sign_test_pvalue(0, 0) == 1.0
    assert 0.0 < sign_test_pvalue(900, 700) < 0.001  # normal-approx branch


def test_benjamini_hochberg():
    assert benjamini_hochberg([]) == []
    # monotone, bounded by 1, and the smallest p gets the largest correction factor
    q = benjamini_hochberg([0.01, 0.02, 0.5])
    assert all(0.0 <= x <= 1.0 for x in q)
    assert math.isclose(q[0], 0.03)  # 0.01 * 3 / 1
    assert q[2] >= q[0]


def test_emission_loglik_prefers_fewer_mismatches_and_weights_n():
    # fewer mismatches -> higher (less negative) emission
    assert emission_loglik(99, 100, 0.95) > emission_loglik(90, 100, 0.95)
    # an empty window carries no information
    assert emission_loglik(0, 0, 0.95) == 0.0
    # value matches the closed form
    assert math.isclose(
        emission_loglik(8, 10, 0.9), 8 * math.log(0.9) + 2 * math.log(0.1)
    )


def _exact_log_tail(a: int, b: int) -> float:
    """log P(X >= a), X ~ Bin(a + b, 0.5), from exact integer binomial coefficients."""
    n = a + b
    return math.log(sum(math.comb(n, i) for i in range(a, n + 1))) - n * math.log(2.0)


def test_sign_test_pvalue_is_exact_above_a_thousand_sites():
    # Above n = 1000 the test used a normal approximation whose 1 - cdf is exactly 0.0
    # beyond z ~ 8.3, so 909:101 reported p = 0 while 900:100 reported 6.7e-162.
    for a, b in ((909, 101), (1500, 300), (600, 500), (5200, 4800)):
        got = sign_test_pvalue(a, b)
        assert got > 0.0
        assert math.isclose(math.log(got), _exact_log_tail(a, b), rel_tol=1e-9)


def test_sign_test_pvalue_is_monotone_across_the_old_cliff():
    # More discordant sites in the donor's favour can only lower the p-value.
    assert 0.0 < sign_test_pvalue(909, 101) < sign_test_pvalue(900, 100)
    assert 0.0 < sign_test_pvalue(1001, 0) < sign_test_pvalue(1000, 0)


def test_sign_test_pvalue_below_the_float_range_is_floored_not_zero():
    from tessera.recomb import stats

    # The true value, about 1e-322, is below the smallest normal double.
    assert stats.sign_test_pvalue(1827, 200) == stats.P_FLOOR > 0.0
    assert stats.sign_test_log_pvalue(1827, 200) < math.log(stats.P_FLOOR)
    assert math.isclose(stats.sign_test_log_pvalue(1827, 200), _exact_log_tail(1827, 200),
                        rel_tol=1e-9)


def test_sign_test_pvalue_one_sided_edges_at_large_n():
    # Review focus: every site against the donor, or every site for it, far above 1000.
    from tessera.recomb import stats

    assert sign_test_pvalue(0, 5000) == 1.0
    assert stats.sign_test_pvalue(5000, 0) == stats.P_FLOOR
    assert math.isclose(sign_test_pvalue(1, 0), 0.5)
    assert math.isclose(math.log(sign_test_pvalue(2500, 2500)), _exact_log_tail(2500, 2500),
                        rel_tol=1e-9)
