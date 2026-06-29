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
