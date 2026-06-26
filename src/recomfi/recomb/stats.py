"""Confidence-aware primitives behind the recombination scoring.

Closed-form and dependency-free (the normal quantile comes from the standard
library's :class:`statistics.NormalDist`, so no SciPy):

- :func:`wilson_ci` -- a Wilson score interval for a per-window identity ``k/n``,
  so a 950/1000 window is distinguished from a 5/5 window;
- :func:`proportion_diff_significant` -- a one-sided two-proportion test for "is
  reference A genuinely closer to the query than reference B?", which replaces the
  arbitrary ``margin = 0`` window rule with a noise-aware criterion;
- :func:`emission_loglik` -- the per-site copying log-likelihood used as the HMM
  emission (the query is modelled as copied from a source with match probability
  ``p``).
"""

from __future__ import annotations

import math
from statistics import NormalDist

_NORMAL = NormalDist()


def _z_two_sided(alpha: float) -> float:
    return _NORMAL.inv_cdf(1.0 - alpha / 2.0)


def _z_one_sided(alpha: float) -> float:
    return _NORMAL.inv_cdf(1.0 - alpha)


def wilson_ci(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Wilson score interval for an identity ``k/n`` at confidence ``1 - alpha``.

    Returns ``(0.0, 1.0)`` for an empty window (``n == 0``) -- no information.
    """
    if n <= 0:
        return 0.0, 1.0
    z = _z_two_sided(alpha)
    phat = k / n
    denom = 1.0 + z * z / n
    center = (phat + z * z / (2.0 * n)) / denom
    half = z * math.sqrt(phat * (1.0 - phat) / n + z * z / (4.0 * n * n)) / denom
    return max(0.0, center - half), min(1.0, center + half)


def proportion_diff_significant(
    k1: int, n1: int, k2: int, n2: int, alpha: float = 0.05
) -> bool:
    """One-sided two-proportion z-test: is ``k1/n1`` significantly above ``k2/n2``?

    Returns ``False`` unless ``p1 > p2`` and the pooled z-test rejects at ``alpha``
    -- the noise-aware replacement for "minor beats major by any positive margin".
    """
    if n1 <= 0 or n2 <= 0:
        return False
    p1, p2 = k1 / n1, k2 / n2
    if p1 <= p2:
        return False
    p_pool = (k1 + k2) / (n1 + n2)
    se = math.sqrt(p_pool * (1.0 - p_pool) * (1.0 / n1 + 1.0 / n2))
    if se == 0.0:
        return False
    return (p1 - p2) / se >= _z_one_sided(alpha)


def sign_test_greater(a: int, b: int, alpha: float = 0.05) -> bool:
    """One-sided exact sign test: is ``a`` significantly more than ``b``?

    Under the null (no preference) ``a ~ Binomial(a+b, 0.5)``; returns ``True``
    when ``a > b`` and ``P(X >= a) <= alpha``. Used on the discordant sites that
    distinguish two candidate parents, so it is unaffected by window overlap.
    """
    n = a + b
    if n == 0 or a <= b:
        return False
    if n <= 1000:  # exact upper tail of Binomial(n, 0.5)
        tail = sum(math.comb(n, i) for i in range(a, n + 1)) / (2.0**n)
        return tail <= alpha
    # normal approximation with a continuity correction for large n
    z = (a - n / 2.0 - 0.5) / math.sqrt(n / 4.0)
    return z >= _z_one_sided(alpha)


def emission_loglik(k: int, n: int, p: float) -> float:
    """Per-site copying log-likelihood ``k*log(p) + (n-k)*log(1-p)``.

    The HMM emission for "the query was copied from this reference here": every
    matching comparable site contributes ``log(p)``, every mismatch ``log(1-p)``.
    The combinatorial term is omitted (the sites are observed, not chosen), so the
    reference with the fewest local mismatches emits highest, weighted by ``n``
    (a thin window carries little evidence). ``n == 0`` emits ``0`` (uninformative).
    """
    if n <= 0:
        return 0.0
    p = min(max(p, 1e-6), 1.0 - 1e-6)
    return k * math.log(p) + (n - k) * math.log(1.0 - p)
