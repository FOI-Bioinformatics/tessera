"""MaxChi / Chimaera caller (Maynard Smith 1992; RDP family).

A second triplet test on the same discriminating-site +/-1 walk the 3SEQ caller uses
(:func:`tessera.recomb.threeseq.triplet_steps`), but with an independent statistic. 3SEQ
asks whether the minor-matching sites form a deep contiguous *run* (maximum drawdown);
MaxChi asks whether the **proportion** of minor matches differs sharply inside a tract
versus outside it -- a 2x2 chi-square. The two are sensitive to different mosaics (a
shallow-but-wide tract has a weak drawdown but a strong proportion shift, and vice
versa), so in an ensemble a region significant under *both* nulls is corroborated by two
independent tests, exactly as RDP treats MaxChi and 3SEQ as separate methods.

The candidate tract is the same maximal minor run 3SEQ localises (reusing
:func:`max_descent`); the significance is the chi-square of the tract's major:minor
counts against the rest, with a **scan-aware permutation** null (each permutation finds
its own best tract) so scanning every breakpoint does not inflate the call. Dependency-
free numpy.
"""

from __future__ import annotations

import numpy as np

from .similarity import WindowSimilarity
from .stats import benjamini_hochberg
from .threeseq import (
    _MIN_DESCENT,
    _MIN_DISCRIMINATING,
    _region_similarity,
    max_descent,
    triplet_steps,
)

_PERMUTATIONS = 20_000  # permutation null resolution (~1/_PERMUTATIONS floor)


def _chi2_2x2(a: int, b: int, c: int, d: int) -> float:
    """Pearson chi-square of the 2x2 table ``[[a, b], [c, d]]`` (0 if any margin is 0)."""
    n = a + b + c + d
    r1, r2, c1, c2 = a + b, c + d, a + c, b + d
    if n == 0 or min(r1, r2, c1, c2) == 0:
        return 0.0
    return n * float(a * d - b * c) ** 2 / (r1 * r2 * c1 * c2)


def maxchi_statistic(steps: np.ndarray, peak: int, trough: int) -> float:
    """Chi-square contrasting major(+1):minor(-1) counts inside ``[peak, trough)`` vs
    outside -- large when the tract's match proportion differs sharply from the rest."""
    inside = steps[peak:trough]
    in_major = int(np.count_nonzero(inside == 1))
    in_minor = int(np.count_nonzero(inside == -1))
    total_major = int(np.count_nonzero(steps == 1))
    total_minor = int(np.count_nonzero(steps == -1))
    return _chi2_2x2(in_major, in_minor, total_major - in_major, total_minor - in_minor)


def maxchi_pvalue(steps: np.ndarray, observed: float, *, seed: int = 0) -> float:
    """Scan-aware permutation ``P(chi-square >= observed)``: each permutation of the walk
    finds its own maximal-drawdown tract and that tract's chi-square (vectorised)."""
    if observed <= 0.0 or steps.size < 2:
        return 1.0
    rng = np.random.default_rng(seed)
    k, length = _PERMUTATIONS, steps.size
    perm = steps[np.argsort(rng.random((k, length)), axis=1)]
    cum = np.hstack([np.zeros((k, 1), dtype=np.int64), np.cumsum(perm, axis=1)])
    runmax = np.maximum.accumulate(cum, axis=1)
    trough = np.argmax(runmax - cum, axis=1)  # (k,) end of each permutation's tract
    cols = np.arange(length + 1)
    # peak = where the running max was reached, on or before the trough
    masked = np.where(cols[None, :] <= trough[:, None], cum, np.iinfo(np.int64).min)
    peak = np.argmax(masked, axis=1)
    cum_major = np.hstack([np.zeros((k, 1), np.int64), np.cumsum(perm == 1, axis=1)])
    cum_minor = np.hstack([np.zeros((k, 1), np.int64), np.cumsum(perm == -1, axis=1)])
    rows = np.arange(k)
    in_major = cum_major[rows, trough] - cum_major[rows, peak]
    in_minor = cum_minor[rows, trough] - cum_minor[rows, peak]
    tot_major, tot_minor = cum_major[:, -1], cum_minor[:, -1]
    out_major, out_minor = tot_major - in_major, tot_minor - in_minor
    n = length
    r1 = in_major + in_minor
    r2 = out_major + out_minor
    det = in_major.astype(np.float64) * out_minor - in_minor.astype(np.float64) * out_major
    denom = r1 * r2 * tot_major * tot_minor
    chi2 = np.where(denom > 0, n * det**2 / np.where(denom > 0, denom, 1), 0.0)
    return float((np.count_nonzero(chi2 >= observed) + 1) / (k + 1))


def call_regions_maxchi(result: WindowSimilarity, analysis, params):
    """Call recombinant regions with the MaxChi chi-square triplet test.

    Mirrors the 3SEQ caller's structure: the major parent wins the most windows; each
    plausible donor is tested over its discriminating sites; a tract whose chi-square is
    significant after Benjamini-Hochberg across donors is a region. Returns
    ``(regions, major, [])`` to match ``call_regions``.
    """
    from .analyze import rank_by_wins
    from .regions import Region, _signif

    labels = list(result.similarities)
    if len(labels) < 2:
        return [], (labels[0] if labels else None), []
    ranked = rank_by_wins(analysis.winners_with_ties, len(labels)) or labels
    major = ranked[0]
    candidates = []
    for minor in ranked:
        if minor == major:
            continue
        steps, cols = triplet_steps(result.rows, result.query, major, minor)
        if steps.size < _MIN_DISCRIMINATING:
            continue
        descent = max_descent(steps)
        if descent.depth < _MIN_DESCENT:
            continue
        chi2 = maxchi_statistic(steps, descent.start_site, descent.end_site)
        p = maxchi_pvalue(steps, chi2)
        candidates.append((minor, steps, cols, descent, chi2, p))
    if not candidates:
        return [], major, []

    qvalues = benjamini_hochberg([c[5] for c in candidates])
    regions: list[Region] = []
    for (minor, steps, cols, descent, _chi2, p), q in sorted(
        zip(candidates, qvalues, strict=True), key=lambda t: t[0][5]
    ):
        if p > params.alpha:
            continue
        peak, trough = descent.start_site, descent.end_site
        if trough <= peak or peak >= cols.size:
            continue
        msa_start = int(cols[peak])
        msa_end = int(cols[min(trough, cols.size) - 1]) + 1
        if msa_end <= msa_start or any(
            msa_start < r.msa_end and r.msa_start < msa_end for r in regions
        ):
            continue  # keep the most significant donor for each tract
        tract = steps[peak:trough]
        support = float(np.count_nonzero(tract == -1) / tract.size) if tract.size else 0.0
        sim_minor = _region_similarity(result.rows, result.query, minor, msa_start, msa_end)
        sim_major = _region_similarity(result.rows, result.query, major, msa_start, msa_end)
        bp_lo = result.column_to_query(int(cols[max(0, peak - 1)]))
        bp_hi = result.column_to_query(int(cols[peak]))
        regions.append(Region(
            minor_parent=minor, major_parent=major,
            msa_start=msa_start, msa_end=msa_end,
            query_start=result.column_to_query(msa_start),
            query_end=result.column_to_query(msa_end),
            length_bp=msa_end - msa_start, n_windows=trough - peak,
            mean_sim_minor=round(sim_minor, 4), mean_sim_major=round(sim_major, 4),
            margin=round(sim_minor - sim_major, 4),
            support=round(support, 3), pvalue=_signif(p), qvalue=_signif(q),
            breakpoint_lo=min(bp_lo, bp_hi), breakpoint_hi=max(bp_lo, bp_hi),
        ))
    regions.sort(key=lambda r: r.msa_start)
    return regions, major, []
