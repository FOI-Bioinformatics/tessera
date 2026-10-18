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
    _MIN_DISCRIMINATING,
    _region_similarity,
    max_descent,
    triplet_steps,
)

_PERMUTATIONS = 20_000  # permutation null resolution (~1/_PERMUTATIONS floor)
# Peak cells in any one permutation batch. Each batch holds several (batch, length)
# int64 arrays, so an uncapped run on a donor with a few thousand discriminating sites
# allocated hundreds of megabytes per array -- per candidate donor, once per reference.
# This bounds the batch, not the permutation count: see maxchi_pvalue.
_MAX_PERM_CELLS = 20_000_000
# Below this descent depth a donor enters the BH family with p = 1 instead of being
# permuted; see call_regions_maxchi.
_MIN_PERMUTED_DEPTH = 4



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


def _exceedances(steps: np.ndarray, observed: float, k: int, rng) -> int:
    """How many of ``k`` permutations reach ``observed``; vectorised over the batch."""
    length = steps.size
    perm = steps[np.argsort(rng.random((k, length)), axis=1)]
    cum = np.hstack([np.zeros((k, 1), dtype=np.int64), np.cumsum(perm, axis=1)])
    runmax = np.maximum.accumulate(cum, axis=1)
    trough = np.argmax(runmax - cum, axis=1)  # (k,) end of each permutation's tract
    cols = np.arange(length + 1)
    # peak = the LAST place the running max was reached, on or before the trough --
    # the same rule threeseq.max_descent applies to the observed walk, so the null
    # and the observed statistic are computed on the same kind of tract.
    masked = np.where(cols[None, :] <= trough[:, None], cum, np.iinfo(np.int64).min)
    peak = length - np.argmax(masked[:, ::-1], axis=1)
    cum_major = np.hstack([np.zeros((k, 1), np.int64), np.cumsum(perm == 1, axis=1)])
    cum_minor = np.hstack([np.zeros((k, 1), np.int64), np.cumsum(perm == -1, axis=1)])
    rows = np.arange(k)
    in_major = cum_major[rows, trough] - cum_major[rows, peak]
    in_minor = cum_minor[rows, trough] - cum_minor[rows, peak]
    chi2 = _chi2_batch(in_major, in_minor, cum_major[:, -1], cum_minor[:, -1])
    return int(np.count_nonzero(chi2 >= observed))


def _chi2_batch(
    in_major: np.ndarray, in_minor: np.ndarray, tot_major: np.ndarray, tot_minor: np.ndarray
) -> np.ndarray:
    """Vectorised :func:`_chi2_2x2` for tract-versus-rest tables.

    Every factor is converted to float64 before it is multiplied. The denominator is a
    product of four counts; in int64 it wraps above about 110 000 discriminating sites
    with balanced counts, and the positivity guard below then set that permutation's
    statistic to zero (or the wrapped value was positive and wrong) -- while the observed
    statistic, computed with Python integers, did not wrap. float64 holds the product
    without wrapping at any size a genome alignment reaches.
    """
    a = in_major.astype(np.float64)
    b = in_minor.astype(np.float64)
    c = tot_major.astype(np.float64) - a
    d = tot_minor.astype(np.float64) - b
    denom = (a + b) * (c + d) * (a + c) * (b + d)
    det = a * d - b * c
    n = a + b + c + d
    return np.where(denom > 0, n * det**2 / np.where(denom > 0, denom, 1.0), 0.0)


def maxchi_pvalue(steps: np.ndarray, observed: float, *, seed: int = 0) -> float:
    """Scan-aware permutation ``P(chi-square >= observed)``: each permutation of the walk
    finds its own maximal-drawdown tract and that tract's chi-square (vectorised).

    The permutations are evaluated in batches sized so no intermediate exceeds
    :data:`_MAX_PERM_CELLS`. Batching bounds peak memory *without* reducing the number
    of permutations: each row is independent, so the exceedance counts simply add.
    Capping the count instead was tried and reverted -- it engaged above a mere 1000
    discriminating sites, which is an ordinary panel, and the coarser null cost real
    detections on the specificity harness's positive control.
    """
    if observed <= 0.0 or steps.size < 2:
        return 1.0
    rng = np.random.default_rng(seed)
    length = steps.size
    batch = max(1, min(_PERMUTATIONS, _MAX_PERM_CELLS // max(1, length)))
    exceed = 0
    remaining = _PERMUTATIONS
    while remaining > 0:
        this = min(batch, remaining)
        exceed += _exceedances(steps, observed, this, rng)
        remaining -= this
    return float((exceed + 1) / (_PERMUTATIONS + 1))


def call_regions_maxchi(result: WindowSimilarity, analysis, params):
    """Call recombinant regions with the MaxChi chi-square triplet test.

    Mirrors the 3SEQ caller's structure: the major parent wins the most windows; each
    plausible donor is tested over its discriminating sites; a tract whose chi-square is
    significant after Benjamini-Hochberg across donors is a region. Returns
    ``(regions, major, [])`` to match ``call_regions``.
    """
    from .analyze import rank_datasets
    from .regions import Region

    labels = list(result.similarities)
    if len(labels) < 2:
        return [], (labels[0] if labels else None), []
    # Keep non-winning references as candidate donors: a sub-window tract's donor may win no
    # window yet carry a strong discriminating-site signal.
    ranked = rank_datasets(analysis, len(labels)) or labels
    major = ranked[0]
    candidates = []
    for minor in ranked:
        if minor == major:
            continue
        steps, cols = triplet_steps(result.rows, result.query, major, minor)
        if steps.size < _MIN_DISCRIMINATING:
            continue
        descent = max_descent(steps)
        chi2 = maxchi_statistic(steps, descent.start_site, descent.end_site)
        # A donor with a shallow descent is still a member of the Benjamini-Hochberg
        # family: it is counted, with p = 1, rather than dropped before the correction,
        # which shrank the family to the survivors and inflated their significance. Its
        # permutation null is not run: 20 000 permutations of a long walk per such donor
        # made a 16-genome scan ten times slower, and 3SEQ, which scores the same tract
        # exactly and cheaply, tests these shallow tracts (see threeseq).
        p = maxchi_pvalue(steps, chi2) if descent.depth >= _MIN_PERMUTED_DEPTH else 1.0
        candidates.append((minor, steps, cols, descent, chi2, p))
    if not candidates:
        return [], major, []

    qvalues = benjamini_hochberg([c[5] for c in candidates])
    regions: list[Region] = []
    for (minor, steps, cols, descent, _chi2, p), q in sorted(
        zip(candidates, qvalues, strict=True), key=lambda t: t[0][5]
    ):
        if q > params.alpha:   # FDR-controlled gate (BH q, not raw p)
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
        # The breakpoint lies between the last site before the tract and its first
        # site. A tract starting at the first discriminating site has no site before
        # it, so the interval runs from the alignment start, not zero-width.
        bp_lo = result.column_to_query(int(cols[peak - 1])) if peak > 0 else 0
        bp_hi = result.column_to_query(int(cols[peak]))
        regions.append(Region(
            minor_parent=minor, major_parent=major,
            msa_start=msa_start, msa_end=msa_end,
            query_start=result.column_to_query(msa_start),
            query_end=result.column_to_query(msa_end),
            n_windows=trough - peak,
            mean_sim_minor=round(sim_minor, 4), mean_sim_major=round(sim_major, 4),
            margin=round(sim_minor - sim_major, 4),
            support=round(support, 3), pvalue=p, qvalue=q,
            test="MaxChi chi-square (scan-aware permutation)",
            statistic="fraction of tract sites matching the minor parent",
            breakpoint_lo=min(bp_lo, bp_hi), breakpoint_hi=max(bp_lo, bp_hi),
        ))
    regions.sort(key=lambda r: r.msa_start)
    return regions, major, []
