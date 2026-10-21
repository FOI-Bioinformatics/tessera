"""3SEQ-style triplet recombination test (Boni, Posada & Feldman 2007).

The HMM caller segments the query and confirms each segment with a discordant-site
sign test -- but that test is applied to a segment chosen *by looking at the data*,
so it is not aware that many breakpoints were scanned. The 3SEQ test fixes this: for
a triplet (query C, major parent P, minor parent Q) it looks only at the
**discriminating sites** where P and Q differ and C matches exactly one of them
(+1 for P, -1 for Q), and asks whether the query's matches to P and Q are *clustered*
(a P-region then a Q-region -- a mosaic) more than random interspersion would give.

The statistic is the **maximum drawdown** of the +/-1 walk over the discriminating
sites in genomic order: a deep, sustained run of Q-matches inside a P background. Its
exact p-value -- the probability that a random arrangement of the same number of P-
and Q-matches reaches that drawdown -- is computed by a dynamic program over the walk
depth (no precomputed table, no dependency), with a vectorised permutation fallback
for inputs too large for the exact DP.

Because it works on discriminating sites alone, the test keeps full power when the
panel is near-identical (intra-species sets, DNA viruses), where base-pair windows
have almost no signal; and being scan-aware it does not over-call.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .similarity import WindowSimilarity, _canonical_mask
from .similarity import region_identity as _region_similarity
from .stats import benjamini_hochberg

# Above this many DP cells (m * n * depth) the exact recursion is replaced by a
# vectorised permutation p-value -- bounds worst-case cost on highly divergent pairs.
_MAX_EXACT_CELLS = 60_000_000
_PERMUTATIONS = 20_000  # permutation-fallback resolution (~1/_PERMUTATIONS floor)


@dataclass
class Descent:
    """A maximal Q-run (the recombinant tract) within a P background."""

    depth: int  # the maximum drawdown, in discriminating sites
    start_site: int  # discriminating-site index where the descent begins (the peak)
    end_site: int  # ... where it ends (the trough)
    m: int  # discriminating sites where the query matches the major parent (P)
    n: int  # ... matches the minor parent (Q)


def triplet_steps(
    rows: dict[str, np.ndarray], query: str, major: str, minor: str
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(steps, columns)`` for the (query, major, minor) triplet.

    A *discriminating* column is one where the major and minor parents carry different
    canonical bases and the query carries a canonical base matching one of them.
    ``steps`` is +1 where the query matches the major and -1 where it matches the
    minor; ``columns`` are the MSA column indices of those sites (genomic order).
    Columns where the query matches both/neither, or any of the three is non-canonical,
    are uninformative and dropped.
    """
    q, p, r = rows[query], rows[major], rows[minor]
    canon = _canonical_mask(q) & _canonical_mask(p) & _canonical_mask(r)
    discriminating = canon & (p != r)
    q_major = discriminating & (q == p) & (q != r)
    q_minor = discriminating & (q == r) & (q != p)
    keep = q_major | q_minor
    columns = np.nonzero(keep)[0]
    steps = np.where(q_major[keep], 1, -1).astype(np.int64)
    return steps, columns


def max_descent(steps: np.ndarray) -> Descent:
    """The maximum drawdown of the +/-1 walk and the site interval that realises it."""
    cumulative = np.concatenate(([0], np.cumsum(steps)))  # S_0 .. S_N
    running_max = np.maximum.accumulate(cumulative)
    drawdown = running_max - cumulative
    trough = int(np.argmax(drawdown))
    depth = int(drawdown[trough])
    # The tract starts at the LAST maximum before the trough. The first one (plain
    # argmax) would pull in any balanced excursion that returns to the same height,
    # reporting the tract as starting earlier than the donor run does.
    peak = trough - int(np.argmax(cumulative[trough::-1])) if trough > 0 else 0
    m = int(np.count_nonzero(steps == 1))
    n = int(steps.size - m)
    return Descent(depth=depth, start_site=peak, end_site=trough, m=m, n=n)


def descent_pvalue_exact(m: int, n: int, depth: int) -> float:
    """Exact ``P(max drawdown >= depth)`` over the uniform arrangements of ``m`` +1 and
    ``n`` -1 steps.

    Dynamic program over the walk's current drawdown ``delta`` in ``[0, depth)``: with
    ``i`` up-steps and ``j`` down-steps still to come, an up-step moves
    ``delta -> max(0, delta - 1)`` with probability ``i / (i + j)`` and a down-step
    ``delta -> delta + 1`` with probability ``j / (i + j)``; a down-step from
    ``delta == depth - 1`` reaches the depth. ``g(i, j, delta)`` is the probability of
    reaching it, and the answer is ``g(m, n, 0)``.

    The reaching probability is propagated directly rather than as one minus the
    probability of never reaching it: every term is a sum of non-negative products, so
    a very small p-value keeps its magnitude instead of cancelling to exactly ``0.0``.
    Probabilities (not path counts) are propagated, so there is no big-integer overflow.

    The table is filled one anti-diagonal ``i + j = t`` at a time: every cell on it
    depends only on the previous diagonal, so each diagonal is one vectorised step. The
    values are those of the cell-by-cell recursion; only the order of evaluation
    changed, which made a long walk with a shallow descent (``m`` in the thousands) ten
    to forty times faster.
    """
    if depth <= 0:
        return 1.0
    if n < depth:  # need at least ``depth`` consecutive down-steps to reach it
        return 0.0
    prev = np.zeros((m + 1, depth))  # diagonal t-1, indexed by i (ups remaining)
    cur = np.zeros((m + 1, depth))
    for t in range(1, m + n + 1):
        lo, hi = max(0, t - n), min(m, t)
        i = np.arange(lo, hi + 1, dtype=np.float64)
        j = t - i
        block = np.zeros((hi - lo + 1, depth))
        # an up-step: (i, j, d) -> (i - 1, j, max(0, d - 1)), from row i - 1 of diagonal t-1
        a = 1 if lo == 0 else 0  # row i = 0 has no up-step
        if hi >= lo + a:
            src = prev[lo + a - 1 : hi]
            up = np.empty_like(src)
            up[:, 0] = src[:, 0]
            up[:, 1:] = src[:, :-1]
            block[a:] += (i[a:] / t)[:, None] * up
        # a down-step: (i, j, d) -> (i, j - 1, d + 1), or the depth is reached
        b = 1 if hi == t else 0  # row j = 0 has no down-step
        if hi - b >= lo:
            src = prev[lo : hi - b + 1]
            down = np.empty_like(src)
            down[:, :-1] = src[:, 1:]
            down[:, -1] = 1.0
            block[: hi - lo + 1 - b] += (j[: hi - lo + 1 - b] / t)[:, None] * down
        cur[lo : hi + 1] = block
        prev, cur = cur, prev
    return float(prev[m][0])


def _descent_pvalue_perm(steps: np.ndarray, depth: int, *, seed: int = 0) -> float:
    """Vectorised permutation ``P(max drawdown >= depth)`` (exact-in-the-limit)."""
    rng = np.random.default_rng(seed)
    n_steps = steps.size
    order = np.argsort(rng.random((_PERMUTATIONS, n_steps)), axis=1)
    permuted = steps[order]
    cumulative = np.cumsum(permuted, axis=1)
    cumulative = np.hstack([np.zeros((_PERMUTATIONS, 1), dtype=cumulative.dtype), cumulative])
    drawdowns = (np.maximum.accumulate(cumulative, axis=1) - cumulative).max(axis=1)
    return float((np.count_nonzero(drawdowns >= depth) + 1) / (_PERMUTATIONS + 1))


def descent_pvalue(steps: np.ndarray, descent: Descent, *, seed: int = 0) -> tuple[float, str]:
    """``(p_value, method)`` for the observed drawdown -- exact DP, or a permutation
    fallback when the exact state space is too large."""
    if descent.depth <= 0:
        return 1.0, "exact"
    if descent.m * descent.n * descent.depth <= _MAX_EXACT_CELLS:
        return descent_pvalue_exact(descent.m, descent.n, descent.depth), "exact"
    return _descent_pvalue_perm(steps, descent.depth, seed=seed), "permutation"


# The triplet test needs enough discriminating sites to be meaningful. There is no floor
# on the descent depth: a shallow descent can still be significant when the donor
# matches are few (three consecutive ones among 189 sites give p = 1.7e-4), so the
# p-value decides, and every donor tested counts towards the Benjamini-Hochberg family.
_MIN_DISCRIMINATING = 20


def call_regions_3seq(result: WindowSimilarity, analysis, params):
    """Call recombinant regions with the 3SEQ triplet test.

    The major parent is the reference winning the most windows. Each plausible donor
    (a reference that wins somewhere) is tested as the minor parent: over the major /
    minor discriminating sites, a significant maximum drawdown (after Benjamini-Hochberg
    across the tested donors) is a recombinant tract. Overlapping tracts keep the most
    significant donor. Returns ``(regions, major, [])`` to match ``call_regions``.
    """
    from .analyze import rank_datasets
    from .regions import Region

    labels = list(result.similarities)
    if len(labels) < 2:
        return [], (labels[0] if labels else None), []
    # Rank by window wins but keep non-winning references as candidate donors too: a
    # sub-window recombinant tract is diluted in the windowed vote, so its true donor may
    # win no window yet still carry a strong discriminating-site signal. Each donor tested
    # is one member of the Benjamini-Hochberg family below, whatever its p-value.
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
        # `kind` is "exact" or "permutation": which one answered matters to a reader,
        # so it is carried through to the region's `test` column rather than discarded.
        p, kind = descent_pvalue(steps, descent)
        candidates.append((minor, steps, cols, descent, p, kind))
    if not candidates:
        return [], major, []

    qvalues = benjamini_hochberg([c[4] for c in candidates])
    regions: list[Region] = []
    for (minor, steps, cols, descent, p, kind), q in sorted(
        zip(candidates, qvalues, strict=True), key=lambda t: t[0][4]
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
            test=f"3SEQ max-descent ({kind})",
            statistic="fraction of tract sites matching the minor parent",
            breakpoint_lo=min(bp_lo, bp_hi), breakpoint_hi=max(bp_lo, bp_hi),
        ))
    regions.sort(key=lambda r: r.msa_start)
    return regions, major, []
