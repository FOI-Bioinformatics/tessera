"""GENECONV-style caller: the longest uninterrupted donor-match run (Sawyer 1989).

3SEQ scores the maximum *drawdown* of the +/-1 walk over discriminating sites -- a deep,
sustained lean toward the minor parent that tolerates a few interspersed major-matches.
GENECONV instead looks for a **clean fragment**: the longest run of *consecutive*
discriminating sites where the query matches the minor parent and not the major, with no
interruption. Its significance is the permutation probability that a random arrangement of
the same numbers of major- and minor-matches produces a minor-run at least that long.

This is the query-vs-panel adaptation of GENECONV's inner-fragment test: it reuses the
triplet discriminating-site series (``triplet_steps``), so like 3SEQ it keeps power on
near-identical panels where base-pair windows carry little signal. It is opt-in (not in
``DEFAULT_METHODS``); the run-based statistic is complementary to 3SEQ's drawdown and
MaxChi's boundary chi-square.
"""

from __future__ import annotations

import numpy as np

from .similarity import WindowSimilarity
from .stats import benjamini_hochberg
from .threeseq import _region_similarity, triplet_steps

# A fragment needs enough discriminating sites to assess, and a run longer than trivial
# noise; the permutation p-value gate does the real work above these floors.
_MIN_DISCRIMINATING = 20
_MIN_RUN = 4
_PERMUTATIONS = 20_000
# Cap the vectorised permutation work (permutations * sites) so a long, dense series does
# not blow up; permutations are scaled down (never below this floor) to stay under it.
_MAX_PERM_CELLS = 40_000_000
_MIN_PERMUTATIONS = 2_000


def longest_minor_run(steps: np.ndarray) -> tuple[int, int, int]:
    """``(length, start_site, end_site)`` of the longest run of consecutive minor-matches
    (``steps == -1``), as inclusive discriminating-site indices; ``(0, -1, -1)`` if none."""
    best = cur = 0
    best_end = end = -1
    for i in range(steps.size):
        if steps[i] == -1:
            cur += 1
            end = i
        else:
            cur = 0
        if cur > best:
            best, best_end = cur, end
    return best, (best_end - best + 1 if best else -1), best_end


def run_pvalue_perm(steps: np.ndarray, run_len: int, *, seed: int = 0) -> float:
    """Permutation ``P(longest minor-run >= run_len)`` over random arrangements of the
    same major/minor counts (exact-in-the-limit; mirrors the 3SEQ permutation fallback)."""
    if run_len <= 0:
        return 1.0
    n_steps = steps.size
    perms = _PERMUTATIONS
    if n_steps * perms > _MAX_PERM_CELLS:
        perms = max(_MIN_PERMUTATIONS, _MAX_PERM_CELLS // max(1, n_steps))
    rng = np.random.default_rng(seed)
    order = np.argsort(rng.random((perms, n_steps)), axis=1)
    is_minor = steps[order] == -1
    longest = np.zeros(perms, dtype=np.int32)
    cur = np.zeros(perms, dtype=np.int32)
    for j in range(n_steps):  # vectorised over permutations; reset-cumulative run length
        cur = np.where(is_minor[:, j], cur + 1, 0)
        longest = np.maximum(longest, cur)
    return float((np.count_nonzero(longest >= run_len) + 1) / (perms + 1))


def call_regions_geneconv(result: WindowSimilarity, analysis, params):
    """Call recombinant regions with the GENECONV-style clean-fragment test.

    The major parent is the reference winning the most windows. Each plausible donor is
    tested as the minor parent: the longest uninterrupted run of minor-matches over the
    major/minor discriminating sites, significant after Benjamini-Hochberg across donors,
    is a recombinant fragment. Overlapping fragments keep the most significant donor.
    Returns ``(regions, major, [])`` to match ``call_regions``.
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
        run_len, start, end = longest_minor_run(steps)
        if run_len < _MIN_RUN:
            continue
        p = run_pvalue_perm(steps, run_len)
        candidates.append((minor, cols, run_len, start, end, p))
    if not candidates:
        return [], major, []

    qvalues = benjamini_hochberg([c[5] for c in candidates])
    regions: list[Region] = []
    for (minor, cols, run_len, start, end, p), q in sorted(
        zip(candidates, qvalues, strict=True), key=lambda t: t[0][5]
    ):
        if p > params.alpha:
            continue
        msa_start = int(cols[start])
        msa_end = int(cols[end]) + 1
        if msa_end <= msa_start or any(
            msa_start < r.msa_end and r.msa_start < msa_end for r in regions
        ):
            continue  # keep the most significant donor for each fragment
        sim_minor = _region_similarity(result.rows, result.query, minor, msa_start, msa_end)
        sim_major = _region_similarity(result.rows, result.query, major, msa_start, msa_end)
        bp_lo = result.column_to_query(int(cols[max(0, start - 1)]))
        bp_hi = result.column_to_query(int(cols[start]))
        regions.append(Region(
            minor_parent=minor, major_parent=major,
            msa_start=msa_start, msa_end=msa_end,
            query_start=result.column_to_query(msa_start),
            query_end=result.column_to_query(msa_end),
            length_bp=msa_end - msa_start, n_windows=run_len,
            mean_sim_minor=round(sim_minor, 4), mean_sim_major=round(sim_major, 4),
            margin=round(sim_minor - sim_major, 4),
            support=1.0, pvalue=_signif(p), qvalue=_signif(q),
            breakpoint_lo=min(bp_lo, bp_hi), breakpoint_hi=max(bp_lo, bp_hi),
        ))
    regions.sort(key=lambda r: r.msa_start)
    return regions, major, []
