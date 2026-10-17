"""Bounded computation: caps that keep a large panel from exhausting memory.

Three computations grew unbounded with the panel: the PHI incompatibility matrix
(quadratic in informative columns, with no guard at all), reference clustering
(all-pairs, gated only by a *lower* bound), and MaxChi's permutation null (several
large arrays per candidate donor). `bootscan.py` already degraded gracefully; these
now follow it.
"""

from __future__ import annotations

import numpy as np

from tessera.recomb.diagnostics import (
    _MAX_MATRIX_CELLS,
    _MIN_INFORMATIVE,
    max_columns,
    subsample_columns,
)
from tessera.recomb.maxchi import _MAX_PERM_CELLS, _PERMUTATIONS, maxchi_pvalue
from tessera.recomb.regions import RegionParams

# --- PHI: the quadratic matrix -------------------------------------------

def test_small_panels_are_not_thinned() -> None:
    # The panels this was written for must behave exactly as before.
    assert subsample_columns(50) is None
    assert subsample_columns(max_columns()) is None


def test_large_panels_are_thinned_to_the_cap() -> None:
    limit = max_columns()
    keep = subsample_columns(limit * 4)
    assert keep is not None
    assert keep.size == limit
    assert limit * limit <= _MAX_MATRIX_CELLS  # the matrix now fits


def test_thinning_keeps_genomic_order() -> None:
    """Sorted, so a PHI window of n ranks still means n nearby informative sites."""
    keep = subsample_columns(max_columns() * 3)
    assert np.all(np.diff(keep) > 0)


def test_thinning_is_deterministic_and_seed_dependent() -> None:
    z = max_columns() * 3
    assert np.array_equal(subsample_columns(z, seed=0), subsample_columns(z, seed=0))
    assert not np.array_equal(subsample_columns(z, seed=0), subsample_columns(z, seed=1))


def test_max_columns_never_drops_below_the_testable_minimum() -> None:
    # A tiny cap must not thin below the point where the statistics mean anything.
    assert max_columns(cap=1) == _MIN_INFORMATIVE


# --- clustering: an upper bound as well as a lower one -------------------

def test_clustering_has_both_bounds() -> None:
    params = RegionParams.with_defaults(window_size=1000)
    assert params.cluster_min_refs < params.cluster_max_refs
    # All-pairs over a 1000-genome panel is ~500k comparisons, each allocating two
    # full-width cumulative-sum arrays.
    assert params.cluster_max_refs <= 1000


# --- MaxChi: bounded permutation null ------------------------------------

def test_maxchi_batches_without_losing_permutations() -> None:
    """Memory is bounded by batching, not by running fewer permutations.

    Capping the permutation *count* was tried first and reverted: it engaged above
    1000 discriminating sites -- an ordinary panel -- and the coarser null dropped
    the specificity harness's positive control from 3/3 to 1/3. Batching keeps the
    full null, because permutation rows are independent and their exceedance counts
    simply add.
    """
    rng = np.random.default_rng(0)
    long_steps = rng.choice([1, -1], size=4000).astype(np.int64)
    assert long_steps.size * _PERMUTATIONS > _MAX_PERM_CELLS  # batching engages
    batch = _MAX_PERM_CELLS // long_steps.size
    assert batch < _PERMUTATIONS  # ... in more than one pass

    # The p-value floor still reflects every permutation, not just one batch.
    p = maxchi_pvalue(long_steps, observed=1e9)  # unreachable -> nothing exceeds
    assert p == 1.0 / (_PERMUTATIONS + 1)


def test_maxchi_floor_is_the_same_batched_or_not() -> None:
    # A short input needs no batching; a long one does. Both resolve to the same floor.
    rng = np.random.default_rng(0)
    short = rng.choice([1, -1], size=200).astype(np.int64)
    assert short.size * _PERMUTATIONS <= _MAX_PERM_CELLS  # single pass
    assert maxchi_pvalue(short, observed=1e9) == 1.0 / (_PERMUTATIONS + 1)


# --- PHI, Hudson-Kaplan and GENECONV: peak memory -------------------------

def _random_gametes(z: int, n: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    allele1 = (rng.random((z, n)) < 0.4).astype(np.float64)
    return allele1, 1.0 - allele1


def _peak_mb(fn, *args):
    import tracemalloc

    tracemalloc.start()
    try:
        out = fn(*args)
        return out, tracemalloc.get_traced_memory()[1] / 1e6
    finally:
        tracemalloc.stop()


def test_incompatibility_matrix_does_not_hold_four_float64_matrices() -> None:
    # z = 2000: one float64 z*z matrix is 32 MB. Four live together peaked at 136 MB;
    # one float32 matrix at a time stays near 25 MB.
    from tessera.recomb.diagnostics import incompatibility_matrix

    a1, a0 = _random_gametes(2000, 40)
    matrix, peak = _peak_mb(incompatibility_matrix, a1, a0)
    assert matrix.dtype == bool and matrix.shape == (2000, 2000)
    assert peak < 60


def _rmin_by_sorted_pairs(incompatible: np.ndarray, positions: np.ndarray):
    """The textbook greedy: every incompatible pair, sorted by right endpoint."""
    pairs = [(int(i), int(j)) for i, j in zip(*np.where(np.triu(incompatible, 1)), strict=True)]
    pairs.sort(key=lambda p: (p[1], p[0]))
    intervals, last = [], -1
    for i, j in pairs:
        if i >= last:
            intervals.append((int(positions[i]), int(positions[j])))
            last = j
    return len(intervals), intervals


def test_hudson_kaplan_matches_the_sorted_pair_greedy() -> None:
    from tessera.recomb.diagnostics import hudson_kaplan_rmin

    rng = np.random.default_rng(1)
    for density in (0.0, 0.02, 0.3, 1.0):
        upper = np.triu(rng.random((60, 60)) < density, 1)
        matrix = upper | upper.T
        positions = np.arange(60) * 7
        assert hudson_kaplan_rmin(matrix, positions) == _rmin_by_sorted_pairs(matrix, positions)


def test_hudson_kaplan_does_not_materialise_the_pairs() -> None:
    # A dense 2000 x 2000 matrix holds ~2 million incompatible pairs; as Python tuples
    # they cost ~380 MB. The column-wise greedy needs almost nothing beyond the matrix.
    from tessera.recomb.diagnostics import hudson_kaplan_rmin

    matrix = np.ones((2000, 2000), dtype=bool)
    np.fill_diagonal(matrix, False)
    (rmin, intervals), peak = _peak_mb(hudson_kaplan_rmin, matrix, np.arange(2000))
    assert rmin == len(intervals) == 1999  # every adjacent pair is its own interval
    assert peak < 20


def test_geneconv_permutations_are_batched() -> None:
    # 10 000 discriminating sites: 4000 permutations, 40 million cells. Unbatched, the
    # random draw, the argsort and the mask peaked at about 680 MB per candidate donor.
    from tessera.recomb import geneconv

    steps = np.random.default_rng(2).choice([1, -1], size=10_000).astype(np.int64)
    p, peak = _peak_mb(geneconv.run_pvalue_perm, steps, 12)
    assert 0.0 < p <= 1.0
    assert peak < 300


def test_geneconv_batching_does_not_change_the_p_value(monkeypatch) -> None:
    # The batches draw the same random stream in the same order as one full draw, so
    # the p-value is identical. Shrink the batch so a small input takes several.
    from tessera.recomb import geneconv

    steps = np.random.default_rng(3).choice([1, -1], size=200).astype(np.int64)
    perms = geneconv._PERMUTATIONS
    order = np.argsort(np.random.default_rng(0).random((perms, steps.size)), axis=1)
    is_minor = steps[order] == -1
    longest = np.zeros(perms, dtype=np.int32)
    cur = np.zeros(perms, dtype=np.int32)
    for j in range(steps.size):
        cur = np.where(is_minor[:, j], cur + 1, 0)
        longest = np.maximum(longest, cur)
    unbatched = float((np.count_nonzero(longest >= 6) + 1) / (perms + 1))

    monkeypatch.setattr(geneconv, "_BATCH_CELLS", 50_000, raising=False)
    assert geneconv.run_pvalue_perm(steps, 6) == unbatched
