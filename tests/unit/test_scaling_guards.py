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
