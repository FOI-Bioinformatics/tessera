"""Parent-free recombination diagnostics: PHI test and Hudson-Kaplan Rmin."""

from __future__ import annotations

import numpy as np

from recomfi.recomb.diagnostics import (
    biallelic_columns,
    hudson_kaplan_rmin,
    incompatibility_matrix,
    phi,
    phi_pvalue,
    recombination_signal,
)


def _rows(seqs: dict[str, str]) -> dict[str, np.ndarray]:
    return {label: np.frombuffer(s.encode(), dtype=np.uint8) for label, s in seqs.items()}


def test_biallelic_columns_selects_only_segregating_biallelic_sites() -> None:
    # col0 monomorphic (all A); col1 biallelic informative (AA/CC); col2 has a
    # singleton third allele (G in one seq) -> three alleles, dropped; col3 minor
    # allele in only one seq (singleton) -> not informative, dropped.
    rows = _rows({
        "s0": "AACA",
        "s1": "ACCA",
        "s2": "ACGC",
        "s3": "AAAA",
    })
    a1, a0, positions = biallelic_columns(rows)
    # Only col1 has exactly two canonical bases each carried by >= 2 sequences.
    assert positions.tolist() == [1]
    # allele indicators partition the four sequences (A:s0,s3 vs C:s1,s2) at col1.
    assert a0[0].tolist() == [1.0, 0.0, 0.0, 1.0]  # major A (count 2)
    assert a1[0].tolist() == [0.0, 1.0, 1.0, 0.0]  # minor C (count 2)


def test_four_gamete_incompatible_vs_compatible_pair() -> None:
    # col0 = AACC, col1 = ACAC -> all four gametes (AA,AC,CA,CC) -> incompatible.
    incompatible = _rows({"s0": "AA", "s1": "AC", "s2": "CA", "s3": "CC"})
    a1, a0, _ = biallelic_columns(incompatible)
    assert bool(incompatibility_matrix(a1, a0)[0, 1]) is True

    # col0 = AACC, col1 = AACC -> only AA and CC gametes -> compatible.
    compatible = _rows({"s0": "AA", "s1": "AA", "s2": "CC", "s3": "CC"})
    a1, a0, _ = biallelic_columns(compatible)
    assert bool(incompatibility_matrix(a1, a0)[0, 1]) is False


def test_incompatibility_ignores_missing_data() -> None:
    # Both columns stay biallelic-informative (each allele in >= 2 sequences), but the
    # one sequence that would supply the fourth (C,C) gamete is gapped at col1, so the
    # gap must not manufacture incompatibility.
    rows = _rows({"s0": "AA", "s1": "AC", "s2": "CA", "s3": "C-", "s4": "AC"})
    a1, a0, positions = biallelic_columns(rows)
    # col0: A in s0,s1,s4 / C in s2,s3; col1: A in s0,s2 / C in s1,s4 (s3 missing).
    assert positions.tolist() == [0, 1]
    # Gametes over the covered sequences: AA, AC, CA, AC -> only three -> compatible.
    assert bool(incompatibility_matrix(a1, a0)[0, 1]) is False


def test_hudson_kaplan_rmin_textbook_cases() -> None:
    positions = np.array([10, 20, 30, 40])
    # Two disjoint incompatible intervals (0,1) and (2,3) -> Rmin 2.
    disjoint = np.zeros((4, 4), dtype=bool)
    disjoint[0, 1] = disjoint[1, 0] = True
    disjoint[2, 3] = disjoint[3, 2] = True
    rmin, intervals = hudson_kaplan_rmin(disjoint, positions)
    assert rmin == 2
    assert intervals == [(10, 20), (30, 40)]

    # Nested incompatible intervals (0,3) and (1,2) -> one disjoint interval -> Rmin 1.
    nested = np.zeros((4, 4), dtype=bool)
    nested[0, 3] = nested[3, 0] = True
    nested[1, 2] = nested[2, 1] = True
    rmin, intervals = hudson_kaplan_rmin(nested, positions)
    assert rmin == 1
    assert intervals == [(20, 30)]  # the tighter (greedy-by-right-endpoint) interval

    # No incompatibility -> Rmin 0.
    assert hudson_kaplan_rmin(np.zeros((4, 4), dtype=bool), positions) == (0, [])


def _block_alignment(*, recombinant: bool, n_seq: int = 6, per_block: int = 12) -> dict:
    """Two genome halves. Clonal: both halves share one partition (compatible).
    Recombinant: the right half follows an incongruent partition, so left and right
    sites are four-gamete incompatible -- the signature of a breakpoint between them."""
    left = [0, 0, 0, 1, 1, 1][:n_seq]
    right = ([1, 0, 1, 0, 1, 0][:n_seq] if recombinant else left)
    cols = [left] * per_block + [right] * per_block
    seqs = {
        f"s{i}": "".join("A" if col[i] == 0 else "C" for col in cols) for i in range(n_seq)
    }
    return _rows(seqs)


def test_phi_clonal_not_significant_recombinant_significant() -> None:
    identity = lambda c: c  # noqa: E731
    clonal = recombination_signal(_block_alignment(recombinant=False), "s0", identity,
                                  window=8, seed=1)
    recomb = recombination_signal(_block_alignment(recombinant=True), "s0", identity,
                                  window=8, seed=1)
    assert clonal is not None and recomb is not None
    assert clonal.phi_p > 0.05  # no cross-block incompatibility -> nothing to detect
    assert recomb.phi_p < 0.05  # local sites compatible, distant ones not -> low PHI
    assert recomb.rmin >= 1  # at least one recombination event is forced
    assert clonal.rmin == 0


def test_phi_pvalue_reproducible_with_fixed_seed() -> None:
    rows = _block_alignment(recombinant=True)
    a1, a0, _ = biallelic_columns(rows)
    incompatible = incompatibility_matrix(a1, a0)
    p1, obs1 = phi_pvalue(incompatible, 8, seed=42)
    p2, obs2 = phi_pvalue(incompatible, 8, seed=42)
    assert p1 == p2 and obs1 == obs2
    # The observed PHI is the banded-mean incompatibility, independent of the seed.
    assert obs1 == phi(incompatible, 8)


def test_recombination_signal_none_below_min_informative() -> None:
    # Three biallelic informative columns is below the _MIN_INFORMATIVE floor (10).
    rows = _rows({"s0": "AAA", "s1": "AAA", "s2": "CCC", "s3": "CCC"})
    assert recombination_signal(rows, "s0", lambda c: c) is None
