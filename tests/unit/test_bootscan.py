"""Bootscan caller: bootstrap support and region calling."""

from __future__ import annotations

import random
from pathlib import Path

import numpy as np

from tessera.recomb.analyze import analyze
from tessera.recomb.bootscan import bootstrap_support
from tessera.recomb.regions import RegionParams, call_regions
from tessera.recomb.similarity import compute_similarity

from ..conftest import recombinant_msa, write_fasta


def test_bootstrap_support_favours_the_closer_reference() -> None:
    # ref 0 matches the query at 18/20 columns, ref 1 at 2/20: support concentrates on 0
    match = np.zeros((2, 20), dtype=bool)
    match[0, :18] = True
    match[1, 18:] = True
    comparable = np.ones((2, 20), dtype=bool)
    sup = bootstrap_support(match, comparable, seed=1)
    assert sup[0] > 0.95 and sup[1] < 0.05
    assert abs(sup.sum() - 1.0) < 1e-9
    # reproducible with a fixed seed
    assert np.array_equal(bootstrap_support(match, comparable, seed=3),
                          bootstrap_support(match, comparable, seed=3))


def test_call_regions_bootscan_finds_the_recombinant(tmp_path: Path) -> None:
    result = compute_similarity(str(recombinant_msa(tmp_path, recombinant=True)),
                                "query", window_size=500, window_step=50)
    regions, major, _ = call_regions(
        result, analyze(result), 500, RegionParams.with_defaults(500, method="bootscan")
    )
    assert major == "A"
    assert len(regions) == 1
    r = regions[0]
    assert r.minor_parent == "B" and r.major_parent == "A"
    assert r.query_start < 2300 and r.query_end > 3700  # overlaps the true insert
    assert r.support is not None and r.support >= 0.70  # bootstrap support
    assert r.methods == ("bootscan",)


def test_call_regions_bootscan_no_false_positive(tmp_path: Path) -> None:
    result = compute_similarity(str(recombinant_msa(tmp_path, recombinant=False)),
                                "query", window_size=500, window_step=50)
    regions, major, _ = call_regions(
        result, analyze(result), 500, RegionParams.with_defaults(500, method="bootscan")
    )
    assert major == "A"
    assert regions == []


def clonal_msa(tmp_path: Path) -> Path:
    """Non-recombinant query against four near-equidistant relatives.

    Each reference is an independent 2 % mutation of one ancestor, and so is the
    query: no recombination exists anywhere, so no region may be called. This is
    the panel shape -- several near-equidistant candidates -- in which a caller
    that merely counts window votes flips winner between windows by chance.
    """
    rng = random.Random(11)
    base = "".join(rng.choice("ACGT") for _ in range(6000))

    def mut(seq: str, frac: float) -> str:
        chars = list(seq)
        for i in range(len(chars)):
            if rng.random() < frac:
                chars[i] = rng.choice("ACGT")
        return "".join(chars)

    records = {"query": mut(base, 0.02)}
    for name in ("A", "B", "C", "D"):
        records[name] = mut(base, 0.02)
    return write_fasta(tmp_path / "clonal.fasta", records)


def test_bootscan_region_carries_a_significance_value(tmp_path: Path) -> None:
    """A called region must expose the p/q-value the FDR gate is applied to."""
    result = compute_similarity(str(recombinant_msa(tmp_path, recombinant=True)),
                                "query", window_size=500, window_step=50)
    regions, _, _ = call_regions(
        result, analyze(result), 500, RegionParams.with_defaults(500, method="bootscan")
    )
    assert len(regions) == 1
    assert regions[0].pvalue is not None
    assert regions[0].qvalue is not None
    assert 0.0 <= regions[0].qvalue <= 1.0


def test_bootscan_does_not_call_a_clonal_panel(tmp_path: Path) -> None:
    """Stochastic window-vote flips among near-identical relatives are not events."""
    result = compute_similarity(str(clonal_msa(tmp_path)), "query",
                                window_size=500, window_step=50)
    regions, _, _ = call_regions(
        result, analyze(result), 500, RegionParams.with_defaults(500, method="bootscan")
    )
    assert regions == []


# --- the run-length null and lineage clusters --------------------------------

def test_block_permutation_keeps_every_win_when_the_last_block_is_partial() -> None:
    # 11 windows, block 10: the second block is one window long. Every window is a win,
    # so every arrangement has a run of 11 and p must be 1. Cutting the shuffled vector
    # back to 11 cells dropped ten wins whenever the short block came first.
    from tessera.recomb.bootscan import run_pvalue

    assert run_pvalue(np.ones(11, dtype=bool), 11, block=10) == 1.0


def test_block_permutations_hold_the_win_count_and_length() -> None:
    from tessera.recomb.bootscan import _block_permutations

    wins = np.zeros(103, dtype=bool)
    wins[[0, 1, 2, 50, 51, 100, 101, 102]] = True
    perms = _block_permutations(wins, 10, 400, np.random.default_rng(0))
    assert perms.shape == (400, 103)
    assert set(perms.sum(axis=1)) == {8}


def test_block_permutations_with_fewer_windows_than_one_block() -> None:
    from tessera.recomb.bootscan import _block_permutations, run_pvalue

    wins = np.array([1, 0, 1, 1, 0], dtype=bool)
    perms = _block_permutations(wins, 5, 50, np.random.default_rng(0))
    assert perms.shape == (50, 5)
    assert (perms == wins).all()  # one block: the only arrangement is the observed one
    assert run_pvalue(wins, 2, block=50) == 1.0  # block is clamped to n_win


def test_block_permutations_unchanged_when_windows_are_a_multiple_of_the_block() -> None:
    from tessera.recomb.bootscan import _block_permutations

    wins = np.zeros(30, dtype=bool)
    wins[[3, 4, 5, 20]] = True
    rng = np.random.default_rng(7)
    got = _block_permutations(wins, 10, 20, rng)
    ref_rng = np.random.default_rng(7)
    order = np.argsort(ref_rng.random((20, 3)), axis=1)
    expected = wins.reshape(3, 10)[order].reshape(20, 30)
    assert (got == expected).all()


def _redundant_donor_msa(tmp_path: Path) -> Path:
    """recombinant_msa with the donor lineage B represented by two genomes ~0.4 % apart,
    each 0.2 % from the B sequence the query carries (neither is closer by design)."""
    from tessera.core.io import read_fasta

    from ..factories import _mutate

    rng = random.Random(3)
    records = dict(read_fasta(recombinant_msa(tmp_path, recombinant=True)))
    original_b = records["B"]
    records["B"] = _mutate(original_b, 0.002, rng)
    records["B2"] = _mutate(original_b, 0.002, rng)
    return write_fasta(tmp_path / "redundant.fasta", records)


def test_bootscan_counts_support_per_lineage_cluster(tmp_path: Path) -> None:
    # Two near-identical genomes of the donor lineage split the bootstrap support between
    # them, so neither reaches 70 % and the tract is lost. Counted per cluster, it is not.
    result = compute_similarity(str(_redundant_donor_msa(tmp_path)), "query",
                                window_size=500, window_step=50)
    split, _, _ = call_regions(
        result, analyze(result), 500, RegionParams.with_defaults(500, method="bootscan")
    )
    assert split == []
    pooled, major, _ = call_regions(
        result, analyze(result), 500,
        RegionParams.with_defaults(500, method="bootscan", reference_clusters=[["B", "B2"]]),
    )
    assert major == "A"
    assert len(pooled) == 1
    assert pooled[0].minor_parent in {"B", "B2"}
    assert pooled[0].query_start < 2300 and pooled[0].query_end > 3700
