"""3SEQ-style triplet recombination test: statistic, exact p-value, and caller."""

from __future__ import annotations

import itertools
import random
from pathlib import Path

import numpy as np

from recomfi.recomb.analyze import analyze
from recomfi.recomb.regions import RegionParams, call_regions
from recomfi.recomb.similarity import compute_similarity
from recomfi.recomb.threeseq import (
    descent_pvalue_exact,
    max_descent,
    triplet_steps,
)

from ..conftest import write_fasta


def _brute_pvalue(m: int, n: int, depth: int) -> float:
    """P(max drawdown >= depth) by enumerating every arrangement of m +1 and n -1."""
    total = hit = 0
    for down_positions in itertools.combinations(range(m + n), n):
        steps = np.ones(m + n, dtype=np.int64)
        for p in down_positions:
            steps[p] = -1
        total += 1
        if max_descent(steps).depth >= depth:
            hit += 1
    return hit / total


def test_max_descent_finds_the_q_run() -> None:
    # P^10 then Q^8 then P^5: the drawdown is the 8-deep Q run at sites [10, 18).
    steps = np.array([1] * 10 + [-1] * 8 + [1] * 5)
    d = max_descent(steps)
    assert d.depth == 8
    assert (d.start_site, d.end_site) == (10, 18)
    assert d.m == 15 and d.n == 8


def test_descent_pvalue_exact_matches_brute_force() -> None:
    for m in range(0, 7):
        for n in range(0, 7):
            for depth in range(1, m + n + 2):
                assert abs(descent_pvalue_exact(m, n, depth) - _brute_pvalue(m, n, depth)) < 1e-9


def test_descent_pvalue_edges() -> None:
    assert descent_pvalue_exact(5, 5, 0) == 1.0   # depth 0 is always reached
    assert descent_pvalue_exact(5, 3, 4) == 0.0   # need >= depth down-steps to reach it
    # a deep clustered run is highly significant; random interspersion is not
    assert descent_pvalue_exact(50, 50, 40) < 1e-6
    assert descent_pvalue_exact(50, 50, 3) > 0.5


def test_triplet_steps_keeps_only_discriminating_sites(tmp_path: Path) -> None:
    # query ACGT; col0 P==Q (not discriminating); col1 P!=Q & query==P (+1);
    # col2 P!=Q & query==Q (-1); col3 P!=Q but query matches neither (dropped).
    msa = write_fasta(tmp_path / "m.fasta", {
        "query": "ACGT",
        "major": "ACAA",
        "minor": "ATGC",
    })
    from recomfi.recomb.similarity import _read_alignment
    rows = _read_alignment(str(msa))
    steps, cols = triplet_steps(rows, "query", "major", "minor")
    assert cols.tolist() == [1, 2]
    assert steps.tolist() == [1, -1]


def _recombinant_msa(tmp_path: Path, *, recombinant: bool) -> Path:
    rng = random.Random(7)
    base = "".join(rng.choice("ACGT") for _ in range(6000))

    def mut(seq: str, frac: float) -> str:
        chars = list(seq)
        for i in range(len(chars)):
            if rng.random() < frac:
                chars[i] = rng.choice("ACGT")
        return "".join(chars)

    a, b, other = mut(base, 0.03), mut(base, 0.03), mut(base, 0.10)
    query = list(a)
    if recombinant:
        query[2000:4000] = list(b[2000:4000])  # A backbone, B insert
    return write_fasta(tmp_path / "m.fasta",
                       {"query": "".join(query), "A": a, "B": b, "other": other})


def test_call_regions_3seq_finds_the_recombinant(tmp_path: Path) -> None:
    result = compute_similarity(str(_recombinant_msa(tmp_path, recombinant=True)),
                                "query", window_size=500, window_step=50)
    regions, major, dropped = call_regions(
        result, analyze(result), 500, RegionParams.with_defaults(500, method="3seq")
    )
    assert major == "A"
    assert dropped == []
    assert len(regions) == 1
    r = regions[0]
    assert r.minor_parent == "B" and r.major_parent == "A"
    assert r.query_start < 2200 and r.query_end > 3800   # spans the true insert
    assert r.pvalue is not None and r.pvalue < 0.01


def test_call_regions_3seq_no_false_positive_on_pure_query(tmp_path: Path) -> None:
    result = compute_similarity(str(_recombinant_msa(tmp_path, recombinant=False)),
                                "query", window_size=500, window_step=50)
    regions, major, _ = call_regions(
        result, analyze(result), 500, RegionParams.with_defaults(500, method="3seq")
    )
    assert major == "A"
    assert regions == []
