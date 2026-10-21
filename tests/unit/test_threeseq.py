"""3SEQ-style triplet recombination test: statistic, exact p-value, and caller."""

from __future__ import annotations

import itertools
import math
from pathlib import Path

import numpy as np

from tessera.recomb.analyze import analyze
from tessera.recomb.regions import RegionParams, call_regions
from tessera.recomb.similarity import compute_similarity
from tessera.recomb.threeseq import (
    descent_pvalue,
    descent_pvalue_exact,
    max_descent,
    triplet_steps,
)

from ..conftest import recombinant_msa, write_fasta


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
    from tessera.recomb.similarity import _read_alignment
    rows = _read_alignment(str(msa))
    steps, cols = triplet_steps(rows, "query", "major", "minor")
    assert cols.tolist() == [1, 2]
    assert steps.tolist() == [1, -1]


def test_call_regions_3seq_finds_the_recombinant(tmp_path: Path) -> None:
    result = compute_similarity(str(recombinant_msa(tmp_path, recombinant=True)),
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
    result = compute_similarity(str(recombinant_msa(tmp_path, recombinant=False)),
                                "query", window_size=500, window_step=50)
    regions, major, _ = call_regions(
        result, analyze(result), 500, RegionParams.with_defaults(500, method="3seq")
    )
    assert major == "A"
    assert regions == []


def test_max_descent_starts_at_the_last_maximum_before_the_trough() -> None:
    # The walk reaches its maximum (2) at site 2, dips, returns to 2 at site 4 and only
    # then descends. The donor run is sites [4, 8); sites 2-3 are a balanced excursion.
    steps = np.array([1, 1, -1, 1, -1, -1, -1, -1])
    d = max_descent(steps)
    assert d.depth == 4
    assert (d.start_site, d.end_site) == (4, 8)
    assert set(steps[d.start_site:d.end_site]) == {-1}


def test_descent_pvalue_exact_keeps_small_values() -> None:
    # depth == n: every down-step must be consecutive, so the count of favourable
    # arrangements is the m + 1 positions of that block among the up-steps.
    for m, n in ((200, 100), (50, 30), (12, 9)):
        expected = (m + 1) / math.comb(m + n, n)
        got = descent_pvalue_exact(m, n, n)
        assert got > 0.0
        assert math.isclose(got, expected, rel_tol=1e-9)


def _two_donor_msa(tmp_path: Path) -> Path:
    """Query = backbone A with a 1 kb tract carried by two donors.

    Both donors share the tract (40 sites that differ from A) and so give the same
    descent depth; they differ in how many extra sites outside the tract differ from A
    (50 versus 150), which sets the number of major-parent matches and so the p-value.
    """
    import random

    rng = random.Random(11)
    n = 3000
    a = [rng.choice("ACGT") for _ in range(n)]

    def other(c: str) -> str:
        return rng.choice([x for x in "ACGT" if x != c])

    tract = range(1000, 2000)
    tract_sites = rng.sample(tract, 40)
    outside = [i for i in range(n) if i not in tract]
    out_sites = rng.sample(outside, 150)
    base = {i: other(a[i]) for i in tract_sites}
    extra = {i: other(a[i]) for i in out_sites}

    def donor(n_extra: int) -> str:
        seq = list(a)
        for i, c in base.items():
            seq[i] = c
        for i in out_sites[:n_extra]:
            seq[i] = extra[i]
        return "".join(seq)

    query = list(a)
    for i, c in base.items():
        query[i] = c
    return write_fasta(tmp_path / "m.fasta", {
        "query": "".join(query), "A": "".join(a), "D1": donor(50), "D2": donor(150),
    })


def test_3seq_reports_the_more_significant_donor_among_tiny_exact_pvalues(
    tmp_path: Path,
) -> None:
    # Candidates whose exact p-values are far below 1e-16 used to tie at 0.0, so the
    # donor kept for a tract was whichever came first in window-rank order. They are now
    # ordered by their p-values.
    from tessera.recomb.analyze import rank_datasets

    result = compute_similarity(str(_two_donor_msa(tmp_path)), "query",
                                window_size=500, window_step=50)
    analysis = analyze(result)
    ranked = rank_datasets(analysis, len(result.similarities))
    assert ranked[0] == "A"
    assert ranked.index("D1") < ranked.index("D2")   # D1 is tried first ...

    pvalues = {}
    for donor in ("D1", "D2"):
        steps, _ = triplet_steps(result.rows, "query", "A", donor)
        d = max_descent(steps)
        p, kind = descent_pvalue(steps, d)
        assert kind == "exact"
        assert 0.0 < p < 1e-16
        pvalues[donor] = p
    assert pvalues["D2"] < pvalues["D1"]             # ... but D2 is more significant

    regions, major, _ = call_regions(
        result, analysis, 500, RegionParams.with_defaults(500, method="3seq")
    )
    assert major == "A"
    assert [r.minor_parent for r in regions] == ["D2"]
