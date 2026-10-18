"""3SEQ-style triplet recombination test: statistic, exact p-value, and caller."""

from __future__ import annotations

import itertools
import math
from pathlib import Path

import numpy as np
import pytest

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


def _call(path: Path, method: str):
    result = compute_similarity(str(path), "query", window_size=500, window_step=50)
    return call_regions(result, analyze(result), 500,
                        RegionParams.with_defaults(500, method=method))


def test_descent_pvalue_switches_to_permutation_above_the_cell_cap(monkeypatch) -> None:
    """Guard: which branch answered is reported, and the two branches agree numerically.

    The drawdown is deliberately shallow. A deep one has an exact p-value far below the
    permutation estimator's resolution of ``1 / (_PERMUTATIONS + 1)``, so the two cannot
    be compared at all and any tolerance wide enough to pass is wide enough to hide a
    broken fallback. Here the exact p-value is about 0.048, where 20 000 permutations
    have a standard error of about 0.0015, so the tolerance below is four standard
    errors: close agreement, not a formality.
    """
    from tessera.recomb import threeseq

    steps = np.array([1] * 5 + [-1] * 2 + [1] * 35, dtype=np.int64)
    d = max_descent(steps)
    exact, kind = descent_pvalue(steps, d)
    assert kind == "exact"
    assert 0.01 < exact < 0.25, exact  # within the permutation estimator's resolution
    monkeypatch.setattr(threeseq, "_MAX_EXACT_CELLS", d.m * d.n * d.depth - 1)
    perm, kind = descent_pvalue(steps, d)
    assert kind == "permutation"
    standard_error = (exact * (1 - exact) / threeseq._PERMUTATIONS) ** 0.5
    assert perm == pytest.approx(exact, abs=4 * standard_error)


def test_a_tract_at_the_first_site_has_a_breakpoint_interval_from_the_start(tmp_path) -> None:
    # The tract starts at the first discriminating site, so there is no site before it;
    # the interval ran from that site to itself (zero width) instead of from column 0.
    import random

    from ..factories import _mutate

    rng = random.Random(5)
    a = "".join(rng.choice("ACGT") for _ in range(6000))
    b, other = _mutate(a, 0.05, rng), _mutate(a, 0.15, rng)
    query = b[:1500] + a[1500:]
    path = write_fasta(tmp_path / "terminal.fasta",
                       {"query": query, "A": a, "B": b, "other": other})
    for method in ("3seq", "maxchi", "geneconv"):
        regions, _, _ = _call(path, method)
        terminal = [r for r in regions if r.minor_parent == "B"]
        assert len(terminal) == 1, method
        assert terminal[0].breakpoint_lo == 0, method
        assert terminal[0].breakpoint_hi > 0, method


# --- pre-test floors, the BH family, the exact/permutation switch ------------

def _three_site_alignment(tmp_path: Path, *, decoy: bool = False, sites: int = 3) -> Path:
    """Major P and minor Q differ at 189 columns; the query is P except at ``sites``
    consecutive ones, where it carries Q. With ``decoy`` a third reference R differs from
    P at 30 other columns, where the query matches P (a donor with no descent at all)."""
    rng = np.random.default_rng(0)
    width = 4000
    p = rng.integers(0, 4, size=width)
    q = p.copy()
    q_cols = np.sort(rng.choice(np.arange(0, width, 2), size=189, replace=False))
    q[q_cols] = (p[q_cols] + 1) % 4
    query = p.copy()
    query[q_cols[100:100 + sites]] = q[q_cols[100:100 + sites]]
    seqs = {"query": query, "P": p, "Q": q}
    if decoy:
        r = p.copy()
        r_cols = rng.choice(np.arange(1, width, 2), size=30, replace=False)  # not Q's
        r[r_cols] = (p[r_cols] + 2) % 4
        seqs["R"] = r
    bases = np.frombuffer(b"ACGT", dtype=np.uint8)
    text = {k: bases[v].tobytes().decode() for k, v in seqs.items()}
    return write_fasta(tmp_path / "three.fasta", text)


def test_a_three_site_tract_is_tested_not_discarded_before_the_p_value(tmp_path) -> None:
    # Three consecutive donor matches among 189 sites: exact p = 187 / C(189, 3) = 1.7e-4.
    # The depth floor of 4 discarded the triplet before any p-value was computed.
    regions, major, _ = _call(_three_site_alignment(tmp_path), "3seq")
    assert major == "P"
    assert [r.minor_parent for r in regions] == ["Q"]
    assert math.isclose(regions[0].pvalue, 187 / math.comb(189, 3), rel_tol=1e-9)


def test_the_bh_family_counts_every_tested_donor(tmp_path) -> None:
    # R has 30 discriminating sites and no descent (p = 1). It was pruned by the depth
    # floor before Benjamini-Hochberg, so the family had one member and q equalled p.
    path = _three_site_alignment(tmp_path, decoy=True, sites=8)
    for method in ("3seq", "maxchi"):
        regions, _, _ = _call(path, method)
        donors = [r for r in regions if r.minor_parent == "Q"]
        assert len(donors) == 1, method
        assert math.isclose(donors[0].qvalue, 2 * donors[0].pvalue, rel_tol=1e-9), method


def test_descent_pvalue_exact_agrees_with_the_permutation_null_on_a_long_shallow_walk() -> None:
    """The anti-diagonal DP against the module's independent permutation estimator.

    ``test_descent_pvalue_exact_matches_brute_force`` pins the DP against enumeration for
    ``m, n <= 6``, which exercises every boundary row of the vectorised fill but only at
    small ``t``. The rewrite exists for the opposite regime -- ``m`` in the thousands with
    a shallow descent -- where an off-by-one in a boundary row would change the magnitude
    and no enumeration can check it. The permutation estimator is a separate
    implementation, so the two agreeing there is a real cross-check; the tolerance is four
    standard errors of its 20 000 draws.
    """
    from tessera.recomb import threeseq

    for m, n, depth in ((1217, 74, 3), (1795, 116, 3), (600, 40, 2)):
        exact = descent_pvalue_exact(m, n, depth)
        steps = np.array([1] * m + [-1] * n, dtype=np.int64)
        perm = threeseq._descent_pvalue_perm(steps, depth)
        standard_error = (perm * (1 - perm) / threeseq._PERMUTATIONS) ** 0.5
        assert exact == pytest.approx(perm, abs=4 * standard_error), (m, n, depth, exact, perm)


def test_descent_pvalue_exact_closed_forms_at_the_boundaries() -> None:
    # Cases with an answer that can be written down, at a scale enumeration cannot reach.
    assert descent_pvalue_exact(5000, 300, 1) == 1.0    # any down-step is a drawdown of 1
    assert descent_pvalue_exact(5000, 0, 1) == 0.0      # no down-step, no drawdown
    assert descent_pvalue_exact(0, 300, 300) == 1.0     # all down-steps: the drawdown is n
    assert descent_pvalue_exact(0, 300, 301) == 0.0     # ... and cannot exceed n
    assert descent_pvalue_exact(5000, 40, 41) == 0.0    # depth > n is unreachable
    # depth == n at a scale the closed-form test above does not cover.
    assert math.isclose(descent_pvalue_exact(3000, 4, 4), 3001 / math.comb(3004, 4),
                        rel_tol=1e-9)
