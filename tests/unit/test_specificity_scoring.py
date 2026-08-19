"""Unit tests for the specificity harness's pure logic (no binaries, no network)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_specificity.py"
_SPEC = importlib.util.spec_from_file_location("run_specificity", _PATH)
rs = importlib.util.module_from_spec(_SPEC)
sys.modules["run_specificity"] = rs
_SPEC.loader.exec_module(rs)


# --- simulation ------------------------------------------------------------

def test_clonal_simulation_is_deterministic():
    a = rs.simulate_clonal("clean", seed=3)
    b = rs.simulate_clonal("clean", seed=3)
    assert a == b
    assert rs.simulate_clonal("clean", seed=4) != a


def test_clonal_simulation_has_a_query_and_aligned_rows():
    seqs = rs.simulate_clonal("clean", seed=1)
    assert rs.QUERY in seqs
    assert len({len(s) for s in seqs.values()}) == 1  # already aligned: equal lengths


def test_every_scenario_simulates():
    for scenario in rs.SCENARIOS:
        seqs = rs.simulate_clonal(scenario, seed=0)
        assert len(seqs) > 2


def test_recombinant_simulation_splices_the_donor_tract():
    seqs, (lo, hi) = rs.simulate_recombinant(seed=1)
    query = seqs[rs.QUERY]
    # inside the tract the query must match a clade-B tip better than a clade-A tip
    def ident(a, b, s, e):
        return sum(1 for i in range(s, e) if a[i] == b[i]) / (e - s)
    assert ident(query, seqs["B0"], lo, hi) > ident(query, seqs["A0"], lo, hi)
    assert ident(query, seqs["A0"], 0, lo) > ident(query, seqs["B0"], 0, lo)


# --- scoring ---------------------------------------------------------------

def _row(**kw):
    row = {"minor_parent": "X", "query_start": "100", "query_end": "200",
           "donor_absent": "no", "methods": "hmm"}
    row.update(kw)
    return row


def test_score_negative_ignores_donor_absent_rows():
    """A coverage gap is a 'no close reference here' flag, not a recombination claim.

    validation/run_hybrids.py::_score_neg_pure excludes them; so must this.
    """
    rows = [_row(), _row(donor_absent="yes"), _row(methods="hmm,maxchi")]
    n_false, per_caller = rs.score_negative(rows)
    assert n_false == 2
    assert per_caller["hmm"] == 2
    assert per_caller["maxchi"] == 1


def test_score_negative_on_a_clean_run():
    assert rs.score_negative([]) == (0, {})


def test_score_positive_detects_and_attributes():
    rows = [_row(minor_parent="B2", query_start="4000", query_end="7000")]
    got = rs.score_positive(rows, (4000, 7000), donor_prefix="B")
    assert got["detected"] and got["donor_ok"]
    assert got["breakpoint_error"] == 0


def test_score_positive_wrong_donor_is_not_attributed():
    rows = [_row(minor_parent="C1", query_start="4000", query_end="7000")]
    got = rs.score_positive(rows, (4000, 7000), donor_prefix="B")
    assert got["detected"] and not got["donor_ok"]


def test_score_positive_missed_tract():
    rows = [_row(minor_parent="B1", query_start="10", query_end="20")]
    got = rs.score_positive(rows, (4000, 7000), donor_prefix="B")
    assert not got["detected"]
    assert got["breakpoint_error"] is None


def test_false_positive_rate_carries_a_confidence_interval():
    """A rate from a handful of replicates is not a point estimate; show the interval."""
    rate, lo, hi = rs.false_positive_rate(2, 40)
    assert rate == 0.05
    assert 0.0 <= lo < rate < hi <= 1.0


def test_false_positive_rate_with_no_runs():
    assert rs.false_positive_rate(0, 0) == (0.0, 0.0, 1.0)
