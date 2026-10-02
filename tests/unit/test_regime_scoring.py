"""Unit tests for the regime harness's simulation and scoring (no binaries, no network)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_regimes.py"
_SPEC = importlib.util.spec_from_file_location("run_regimes", _PATH)
rg = importlib.util.module_from_spec(_SPEC)
sys.modules["run_regimes"] = rg
_SPEC.loader.exec_module(rg)

SMALL = {"length": 6000, "distance": 0.004, "tract": (2000, 3000)}


def _identity(a: str, b: str, lo: int, hi: int) -> float:
    return sum(1 for i in range(lo, hi) if a[i] == b[i]) / (hi - lo)


# --- simulation ------------------------------------------------------------

def test_near_identical_simulation_is_deterministic_and_aligned():
    a, tract = rg.simulate_near_identical(3, recombinant=True, **SMALL)
    b, _ = rg.simulate_near_identical(3, recombinant=True, **SMALL)
    assert a == b and tract == SMALL["tract"]
    assert rg.simulate_near_identical(4, recombinant=True, **SMALL)[0] != a
    assert rg.QUERY in a and set(rg.NEAR_REFS) <= set(a)
    assert len({len(s) for s in a.values()}) == 1


def test_near_identical_panel_is_near_identical():
    seqs, _ = rg.simulate_near_identical(1, recombinant=False, **SMALL)
    # two references each 0.4 % from the root are ~0.8 % apart
    assert 0.985 < _identity(seqs["R0"], seqs["R1"], 0, SMALL["length"]) < 0.998


def test_near_identical_recombinant_carries_the_donor_over_the_tract_only():
    seqs, tract = rg.simulate_near_identical(1, recombinant=True, **SMALL)
    lo, hi = tract
    query = seqs[rg.QUERY]
    assert query[lo:hi] == seqs["R1"][lo:hi]
    assert _identity(query, seqs["R0"], 0, lo) > _identity(query, seqs["R1"], 0, lo)


def test_near_identical_clonal_has_no_tract():
    seqs, tract = rg.simulate_near_identical(1, recombinant=False, **SMALL)
    assert tract is None
    assert _identity(seqs[rg.QUERY], seqs["R0"], 0, SMALL["length"]) > 0.998


def test_sibling_tract_differs_from_the_query_only_over_the_span():
    seqs, (lo, hi) = rg.simulate_sibling_tract(
        2, length=6000, distance=0.05, tract=(2000, 3000), extra=300
    )
    query, sibling = seqs[rg.QUERY], seqs["sibling"]
    assert (lo, hi) == (3000, 3300)
    assert query[:lo] == sibling[:lo] and query[hi:] == sibling[hi:]
    assert query[lo:hi] != sibling[lo:hi]
    assert query[lo:hi] == seqs["parent_A"][lo:hi]
    assert sibling[hi:] == seqs["parent_A"][hi:]  # identical downstream of the span


# --- scoring ---------------------------------------------------------------

def _row(**kw):
    row = {"minor_parent": "R1", "query_start": "100", "query_end": "200",
           "donor_absent": "no", "methods": "3seq,maxchi"}
    row.update(kw)
    return row


def test_score_tract_reports_detection_donor_and_error():
    got = rg.score_tract([_row(query_start="90", query_end="230")], (100, 200), donor="R1")
    assert got == {"detected": True, "hmm": False, "donor_ok": True, "breakpoint_error": 20}


def test_score_tract_notes_when_the_hmm_was_among_the_callers():
    rows = [_row(methods="hmm,3seq")]
    assert rg.score_tract(rows, (100, 200), donor="R1")["hmm"] is True


def test_score_tract_miss_and_wrong_donor():
    assert rg.score_tract([_row(query_start="500", query_end="600")], (100, 200),
                          donor="R1")["detected"] is False
    assert rg.score_tract([_row(minor_parent="R3")], (100, 200),
                          donor="R1")["donor_ok"] is False


def test_score_tract_ignores_donor_absent_rows():
    rows = [_row(donor_absent="yes", methods="")]
    assert rg.score_tract(rows, (100, 200), donor="R1")["detected"] is False


def test_score_clonal_counts_present_regions_per_caller():
    rows = [_row(methods="hmm"), _row(methods="hmm,3seq"), _row(donor_absent="yes")]
    assert rg.score_clonal(rows) == (2, {"hmm": 2, "3seq": 1})
    assert rg.score_clonal([]) == (0, {})


def test_score_span_reports_the_excess_over_the_true_difference():
    rows = [_row(query_start="17900", query_end="29600"),
            _row(query_start="18007", query_end="18591")]
    got = rg.score_span(rows, (18000, 18600))
    assert got == {"detected": True, "reported_bp": 11700, "excess_bp": 11100}
    tight = rg.score_span(rows[1:], (18000, 18600))
    assert tight["excess_bp"] == 0
    assert rg.score_span([], (18000, 18600))["detected"] is False
