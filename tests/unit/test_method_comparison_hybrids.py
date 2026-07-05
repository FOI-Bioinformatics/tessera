"""Unit tests for the GENECONV gap-probe's pure aggregation (no binaries/network)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_method_comparison_hybrids.py"
_SPEC = importlib.util.spec_from_file_location("run_method_comparison_hybrids", _PATH)
gp = importlib.util.module_from_spec(_SPEC)
sys.modules["run_method_comparison_hybrids"] = gp
_SPEC.loader.exec_module(gp)

_METHODS = ("geneconv", "maxchi")


def _pos(tessera_pass, openrdp):
    return {"is_negative": False, "tessera_pass": tessera_pass, "openrdp": openrdp}


def _neg(tessera_pass, openrdp):
    return {"is_negative": True, "tessera_pass": tessera_pass, "openrdp": openrdp}


def test_gap_summary_counts_rescues_and_false_calls():
    rows = [
        # positive Tessera missed, GENECONV detected -> a rescue for geneconv, not maxchi
        _pos(False, {"geneconv": True, "maxchi": False}),
        # positive Tessera passed -> not a rescue even though geneconv fired
        _pos(True, {"geneconv": True, "maxchi": True}),
        # positive Tessera missed and geneconv also missed -> no rescue
        _pos(False, {"geneconv": False, "maxchi": False}),
        # negative: geneconv false-calls, maxchi clean
        _neg(True, {"geneconv": True, "maxchi": False}),
    ]
    summ = gp.gap_summary(rows, _METHODS)
    assert summ["n_positive"] == 3 and summ["n_negative"] == 1
    assert summ["tessera_missed"] == 2
    assert summ["per_method"]["geneconv"] == {"rescues": 1, "false_calls": 1}
    assert summ["per_method"]["maxchi"] == {"rescues": 0, "false_calls": 0}


def test_gap_summary_ignores_rows_where_openrdp_did_not_run():
    rows = [_pos(False, None), _neg(True, None)]
    summ = gp.gap_summary(rows, _METHODS)
    assert summ["n_positive"] == 0 and summ["n_negative"] == 0
    assert summ["per_method"]["geneconv"] == {"rescues": 0, "false_calls": 0}


def test_gap_summary_tessera_false_calls_on_negatives():
    rows = [_neg(False, {"geneconv": False, "maxchi": False})]
    summ = gp.gap_summary(rows, _METHODS)
    assert summ["tessera_false_calls"] == 1


def test_tessera_callers_positive_uses_span_overlap():
    regions = [
        {"query_start": "100", "query_end": "200", "donor_absent": "no", "methods": "hmm,3seq"},
        {"query_start": "900", "query_end": "950", "donor_absent": "no", "methods": "maxchi"},
    ]
    # true span 120-180 overlaps only the first region
    fired = gp.tessera_callers(regions, [(120, 180, "B")], is_negative=False)
    assert fired == {"hmm", "3seq"}


def test_tessera_callers_negative_counts_any_present_region():
    regions = [{"query_start": "10", "query_end": "40", "donor_absent": "no",
                "methods": "bootscan"}]
    fired = gp.tessera_callers(regions, [], is_negative=True)
    assert fired == {"bootscan"}


def test_tessera_callers_ignores_donor_absent_regions():
    regions = [{"query_start": "10", "query_end": "40", "donor_absent": "yes", "methods": "hmm"}]
    assert gp.tessera_callers(regions, [(0, 100, "B")], is_negative=False) == set()


def test_select_adversarial_returns_only_adversarial_tier():
    cases = gp._select(names=[], include_core=False, adversarial=True)
    assert cases, "expected adversarial-tier cases to be defined"
    assert all(c.get("tier") == "adversarial" for c in cases)
    assert all(c.get("pattern") == "AB_micro" for c in cases)


def test_select_default_excludes_adversarial_tier():
    cases = gp._select(names=[], include_core=False, adversarial=False)
    assert all(c.get("tier", "must_pass") == "must_pass" for c in cases)
