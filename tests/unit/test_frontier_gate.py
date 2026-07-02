"""The tier/--frontier gate: frontier cases run only under --frontier."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "run_hybrids", Path(__file__).resolve().parents[2] / "validation" / "run_hybrids.py")
rh = importlib.util.module_from_spec(_SPEC)
sys.modules["run_hybrids"] = rh
_SPEC.loader.exec_module(rh)


def test_select_cases_hides_frontier_without_flag():
    cases = [{"name": "a", "tier": "must_pass"}, {"name": "f", "tier": "frontier"}]
    assert [c["name"] for c in rh._select_cases(cases, names=[], frontier=False)] == ["a"]
    assert [c["name"] for c in rh._select_cases(cases, names=[], frontier=True)] == ["f"]


def test_select_cases_default_tier_is_must_pass():
    cases = [{"name": "a"}]  # no tier -> must_pass
    assert [c["name"] for c in rh._select_cases(cases, names=[], frontier=False)] == ["a"]


def test_make_cross_hybrid_splices_by_fraction():
    a = "A" * 100
    b = "B" * 100
    q, s, e = rh.make_cross_hybrid(a, b, insert=(0.35, 0.65))
    assert q == "A" * 35 + "B" * 30 + "A" * 35 and (s, e) == (35, 65)


def test_envelope_known_limit_below_floor():
    v, _d = rh._score_frontier_envelope(0.74, detected=True, backbone_ok=False, donor_ok=False)
    assert v == "KNOWN-LIMIT"


def test_envelope_xpass_in_envelope_detected_attributed():
    v, _d = rh._score_frontier_envelope(0.83, detected=True, backbone_ok=True, donor_ok=True)
    assert v == "XPASS"


def test_envelope_xfail_in_envelope_not_attributed():
    v, _d = rh._score_frontier_envelope(0.83, detected=True, backbone_ok=False, donor_ok=False)
    assert v == "XFAIL"


def _reassort_result(verdict, segs):
    from tessera.reassort.assign import ReassortmentResult, SegmentAssignment
    segments = [SegmentAssignment(seg, "ds", strain, clade, ani, status)
                for seg, strain, clade, ani, status in segs]
    return ReassortmentResult(segments=segments, verdict=verdict)


def test_reassortant_verdict_xpasses(tmp_path):
    result = _reassort_result("reassortant", [
        ("HA", "A/X", "J.2.2", 100.0, "assigned"),
        ("NA", "A/Y", "B.4.2", 100.0, "assigned")])
    v, _d = rh._score_frontier_reassortant(result)
    assert v == "XPASS"


def test_clonal_verdict_xfails(tmp_path):
    result = _reassort_result("clonal", [
        ("HA", "A/X", "J.2.2", 100.0, "assigned"),
        ("NA", "A/X", "J.2.2", 100.0, "assigned")])
    v, _d = rh._score_frontier_reassortant(result)
    assert v == "XFAIL"


def test_undetermined_verdict_is_known_limit(tmp_path):
    result = _reassort_result("undetermined", [
        ("HA", "A/X", "J.2.2", 100.0, "assigned"),
        ("NA", None, None, 0.0, "unassigned")])
    v, _d = rh._score_frontier_reassortant(result)
    assert v == "KNOWN-LIMIT"
