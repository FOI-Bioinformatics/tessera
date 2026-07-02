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


def _reassort_regions(tmp_path, start, end):
    cols = ["minor_parent", "major_parent", "query_start", "query_end", "methods",
            "donor_absent", "donor_undercovered"]
    row = ["m", "M", str(start), str(end), "hmm", "no", "no"]
    (tmp_path / "recombination_regions.tsv").write_text(
        "\t".join(cols) + "\n" + "\t".join(row) + "\n")
    (tmp_path / "recombination_profile.tsv").write_text("phi_p\trmin\n-\t0\n")


def test_reassortant_junction_hit(tmp_path):
    _reassort_regions(tmp_path, 1400, 1700)
    v, _d = rh._score_frontier_reassortant(tmp_path, junction=1500)
    assert v == "XPASS"


def test_reassortant_junction_miss(tmp_path):
    _reassort_regions(tmp_path, 100, 300)
    v, _d = rh._score_frontier_reassortant(tmp_path, junction=1500)
    assert v == "XFAIL"
