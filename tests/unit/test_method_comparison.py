"""Unit tests for the method-comparison harness's pure helpers (no binaries)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_method_comparison.py"
_SPEC = importlib.util.spec_from_file_location("run_method_comparison", _PATH)
mc = importlib.util.module_from_spec(_SPEC)
sys.modules["run_method_comparison"] = mc
_SPEC.loader.exec_module(mc)


# -- Tessera-side detection ------------------------------------------------------------------

def test_no_regions_is_not_detected():
    assert mc.any_region_detected([]) is False


def test_present_region_is_detected():
    assert mc.any_region_detected([{"minor_parent": "X", "donor_absent": "no"}]) is True


def test_only_donor_absent_regions_is_not_detected():
    regions = [{"minor_parent": "X", "donor_absent": "yes"},
               {"minor_parent": "Y", "donor_absent": "yes"}]
    assert mc.any_region_detected(regions) is False


def test_mixed_regions_is_detected():
    regions = [{"minor_parent": "X", "donor_absent": "yes"},
               {"minor_parent": "Y", "donor_absent": "no"}]
    assert mc.any_region_detected(regions) is True


# -- alignment sanitising --------------------------------------------------------------------

def test_sanitize_maps_ambiguity_to_n_and_preserves_headers():
    src = ">seq1\nACGTRYSWN-\n>seq2\nacgtn\n"
    out = mc.sanitize_alignment(src)
    lines = out.splitlines()
    assert lines[0] == ">seq1"
    assert lines[1] == "ACGTNNNNN-"   # R,Y,S,W -> N; existing N/- kept; lower-cased upper
    assert lines[2] == ">seq2"
    assert lines[3] == "ACGTN"


# -- OpenRDP-side parsing + detection --------------------------------------------------------

_CSV = (
    "Method,Start,End,Recombinant,Parent1,Parent2,Pvalue\n"
    "Threeseq,3186,9646,GT2a_JFH1,GT1b_Con1,RF1_2k1b,3.3e-162\n"
    "Bootscan,2280,3160,RF1_2k1b,GT1b_Con1,GT2a_JFH1,6.1e-05\n"
    "Maxchi,521,569,RF1_2k1b,GT1a_H77,GT2a_JFH1,0.0033\n"
    "Geneconv,3187,9410,GT1b_Con1,RF1_2k1b,-,0.0\n"
)


def test_parse_openrdp_csv_reads_all_rows():
    rows = mc.parse_openrdp_csv(_CSV)
    assert len(rows) == 4
    assert rows[0]["Method"] == "Threeseq"
    assert rows[0]["Recombinant"] == "GT2a_JFH1"


def test_openrdp_detection_uses_full_triplet():
    det = mc.openrdp_detections(mc.parse_openrdp_csv(_CSV), "RF1_2k1b")
    # Query appears as a parent (Threeseq/Geneconv) or the recombinant (Bootscan/Maxchi).
    assert det["threeseq"] is True
    assert det["bootscan"] is True
    assert det["maxchi"] is True
    assert det["geneconv"] is True
    # Methods that produced no rows do not detect.
    assert det["siscan"] is False
    assert det["chimaera"] is False
    assert det["rdp"] is False


def test_openrdp_detection_absent_query_is_all_false():
    det = mc.openrdp_detections(mc.parse_openrdp_csv(_CSV), "NOT_PRESENT")
    assert not any(det.values())


def test_openrdp_command_prefers_env_override(monkeypatch):
    monkeypatch.setenv("OPENRDP_CMD", "conda run -n foo openrdp")
    assert mc.openrdp_command() == ["conda", "run", "-n", "foo", "openrdp"]
