"""Unit tests for the method-comparison harness's pure detection helper (no binaries)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_method_comparison.py"
_SPEC = importlib.util.spec_from_file_location("run_method_comparison", _PATH)
mc = importlib.util.module_from_spec(_SPEC)
sys.modules["run_method_comparison"] = mc
_SPEC.loader.exec_module(mc)


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
