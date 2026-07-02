"""Unit test for the reassort-scan probe's localization scorer (no binaries)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_reassort_scan.py"
_SPEC = importlib.util.spec_from_file_location("run_reassort_scan", _PATH)
rrs = importlib.util.module_from_spec(_SPEC)
sys.modules["run_reassort_scan"] = rrs
_SPEC.loader.exec_module(rrs)


def _rows(*triples):
    return [{"query_start": str(a), "query_end": str(b), "donor_absent": d}
            for a, b, d in triples]


def test_overlapping_region_is_localized():
    assert rrs.region_overlaps_span(_rows((600, 1200, "no")), 601, 1122) is True


def test_disjoint_region_is_not_localized():
    assert rrs.region_overlaps_span(_rows((0, 100, "no")), 601, 1122) is False


def test_touching_edge_counts_as_overlap():
    # region starts exactly at the span end -> touching interval counts.
    assert rrs.region_overlaps_span(_rows((1122, 1300, "no")), 601, 1122) is True


def test_donor_absent_region_is_excluded():
    assert rrs.region_overlaps_span(_rows((600, 1200, "yes")), 601, 1122) is False


def test_no_regions_is_not_localized():
    assert rrs.region_overlaps_span([], 601, 1122) is False
