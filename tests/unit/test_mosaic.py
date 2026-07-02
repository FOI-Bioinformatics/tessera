"""make_mosaic span math and pick_parents_n selection (no aligner)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "run_hybrids", Path(__file__).resolve().parents[2] / "validation" / "run_hybrids.py")
rh = importlib.util.module_from_spec(_SPEC)
sys.modules["run_hybrids"] = rh
_SPEC.loader.exec_module(rh)


def test_make_mosaic_emits_non_backbone_spans_in_order():
    ref = "A" * 100
    # No mutations -> each segment reconstructs to the reference slice; boundaries at
    # 0/30/45/70/100. Non-backbone segments (clade set) become true_spans.
    segs = [([], 0.0, 0.30, None), ([], 0.30, 0.45, "B"),
            ([], 0.45, 0.70, None), ([], 0.70, 1.0, "C")]
    query, spans = rh.make_mosaic(ref, segs)
    assert len(query) == 100
    assert [c for _s, _e, c in spans] == ["B", "C"]
    assert spans[0] == (30, 45, "B") and spans[1] == (70, 100, "C")


def test_make_mosaic_single_insert_matches_make_hybrid_span():
    ref = "ACGT" * 25  # 100 bp
    segs = [([], 0.0, 0.35, None), ([], 0.35, 0.65, "B"), ([], 0.65, 1.0, None)]
    _q, spans = rh.make_mosaic(ref, segs)
    assert spans == [(35, 65, "B")]
