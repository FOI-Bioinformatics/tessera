"""Unit tests for the reassortment benchmark's pure scoring (no binaries)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_reassort_benchmark.py"
_SPEC = importlib.util.spec_from_file_location("run_reassort_benchmark", _PATH)
rb = importlib.util.module_from_spec(_SPEC)
sys.modules["run_reassort_benchmark"] = rb
_SPEC.loader.exec_module(rb)


def test_prf1_basic():
    p, r, f = rb.prf1(tp=8, fp=2, fn=2)
    assert p == 0.8 and r == 0.8 and abs(f - 0.8) < 1e-9


def test_prf1_zero_when_undefined():
    assert rb.prf1(0, 0, 0) == (0.0, 0.0, 0.0)


def test_score_confusion_and_f1():
    labeled = [
        (True, "reassortant"), (True, "reassortant"),   # 2 TP
        (True, "clonal"),                                # 1 FN
        (True, "undetermined"),                          # FN (miss on a reassortant)
        (False, "clonal"), (False, "clonal"),            # 2 TN
        (False, "reassortant"),                          # 1 FP
        (False, "undetermined"),                         # undetermined on clonal, not a FP
    ]
    res = rb.score(labeled)
    assert res["tp"] == 2 and res["fn"] == 2 and res["fp"] == 1 and res["tn"] == 2
    assert res["undetermined_clonal"] == 1
    assert res["precision"] == 2 / 3
    assert res["recall"] == 2 / 4
    assert abs(res["f1"] - (2 * (2 / 3) * 0.5) / ((2 / 3) + 0.5)) < 1e-9


def test_score_perfect():
    labeled = [(True, "reassortant"), (False, "clonal")]
    res = rb.score(labeled)
    assert res["precision"] == 1.0 and res["recall"] == 1.0 and res["f1"] == 1.0
