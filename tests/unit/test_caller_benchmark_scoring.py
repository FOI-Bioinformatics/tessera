"""Unit tests for the region-caller benchmark's pure logic (no data, no binaries)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_caller_benchmark.py"
_SPEC = importlib.util.spec_from_file_location("run_caller_benchmark", _PATH)
cb = importlib.util.module_from_spec(_SPEC)
sys.modules["run_caller_benchmark"] = cb
_SPEC.loader.exec_module(cb)


def test_pick_queries_is_deterministic_and_bounded():
    labels = [f"seq{i}" for i in range(100)]
    a = cb.pick_queries(labels, 10, seed=0)
    assert a == cb.pick_queries(labels, 10, seed=0)
    assert len(a) == 10
    assert set(a) <= set(labels)
    assert cb.pick_queries(labels, 10, seed=1) != a


def test_pick_queries_when_fewer_labels_than_asked():
    labels = ["a", "b", "c"]
    assert sorted(cb.pick_queries(labels, 10, seed=0)) == ["a", "b", "c"]


def test_pick_queries_empty():
    assert cb.pick_queries([], 5, seed=0) == []


def test_alignment_is_detected_when_any_query_calls_a_region():
    assert cb.alignment_detected([0, 0, 2, 0]) is True
    assert cb.alignment_detected([0, 0, 0]) is False
    assert cb.alignment_detected([]) is False


def test_reuses_the_phi_benchmark_scorer_so_the_two_are_comparable():
    """Both benchmarks must aggregate identically or their numbers cannot be compared."""
    results = [(0.0, 1.0), (0.0, 1.0), (0.05, 0.0), (0.05, 1.0)]
    agg = cb.power_specificity(results, 0.5)
    assert agg["n_clonal"] == 2
    assert agg["n_recombining"] == 2
    assert agg["power"] == 0.5          # one of two recombining scored below alpha
    assert agg["specificity"] == 1.0    # both clonal scored at/above alpha
