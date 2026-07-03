"""Unit tests for the coalescent (Posada & Crandall) benchmark harness (no msprime)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_coalescent_benchmark.py"
_SPEC = importlib.util.spec_from_file_location("run_coalescent_benchmark", _PATH)
rcb = importlib.util.module_from_spec(_SPEC)
sys.modules["run_coalescent_benchmark"] = rcb
_SPEC.loader.exec_module(rcb)


def test_grid_includes_a_clonal_control():
    # The design needs recombination_rate == 0 as the specificity (false-positive) control.
    assert 0.0 in rcb.REC_RATES
    assert any(r > 0 for r in rcb.REC_RATES)
    assert rcb.MUT_RATES  # at least one diversity level


def test_reuses_shared_power_specificity_scorer():
    # It imports the same scorer that run_benchmark.py exposes (CI-tested there).
    agg = rcb.power_specificity([(0.0, 0.9), (1e-6, 0.01)], alpha=0.05)
    assert agg["specificity"] == 1.0 and agg["power"] == 1.0


def test_main_skips_without_msprime(monkeypatch):
    # With msprime unavailable, main() prints a SKIP and returns 0 (opt-in, like CI).
    import builtins
    real_import = builtins.__import__

    def no_msprime(name, *a, **k):
        if name == "msprime":
            raise ImportError("no msprime")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", no_msprime)
    assert rcb.main([]) == 0
