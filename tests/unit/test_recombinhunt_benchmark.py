"""Unit tests for the RecombinHunt-style noise harness's pure logic (no binaries)."""

from __future__ import annotations

import importlib.util
import random
import sys
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_recombinhunt_benchmark.py"
_SPEC = importlib.util.spec_from_file_location("run_recombinhunt_benchmark", _PATH)
rh = importlib.util.module_from_spec(_SPEC)
sys.modules["run_recombinhunt_benchmark"] = rh
_SPEC.loader.exec_module(rh)


def test_inject_noise_changes_exactly_k_positions():
    seq = "ACGT" * 25  # 100 bp
    out = rh.inject_noise(seq, 10, random.Random(0))
    diffs = sum(a != b for a, b in zip(seq, out, strict=True))
    assert diffs == 10
    assert len(out) == len(seq)


def test_inject_noise_zero_is_identity():
    seq = "ACGTACGT"
    assert rh.inject_noise(seq, 0, random.Random(0)) == seq


def test_inject_noise_substitutes_a_different_base():
    seq = "AAAA"
    out = rh.inject_noise(seq, 4, random.Random(1))
    assert all(c != "A" for c in out)  # every position changed to a different base


def test_sensitivity_by_noise_aggregates_per_level():
    results = [(0, True), (0, True), (5, True), (5, False), (30, False)]
    agg = rh.sensitivity_by_noise(results)
    assert agg == {0: (2, 2), 5: (1, 2), 30: (0, 1)}


def test_noise_grid_starts_at_zero_and_increases():
    assert rh.NOISE_GRID[0] == 0
    assert list(rh.NOISE_GRID) == sorted(rh.NOISE_GRID)
