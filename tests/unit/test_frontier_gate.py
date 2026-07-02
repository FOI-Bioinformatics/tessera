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
