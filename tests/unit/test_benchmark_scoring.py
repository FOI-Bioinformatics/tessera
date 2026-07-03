"""Unit tests for the PHI benchmark harness's pure logic (no binaries)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_benchmark.py"
_SPEC = importlib.util.spec_from_file_location("run_benchmark", _PATH)
rb = importlib.util.module_from_spec(_SPEC)
sys.modules["run_benchmark"] = rb
_SPEC.loader.exec_module(rb)


def test_parse_msa_name_extracts_parameters():
    p = rb.parse_msa_name("msa_m0.001_rc0.05_n100_dual1_rep2.fasta")
    assert p == {"mut": 0.001, "rec": 0.05, "n": 100, "dual": 1.0, "rep": 2}


def test_parse_msa_name_zero_recombination():
    assert rb.parse_msa_name("msa_m0.03_rc0_n24_dual0.1_rep1.fasta")["rec"] == 0.0


def test_parse_msa_name_rejects_unrelated():
    assert rb.parse_msa_name("README.md") is None
    assert rb.parse_msa_name("something.fasta") is None


def test_power_specificity_counts_correctly():
    # rec>0: two of three testable are significant -> power 2/3; rec=0: one testable,
    # non-significant -> specificity 1/1; one untestable (None) excluded from the rates.
    results = [
        (0.05, 0.01), (0.05, 0.2), (0.05, 0.001),
        (0.0, 0.9), (0.0, None),
    ]
    agg = rb.power_specificity(results, alpha=0.05)
    assert agg["power"] == 2 / 3
    assert agg["specificity"] == 1.0
    assert agg["n_recombining"] == 3
    assert agg["n_clonal"] == 2
    assert agg["n_untestable"] == 1


def test_power_specificity_handles_empty_groups():
    agg = rb.power_specificity([(0.0, 0.9)], alpha=0.05)
    assert agg["power"] is None  # no recombining alignments
    assert agg["specificity"] == 1.0
