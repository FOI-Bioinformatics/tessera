"""Unit tests for run_hybrids scoring branches (no aligner; synthetic regions TSV)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "run_hybrids",
    Path(__file__).resolve().parents[2] / "validation" / "run_hybrids.py",
)
rh = importlib.util.module_from_spec(_SPEC)
sys.modules["run_hybrids"] = rh  # required for @dataclass to resolve its own module
_SPEC.loader.exec_module(rh)


def _setup(**kw):
    d = dict(name="t", out=Path("/tmp"), clade_a="A", clade_b="B", divergence=10.0,
             query=Path("/tmp/q"), query_label="q", q_start=100, q_end=200,
             window=300, step=30, sel_window=1500, aligner="mafft", reference="",
             tips={}, pool=[], members_by_clade={}, case_type="single_insert")
    d.update(kw)
    return rh.CaseSetup(**d)


def _write_regions(tmp_path, rows):
    """rows: list of dicts with query_start/query_end/major_parent/minor_parent/methods/
    donor_absent. Writes recombination_regions.tsv + an empty profile."""
    cols = ["minor_parent", "major_parent", "query_start", "query_end", "methods",
            "donor_absent", "donor_undercovered"]
    lines = ["\t".join(cols)]
    for r in rows:
        lines.append("\t".join(str(r.get(c, "")) for c in cols))
    (tmp_path / "recombination_regions.tsv").write_text("\n".join(lines) + "\n")
    (tmp_path / "recombination_profile.tsv").write_text("phi_p\trmin\n-\t0\n")


def test_neg_pure_pass_when_no_regions(tmp_path):
    _write_regions(tmp_path, [])  # no calls
    setup = _setup(out=tmp_path, case_type="neg_pure", clade_b="", q_start=0, q_end=0)
    res = rh._score_regions(tmp_path, (lambda x: "A"), setup, 5, "tip", 1.0)
    assert res["pass"] is True and res["n_false_regions"] == 0


def test_neg_pure_fail_on_false_call(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 10,
         "query_end": 90, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="neg_pure", clade_b="", q_start=0, q_end=0)
    res = rh._score_regions(tmp_path, {"gA": "A", "gB": "B"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is False and res["n_false_regions"] == 1


def test_neg_within_passes_same_clade_call(tmp_path):
    # A region whose minor and major are the SAME top-level clade is not a cross-clade FP.
    _write_regions(tmp_path, [
        {"minor_parent": "gA2", "major_parent": "gA1", "query_start": 10,
         "query_end": 90, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="neg_within", clade_a="A", clade_b="A")
    res = rh._score_regions(tmp_path, {"gA1": "A.1", "gA2": "A.2"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is True


def test_neg_within_fails_cross_clade_call(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA1", "query_start": 10,
         "query_end": 90, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="neg_within", clade_a="A", clade_b="A")
    res = rh._score_regions(tmp_path, {"gA1": "A.1", "gB": "B"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is False


def test_single_insert_pass(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "hmm,3seq", "donor_absent": "no"},
    ])
    clade_of = {"gA": "A", "gB": "B"}.get
    setup = _setup(out=tmp_path)
    res = rh._score_regions(tmp_path, clade_of, setup, 5, "tip", 1.0)
    assert res["pass"] is True and res["backbone_ok"] and res["donor_ok"]


def test_low_div_requires_backbone_top_level(tmp_path):
    # donor correct top-level (B), backbone WRONG top-level (C) -> FAIL even at low div.
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gC", "query_start": 120,
         "query_end": 180, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="low_div", divergence=2.0,
                   clade_a="A", clade_b="B")
    res = rh._score_regions(tmp_path, {"gB": "B.1", "gC": "C.1"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is False  # backbone C != A, not relaxed away


def test_low_div_pass_top_level(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="low_div", divergence=2.0,
                   clade_a="A", clade_b="B")
    res = rh._score_regions(tmp_path, {"gA": "A.2", "gB": "B.1"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is True


def test_donor_absent_pass_on_gap(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "coverage", "donor_absent": "yes"}])
    setup = _setup(out=tmp_path, case_type="panel_donor_absent", clade_a="A", clade_b="B")
    res = rh._score_regions(tmp_path, {"gA": "A"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is True


def test_donor_absent_fail_on_misattribution(tmp_path):
    # No donor_absent region; a confident cross-clade present region mis-attributes the span.
    _write_regions(tmp_path, [
        {"minor_parent": "gC", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "hmm,3seq", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="panel_donor_absent", clade_a="A", clade_b="B")
    res = rh._score_regions(tmp_path, {"gA": "A", "gC": "C"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is False


def test_donor_absent_pass_on_undercovered_caveat(tmp_path):
    # No donor_absent region, but the span-overlapping region is caveated
    # donor_undercovered -> the tool acknowledged the missing donor -> PASS.
    _write_regions(tmp_path, [
        {"minor_parent": "gC", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "hmm,3seq", "donor_absent": "no",
         "donor_undercovered": "yes"}])
    setup = _setup(out=tmp_path, case_type="panel_donor_absent", clade_a="A", clade_b="B")
    res = rh._score_regions(tmp_path, {"gA": "A", "gC": "C"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is True


def test_donor_absent_fail_on_clean_confident_misattribution(tmp_path):
    # Confident, cross-clade, NOT undercovered -> a silent mis-attribution -> FAIL.
    _write_regions(tmp_path, [
        {"minor_parent": "gC", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "hmm,3seq", "donor_absent": "no",
         "donor_undercovered": "no"}])
    setup = _setup(out=tmp_path, case_type="panel_donor_absent", clade_a="A", clade_b="B")
    res = rh._score_regions(tmp_path, {"gA": "A", "gC": "C"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is False


def test_equidistant_pass_when_B_wins(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="panel_equidistant", clade_a="A", clade_b="B")
    res = rh._score_regions(tmp_path, {"gA": "A", "gB": "B", "gC": "C"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is True


def test_mosaic_pass_all_spans_recovered(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 300,
         "query_end": 450, "methods": "hmm", "donor_absent": "no"},
        {"minor_parent": "gC", "major_parent": "gA", "query_start": 700,
         "query_end": 1000, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="mosaic", clade_a="A",
                   true_spans=[(300, 450, "B"), (700, 1000, "C")], pattern="ABAC")
    res = rh._score_regions(tmp_path, {"gA": "A", "gB": "B", "gC": "C"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is True and res["spans_hit"] == 2


def test_mosaic_fail_when_a_span_missed(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 300,
         "query_end": 450, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="mosaic", clade_a="A",
                   true_spans=[(300, 450, "B"), (700, 1000, "C")], pattern="ABAC")
    res = rh._score_regions(tmp_path, {"gA": "A", "gB": "B", "gC": "C"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is False and res["spans_hit"] == 1


def test_mosaic_short_is_detection_gated(tmp_path):
    _write_regions(tmp_path, [  # a region exists but does not overlap the true span
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 10,
         "query_end": 40, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="mosaic", clade_a="A",
                   true_spans=[(500, 560, "B")], pattern="AB_short")
    res = rh._score_regions(tmp_path, {"gA": "A", "gB": "B"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is True  # detection-gated: a call exists, span recovery only reported
