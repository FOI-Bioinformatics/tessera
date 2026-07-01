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
            "donor_absent"]
    lines = ["\t".join(cols)]
    for r in rows:
        lines.append("\t".join(str(r.get(c, "")) for c in cols))
    (tmp_path / "recombination_regions.tsv").write_text("\n".join(lines) + "\n")
    (tmp_path / "recombination_profile.tsv").write_text("phi_p\trmin\n-\t0\n")


def test_single_insert_pass(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "hmm,3seq", "donor_absent": "no"},
    ])
    clade_of = {"gA": "A", "gB": "B"}.get
    setup = _setup(out=tmp_path)
    res = rh._score_regions(tmp_path, clade_of, setup, 5, "tip", 1.0)
    assert res["pass"] is True and res["backbone_ok"] and res["donor_ok"]
