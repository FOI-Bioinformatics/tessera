"""run_recomb helpers: lineage-map auto-discovery beside the output or the MSA."""

from __future__ import annotations

from pathlib import Path

from recomfi.recomb.run import _discover_lineage_tsv


def test_discover_lineage_tsv_in_output_dir(tmp_path: Path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "lineages.tsv").write_text("label\tgenotype\tsource\n")
    assert _discover_lineage_tsv(out, tmp_path / "panel.msa.fasta") == out / "lineages.tsv"


def test_discover_lineage_tsv_beside_msa(tmp_path: Path):
    msa = tmp_path / "panel.msa.fasta"
    msa.write_text(">q\nACGT\n")
    (tmp_path / "lineages.tsv").write_text("label\tgenotype\tsource\n")
    assert _discover_lineage_tsv(tmp_path / "elsewhere", msa) == tmp_path / "lineages.tsv"


def test_discover_lineage_tsv_none_when_absent(tmp_path: Path):
    assert _discover_lineage_tsv(tmp_path / "out", tmp_path / "panel.msa.fasta") is None
