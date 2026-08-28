"""run_recomb helpers: lineage-map auto-discovery beside the output or the MSA."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from tessera.core.errors import UserInputError
from tessera.recomb.run import RecombParams, _discover_lineage_tsv, run_recomb

from ..conftest import write_fasta


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


@pytest.mark.parametrize("refs", [{}, {"refA": "ACGTACGTAC"}])
def test_undersized_panel_refuses_rather_than_reporting_a_negative(
    tmp_path: Path, refs: dict[str, str]
) -> None:
    """A panel with fewer than two references cannot express a donor switch.

    Every caller early-returns below two labels, so the scan would otherwise finish
    and report no regions -- a result a reader cannot tell apart from a genuine
    clean negative.
    """
    records = {"query": "ACGTACGTAC", **refs}
    msa = write_fasta(tmp_path / "panel.msa.fasta", records)
    params = RecombParams(
        msa=msa, query="query", output=tmp_path / "out", window_size=5, window_step=5
    )
    with pytest.raises(UserInputError, match="at least 2"):
        run_recomb(params, logging.getLogger("tessera"))
