"""Provenance: the aligner sidecar written beside an MSA, and its use by the scan.

An alignment is an interpretation of the input, not a fact about it, so a result
that cannot name the aligner behind it cannot be reproduced by a reader.
"""

from __future__ import annotations

import json
from pathlib import Path

from tessera.msa.build import provenance_path, read_provenance
from tessera.recomb.run import _msa_provenance_lines, _write_run_provenance


def test_provenance_path_sits_beside_the_msa(tmp_path: Path) -> None:
    assert provenance_path(tmp_path / "panel.msa.fasta").name == "panel.msa.provenance.json"


def test_read_provenance_absent_is_none(tmp_path: Path) -> None:
    # A user-supplied alignment has no sidecar. That is a normal input, not an error.
    assert read_provenance(tmp_path / "hand_made.fasta") is None


def test_read_provenance_unreadable_is_none(tmp_path: Path) -> None:
    msa = tmp_path / "panel.msa.fasta"
    msa.write_text(">q\nACGT\n")
    provenance_path(msa).write_text("{ truncated")
    assert read_provenance(msa) is None


def test_read_provenance_roundtrip(tmp_path: Path) -> None:
    msa = tmp_path / "panel.msa.fasta"
    provenance_path(msa).write_text(json.dumps({"aligner": "mafft", "backbone": "ref"}))
    assert read_provenance(msa)["aligner"] == "mafft"


def test_msa_provenance_lines_name_the_aligner_and_version() -> None:
    lines = _msa_provenance_lines({
        "aligner": "mafft",
        "aligner versions": {"mafft": "7.526.0"},
        "aligner args": {"maxiterate": "2"},
        "backbone": "outgroup",
    })
    assert lines["aligner"] == "mafft (mafft 7.526.0)"
    assert lines["aligner args"] == "maxiterate=2"
    assert lines["backbone"] == "outgroup"


def test_msa_provenance_lines_claim_nothing_without_a_sidecar() -> None:
    # Better to say nothing than to report an aligner the run cannot actually know.
    assert _msa_provenance_lines(None) == {}
    assert _msa_provenance_lines({}) == {}


def test_msa_provenance_lines_tolerate_a_sparse_record() -> None:
    assert _msa_provenance_lines({"aligner": "minimap2"}) == {"aligner": "minimap2"}


def test_run_provenance_json_carries_run_and_msa(tmp_path: Path, logger) -> None:
    out = tmp_path / "out"
    written = _write_run_provenance(
        out, {"tessera version": "1.0.0", "query": "q"}, {"aligner": "mafft"}, logger
    )
    assert written == out / "run_provenance.json"
    record = json.loads(written.read_text())
    assert record["run"]["query"] == "q"
    assert record["msa"]["aligner"] == "mafft"


def test_run_provenance_json_records_a_missing_sidecar_as_null(tmp_path: Path, logger) -> None:
    written = _write_run_provenance(tmp_path / "out", {"query": "q"}, None, logger)
    assert json.loads(written.read_text())["msa"] is None
