"""``msa`` warns when the collection is too small for ``recomb`` to use."""

from __future__ import annotations

import logging
from pathlib import Path

from tessera.aligners.base import AlignResult
from tessera.msa import build as build_module
from tessera.msa.build import MsaParams, build_msa


class _StubAligner:
    def preflight(self) -> dict[str, str]:
        return {"stub": "0"}

    def align(self, genomes, reference, out_dir, params, logger) -> AlignResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        msa = out_dir / "msa.fasta"
        msa.write_text("".join(f">{g.stem}\nACGT\n" for g in genomes))
        return AlignResult(msa_fasta=msa)


def _build(tmp_path: Path, n_references: int, monkeypatch, caplog) -> None:
    monkeypatch.setattr(build_module.aligner_registry, "create", lambda name: _StubAligner())
    collection = tmp_path / "collection"
    collection.mkdir()
    for i in range(n_references):
        (collection / f"ref{i}.fasta").write_text(f">r{i}\nACGT\n")
    query = tmp_path / "q.fasta"
    query.write_text(">q\nACGT\n")
    with caplog.at_level(logging.WARNING, logger="tessera.test"):
        build_msa(
            MsaParams(query=query, collection=collection, output=tmp_path / "m.fasta"),
            logging.getLogger("tessera.test"),
        )


def test_a_single_reference_collection_is_warned_about(tmp_path: Path, monkeypatch, caplog) -> None:
    _build(tmp_path, 1, monkeypatch, caplog)
    assert (tmp_path / "m.fasta").exists()  # still built: the alignment itself is valid
    text = " ".join(caplog.messages)
    assert "holds 1 reference(s) besides the query" in text
    assert "needs at least 2" in text


def test_two_references_are_not_warned_about(tmp_path: Path, monkeypatch, caplog) -> None:
    _build(tmp_path, 2, monkeypatch, caplog)
    assert "besides the query" not in " ".join(caplog.messages)
