"""Genome staging and reference selection."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import pytest

from recomfi.core.errors import UserInputError
from recomfi.core.io import select_reference, stage_genomes, strip_sequence_extension

_LOG = logging.getLogger("recomfi.test")


def test_strip_sequence_extension() -> None:
    assert strip_sequence_extension("cowpox.fasta.gz") == "cowpox"
    assert strip_sequence_extension("variola.fa") == "variola"
    assert strip_sequence_extension("sample.fna") == "sample"
    assert strip_sequence_extension("noext") == "noext"


def test_stage_genomes_handles_gz_and_plain(tmp_path: Path) -> None:
    collection = tmp_path / "collection"
    collection.mkdir()
    # one gzipped, one plain reference
    with gzip.open(collection / "refA.fasta.gz", "wt") as fo:
        fo.write(">refA\nACGTACGT\n")
    (collection / "refB.fasta").write_text(">refB\nACGTACGT\n")
    query = tmp_path / "query.fasta.gz"
    with gzip.open(query, "wt") as fo:
        fo.write(">query\nACGTACGT\n")

    target = tmp_path / "genomes"
    staged, query_staged = stage_genomes(query, collection, target, _LOG)

    names = sorted(p.name for p in staged)
    assert names == ["query.fasta", "refA.fasta", "refB.fasta"]
    # decompressed content is readable plain text
    assert (target / "refA.fasta").read_text().startswith(">refA")
    assert query_staged.name == "query.fasta"


def test_stage_genomes_missing_inputs(tmp_path: Path) -> None:
    with pytest.raises(UserInputError):
        stage_genomes(tmp_path / "nope.fasta", tmp_path, tmp_path / "g", _LOG)


def test_select_reference_default_is_first_non_query(tmp_path: Path) -> None:
    genomes = [tmp_path / "a.fasta", tmp_path / "q.fasta"]
    query = tmp_path / "q.fasta"
    ref = select_reference(genomes, query, query_as_backbone=False)
    assert ref.name == "a.fasta"


def test_select_reference_query_as_backbone(tmp_path: Path) -> None:
    genomes = [tmp_path / "a.fasta", tmp_path / "q.fasta"]
    query = tmp_path / "q.fasta"
    ref = select_reference(genomes, query, query_as_backbone=True)
    assert ref == query


def test_select_reference_explicit_match_and_miss(tmp_path: Path) -> None:
    genomes = [tmp_path / "a.fasta", tmp_path / "b.fasta", tmp_path / "q.fasta"]
    query = tmp_path / "q.fasta"
    assert select_reference(genomes, query, False, "b").name == "b.fasta"
    assert select_reference(genomes, query, False, "b.fasta").name == "b.fasta"
    with pytest.raises(UserInputError):
        select_reference(genomes, query, False, "missing")
