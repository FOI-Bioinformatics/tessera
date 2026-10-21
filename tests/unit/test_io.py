"""Genome staging and reference selection."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import pytest

from tessera.core.errors import UserInputError
from tessera.core.io import select_reference, stage_genomes, strip_sequence_extension

_LOG = logging.getLogger("tessera.test")


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


def test_colliding_labels_are_rejected(tmp_path: Path) -> None:
    """Two files reducing to one label cannot be told apart in the alignment.

    Staging would otherwise fail on the symlink with a bare FileExistsError, which
    the CLI reports as "Unexpected error".
    """
    coll = tmp_path / "coll"
    coll.mkdir()
    (coll / "refA.fa").write_text(">a\nACGT\n")
    (coll / "refA.fasta").write_text(">a\nACGT\n")
    query = tmp_path / "q.fasta"
    query.write_text(">q\nACGT\n")

    with pytest.raises(UserInputError, match="share the label 'refA'"):
        stage_genomes(query, coll, tmp_path / "staged", _LOG)


def test_query_colliding_with_a_collection_member_is_rejected(tmp_path: Path) -> None:
    coll = tmp_path / "coll"
    coll.mkdir()
    (coll / "sample.fa").write_text(">a\nACGT\n")
    query = tmp_path / "sample.fasta"
    query.write_text(">q\nACGT\n")

    with pytest.raises(UserInputError, match="share the label 'sample'"):
        stage_genomes(query, coll, tmp_path / "staged", _LOG)


def test_stage_genomes_cleans_whitespace_out_of_sequence_lines(tmp_path: Path) -> None:
    """The aligner reads the staged file; Tessera reads the same genome through
    read_fasta. If only one of them ignores a trailing space the two disagree about the
    genome's length -- minimap2 counted the spaces and every later column was shifted."""
    coll = tmp_path / "coll"
    coll.mkdir()
    (coll / "clean.fasta").write_text(">c desc\nACGT\nTTAA\n")
    (coll / "spaces.fasta").write_text(">s desc\nACGT \nTT AA\t\n")
    (coll / "crlf.fasta").write_bytes(b">w desc\r\nACGT\r\nTTAA\r\n")
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGT\n")

    staged, _ = stage_genomes(query, coll, tmp_path / "stage", _LOG)
    by_stem = {p.stem: p for p in staged}

    assert by_stem["clean"].is_symlink()  # untouched input is still linked, not copied
    for stem, header in (("spaces", ">s desc"), ("crlf", ">w desc")):
        assert not by_stem[stem].is_symlink()
        assert by_stem[stem].read_text() == f"{header}\nACGT\nTTAA\n"
    # The user's own files are not modified.
    assert (coll / "spaces.fasta").read_text() == ">s desc\nACGT \nTT AA\t\n"


def test_stage_genomes_cleans_a_gzipped_genome(tmp_path: Path) -> None:
    coll = tmp_path / "coll"
    coll.mkdir()
    with gzip.open(coll / "z.fasta.gz", "wt") as fo:
        fo.write(">z\nACGT \n\nTTAA\n")
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGT\n")
    staged, _ = stage_genomes(query, coll, tmp_path / "stage", _LOG)
    assert {p.stem: p for p in staged}["z"].read_text() == ">z\nACGT\nTTAA\n"
