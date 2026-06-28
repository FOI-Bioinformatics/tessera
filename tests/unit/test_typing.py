"""Typed genotype/lineage names: title mining, source merge, sidecar IO, rendering."""

from __future__ import annotations

from pathlib import Path

from recomfi.recomb.typing import (
    build_lineage_map,
    dominant_lineage_token,
    genotype_from_title,
    lineage_map_from_rows,
    lineage_of,
    load_lineage_map,
    read_lineage_rows,
    titles_from_collection,
    typed,
    write_lineage_map,
)

# --- genotype_from_title ---------------------------------------------------

def test_genotype_from_norovirus_title():
    title = "MK573073.1 Norovirus GII.P16-GII.4 isolate Hu/2015, complete genome"
    assert genotype_from_title(title, "Norovirus GII") == "GII.P16-GII.4"


def test_genotype_from_hiv_title_prefers_crf_over_organism_abbrev():
    title = "U54771.1 Human immunodeficiency virus 1 (HIV-1) CRF01_AE, complete genome"
    assert genotype_from_title(title, "Human immunodeficiency virus 1") == "CRF01_AE"


def test_genotype_from_pango_note():
    # The NCBI datasets lineage note appended after the accession.
    assert genotype_from_title("OM958567.1 B.1.1.7") == "B.1.1.7"


def test_genotype_none_when_only_organism():
    title = "X12345.1 Human immunodeficiency virus 1, complete genome"
    assert genotype_from_title(title, "Human immunodeficiency virus 1") is None


def test_genotype_ignores_accession_token():
    # The leading accession (carries a digit) must not be returned as a genotype.
    assert genotype_from_title("MK573073.1, complete genome") is None


def test_dominant_lineage_token_across_titles():
    titles = [
        "Norovirus GII.P16-GII.1 isolate A",
        "Norovirus GII.P16-GII.1 isolate B",
        "Norovirus GII.4 isolate C",
    ]
    assert dominant_lineage_token(titles) == "GII.P16-GII.1"


# --- build_lineage_map: priority + normalization ---------------------------

def test_build_lineage_map_priority_user_overrides_title(tmp_path: Path):
    user = tmp_path / "user.tsv"
    write_lineage_map(user, [("MK573073", "GII.P16-GII.4 (curated)", "user")])
    rows = build_lineage_map(
        user_tsv=user,
        title_by_label={"MK573073.1": "MK573073.1 Norovirus GII.P16-GII.4 isolate Hu"},
        organism="Norovirus GII",
    )
    out = dict((label, g) for label, g, _ in rows)
    # Keys are normalized and the user value wins over the mined token.
    assert out["MK573073"] == "GII.P16-GII.4 (curated)"


def test_build_lineage_map_datasets_over_title():
    rows = build_lineage_map(
        datasets_rows=[("OM958567", "B.1.1.7")],
        title_by_label={"OM958567": "OM958567.1 Severe acute respiratory syndrome coronavirus 2"},
    )
    assert dict((label, g) for label, g, _ in rows)["OM958567"] == "B.1.1.7"


def test_build_lineage_map_skips_untyped_titles():
    rows = build_lineage_map(title_by_label={"X1": "X1 Norovirus, complete genome"})
    assert rows == []


# --- sidecar IO round-trip -------------------------------------------------

def test_lineage_tsv_round_trip(tmp_path: Path):
    path = tmp_path / "lineages.tsv"
    write_lineage_map(path, [("acc1", "GII.4", "title"), ("acc2", "CRF01_AE", "user")])
    rows = read_lineage_rows(path)
    assert rows == [("acc1", "GII.4", "title"), ("acc2", "CRF01_AE", "user")]
    assert load_lineage_map(path) == {"acc1": "GII.4", "acc2": "CRF01_AE"}


def test_read_lineage_rows_skips_garbled(tmp_path: Path):
    path = tmp_path / "lineages.tsv"
    path.write_text("label\tgenotype\tsource\n\nacc1\tGII.4\ttitle\nbroken-line\n\t\t\n")
    assert read_lineage_rows(path) == [("acc1", "GII.4", "title")]


def test_load_lineage_map_missing_file_is_empty(tmp_path: Path):
    assert load_lineage_map(tmp_path / "nope.tsv") == {}
    assert load_lineage_map(None) == {}


# --- titles_from_collection ------------------------------------------------

def test_titles_from_collection_reads_raw_defline(tmp_path: Path):
    (tmp_path / "MK573073.fasta").write_text(">MK573073.1 Norovirus GII.P16-GII.4\nACGT\n")
    titles = titles_from_collection([tmp_path / "MK573073.fasta"])
    assert titles == {"MK573073": ">MK573073.1 Norovirus GII.P16-GII.4"[1:]}


# --- typed() / lineage_of() ------------------------------------------------

def test_typed_mapped_and_unmapped():
    lm = {"MK573073": "GII.P16-GII.4"}
    assert typed("MK573073", lm) == "GII.P16-GII.4 (MK573073)"
    assert typed("OTHER", lm) == "OTHER"  # unmapped -> bare label


def test_typed_no_map_is_bare_label():
    assert typed("MK573073", None) == "MK573073"
    assert typed("MK573073", {}) == "MK573073"


def test_typed_version_insensitive_fallback():
    # A versioned label resolves against an unversioned map key (and vice versa).
    lm = {"MK573073": "GII.4"}
    assert typed("MK573073.1", lm) == "GII.4 (MK573073.1)"
    assert lineage_of("MK573073.2", lm) == "GII.4"


def test_lineage_map_from_rows_strips_file_ext_keeps_version():
    # Keys lose a sequence extension but keep the version; version skew is handled at
    # lookup time (test_typed_version_insensitive_fallback).
    rows = [("MK573073.1.fasta", "GII.4", "title")]
    assert lineage_map_from_rows(rows) == {"MK573073.1": "GII.4"}
