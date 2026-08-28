"""Native sliding-window similarity engine and query-coordinate mapping."""

from __future__ import annotations

from math import isnan
from pathlib import Path

import pytest

from tessera.core.errors import UserInputError
from tessera.recomb.similarity import (
    _informative_column_mask,
    _read_alignment,
    compute_similarity,
    compute_similarity_informative,
    informative_site_count,
    low_canonical_records,
)

from ..conftest import write_fasta


def test_identity_and_mismatch_windows(tmp_path: Path) -> None:
    msa = write_fasta(
        tmp_path / "m.fasta",
        {
            "query": "AAAAAAAAAAGGGGGGGGGG",
            "refA": "AAAAAAAAAAGGGGGGGGGG",  # identical
            "refB": "AAAAAAAAAATTTTTTTTTT",  # second half differs
        },
    )
    result = compute_similarity(str(msa), "query", window_size=10, window_step=10)
    assert result.positions == [5, 15]
    assert result.similarities["refA"] == [1.0, 1.0]
    assert result.similarities["refB"] == [1.0, 0.0]


def test_query_excluded_from_datasets(tmp_path: Path) -> None:
    msa = write_fasta(tmp_path / "m.fasta", {"query": "ACGTACGTAC", "refA": "ACGTACGTAC"})
    result = compute_similarity(str(msa), "query", window_size=5, window_step=5)
    assert "query" not in result.similarities
    assert set(result.similarities) == {"refA"}


def test_gap_columns_excluded(tmp_path: Path) -> None:
    # The gapped column (position 2) is ignored; the remaining columns match.
    msa = write_fasta(
        tmp_path / "m.fasta",
        {"query": "AC-GT", "refA": "AC-GT", "refB": "A--GT"},
    )
    result = compute_similarity(str(msa), "query", window_size=5, window_step=5)
    assert result.similarities["refA"] == [1.0]
    assert result.similarities["refB"] == [1.0]


def test_query_coordinate_mapping_with_gaps(tmp_path: Path) -> None:
    # The query has gaps in columns 2-3, so MSA columns and query bases diverge.
    #   col:   0 1 2 3 4 5 6 7
    #   query: A C - - G T A C   -> query bases before col 6 == 4 (A,C,G,T)
    msa = write_fasta(
        tmp_path / "m.fasta",
        {"query": "AC--GTAC", "refA": "ACGTGTAC"},
    )
    result = compute_similarity(str(msa), "query", window_size=4, window_step=2)
    # column_to_query counts non-gap query bases before a column
    assert result.column_to_query(0) == 0
    assert result.column_to_query(4) == 2   # only A,C precede col 4
    assert result.column_to_query(8) == 6   # all six query bases
    # query_positions parallels the window centers
    assert len(result.query_positions) == len(result.positions)


def test_ambiguous_bases_not_comparable(tmp_path: Path) -> None:
    # N must not count as identity. Columns 0,1 are canonical and match; columns
    # 2,3 are N in the query (excluded), so similarity is 1.0 (2/2), not 0.5.
    msa = write_fasta(tmp_path / "m.fasta", {"query": "ACNN", "refA": "ACAA"})
    result = compute_similarity(str(msa), "query", window_size=4, window_step=4)
    assert result.similarities["refA"] == [1.0]


def test_zero_comparable_window_is_nan(tmp_path: Path) -> None:
    # A window where the query has no canonical base is uninformative -> nan.
    msa = write_fasta(tmp_path / "m.fasta", {"query": "NNNN", "refA": "ACGT"})
    result = compute_similarity(str(msa), "query", window_size=4, window_step=4)
    assert isnan(result.similarities["refA"][0])


def test_unknown_query_raises(tmp_path: Path) -> None:
    msa = write_fasta(tmp_path / "m.fasta", {"a": "ACGTACGTAC", "b": "ACGTACGTAC"})
    with pytest.raises(UserInputError):
        compute_similarity(str(msa), "missing", window_size=5, window_step=5)


def test_unknown_metric_raises(tmp_path: Path) -> None:
    msa = write_fasta(tmp_path / "m.fasta", {"q": "ACGTAC", "r": "ACGTAC"})
    with pytest.raises(UserInputError):
        compute_similarity(str(msa), "q", window_size=3, window_step=3, metric="jc69")


def test_window_larger_than_alignment_raises(tmp_path: Path) -> None:
    msa = write_fasta(tmp_path / "m.fasta", {"q": "ACGT", "r": "ACGT"})
    with pytest.raises(UserInputError):
        compute_similarity(str(msa), "q", window_size=100, window_step=10)


def test_duplicate_identifiers_raise(tmp_path: Path) -> None:
    # A repeated header would otherwise overwrite the earlier row in the label ->
    # sequence mapping, shrinking the panel without saying so.
    msa = tmp_path / "dup.fasta"
    msa.write_text(">query\nACGTACGT\n>refA\nACGTACGT\n>refA\nTTTTTTTT\n")
    with pytest.raises(UserInputError, match="Duplicate sequence identifier"):
        _read_alignment(str(msa))


def test_ragged_rows_name_the_offending_record(tmp_path: Path) -> None:
    # Biopython's own alignment reader raises first and names no record; the error
    # the user sees should identify which sequence is the wrong length.
    msa = write_fasta(tmp_path / "m.fasta", {"query": "ACGTACGT", "refA": "ACGT"})
    with pytest.raises(UserInputError, match="refA"):
        _read_alignment(str(msa))


# --- informative-site windowing (low-divergence mode) ----------------------

def test_informative_column_mask_flags_only_variable_columns(tmp_path: Path) -> None:
    msa = write_fasta(tmp_path / "m.fasta", {
        "query": "ACGT",
        "refA": "ACGT",   # matches everywhere
        "refB": "ATGT",   # differs at col 1
        "refC": "ACGA",   # differs at col 3
    })
    rows = _read_alignment(str(msa))
    mask = _informative_column_mask(rows, ["refA", "refB", "refC"])
    assert mask.tolist() == [False, True, False, True]
    assert informative_site_count(rows, "query") == 2


def test_informative_mask_ignores_gaps_and_single_ref(tmp_path: Path) -> None:
    # col0: A vs gap -> only one canonical ref -> not informative.
    # col1: C vs C -> agree -> not informative. col2: G vs T -> informative.
    msa = write_fasta(tmp_path / "m.fasta", {
        "query": "ACG", "refA": "ACG", "refB": "-CT",
    })
    rows = _read_alignment(str(msa))
    assert _informative_column_mask(rows, ["refA", "refB"]).tolist() == [False, False, True]


def test_informative_windowing_recovers_low_divergence_recombinant(tmp_path: Path) -> None:
    # 600 bp, only 5 SNPs (<1% divergence): bp-windowing has almost no contrast,
    # informative-site windowing recovers the A-backbone / B-insert mosaic.
    pa = bytearray(b"A" * 600)
    pb = bytearray(b"A" * 600)
    pa[100] = pa[200] = ord("G")   # parentA-specific (first half)
    pb[400] = pb[500] = ord("C")   # parentB-specific (second half)
    pa[300] = ord("G")
    pb[300] = ord("C")             # col 300 distinguishes the two parents
    query = bytes(pa[:300]) + bytes(pb[300:])  # backbone A, insert B
    msa = write_fasta(tmp_path / "msa.fasta", {
        "query": query.decode(),
        "parentA": bytes(pa).decode(),
        "parentB": bytes(pb).decode(),
    })

    result = compute_similarity_informative(str(msa), "query", info_window=2, info_step=1)
    # windows span the informative columns, not fixed bp
    assert result.window_spans
    assert result.numerators["parentA"][0] > result.numerators["parentB"][0]   # A early
    assert result.numerators["parentB"][-1] > result.numerators["parentA"][-1]  # B late

    from tessera.recomb.hmm import segment_query
    segments, major = segment_query(result, jump_rate=1e-3)
    states = [s.state for s in segments]
    assert states[0] == "parentA" and "parentB" in states
    switch = next(s for s in segments if s.state == "parentB")
    assert 200 <= switch.query_start <= 400  # breakpoint near the true splice at 300


def test_informative_too_few_sites_raises(tmp_path: Path) -> None:
    msa = write_fasta(tmp_path / "m.fasta", {"query": "ACGTACGT", "refA": "ACGTACGT"})
    with pytest.raises(UserInputError, match="informative site"):
        compute_similarity_informative(str(msa), "query", info_window=5, info_step=1)


def test_rna_alphabet_reads_identically_to_dna(tmp_path: Path) -> None:
    """U is the same base as T for identity; an RNA record must not go dark.

    Tessera's whole target domain is RNA viruses, and NCBI/Nextclade will sometimes
    serve U-alphabet FASTA. Treating U as non-comparable silently drops every such
    column instead of failing loudly.
    """
    records = {
        "query": "ACGTACGTACGTACGTACGT",
        "refA": "ACGTACGTACGTACGTACGT",
        "refB": "ACGTACGTACGTTTTTTTTT",
    }
    dna = compute_similarity(
        str(write_fasta(tmp_path / "dna.fasta", records)), "query",
        window_size=10, window_step=10,
    )
    rna = compute_similarity(
        str(write_fasta(tmp_path / "rna.fasta",
                        {k: v.replace("T", "U") for k, v in records.items()})),
        "query", window_size=10, window_step=10,
    )
    assert rna.similarities == dna.similarities
    assert rna.denominators == dna.denominators


def test_read_alignment_normalises_u_to_t(tmp_path: Path) -> None:
    msa = write_fasta(tmp_path / "m.fasta", {"a": "ACGU", "b": "acgu"})
    rows = _read_alignment(str(msa))
    assert rows["a"].tobytes() == b"ACGT"
    assert rows["b"].tobytes() == b"ACGT"


def test_low_canonical_records_flags_unreadable_rows(tmp_path: Path) -> None:
    """A record that is mostly not A/C/G/T contributes almost nothing; say so."""
    msa = write_fasta(tmp_path / "m.fasta", {
        "good": "ACGTACGTAC",
        "gappy": "AC--------",   # gaps are expected in an MSA -- not flagged
        "junk": "ACRYKMSWBV",    # IUPAC ambiguity throughout -- flagged
    })
    rows = _read_alignment(str(msa))
    flagged = dict(low_canonical_records(rows, min_fraction=0.5))
    assert "junk" in flagged
    assert "good" not in flagged
    assert "gappy" not in flagged   # judged over non-gap positions only
