"""Native sliding-window similarity engine and query-coordinate mapping."""

from __future__ import annotations

from math import isnan
from pathlib import Path

import pytest

from recomfi.core.errors import UserInputError
from recomfi.recomb.similarity import compute_similarity

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
