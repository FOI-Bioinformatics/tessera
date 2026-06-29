"""Per-lineage consensus sequence builder."""

from __future__ import annotations

import pytest

from tessera.recomb.consensus import consensus_sequence


def test_majority_base_per_column() -> None:
    # col0 all A; col1 A,A,C -> A; col2 G,T,T -> T; col3 C,C,C
    seqs = ["AAGC", "AATC", "ACTC"]
    assert consensus_sequence(seqs) == "AATC"


def test_single_sequence_returned_as_is() -> None:
    assert consensus_sequence(["acgt"]) == "ACGT"


def test_non_canonical_column_keeps_first_row() -> None:
    # col1 has no canonical base in any row (gap / N) -> first row's char kept
    seqs = ["A-CC", "A-CC", "ANCC"]  # col1 = '-','-','N'
    result = consensus_sequence(seqs)
    assert result[0] == "A" and result[2:] == "CC"
    assert result[1] == "-"  # the first row's character at that column


def test_unequal_lengths_rejected() -> None:
    with pytest.raises(ValueError):
        consensus_sequence(["ACGT", "ACG"])
