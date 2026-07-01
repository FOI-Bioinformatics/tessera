"""The shared recombinant-lineage name test."""

from __future__ import annotations

from tessera.recomb.typing import is_recombinant_lineage


def test_recombinant_names_detected():
    assert is_recombinant_lineage("CRF01_AE")
    assert is_recombinant_lineage("URF1")
    assert is_recombinant_lineage("XBB")          # Pango recombinant (X + letter)
    assert is_recombinant_lineage("recombinant")  # literal, case-insensitive
    assert is_recombinant_lineage("RECOMBINANT")


def test_non_recombinant_names_pass():
    assert not is_recombinant_lineage("A.1")
    assert not is_recombinant_lineage("B.1.1.7")
    assert not is_recombinant_lineage("A.D.1.8")
    assert not is_recombinant_lineage("clade 9")
    assert not is_recombinant_lineage("")
