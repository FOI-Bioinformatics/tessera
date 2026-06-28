"""Nextclade pool source: reconstruction, labeling, resolution, caching."""

from __future__ import annotations

from recomfi.discover import nextclade as nc


def test_reconstruct_applies_subs_and_deletions_and_strips_gaps():
    # reference ACGT; A1G sets pos1 -> G; T4- deletes pos4 -> stripped.
    assert nc._reconstruct_sequence("ACGT", ["A1G", "T4-"]) == "GCG"


def test_reconstruct_ignores_out_of_range_and_malformed():
    assert nc._reconstruct_sequence("ACGT", ["A9G", "garbage", "C2T"]) == "ATGT"


def test_clade_priority_prefers_pango_then_universal_membership():
    assert nc._clade_of({"Nextclade_pango": {"value": "BA.2"},
                         "clade_membership": {"value": "21L"}}) == "BA.2"
    assert nc._clade_of({"clade_membership": "H2"}) == "H2"  # bare scalar
    assert nc._clade_of({"author": {"value": "x"}}) == "NA"  # no clade key


def test_accession_prefers_attr_then_name_and_is_sanitised():
    assert nc._accession_of({"node_attrs": {"accession": {"value": "AY/03 7026"}},
                             "name": "ignored"}) == "AY_03_7026"
    assert nc._accession_of({"node_attrs": {}, "name": "JN417241"}) == "JN417241"
