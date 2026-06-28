"""Nextclade pool source: reconstruction, labeling, resolution, caching."""

from __future__ import annotations

import json as _json
import logging

import pytest

from recomfi.core.errors import UserInputError
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


_FAKE_INDEX = {
    "collections": [{
        "datasets": [
            {"path": "community/neherlab/hiv-1/hxb2",
             "versions": [{"tag": "2026-04-14--11-55-23Z"}, {"tag": "2024-01-01--00-00-00Z"}],
             "files": {"reference": "reference.fasta", "treeJson": "tree.json",
                       "examples": "sequences.fasta"}},
            {"path": "nextstrain/sars-cov-2/XBB",
             "versions": [{"tag": "2025-09-09--12-13-13Z"}],
             "files": {"reference": "reference.fasta", "treeJson": "tree.json",
                       "examples": "sequences.fasta"}},
        ],
    }],
}


class _Resp:
    def __init__(self, data: bytes):
        self._data = data

    def read(self):
        return self._data

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def _patch_index(monkeypatch):
    monkeypatch.setattr(nc, "urlopen",
                        lambda url, timeout=0: _Resp(_json.dumps(_FAKE_INDEX).encode()))


@pytest.fixture
def logger():
    return logging.getLogger("test")


def test_fetch_index_flattens_and_keeps_newest_tag(monkeypatch):
    _patch_index(monkeypatch)
    idx = nc.fetch_index()
    assert set(idx) == {"community/neherlab/hiv-1/hxb2", "nextstrain/sars-cov-2/XBB"}
    assert idx["community/neherlab/hiv-1/hxb2"]["versions"][0]["tag"] == "2026-04-14--11-55-23Z"


def test_resolve_override_validates_against_index(monkeypatch, tmp_path, logger):
    _patch_index(monkeypatch)
    q = tmp_path / "q.fasta"
    q.write_text(">q\nACGT\n")
    ds = nc.resolve_dataset(q, "nextstrain/sars-cov-2/XBB", email=None, logger=logger)
    assert ds.path == "nextstrain/sars-cov-2/XBB"
    assert ds.tag == "2025-09-09--12-13-13Z"
    assert ds.files["treeJson"] == "tree.json"


def test_resolve_override_unknown_path_raises(monkeypatch, tmp_path, logger):
    _patch_index(monkeypatch)
    q = tmp_path / "q.fasta"
    q.write_text(">q\nACGT\n")
    with pytest.raises(UserInputError, match="not a Nextclade dataset"):
        nc.resolve_dataset(q, "nextstrain/nope", email=None, logger=logger)


def test_resolve_blast_fallback_maps_organism_to_dataset(monkeypatch, tmp_path, logger):
    _patch_index(monkeypatch)
    monkeypatch.setattr(nc.shutil, "which", lambda name: None)  # no nextclade binary
    monkeypatch.setattr(nc, "detect_taxon",
                        lambda q, *, email, logger: "Human immunodeficiency virus 1")
    q = tmp_path / "q.fasta"
    q.write_text(">q\nACGT\n")
    ds = nc.resolve_dataset(q, None, email=None, logger=logger)
    assert ds.path == "community/neherlab/hiv-1/hxb2"


def test_resolve_blast_fallback_unmapped_raises(monkeypatch, tmp_path, logger):
    _patch_index(monkeypatch)
    monkeypatch.setattr(nc.shutil, "which", lambda name: None)
    monkeypatch.setattr(nc, "detect_taxon", lambda q, *, email, logger: "Unknown mystery virus")
    q = tmp_path / "q.fasta"
    q.write_text(">q\nACGT\n")
    with pytest.raises(UserInputError, match="--nextclade-dataset"):
        nc.resolve_dataset(q, None, email=None, logger=logger)
