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


_REF = ">NC ref\n" + "ACGTACGTAC" * 1 + "\n"  # 10 bp reference

_TREE = {
    "tree": {
        "name": "root",
        "branch_attrs": {"mutations": {"nuc": ["A1G"]}},
        "children": [
            {"name": "TIP1",
             "branch_attrs": {"mutations": {"nuc": ["C2T"]}},
             "node_attrs": {"accession": "ACC1", "clade_membership": {"value": "X"}},
             "children": []},
            {"name": "TIP2",
             "branch_attrs": {"mutations": {"nuc": ["T4-"]}},
             "node_attrs": {"accession": "ACC2", "clade_membership": "Y"},
             "children": []},
        ],
    },
}

_EXAMPLES = ">EX1 example\n" + "ACGTACGTAC" + "\n"


def _patch_files(monkeypatch):
    payloads = {
        "reference.fasta": _REF.encode(),
        "tree.json": _json.dumps(_TREE).encode(),
        "sequences.fasta": _EXAMPLES.encode(),
    }

    def fake_urlopen(url, timeout=0):
        for name, data in payloads.items():
            if url.endswith(name):
                return _Resp(data)
        raise AssertionError(f"unexpected url {url}")

    monkeypatch.setattr(nc, "urlopen", fake_urlopen)


def _dataset():
    return nc.NextcladeDataset(
        path="x/y", tag="t",
        files={"reference": "reference.fasta", "treeJson": "tree.json",
               "examples": "sequences.fasta"},
    )


def test_build_pool_reconstructs_tips_examples_and_labels(monkeypatch, tmp_path, logger):
    _patch_files(monkeypatch)
    genomes = nc.build_pool(_dataset(), cache_dir=tmp_path / "c", logger=logger)
    by_name = {g.stem: g.read_text() for g in genomes}
    assert {"ACC1", "ACC2", "EX1"} <= set(by_name)
    # TIP1: A1G then C2T on ACGTACGTAC -> GTGTACGTAC, header carries clade X
    assert ">ACC1 X" in by_name["ACC1"]
    assert "GTGTACGTAC" in by_name["ACC1"].replace("\n", "")
    # TIP2: A1G then T4 deleted -> GCGACGTAC (gap stripped), clade Y
    assert ">ACC2 Y" in by_name["ACC2"]


def test_build_pool_uses_cache_on_second_call(monkeypatch, tmp_path, logger):
    _patch_files(monkeypatch)
    cache = tmp_path / "c"
    nc.build_pool(_dataset(), cache_dir=cache, logger=logger)

    def boom(url, timeout=0):
        raise AssertionError("network must not be touched on a cache hit")

    monkeypatch.setattr(nc, "urlopen", boom)
    again = nc.build_pool(_dataset(), cache_dir=cache, logger=logger)
    assert {g.stem for g in again} >= {"ACC1", "ACC2", "EX1"}


def test_build_pool_examples_failure_is_nonfatal(monkeypatch, tmp_path, logger):
    """A URLError on the examples file must not abort build_pool; tips are still returned."""
    from urllib.error import URLError as _URLError

    payloads = {
        "reference.fasta": _REF.encode(),
        "tree.json": _json.dumps(_TREE).encode(),
    }

    def fake_urlopen(url, timeout=0):
        for name, data in payloads.items():
            if url.endswith(name):
                return _Resp(data)
        if url.endswith("sequences.fasta"):
            raise _URLError("boom")
        raise AssertionError(f"unexpected url: {url}")

    monkeypatch.setattr(nc, "urlopen", fake_urlopen)
    genomes = nc.build_pool(_dataset(), cache_dir=tmp_path / "c", logger=logger)
    names = {g.stem for g in genomes}
    assert "ACC1" in names
    assert "ACC2" in names


def test_build_pool_handles_deep_tree(monkeypatch, tmp_path, logger):
    """Iterative tree walk must not raise RecursionError on a 1200-node linear chain.

    Both json.dumps and json.loads also recurse on deeply nested dicts, so we
    build the Python dict iteratively and inject it via patched helpers to keep
    the test self-contained and fast.
    """
    import json as _stdlib_json

    depth = 1200
    # Build the deep Python dict iteratively -- Python dict construction never recurses.
    tip: dict = {
        "name": "TIP_DEEP",
        "branch_attrs": {"mutations": {"nuc": []}},
        "node_attrs": {"accession": "DEEP1", "clade_membership": "deep"},
        "children": [],
    }
    node: dict = tip
    for _ in range(depth):
        node = {
            "name": "internal",
            "branch_attrs": {"mutations": {"nuc": []}},
            "node_attrs": {},
            "children": [node],
        }
    deep_tree_obj = {"tree": node}

    # Bypass _download_text so we never serialise/parse the deep dict via JSON.
    def fake_download_text(dataset, role, logger_):  # type: ignore[override]
        if role == "reference":
            return _REF
        return "__DEEP_SENTINEL__"  # intercepted by patched_loads below

    _real_loads = _stdlib_json.loads

    def patched_loads(s, **kwargs):
        if s == "__DEEP_SENTINEL__":
            return deep_tree_obj
        return _real_loads(s, **kwargs)

    monkeypatch.setattr(nc, "_download_text", fake_download_text)
    monkeypatch.setattr(_stdlib_json, "loads", patched_loads)

    ds = nc.NextcladeDataset(
        path="x/y", tag="t",
        files={"reference": "reference.fasta", "treeJson": "tree.json"},
    )
    genomes = nc.build_pool(ds, cache_dir=tmp_path / "c", logger=logger)
    assert any(g.stem == "DEEP1" for g in genomes)
