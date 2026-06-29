"""Pango recombinant cross-check: parent expansion, cached load, report block."""

from __future__ import annotations

import json
from pathlib import Path

from tessera.core.cache import pango_alias_path
from tessera.recomb.pango import crosscheck_html, expand_recombinant, load_alias_key

# A trimmed alias_key.json: recombinants map to a list of parents, aliases to a string.
ALIAS = {
    "XBB": ["BJ.1", "CJ.1"],
    "XD": ["B.1.617.2", "BA.1"],
    "AY": "B.1.617.2",  # an ordinary alias (single string)
    "XEMPTY": [],
}


def test_expand_recombinant_returns_parents():
    assert expand_recombinant("XBB", ALIAS) == ["BJ.1", "CJ.1"]


def test_expand_recombinant_resolves_sublineage_via_root():
    assert expand_recombinant("XBB.1.5", ALIAS) == ["BJ.1", "CJ.1"]


def test_expand_recombinant_non_recombinant_is_empty():
    assert expand_recombinant("B.1.1.7", ALIAS) == []
    assert expand_recombinant("AY.4", ALIAS) == []  # an alias, not a recombinant
    assert expand_recombinant("", ALIAS) == []


def test_expand_recombinant_unknown_or_empty_list():
    assert expand_recombinant("XZZ", ALIAS) == []
    assert expand_recombinant("XEMPTY", ALIAS) == []


def test_load_alias_key_reads_cache_without_network(tmp_path: Path, logger):
    path = pango_alias_path(tmp_path)
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(ALIAS))
    # offline=True guarantees no fetch is attempted; the cached file is returned.
    assert load_alias_key(cache_override=tmp_path, logger=logger, offline=True) == ALIAS


def test_load_alias_key_missing_offline_is_empty(tmp_path: Path, logger):
    assert load_alias_key(cache_override=tmp_path, logger=logger, offline=True) == {}


def test_crosscheck_html_lists_query_and_parents():
    html = crosscheck_html("XBB.1.5", ["BJ.1", "CJ.1"])
    assert "XBB.1.5" in html and "BJ.1" in html and "CJ.1" in html
