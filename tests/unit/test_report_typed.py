"""Report rendering of typed names: the query-lineage verdict clause."""

from __future__ import annotations

from tessera.recomb.report import _verdict_html

# A minimal summary dict as produced by _summary for a no-recombination call.
_S = {
    "query_len": 1000, "major": "MK573073", "n_regions": 0, "n_absent": 0,
    "recomb_bp": 0, "pct": 0.0, "minors": [], "confidence": None,
}


def test_verdict_shows_query_lineage_when_it_differs_from_label():
    html = _verdict_html(_S, "myquery", {}, None, "GII.P16-GII.4")
    assert "typed as" in html and "GII.P16-GII.4" in html


def test_verdict_omits_clause_when_lineage_equals_label():
    html = _verdict_html(_S, "GII.P16-GII.4", {}, None, "GII.P16-GII.4")
    assert "typed as" not in html


def test_verdict_omits_clause_when_no_query_lineage():
    assert "typed as" not in _verdict_html(_S, "myquery", {}, None, None)
