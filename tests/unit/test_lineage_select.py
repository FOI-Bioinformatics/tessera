"""Unit tests for the lineage-aware panel reduction in discover/pool.py.

``_lineage_select`` is exercised directly with a stubbed ANI map and a stubbed
``dereplicate`` so no skani/skDER binary is needed.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from tessera.discover import pool as pool_mod

_LOG = logging.getLogger("tessera.test")


@pytest.fixture
def stub_dereplicate(monkeypatch):
    """Stub ``dereplicate`` to keep the first genome of the input (a stand-in for
    ANI collapse) and record what it was called with."""
    calls = {}

    def fake(genomes, *, ani, logger):
        calls["genomes"] = list(genomes)
        calls["ani"] = ani
        return (genomes[:1], genomes[1:])

    monkeypatch.setattr(pool_mod, "dereplicate", fake)
    return calls


def _p(name: str) -> Path:
    return Path(f"/pool/{name}.fasta")


def test_keeps_one_query_closest_rep_per_lineage(stub_dereplicate):
    a1, a2, b1 = _p("a1"), _p("a2"), _p("b1")
    lineage_of = {"a1": "A.1", "a2": "A.1", "b1": "B.1"}
    query_ani = {a1: (95.0, 90.0), a2: (98.0, 90.0), b1: (80.0, 90.0)}
    kept = pool_mod._lineage_select(
        [a1, a2, b1], lineage_of, query_ani,
        keep_recombinant=False, derep_ani=99.0, logger=_LOG,
    )
    assert set(kept) == {a2, b1}  # a2 beats a1 within A.1; b1 is B.1's only member


def test_excludes_recombinant_lineages_by_default(stub_dereplicate):
    a1, x1 = _p("a1"), _p("x1")
    lineage_of = {"a1": "A.1", "x1": "CRF01_AE"}
    query_ani = {a1: (95.0, 90.0), x1: (99.0, 90.0)}
    kept = pool_mod._lineage_select(
        [a1, x1], lineage_of, query_ani,
        keep_recombinant=False, derep_ani=99.0, logger=_LOG,
    )
    assert kept == [a1]  # the CRF is dropped even though it is the closest match


def test_keeps_recombinant_lineages_when_toggled(stub_dereplicate):
    a1, x1 = _p("a1"), _p("x1")
    lineage_of = {"a1": "A.1", "x1": "CRF01_AE"}
    query_ani = {a1: (95.0, 90.0), x1: (99.0, 90.0)}
    kept = pool_mod._lineage_select(
        [a1, x1], lineage_of, query_ani,
        keep_recombinant=True, derep_ani=99.0, logger=_LOG,
    )
    assert set(kept) == {a1, x1}  # the CRF competes as its own lineage


def test_untyped_genomes_route_to_ani_fallback(stub_dereplicate):
    a1, u1, u2, u3 = _p("a1"), _p("u1"), _p("u2"), _p("u3")
    lineage_of = {"a1": "A.1"}  # u1..u3 untyped
    query_ani = {a1: (95.0, 90.0)}
    kept = pool_mod._lineage_select(
        [a1, u1, u2, u3], lineage_of, query_ani,
        keep_recombinant=False, derep_ani=99.0, logger=_LOG,
    )
    # a1 kept as its lineage rep; the 3 untyped go through dereplicate (stub keeps u1).
    assert set(kept) == {a1, u1}
    assert stub_dereplicate["genomes"] == [u1, u2, u3]
    assert stub_dereplicate["ani"] == 99.0


def test_empty_map_sends_everything_to_fallback(stub_dereplicate):
    u1, u2, u3 = _p("u1"), _p("u2"), _p("u3")
    kept = pool_mod._lineage_select(
        [u1, u2, u3], {}, {},
        keep_recombinant=False, derep_ani=99.0, logger=_LOG,
    )
    assert kept == [u1]  # all untyped -> the ANI fallback (stub keeps the first)
