"""Lineage assignment: de-novo clustering (pure) and nearest-reference (stubbed skani)."""

from __future__ import annotations

import logging
from pathlib import Path

from tessera.discover import lineage_assign as la


def _g(name: str) -> Path:
    return Path(f"/pool/{name}.fasta")


def test_cluster_by_ani_connected_components():
    a, b, c, d = _g("a"), _g("b"), _g("c"), _g("d")
    # a-b-c linked (transitively one lineage); d isolated (its own)
    clusters = la.cluster_by_ani([a, b, c, d], [(a, b), (b, c)])
    groups = sorted(sorted(p.name for p in members) for members in clusters.values())
    assert groups == [["a.fasta", "b.fasta", "c.fasta"], ["d.fasta"]]
    assert all(k.startswith("denovo_") for k in clusters)


def test_cluster_by_ani_all_singletons_when_no_edges():
    a, b = _g("a"), _g("b")
    clusters = la.cluster_by_ani([a, b], [])
    assert sorted(len(v) for v in clusters.values()) == [1, 1]


def test_assign_by_nearest_reference(monkeypatch):
    g1, g2 = _g("q1"), _g("q2")
    tipA, tipB = _g("refA"), _g("refB")
    labeled = {tipA: "A.1", tipB: "B.2"}

    # q1 is closest to tipA (95%) and clears the floor; q2's best (85%) is below it.
    fake = {
        g1: {tipA: (95.0, 90.0), tipB: (80.0, 88.0)},
        g2: {tipA: (85.0, 90.0), tipB: (82.0, 80.0)},
    }
    monkeypatch.setattr(la, "skani_query_ani", lambda q, refs, logger: fake[q])

    out = la.assign_by_nearest_reference([g1, g2], labeled, ani_floor=90.0,
                                         logger=logging.getLogger("t"))
    assert out == {"q1": "A.1"}          # q1 assigned; q2 below floor -> omitted
