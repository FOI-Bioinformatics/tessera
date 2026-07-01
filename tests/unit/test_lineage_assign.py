"""Lineage assignment: de-novo clustering (pure) and nearest-reference (stubbed skani)."""

from __future__ import annotations

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
