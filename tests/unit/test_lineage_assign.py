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


def test_assign_lineages_ladder_precedence(monkeypatch):
    m1, ref1, dn1, dn2 = _g("m1"), _g("ref1"), _g("dn1"), _g("dn2")
    genomes = [m1, ref1, dn1, dn2]

    # (1) metadata types only m1
    monkeypatch.setattr(la, "build_lineage_map",
                        lambda **kw: [("m1", "B.1", "title")])
    monkeypatch.setattr(la, "titles_from_collection", lambda files: {})
    # (2) a reference exists; nearest-neighbour types ref1
    monkeypatch.setattr(la, "_reference_tips", lambda **kw: {_g("tip"): "A.1"})
    monkeypatch.setattr(la, "assign_by_nearest_reference",
                        lambda genomes, tips, *, ani_floor, logger: {"ref1": "A.1"})
    # (3) de-novo clusters the remaining dn1, dn2 together
    monkeypatch.setattr(la, "_ani_edges", lambda genomes, *, threshold, logger: [(dn1, dn2)])

    rows = la.assign_lineages(genomes, taxon="Test virus", logger=logging.getLogger("t"))
    by_label = {label: (lin, src) for label, lin, src in rows}
    assert by_label["m1"] == ("B.1", "title")
    assert by_label["ref1"] == ("A.1", "nextclade-nn")
    assert by_label["dn1"][1] == "denovo" and by_label["dn2"][1] == "denovo"
    assert by_label["dn1"][0] == by_label["dn2"][0]      # same de-novo lineage
