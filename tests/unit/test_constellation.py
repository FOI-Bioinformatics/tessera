"""Unit tests for the pure reassortment call (no binaries, no I/O)."""

from __future__ import annotations

from tessera.reassort.constellation import call_constellation


def _uni(*strains):
    return set(strains)


def test_two_concordant_segments_are_clonal():
    # HA and NA both have strain S1 within-margin of their best -> one parent -> clonal.
    candidates = {"HA": [("S1", 99.0), ("S2", 98.9)], "NA": [("S1", 99.0), ("S3", 98.8)]}
    universes = {"HA": _uni("S1", "S2", "S3"), "NA": _uni("S1", "S2", "S3")}
    c = call_constellation(candidates, universes, margin=0.5)
    assert c.verdict == "clonal"
    assert len(c.groups) == 1
    assert set(c.groups[0].segments) == {"HA", "NA"}
    assert "S1" in c.groups[0].parent_strains


def test_one_discordant_pair_is_reassortant():
    # HA best = S1 (also typed in NA's tree but far there); NA best = S3 (typed in HA's tree,
    # far there) -> cross-typed and disagree -> reassortant, two groups.
    candidates = {"HA": [("S1", 99.0)], "NA": [("S3", 99.0)]}
    universes = {"HA": _uni("S1", "S3"), "NA": _uni("S1", "S3")}
    c = call_constellation(candidates, universes, margin=0.5)
    assert c.verdict == "reassortant"
    assert len(c.groups) == 2


def test_sparse_cross_typing_clonal_isolate_is_not_reassortant():
    # Eight segments of one clonal isolate. Each names a different best strain, and those
    # strains are NOT cross-typed into the other segments' trees -> every pair uninformative
    # -> must NOT be called reassortant (this is the v1 #1 regression).
    segs = [f"seg{i}" for i in range(8)]
    candidates = {s: [(f"strain_{s}", 99.0)] for s in segs}
    universes = {s: {f"strain_{s}"} for s in segs}  # no strain appears in another segment's tree
    c = call_constellation(candidates, universes, margin=0.5)
    assert c.verdict != "reassortant"
    assert c.verdict == "undetermined"


def test_seven_plus_one_reassortant_groups_the_seven():
    # Seven segments share strain A within margin; the eighth prefers B (typed in the others,
    # far there) -> reassortant; the seven in one group, the eighth alone.
    seven = [f"seg{i}" for i in range(7)]
    candidates = {s: [("A", 99.0)] for s in seven}
    candidates["odd"] = [("B", 99.0)]
    universe = {"A", "B"}
    universes = {s: set(universe) for s in [*seven, "odd"]}
    c = call_constellation(candidates, universes, margin=0.5)
    assert c.verdict == "reassortant"
    big = max(c.groups, key=lambda g: len(g.segments))
    assert set(big.segments) == set(seven)
    assert any(g.segments == ["odd"] for g in c.groups)


def test_partial_linkage_without_disagreement_is_undetermined():
    # HA-NA concordant on S1; PB2 best is S9, not typed in HA/NA trees, and HA/NA best not in
    # PB2's tree -> PB2 pairs uninformative -> not all linked, no disagreement -> undetermined.
    candidates = {"HA": [("S1", 99.0)], "NA": [("S1", 99.0)], "PB2": [("S9", 99.0)]}
    universes = {"HA": {"S1"}, "NA": {"S1"}, "PB2": {"S9"}}
    c = call_constellation(candidates, universes, margin=0.5)
    assert c.verdict == "undetermined"


def test_fewer_than_two_segments_is_undetermined():
    c = call_constellation({"HA": [("S1", 99.0)]}, {"HA": {"S1"}}, margin=0.5)
    assert c.verdict == "undetermined"
