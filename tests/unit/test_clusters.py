"""Lineage clustering: near-duplicates merge, distinct parents stay separate."""

from __future__ import annotations

import random
from pathlib import Path

from tessera.recomb.analyze import analyze
from tessera.recomb.clusters import all_singletons, cluster_references, clustered_view
from tessera.recomb.regions import RegionParams, call_regions
from tessera.recomb.similarity import compute_similarity

from ..conftest import write_fasta
from ..factories import _mutate

_BASES = "ACGT"
WIN, STEP, L = 300, 75, 6000
R1 = 4000  # region 1 = lineage A, region 2 = lineage B


def _lanes(seed: int = 1, divergence: float = 0.20) -> tuple[str, str]:
    rng = random.Random(seed)
    lane_a = "".join(rng.choice(_BASES) for _ in range(L))
    return lane_a, _mutate(lane_a, divergence, rng)


def _result(records: dict[str, str], tmp_path: Path):
    msa = tmp_path / "msa.fasta"
    write_fasta(msa, records)
    return compute_similarity(str(msa), "query", window_size=WIN, window_step=STEP)


def _params() -> RegionParams:
    return RegionParams.with_defaults(WIN, min_region=WIN)


def _cluster_labels(result, params) -> list[set[str]]:
    return [set(c) for c in cluster_references(result, WIN, params)]


def test_near_duplicates_merge_distinct_parents_separate(tmp_path) -> None:
    lane_a, lane_b = _lanes()
    query = lane_a[:R1] + lane_b[R1:]
    rng = random.Random(7)
    records = {"query": query}
    for i in range(4):  # a recombinant sibling lineage (near-identical to the query)
        records[f"sib{i}"] = _mutate(query, 0.001, rng)
    for i in range(3):  # parent-A lineage cloud
        records[f"pa{i}"] = _mutate(lane_a, 0.001, rng)
    for i in range(3):  # parent-B lineage cloud
        records[f"pb{i}"] = _mutate(lane_b, 0.001, rng)

    clusters = _cluster_labels(_result(records, tmp_path), _params())
    assert len(clusters) == 3
    assert {"sib0", "sib1", "sib2", "sib3"} in clusters
    assert {"pa0", "pa1", "pa2"} in clusters
    assert {"pb0", "pb1", "pb2"} in clusters


def test_diffuse_distinct_parents_stay_separate(tmp_path) -> None:
    # Two references ~2% apart everywhere (no structural tract) -- only an absolute
    # floor separates them; they must NOT merge (distinct orthopox-like species).
    lane_a, _ = _lanes()
    rng = random.Random(8)
    records = {
        "query": lane_a,
        "variola": _mutate(lane_a, 0.02, rng),
        "cowpox": _mutate(lane_a, 0.001, rng),
    }
    clusters = _cluster_labels(_result(records, tmp_path), _params())
    assert all_singletons(clusters)


def test_structural_close_parents_stay_separate(tmp_path) -> None:
    # Two references ~99% genome-wide but differing over one contiguous tract (a
    # SARS-sublineage-like pair): the contiguous-run test must keep them separate.
    lane_a, lane_b = _lanes()
    rng = random.Random(9)
    base = _mutate(lane_a, 0.001, rng)
    other = list(_mutate(lane_a, 0.001, rng))
    other[R1:] = list(_mutate(lane_b[R1:], 0.003, rng))  # a divergent tract in region 2
    records = {"query": lane_a, "parent_x": base, "parent_y": "".join(other)}
    clusters = _cluster_labels(_result(records, tmp_path), _params())
    assert all_singletons(clusters)


def test_clustered_view_pools_counts_and_is_identity_for_singletons(tmp_path) -> None:
    lane_a, lane_b = _lanes()
    query = lane_a[:R1] + lane_b[R1:]
    rng = random.Random(10)
    records = {"query": query, "sib0": _mutate(query, 0.001, rng),
               "sib1": _mutate(query, 0.001, rng), "pa0": _mutate(lane_a, 0.001, rng)}
    result = _result(records, tmp_path)

    # Pooling two members sums their per-window counts.
    view, sizes = clustered_view(result, [["sib0", "sib1"], ["pa0"]])
    rep = next(r for r, n in sizes.items() if n == 2)
    for i in range(len(result.positions)):
        assert view.numerators[rep][i] == (
            result.numerators["sib0"][i] + result.numerators["sib1"][i]
        )
        assert view.denominators[rep][i] == (
            result.denominators["sib0"][i] + result.denominators["sib1"][i]
        )

    # All-singleton view reproduces the input arrays exactly (backward compatibility).
    singleton_view, sizes2 = clustered_view(result, [["sib0"], ["sib1"], ["pa0"]])
    assert set(sizes2.values()) == {1}
    assert singleton_view.numerators["pa0"] == result.numerators["pa0"]


def test_cluster_caller_excludes_sibling_lineage_and_calls_one_region(tmp_path) -> None:
    # End-to-end: a duplicate-heavy panel. With clustering the sibling lineage is
    # pooled and excluded, and the A/B mosaic is one region with a pooled minor.
    lane_a, lane_b = _lanes()
    query = lane_a[:R1] + lane_b[R1:]
    rng = random.Random(11)
    records = {"query": query}
    for i in range(4):
        records[f"sib{i}"] = _mutate(query, 0.001, rng)  # closest everywhere -> would mask
    for i in range(3):
        records[f"pa{i}"] = _mutate(lane_a, 0.001, rng)  # major lineage cloud
    for i in range(3):
        records[f"pb{i}"] = _mutate(lane_b, 0.001, rng)  # minor lineage cloud (region 2)
    result = _result(records, tmp_path)

    regions, major, dropped = call_regions(result, analyze(result), WIN, _params())
    assert major.startswith("pa")  # a parent-A representative, not a sibling
    assert {ev.label for ev in dropped}  # the sibling lineage was excluded
    minor_regions = [r for r in regions if r.minor_parent.startswith("pb")]
    assert len(minor_regions) == 1  # one region, not fragments
    assert minor_regions[0].minor_cluster_size == 3  # reported as a pooled lineage
