"""Sibling-aware caller: a whole-genome twin is set aside, parents are kept."""

from __future__ import annotations

import random
from pathlib import Path

from recomfi.recomb.analyze import analyze
from recomfi.recomb.regions import RegionParams, call_regions
from recomfi.recomb.similarity import compute_similarity

from ..conftest import write_fasta
from ..factories import _mutate

_BASES = "ACGT"
WIN, STEP, L = 300, 75, 6000  # ~77 windows -> stable lead/beaten fractions
R1, R2 = 4000, 6000  # region 1 = [0,4000) lineage A, region 2 = [4000,6000) lineage B


def _lineages(seed: int = 1, divergence: float = 0.20) -> tuple[str, str]:
    rng = random.Random(seed)
    lane_a = "".join(rng.choice(_BASES) for _ in range(L))
    lane_b = _mutate(lane_a, divergence, rng)  # a clearly distinct lineage
    return lane_a, lane_b


def _recombinant_query(lane_a: str, lane_b: str) -> str:
    return lane_a[:R1] + lane_b[R1:]  # major lineage A, minor tract from lineage B


def _result_and_regions(records: dict[str, str], tmp_path: Path, *, exclude: bool = True):
    msa = tmp_path / "msa.fasta"
    write_fasta(msa, records)
    result = compute_similarity(str(msa), "query", window_size=WIN, window_step=STEP)
    params = RegionParams.with_defaults(
        WIN, min_region=WIN, method="hmm", exclude_siblings=exclude
    )
    return call_regions(result, analyze(result), WIN, params)


def test_excludes_sibling_and_recovers_recombination(tmp_path) -> None:
    lane_a, lane_b = _lineages()
    query = _recombinant_query(lane_a, lane_b)
    rng = random.Random(2)
    records = {
        "query": query,
        "sibling": _mutate(query, 0.005, rng),  # same mosaic, ~99% everywhere
        "parent_A": _mutate(lane_a, 0.10, rng),  # ~95% in region 1, diverges in region 2
        "parent_B": _mutate(lane_b, 0.10, rng),  # ~95% in region 2, diverges in region 1
    }

    # Without exclusion the sibling wins every window -> masked (no region).
    regions_off, major_off, dropped_off = _result_and_regions(records, tmp_path, exclude=False)
    assert major_off == "sibling"
    assert regions_off == []
    assert dropped_off == []

    # With exclusion the sibling is set aside -> the A/B mosaic is recovered.
    regions_on, major_on, dropped_on = _result_and_regions(records, tmp_path, exclude=True)
    assert [ev.label for ev in dropped_on] == ["sibling"]
    assert major_on == "parent_A"
    assert any(r.minor_parent == "parent_B" for r in regions_on)


def test_partial_coverage_parent_is_not_flagged_as_sibling(tmp_path) -> None:
    # A region-only reference (gaps outside region 2) leads its own windows but does
    # not cover the whole query, so it must NOT be taken for a whole-genome sibling.
    lane_a, lane_b = _lineages()
    query = _recombinant_query(lane_a, lane_b)
    rng = random.Random(3)
    env_only = "-" * R1 + _mutate(lane_b[R1:], 0.10, rng)  # a region-2 parent, covers only it
    records = {
        "query": query,
        "sibling": _mutate(query, 0.005, rng),
        "parent_A": _mutate(lane_a, 0.10, rng),
        "env_only": env_only,
    }
    _, major, dropped = _result_and_regions(records, tmp_path, exclude=True)
    labels = {ev.label for ev in dropped}
    assert "sibling" in labels  # the whole-genome twin is dropped
    assert "env_only" not in labels  # the partial-coverage parent is kept


def test_no_sibling_present_nothing_dropped_region_still_called(tmp_path) -> None:
    lane_a, lane_b = _lineages()
    query = _recombinant_query(lane_a, lane_b)
    rng = random.Random(4)
    records = {
        "query": query,
        "parent_A": _mutate(lane_a, 0.10, rng),
        "parent_B": _mutate(lane_b, 0.10, rng),
        "outgroup": _mutate(lane_a, 0.30, rng),  # distant, beaten everywhere
    }
    regions, major, dropped = _result_and_regions(records, tmp_path, exclude=True)
    assert dropped == []  # no whole-genome explainer to remove
    assert major == "parent_A"
    assert any(r.minor_parent == "parent_B" for r in regions)


def test_non_recombinant_query_gains_no_false_region(tmp_path) -> None:
    # Query is pure lineage A; a sibling and a parent both explain it. Removing the
    # sibling must not invent a region (the significance gate reverts).
    lane_a, _ = _lineages()
    rng = random.Random(5)
    records = {
        "query": lane_a,
        "sibling": _mutate(lane_a, 0.01, rng),
        "parent_A": _mutate(lane_a, 0.10, rng),
        "distractor": _mutate(lane_a, 0.30, rng),
    }
    regions_on, _, _ = _result_and_regions(records, tmp_path, exclude=True)
    regions_off, _, _ = _result_and_regions(records, tmp_path, exclude=False)
    assert regions_on == []
    assert regions_off == []
