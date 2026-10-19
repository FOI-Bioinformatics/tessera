"""Lineage clustering: near-duplicates merge, distinct parents stay separate."""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from tessera.recomb.analyze import analyze
from tessera.recomb.clusters import (
    all_singletons,
    cluster_references,
    clustered_view,
    median_window_step,
)
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

    # Pooling two members averages their per-window counts (rounded to whole sites), so
    # the pooled lineage carries one genome's worth of evidence; the similarity is the
    # pooled ratio.
    view, sizes = clustered_view(result, [["sib0", "sib1"], ["pa0"]])
    rep = next(r for r, n in sizes.items() if n == 2)
    for i in range(len(result.positions)):
        num = result.numerators["sib0"][i] + result.numerators["sib1"][i]
        den = result.denominators["sib0"][i] + result.denominators["sib1"][i]
        assert view.numerators[rep][i] == round(num / 2)
        assert view.denominators[rep][i] == round(den / 2)
        assert view.similarities[rep][i] == pytest.approx(num / den)

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


# --- window spacing --------------------------------------------------------

def test_median_window_step_is_robust_to_uneven_spacing():
    """Window spacing must not be read off a single sampled gap.

    Under base-pair windowing every window is `window_step` apart, so the first gap
    is the spacing. Under *informative-site* windowing -- which engages automatically
    below ~8% divergence, exactly where the low-divergence pathogens live -- positions
    are midpoints of informative-site windows and are not uniformly spaced. Taking
    positions[1] - positions[0] then generalises one local gap to the whole genome.
    """
    uniform = _Positions([0, 100, 200, 300, 400])
    assert median_window_step(uniform) == 100

    # dense at the start, sparse later: the first gap is unrepresentative
    uneven = _Positions([0, 5, 10, 20, 400, 900, 1600])
    assert uneven.positions[1] - uneven.positions[0] == 5   # what the old code used
    assert median_window_step(uneven) == 195                # the honest central value


def test_median_window_step_degenerate_inputs():
    assert median_window_step(_Positions([])) == 1
    assert median_window_step(_Positions([42])) == 1
    assert median_window_step(_Positions([7, 7, 7])) == 1   # never returns 0


class _Positions:
    """Minimal stand-in: median_window_step reads only `positions`."""

    def __init__(self, positions):
        self.positions = positions


def test_a_pooled_cluster_carries_one_genomes_evidence(tmp_path) -> None:
    # Four exact copies of one genome pooled into a lineage must emit as that genome
    # does. Summed, the copies carried four times its evidence and every singleton state
    # was charged log(0.5) for each site of the difference.
    from tessera.recomb.hmm import _emission_matrix

    lane_a, lane_b = _lanes()
    query = lane_a[:R1] + lane_b[R1:]
    rng = random.Random(12)
    a = _mutate(lane_a, 0.02, rng)
    records = {"query": query, "A": a, "A1": a, "A2": a, "A3": a,
               "D": _mutate(lane_b, 0.02, rng)}
    result = _result(records, tmp_path)
    view, sizes = clustered_view(result, [["A", "A1", "A2", "A3"], ["D"]])
    assert sizes == {"A": 4, "D": 1}
    assert view.numerators["A"] == result.numerators["A"]
    assert view.denominators["A"] == result.denominators["A"]
    pooled = _emission_matrix(view, ["A", "D"], 0.95)
    single = _emission_matrix(result, ["A", "D"], 0.95)
    assert pooled == pytest.approx(single)


# --- a cluster covering the whole panel is not shared ----------------------

def _whole_panel_pool(tmp_path: Path, *, refs: int, rate: float, donor: str,
                      tract: tuple[int, int], length: int = 12_000,
                      seed: int = 15) -> Path:
    """A panel so uniform that clustering pools every reference into one lineage.

    Under base-pair windowing every pair stays above the 98.5 % floor in every window,
    which is the case documented in docs/detection-methods.md: it leaves the HMM a
    single state. The query carries a real tract from ``donor``.
    """
    rng = random.Random(seed)
    anc = [rng.choice(_BASES) for _ in range(length)]

    def mut(seq: list[str], r: float) -> list[str]:
        return [rng.choice([b for b in _BASES if b != c]) if rng.random() < r else c
                for c in seq]

    records = {f"ref_{i}": "".join(mut(anc, rate)) for i in range(refs)}
    lo, hi = tract
    query = list(records["ref_0"])
    query[lo:hi] = list(records[donor][lo:hi])
    msa = tmp_path / "pool.fasta"
    write_fasta(msa, {"query": "".join(query), **records})
    return msa


def _pipeline_rows(msa: Path, out: Path, logger, **kwargs) -> list[dict[str, str]]:
    import csv

    from tessera.recomb.run import RecombParams, run_recomb

    run_recomb(
        RecombParams(msa=msa, output=out, query="query", plot_format="png",
                     informative_sites=False, **kwargs),
        logger,
    )
    with (out / "recombination_regions.tsv").open() as fo:
        return list(csv.DictReader(fo, delimiter="\t"))


def test_a_whole_panel_cluster_is_not_shared_with_the_merge(tmp_path, logger) -> None:
    # Pooling the whole panel leaves the HMM one state (reported as "not run"), but it
    # must not reach the merge: there it would make every pair of references the same
    # donor, so two callers that named *different* donors would be reported as one row
    # listing both -- agreement manufactured out of a disagreement, which also passes
    # --min-methods 2. On this panel Bootscan names ref_3 and the site callers ref_6.
    msa = _whole_panel_pool(tmp_path, refs=8, rate=0.001, donor="ref_3",
                            tract=(4000, 8000))
    rows = _pipeline_rows(
        msa, tmp_path / "out", logger,
        methods=("hmm", "3seq", "maxchi", "bootscan", "geneconv"),
    )
    by_donor = {r["minor_parent"]: set(r["methods"].split(",")) for r in rows}
    assert set(by_donor) == {"ref_3", "ref_6"}, "the two donors must stay two rows"
    assert by_donor["ref_3"] == {"bootscan"}
    assert by_donor["ref_6"] == {"3seq", "maxchi", "geneconv"}


def test_a_whole_panel_cluster_does_not_silence_bootscan(tmp_path, logger) -> None:
    # Every candidate in the major parent's own group: no group can beat the major, so
    # Bootscan could never claim a window. It must keep competing genomes here.
    msa = _whole_panel_pool(tmp_path, refs=5, rate=0.0025, donor="ref_3",
                            tract=(4000, 8000))
    rows = _pipeline_rows(msa, tmp_path / "out", logger,
                          methods=("3seq", "maxchi", "bootscan"))
    assert rows
    assert any("bootscan" in r["methods"].split(",") for r in rows)


def test_a_whole_panel_cluster_still_leaves_the_hmm_one_state(tmp_path, logger) -> None:
    # Invariant the other way round: the HMM keeps seeing the whole-panel cluster, so it
    # is still reported as "not run" rather than competing every reference individually.
    import json

    msa = _whole_panel_pool(tmp_path, refs=8, rate=0.001, donor="ref_3",
                            tract=(4000, 8000))
    out = tmp_path / "out"
    _pipeline_rows(msa, out, logger,
                   methods=("hmm", "3seq", "maxchi", "bootscan", "geneconv"))
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["callers not run"].startswith(
        "hmm (lineage clustering pooled every reference"
    )
