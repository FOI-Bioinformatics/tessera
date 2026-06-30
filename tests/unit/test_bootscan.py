"""Bootscan caller: bootstrap support and region calling."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from tessera.recomb.analyze import analyze
from tessera.recomb.bootscan import bootstrap_support
from tessera.recomb.regions import RegionParams, call_regions
from tessera.recomb.similarity import compute_similarity

from ..conftest import recombinant_msa


def test_bootstrap_support_favours_the_closer_reference() -> None:
    # ref 0 matches the query at 18/20 columns, ref 1 at 2/20: support concentrates on 0
    match = np.zeros((2, 20), dtype=bool)
    match[0, :18] = True
    match[1, 18:] = True
    comparable = np.ones((2, 20), dtype=bool)
    sup = bootstrap_support(match, comparable, seed=1)
    assert sup[0] > 0.95 and sup[1] < 0.05
    assert abs(sup.sum() - 1.0) < 1e-9
    # reproducible with a fixed seed
    assert np.array_equal(bootstrap_support(match, comparable, seed=3),
                          bootstrap_support(match, comparable, seed=3))


def test_call_regions_bootscan_finds_the_recombinant(tmp_path: Path) -> None:
    result = compute_similarity(str(recombinant_msa(tmp_path, recombinant=True)),
                                "query", window_size=500, window_step=50)
    regions, major, _ = call_regions(
        result, analyze(result), 500, RegionParams.with_defaults(500, method="bootscan")
    )
    assert major == "A"
    assert len(regions) == 1
    r = regions[0]
    assert r.minor_parent == "B" and r.major_parent == "A"
    assert r.query_start < 2300 and r.query_end > 3700  # overlaps the true insert
    assert r.support is not None and r.support >= 0.70  # bootstrap support
    assert r.methods == ("bootscan",)


def test_call_regions_bootscan_no_false_positive(tmp_path: Path) -> None:
    result = compute_similarity(str(recombinant_msa(tmp_path, recombinant=False)),
                                "query", window_size=500, window_step=50)
    regions, major, _ = call_regions(
        result, analyze(result), 500, RegionParams.with_defaults(500, method="bootscan")
    )
    assert major == "A"
    assert regions == []
