"""MaxChi chi-square triplet caller: statistic, scan-aware p-value, and caller."""

from __future__ import annotations

import random
from pathlib import Path

import numpy as np

from tessera.recomb.analyze import analyze
from tessera.recomb.maxchi import (
    _chi2_2x2,
    maxchi_pvalue,
    maxchi_statistic,
)
from tessera.recomb.regions import RegionParams, call_regions
from tessera.recomb.similarity import compute_similarity
from tessera.recomb.threeseq import max_descent

from ..conftest import write_fasta


def test_chi2_2x2_matches_closed_form() -> None:
    # 2x2 = [[10,0],[0,10]] -> perfect association -> chi2 = n = 20
    assert _chi2_2x2(10, 0, 0, 10) == 20.0
    # independent table -> 0
    assert _chi2_2x2(5, 5, 5, 5) == 0.0
    # a zero margin is undefined -> 0
    assert _chi2_2x2(0, 0, 4, 6) == 0.0


def test_maxchi_statistic_high_for_clean_tract() -> None:
    # +1 background with a pure -1 tract in the middle: strong proportion shift
    steps = np.array([1] * 20 + [-1] * 20 + [1] * 20)
    d = max_descent(steps)
    chi2 = maxchi_statistic(steps, d.start_site, d.end_site)
    assert chi2 > 30  # near-perfect 2x2 association
    # a walk with no tract (alternating) -> small statistic over its (shallow) descent
    flat = np.array([1, -1] * 30)
    assert maxchi_statistic(flat, max_descent(flat).start_site, max_descent(flat).end_site) < chi2


def test_maxchi_pvalue_significant_vs_calibrated() -> None:
    rng = np.random.default_rng(0)
    # clear mosaic -> tiny p
    steps = np.array([1] * 30 + [-1] * 30 + [1] * 30)
    chi2 = maxchi_statistic(steps, max_descent(steps).start_site, max_descent(steps).end_site)
    assert maxchi_pvalue(steps, chi2, seed=1) < 0.01
    # random interspersion -> not significant (well above alpha)
    noise = rng.choice([1, -1], size=80)
    nchi = maxchi_statistic(noise, max_descent(noise).start_site, max_descent(noise).end_site)
    assert maxchi_pvalue(noise, nchi, seed=1) > 0.2
    # reproducible with a fixed seed
    assert maxchi_pvalue(steps, chi2, seed=7) == maxchi_pvalue(steps, chi2, seed=7)


def _recombinant_msa(tmp_path: Path, *, recombinant: bool) -> Path:
    rng = random.Random(7)
    base = "".join(rng.choice("ACGT") for _ in range(6000))

    def mut(seq: str, frac: float) -> str:
        chars = list(seq)
        for i in range(len(chars)):
            if rng.random() < frac:
                chars[i] = rng.choice("ACGT")
        return "".join(chars)

    a, b, other = mut(base, 0.03), mut(base, 0.03), mut(base, 0.10)
    query = list(a)
    if recombinant:
        query[2000:4000] = list(b[2000:4000])
    return write_fasta(tmp_path / "m.fasta",
                       {"query": "".join(query), "A": a, "B": b, "other": other})


def test_call_regions_maxchi_finds_the_recombinant(tmp_path: Path) -> None:
    result = compute_similarity(str(_recombinant_msa(tmp_path, recombinant=True)),
                                "query", window_size=500, window_step=50)
    regions, major, _ = call_regions(
        result, analyze(result), 500, RegionParams.with_defaults(500, method="maxchi")
    )
    assert major == "A"
    assert len(regions) == 1
    r = regions[0]
    assert r.minor_parent == "B" and r.major_parent == "A"
    assert r.query_start < 2200 and r.query_end > 3800
    assert r.qvalue is not None and r.qvalue < 0.01
    assert r.methods == ("maxchi",)


def test_call_regions_maxchi_no_false_positive(tmp_path: Path) -> None:
    result = compute_similarity(str(_recombinant_msa(tmp_path, recombinant=False)),
                                "query", window_size=500, window_step=50)
    regions, major, _ = call_regions(
        result, analyze(result), 500, RegionParams.with_defaults(500, method="maxchi")
    )
    assert major == "A"
    assert regions == []
