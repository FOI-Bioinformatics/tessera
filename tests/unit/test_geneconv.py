"""GENECONV-style clean-fragment caller: run statistic, permutation p-value, and caller."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from tessera.recomb.analyze import analyze
from tessera.recomb.geneconv import longest_minor_run, run_pvalue_perm
from tessera.recomb.regions import RegionParams, call_regions
from tessera.recomb.similarity import compute_similarity

from ..conftest import recombinant_msa


def test_longest_minor_run_finds_the_clean_fragment() -> None:
    # P^10 then Q^8 then P^5: the longest uninterrupted minor-run is the 8 at sites [10, 18).
    steps = np.array([1] * 10 + [-1] * 8 + [1] * 5)
    length, start, end = longest_minor_run(steps)
    assert length == 8
    assert (start, end) == (10, 17)   # inclusive site indices


def test_longest_minor_run_breaks_on_interruption() -> None:
    # -1 -1 +1 -1 -1 -1 : an interrupting +1 breaks the run; the longest clean run is 3.
    steps = np.array([-1, -1, 1, -1, -1, -1])
    length, start, end = longest_minor_run(steps)
    assert length == 3 and (start, end) == (3, 5)


def test_longest_minor_run_none() -> None:
    assert longest_minor_run(np.array([1, 1, 1])) == (0, -1, -1)


def test_run_pvalue_clustered_is_significant_random_is_not() -> None:
    # a single long clean run among few interruptions is unlikely by chance
    clustered = np.array([1] * 40 + [-1] * 40)
    length, _, _ = longest_minor_run(clustered)
    assert run_pvalue_perm(clustered, length, seed=0) < 0.01
    # an alternating series has no long clean run; asking for a short one is not significant
    alternating = np.array([1, -1] * 40)
    assert run_pvalue_perm(alternating, 2, seed=0) > 0.5


def test_run_pvalue_edges() -> None:
    assert run_pvalue_perm(np.array([-1, -1, -1]), 0) == 1.0


def test_call_regions_geneconv_finds_the_recombinant(tmp_path: Path) -> None:
    result = compute_similarity(str(recombinant_msa(tmp_path, recombinant=True)),
                                "query", window_size=500, window_step=50)
    regions, major, dropped = call_regions(
        result, analyze(result), 500, RegionParams.with_defaults(500, method="geneconv")
    )
    assert major == "A"
    assert dropped == []
    assert len(regions) == 1
    r = regions[0]
    assert r.minor_parent == "B" and r.major_parent == "A"
    assert r.query_start < 2200 and r.query_end > 3800   # spans the true insert
    assert r.pvalue is not None and r.pvalue < 0.01
    assert r.methods == ("geneconv",)


def test_call_regions_geneconv_no_false_positive_on_pure_query(tmp_path: Path) -> None:
    result = compute_similarity(str(recombinant_msa(tmp_path, recombinant=False)),
                                "query", window_size=500, window_step=50)
    regions, major, _ = call_regions(
        result, analyze(result), 500, RegionParams.with_defaults(500, method="geneconv")
    )
    assert major == "A"
    assert regions == []
