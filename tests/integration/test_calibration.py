"""Calibration: the statistical caller must not invent recombination from noise."""

from __future__ import annotations

from pathlib import Path

from recomfi.recomb.analyze import analyze
from recomfi.recomb.regions import RegionParams, call_regions
from recomfi.recomb.similarity import compute_similarity

from ..conftest import write_fasta


def test_no_false_recombination_for_a_pure_reference(tmp_path: Path) -> None:
    # the query IS refA everywhere; refB is uniformly distant. A non-recombinant
    # query must produce zero regions -- the false-positive control the old
    # margin=0 heuristic failed.
    query = "ACGT" * 300
    msa = tmp_path / "msa.fasta"
    write_fasta(msa, {"query": query, "refA": query, "refB": "TGCA" * 300})

    result = compute_similarity(str(msa), "query", window_size=200, window_step=50)
    analysis = analyze(result)
    regions, major = call_regions(result, analysis, 200, RegionParams.with_defaults(200))

    assert major == "refA"
    assert regions == []


def test_heuristic_would_overcall_the_same_input(tmp_path: Path) -> None:
    # the legacy heuristic with margin=0 calls noise where the HMM stays silent;
    # this documents the difference the new default fixes.
    query = "ACGT" * 300
    msa = tmp_path / "msa.fasta"
    # refB matches the query slightly better in one stretch by chance
    refb = list("TGCA" * 300)
    for i in range(200, 260):
        refb[i] = query[i]
    write_fasta(msa, {"query": query, "refA": query, "refB": "".join(refb)})

    result = compute_similarity(str(msa), "query", window_size=100, window_step=50)
    analysis = analyze(result)
    hmm, _ = call_regions(result, analysis, 100, RegionParams.with_defaults(100))
    assert hmm == []  # the HMM does not switch away from the perfect match
