"""find-references: coverage-gap subsequence extraction, candidate reporting, exclusion."""

from __future__ import annotations

from pathlib import Path

import pytest

from recomfi.discover import run as discover_run
from recomfi.discover.blast import Hit, _organism
from recomfi.discover.run import FindRefParams, find_references

from ..conftest import write_fasta

# Three 100 bp segments. The query matches refA across seg1+seg2; in seg3 it
# differs from BOTH references -> a divergent coverage gap with no good reference.
_Q1 = "ACGTACGTAC" * 10
_Q2 = "TGCATGCATG" * 10
_Q3 = "AAAAACCCCC" * 10
_A3 = "GGGGGTTTTT" * 10
_B2 = "CCCCCGGGGG" * 10
_B3 = "TTTTTGGGGG" * 10

QUERY = _Q1 + _Q2 + _Q3
REF_A = _Q1 + _Q2 + _A3
REF_B = _Q1 + _B2 + _B3


def _msa(tmp_path: Path) -> Path:
    path = tmp_path / "msa.fasta"
    write_fasta(path, {"q": QUERY, "refA": REF_A, "refB": REF_B})
    return path


def test_find_references_blasts_the_gap_and_marks_existing(monkeypatch, tmp_path, logger):
    captured: dict[str, str] = {}

    def fake_blast(seq, *, max_hits, logger, email=None):
        captured["seq"] = seq
        return [
            Hit("refA", "cowpox-like [refA]", 95.0, 90.0, 1e-50),  # already present
            Hit("NEW123", "novel donor virus", 88.0, 92.0, 1e-40),  # new
        ]

    monkeypatch.setattr(discover_run, "blast_subsequence", fake_blast)

    params = FindRefParams(
        msa=_msa(tmp_path), query="q", output=tmp_path / "out",
        window_size=60, window_step=30, top_gaps=1, max_hits=2,
    )
    candidates = find_references(params, logger)

    # the subsequence BLASTed is the divergent seg3, de-gapped
    assert _Q3[:50] in captured["seq"]
    accessions = {c.hit.accession for c in candidates}
    assert accessions == {"refA", "NEW123"}
    by_acc = {c.hit.accession: c for c in candidates}
    assert by_acc["refA"].in_collection is True  # refA is an existing label
    assert by_acc["NEW123"].in_collection is False

    tsv = (tmp_path / "out" / "candidate_references.tsv").read_text().splitlines()
    assert tsv[0].split("\t")[3] == "accession"
    assert any("NEW123" in line for line in tsv[1:])


def test_self_hit_is_auto_skipped(monkeypatch, tmp_path, logger):
    # A near-identical, near-full-coverage hit is the query's own record -> skipped.
    def fake_blast(seq, *, max_hits, logger, email=None):
        return [
            Hit("SELF999", "the query itself", 99.8, 100.0, 0.0),  # self-hit
            Hit("NEW123", "novel donor", 90.0, 95.0, 1e-40),       # genuine candidate
        ]

    monkeypatch.setattr(discover_run, "blast_subsequence", fake_blast)
    params = FindRefParams(
        msa=_msa(tmp_path), query="q", output=tmp_path / "out",
        window_size=60, window_step=30, top_gaps=1,
    )
    candidates = find_references(params, logger)
    assert {c.hit.accession for c in candidates} == {"NEW123"}  # self-hit dropped

    # ...unless explicitly kept
    kept = find_references(
        FindRefParams(
            msa=_msa(tmp_path), query="q", output=tmp_path / "out2",
            window_size=60, window_step=30, top_gaps=1, keep_self_hits=True,
        ),
        logger,
    )
    assert "SELF999" in {c.hit.accession for c in kept}


def test_exclude_drops_accession_version_insensitively(monkeypatch, tmp_path, logger):
    def fake_blast(seq, *, max_hits, logger, email=None):
        return [
            Hit("U54771.1", "query own strain", 95.0, 90.0, 0.0),
            Hit("NEW123", "novel donor", 90.0, 95.0, 1e-40),
        ]

    monkeypatch.setattr(discover_run, "blast_subsequence", fake_blast)
    params = FindRefParams(
        msa=_msa(tmp_path), query="q", output=tmp_path / "out",
        window_size=60, window_step=30, top_gaps=1,
        exclude=("U54771",),  # bare accession matches the versioned hit
    )
    candidates = find_references(params, logger)
    assert {c.hit.accession for c in candidates} == {"NEW123"}


def test_no_gap_means_no_blast(monkeypatch, tmp_path, logger):
    # query identical to a reference everywhere -> no coverage gap -> no BLAST
    path = tmp_path / "msa.fasta"
    write_fasta(path, {"q": QUERY, "refA": QUERY, "refB": QUERY})

    def boom(*a, **k):
        raise AssertionError("BLAST should not be called when there are no gaps")

    monkeypatch.setattr(discover_run, "blast_subsequence", boom)
    params = FindRefParams(
        msa=path, query="q", output=tmp_path / "out", window_size=60, window_step=30,
    )
    assert find_references(params, logger) == []
    assert (tmp_path / "out" / "candidate_references.tsv").exists()


def test_download_writes_a_manifest_of_added_references(monkeypatch, tmp_path, logger):
    def fake_blast(seq, *, max_hits, logger, email=None):
        return [Hit("NEW123", "novel donor virus", 90.0, 95.0, 1e-40)]

    def fake_efetch(accession, collection_dir, logger):
        path = collection_dir / f"{accession}.fasta"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f">{accession}\nACGT\n")
        return path

    monkeypatch.setattr(discover_run, "blast_subsequence", fake_blast)
    monkeypatch.setattr(discover_run, "efetch_available", lambda: True)
    monkeypatch.setattr(discover_run, "efetch_fasta", fake_efetch)

    coll = tmp_path / "coll"
    params = FindRefParams(
        msa=_msa(tmp_path), query="q", output=tmp_path / "out",
        window_size=60, window_step=30, top_gaps=1, download=coll,
    )
    find_references(params, logger)

    assert (coll / "NEW123.fasta").exists()  # the reference was fetched
    manifest = (tmp_path / "out" / "downloaded_references.tsv").read_text().splitlines()
    assert manifest[0].split("\t")[0] == "accession"
    assert manifest[1].split("\t")[0] == "NEW123"


def test_download_without_efetch_is_a_clear_error(monkeypatch, tmp_path, logger):
    monkeypatch.setattr(discover_run, "efetch_available", lambda: False)
    from recomfi.core.errors import UserInputError

    params = FindRefParams(
        msa=_msa(tmp_path), query="q", output=tmp_path / "out",
        download=tmp_path / "coll",
    )
    with pytest.raises(UserInputError, match="entrez-direct"):
        find_references(params, logger)


def test_organism_extraction():
    assert _organism("Some virus strain X, complete genome [Variola virus]") == "Variola virus"
    assert _organism("Cowpox virus, complete genome") == "Cowpox virus"


def test_tiles_geometry():
    from recomfi.discover.run import _tiles

    assert _tiles(300, 400) == [(0, 300)]  # shorter than a tile -> one search
    tiles = _tiles(800, 400)  # overlapping, half-step, last aligned to the end
    assert tiles == [(0, 400), (200, 600), (400, 800)]


def test_collect_candidates_subtiling_surfaces_per_region_donors(monkeypatch, logger):
    from recomfi.discover.run import collect_candidates
    from recomfi.recomb.coverage import CoverageGap

    # An 800 bp gap: an 'A' flank, a 'C' flank, and (in the middle tile) a mix.
    query_row = "A" * 400 + "C" * 400
    gap = CoverageGap(
        msa_start=0, msa_end=800, query_start=0, query_end=800,
        length_bp=800, n_windows=8, best_label="ref", mean_best=0.9, kind="divergent",
    )

    def fake_blast(seq, *, max_hits, logger, email=None):
        if "C" not in seq:
            donor = Hit("FLANK_A", "a", 90.0, 95.0, 1e-9)
        elif "A" not in seq:
            donor = Hit("FLANK_C", "c", 90.0, 95.0, 1e-9)
        else:
            donor = Hit("CORE", "core donor", 88.0, 95.0, 1e-9)
        return [donor, Hit("SHARED", "in every tile", 92.0, 95.0, 1e-9)]

    monkeypatch.setattr(discover_run, "blast_subsequence", fake_blast)
    common = dict(max_hits=5, email=None, exclude=set(), keep_self_hits=True, logger=logger)

    whole = collect_candidates([gap], query_row, set(), subtile=0, **common)
    assert {c.hit.accession for c in whole} == {"CORE", "SHARED"}  # one diluted search

    tiled = collect_candidates([gap], query_row, set(), subtile=400, **common)
    # each region's donor surfaces; SHARED appears in all tiles but is kept once
    assert {c.hit.accession for c in tiled} == {"FLANK_A", "CORE", "FLANK_C", "SHARED"}
    assert sum(c.hit.accession == "SHARED" for c in tiled) == 1
