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
