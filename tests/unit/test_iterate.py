"""Iterative fill-references loop: rounds, convergence, and the run summary."""

from __future__ import annotations

from pathlib import Path

from recomfi.discover import iterate
from recomfi.discover.blast import Hit
from recomfi.discover.iterate import FillParams, fill_references
from recomfi.discover.run import Candidate
from recomfi.recomb.coverage import CoverageGap


class _StubResult:
    similarities = {"refA": [1.0]}


def _gap(mean_best: float) -> CoverageGap:
    return CoverageGap(
        msa_start=100, msa_end=600, query_start=100, query_end=600,
        length_bp=500, n_windows=5, best_label="refA", mean_best=mean_best, kind="divergent",
    )


def _setup(tmp_path: Path) -> tuple[Path, Path, Path]:
    coll = tmp_path / "coll"
    coll.mkdir()
    (coll / "refA.fasta").write_text(">refA\nACGT\n")
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "ACGT" * 100 + "\n")
    return query, coll, tmp_path / "out"


def _common_mocks(monkeypatch, coverage_returns):
    monkeypatch.setattr(iterate, "efetch_available", lambda: True)
    monkeypatch.setattr(iterate, "build_msa", lambda p, logger: p.output.write_text(">q\nACGT\n"))
    monkeypatch.setattr(iterate, "compute_similarity", lambda *a, **k: _StubResult())
    monkeypatch.setattr(iterate, "read_fasta", lambda path: [("q", "ACGT" * 100)])
    monkeypatch.setattr(iterate, "run_recomb", lambda p, logger: None)
    returns = list(coverage_returns)
    monkeypatch.setattr(iterate, "call_coverage_gaps", lambda *a, **k: returns.pop(0))


def test_loop_converges_when_gaps_close(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    # round 1 has a gap; after adding a reference, round 2 has none.
    _common_mocks(monkeypatch, [([_gap(0.84)], 0.94), ([], 0.95)])
    monkeypatch.setattr(
        iterate, "collect_candidates",
        lambda *a, **k: [Candidate(_gap(0.84), Hit("NEW1", "new", 90.0, 95.0, 1e-9), False)],
    )

    def fake_download(cands, dest, logger):
        (dest / "NEW1.fasta").write_text(">NEW1\nACGT\n")
        return cands

    monkeypatch.setattr(iterate, "_download", fake_download)

    trace = fill_references(FillParams(query=query, collection=coll, output=out), logger)

    assert [r.round for r in trace] == [1, 2]
    assert trace[0].added == ["NEW1"]
    assert trace[1].n_gaps == 0  # converged
    assert (out / "collection" / "NEW1.fasta").exists()  # grown copy, original untouched
    assert not (coll / "NEW1.fasta").exists()
    summary = (out / "fill_summary.tsv").read_text().splitlines()
    assert summary[0].split("\t")[0] == "round"
    assert len(summary) == 3  # header + 2 rounds


def test_query_own_accession_is_auto_excluded(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    # query FASTA header carries its own accession
    query.write_text(">MG572182.1 Norovirus GII\n" + "ACGT" * 100 + "\n")
    _common_mocks(monkeypatch, [([_gap(0.84)], 0.94), ([], 0.95)])

    # the query FASTA reports its accession header; the MSA reports the query label
    def fake_read(path):
        if str(path).endswith("q.fasta"):
            return [("MG572182.1", "ACGT" * 100)]
        return [("q", "ACGT" * 100)]

    monkeypatch.setattr(iterate, "read_fasta", fake_read)

    seen: dict[str, set] = {}

    def capture(targets, query_row, existing, *, exclude, **kw):
        seen["exclude"] = exclude
        return []  # no candidates -> loop stops after round 1

    monkeypatch.setattr(iterate, "collect_candidates", capture)
    monkeypatch.setattr(iterate, "_download", lambda c, d, logger: [])

    fill_references(FillParams(query=query, collection=coll, output=out), logger)
    assert "MG572182" in seen["exclude"]  # auto-excluded from its header


def test_loop_stops_when_coverage_stalls(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    # the best reference stays at 0.84 both rounds -> no improvement -> stop
    _common_mocks(monkeypatch, [([_gap(0.84)], 0.94), ([_gap(0.84)], 0.94)])
    monkeypatch.setattr(
        iterate, "collect_candidates",
        lambda *a, **k: [Candidate(_gap(0.84), Hit("NEW1", "new", 90.0, 95.0, 1e-9), False)],
    )

    def fake_download(cands, dest, logger):
        (dest / "NEW1.fasta").write_text(">NEW1\nACGT\n")
        return cands

    monkeypatch.setattr(iterate, "_download", fake_download)

    trace = fill_references(
        FillParams(query=query, collection=coll, output=out, min_improvement=0.01), logger
    )
    assert [r.round for r in trace] == [1, 2]
    assert trace[1].added == []  # round 2 stalled, no second download
