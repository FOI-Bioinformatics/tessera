"""Unit tests for the per-segment intragenic scan (stubbed; no aligner/binaries)."""

from __future__ import annotations

import logging
from pathlib import Path

from tessera.reassort import scan
from tessera.reassort.scan import SegmentScan, _summarize_regions, scan_segment, window_params

LOG = logging.getLogger("test")


class _DS:
    path = "flu/ha"
    tag = "tag"


def _stub_pool(monkeypatch, tips):
    monkeypatch.setattr(scan, "build_pool",
                        lambda ds, *, cache_dir, logger, per_clade_consensus: tips)
    monkeypatch.setattr(scan, "nextclade_cache", lambda p, t, override=None: Path("/x"))


def _regions_tsv(tmp_path, rows):
    p = tmp_path / "recombination_regions.tsv"
    header = "query_start\tquery_end\tdonor_absent"
    p.write_text(header + "\n" + "\n".join(rows) + ("\n" if rows else ""))
    return p


def test_summarize_counts_present_regions(tmp_path):
    p = _regions_tsv(tmp_path, ["100\t200\tno", "300\t400\tno"])
    assert _summarize_regions(p) == (2, True)


def test_summarize_excludes_donor_absent(tmp_path):
    p = _regions_tsv(tmp_path, ["100\t200\tno", "300\t400\tyes"])
    assert _summarize_regions(p) == (1, True)


def test_summarize_empty_file_is_not_recombinant(tmp_path):
    p = _regions_tsv(tmp_path, [])
    assert _summarize_regions(p) == (0, False)


def test_summarize_missing_file_is_not_recombinant(tmp_path):
    assert _summarize_regions(tmp_path / "nope.tsv") == (0, False)


def test_window_params_scales_to_short_segment():
    window, step, sel = window_params(1700)
    assert 120 <= window <= 500
    assert step == max(20, window // 10)
    assert sel >= window


def test_scan_segment_single_clade_panel_is_skipped(tmp_path, monkeypatch):
    # A one-clade panel has nothing to recombine between -> scanned=False, note explains why.
    one = tmp_path / "pool" / "A_consensus.fasta"
    one.parent.mkdir(parents=True)
    one.write_text(">A_consensus A\nACGT\n")
    _stub_pool(monkeypatch, [one])
    result = scan_segment("HA", "ACGTACGT", _DS(), tmp_path / "out", aligner="mafft",
                          cache_dir=None, logger=LOG)
    assert result == SegmentScan("HA", False, False, 0, "single-clade panel")


def test_scan_segment_runs_recomb_and_summarizes(tmp_path, monkeypatch):
    # Two-clade panel: build_msa/run_recomb are stubbed; run_recomb writes a regions file that
    # _summarize_regions reads back -> scanned + recombinant.
    pool = tmp_path / "pool"
    pool.mkdir()
    tips = []
    for c in ("A", "B"):
        t = pool / f"{c}_consensus.fasta"
        t.write_text(f">{c}_consensus {c}\nACGTACGT\n")
        tips.append(t)
    _stub_pool(monkeypatch, tips)
    monkeypatch.setattr(scan, "build_msa", lambda params, logger: params.output)

    def fake_run_recomb(params, logger):
        (params.output / "recombination_regions.tsv").write_text(
            "query_start\tquery_end\tdonor_absent\n100\t200\tno\n")
        return "bp"
    monkeypatch.setattr(scan, "run_recomb", fake_run_recomb)

    result = scan_segment("HA", "ACGTACGT", _DS(), tmp_path / "out", aligner="mafft",
                          cache_dir=None, logger=LOG)
    assert result.segment == "HA" and result.scanned and result.recombinant
    assert result.n_regions == 1
    assert (tmp_path / "out" / "HA" / "recombination_regions.tsv").exists()


def test_scan_segment_keeps_only_consensus_genomes(tmp_path, monkeypatch):
    # build_pool also returns raw example strains; the scan panel must keep only *_consensus.
    # Here one consensus + one example -> only one consensus -> single-clade panel skip.
    pool = tmp_path / "pool"
    pool.mkdir()
    cons = pool / "A_consensus.fasta"
    cons.write_text(">A_consensus A\nACGT\n")
    example = pool / "A_Texas_2_2021.fasta"
    example.write_text(">A_Texas_2_2021 example\nACGT\n")
    _stub_pool(monkeypatch, [cons, example])
    result = scan_segment("HA", "ACGT", _DS(), tmp_path / "out", aligner="mafft",
                          cache_dir=None, logger=LOG)
    assert result == SegmentScan("HA", False, False, 0, "single-clade panel")


def test_scan_segment_failure_is_non_fatal(tmp_path, monkeypatch):
    pool = tmp_path / "pool"
    pool.mkdir()
    tips = []
    for c in ("A", "B"):
        t = pool / f"{c}_consensus.fasta"
        t.write_text(f">{c}_consensus {c}\nACGT\n")
        tips.append(t)
    _stub_pool(monkeypatch, tips)

    def boom(params, logger):
        raise RuntimeError("aligner blew up")
    monkeypatch.setattr(scan, "build_msa", boom)

    result = scan_segment("HA", "ACGT", _DS(), tmp_path / "out", aligner="mafft",
                          cache_dir=None, logger=LOG)
    assert result.scanned is False
    assert "scan failed" in result.note
