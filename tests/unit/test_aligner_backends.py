"""MAFFT and minimap2 backends: converters and command construction."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from recomfi.aligners import mafft as mafft_mod
from recomfi.aligners import minimap2 as mm_mod
from recomfi.aligners.base import AlignParams
from recomfi.converters.mafft_merge import merge_added_fragments
from recomfi.converters.sam_to_fasta import project_sam_to_row

_LOG = logging.getLogger("recomfi.test")


# --- MAFFT keeplength merge ------------------------------------------------
def test_merge_added_fragments_picks_non_gap(tmp_path: Path) -> None:
    # reference first, then two contig rows that each cover part of the reference
    aligned = tmp_path / "a.fasta"
    aligned.write_text(
        ">ref\nACGTACGT\n"
        ">contig1\nAC------\n"
        ">contig2\n----ACGT\n"
    )
    ref_row, merged = merge_added_fragments(aligned)
    assert ref_row == "ACGTACGT"
    assert merged == "AC--ACGT"  # contig1 covers 0-1, contig2 covers 4-7


# --- SAM -> reference projection -------------------------------------------
def _sam(*records: str) -> str:
    header = "@HD\tVN:1.6\n@SQ\tSN:ref\tLN:10\n"
    return header + "".join(r + "\n" for r in records)


def _record(qname, flag, rname, pos, cigar, seq) -> str:
    return f"{qname}\t{flag}\t{rname}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t*"


def test_project_sam_match_insertion_deletion(tmp_path: Path) -> None:
    sam = tmp_path / "a.sam"
    # 2M1I2M: AA at ref0-1, C inserted (dropped), GT at ref2-3
    sam.write_text(_sam(_record("q", 0, "ref", 1, "2M1I2M", "AACGT")))
    row = project_sam_to_row(str(sam), {"ref": 0}, 10)
    assert row == "AAGT------"

    sam.write_text(_sam(_record("q", 0, "ref", 1, "2M1D2M", "AAGT")))
    row = project_sam_to_row(str(sam), {"ref": 0}, 10)
    assert row == "AA-GT-----"  # 1D leaves a gap at ref2


def test_project_sam_skips_secondary_and_unmapped(tmp_path: Path) -> None:
    sam = tmp_path / "a.sam"
    sam.write_text(_sam(
        _record("q", 0x4, "ref", 1, "4M", "AAAA"),     # unmapped: skipped
        _record("q", 0x100, "ref", 1, "4M", "TTTT"),   # secondary: skipped
        _record("q", 0, "ref", 5, "3M", "GGG"),        # primary: placed at 4-6
    ))
    row = project_sam_to_row(str(sam), {"ref": 0}, 10)
    assert row == "----GGG---"


def test_project_sam_multicontig_reference_offsets(tmp_path: Path) -> None:
    sam = tmp_path / "a.sam"
    sam.write_text(
        "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:5\n@SQ\tSN:c2\tLN:5\n"
        + _record("q", 0, "c2", 1, "3M", "TTT") + "\n"
    )
    row = project_sam_to_row(str(sam), {"c1": 0, "c2": 5}, 10)
    assert row == "-----TTT--"  # c2 starts at global offset 5


# --- command construction (monkeypatched run_tool) -------------------------
def _two_genomes(tmp_path: Path) -> list[Path]:
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref\nACGTACGTACGT\n")
    qry = tmp_path / "qry.fasta"
    qry.write_text(">qry\nACGTACGTACGT\n")
    return [ref, qry]


def test_mafft_forwards_keeplength_and_tuning(monkeypatch, tmp_path: Path) -> None:
    captured: dict[str, list[str]] = {}

    def fake_run(caps, cmd, **kw):
        captured["cmd"] = [str(c) for c in cmd]
        raise RuntimeError("stop after capture")

    monkeypatch.setattr(mafft_mod, "run_tool", fake_run)
    ref, qry = _two_genomes(tmp_path)
    params = AlignParams(threads=2, extra={"maxiterate": 1000})
    with pytest.raises(RuntimeError):
        mafft_mod.MafftAligner().align([ref, qry], ref, tmp_path / "out", params, _LOG)

    cmd = captured["cmd"]
    assert "--keeplength" in cmd
    assert "--addfragments" in cmd
    assert cmd[cmd.index("--thread") + 1] == "2"
    assert cmd[cmd.index("--maxiterate") + 1] == "1000"


def test_minimap2_forwards_preset_and_threads(monkeypatch, tmp_path: Path) -> None:
    captured: dict[str, list[str]] = {}

    def fake_run(caps, cmd, **kw):
        captured["cmd"] = [str(c) for c in cmd]
        raise RuntimeError("stop after capture")

    monkeypatch.setattr(mm_mod, "run_tool", fake_run)
    ref, qry = _two_genomes(tmp_path)
    params = AlignParams(threads=3, extra={"preset": "asm10"})
    with pytest.raises(RuntimeError):
        mm_mod.Minimap2Aligner().align([ref, qry], ref, tmp_path / "out", params, _LOG)

    cmd = captured["cmd"]
    assert "-a" in cmd
    assert cmd[cmd.index("-x") + 1] == "asm10"
    assert cmd[cmd.index("-t") + 1] == "3"


def _read_msa(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            out[name] = ""
        elif name is not None:
            out[name] += line
    return out


def test_mafft_assembles_reference_anchored_msa(monkeypatch, tmp_path: Path) -> None:
    def fake_run(caps, cmd, **kw):
        # emulate `mafft --keeplength`: reference row + one aligned contig row
        kw["stdout_path"].write_text(">ref\nACGTACGTACGT\n>qry_contig\nACGTACGTACGT\n")
        return ""

    monkeypatch.setattr(mafft_mod, "run_tool", fake_run)
    ref, qry = _two_genomes(tmp_path)
    result = mafft_mod.MafftAligner().align([ref, qry], ref, tmp_path / "out",
                                            AlignParams(threads=1), _LOG)
    msa = _read_msa(result.msa_fasta)
    assert set(msa) == {"ref", "qry"}
    assert len(msa["ref"]) == len(msa["qry"]) == 12


def test_minimap2_assembles_reference_anchored_msa(monkeypatch, tmp_path: Path) -> None:
    def fake_run(caps, cmd, **kw):
        kw["stdout_path"].write_text(
            "@SQ\tSN:ref\tLN:12\n"
            "q\t0\tref\t1\t60\t12M\t*\t0\t0\tACGTACGTACGT\t*\n"
        )
        return ""

    monkeypatch.setattr(mm_mod, "run_tool", fake_run)
    ref, qry = _two_genomes(tmp_path)
    result = mm_mod.Minimap2Aligner().align([ref, qry], ref, tmp_path / "out",
                                            AlignParams(threads=1), _LOG)
    msa = _read_msa(result.msa_fasta)
    assert set(msa) == {"ref", "qry"}
    assert msa["ref"] == "ACGTACGTACGT"
    assert msa["qry"] == "ACGTACGTACGT"
