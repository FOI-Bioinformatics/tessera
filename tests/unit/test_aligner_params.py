"""Aligner adapters forward tuning options into the tool command line."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from recomfi.aligners import progressivemauve as pm
from recomfi.aligners import sibeliaz as sz
from recomfi.aligners.base import AlignParams

_LOG = logging.getLogger("recomfi.test")


def _genomes(tmp_path: Path, n: int) -> list[Path]:
    paths = []
    for i in range(n):
        p = tmp_path / f"g{i}.fasta"
        p.write_text(f">g{i}\nACGT\n")
        paths.append(p)
    return paths


def test_progressivemauve_forwards_seed_weight(monkeypatch, tmp_path: Path) -> None:
    captured: dict[str, list[str]] = {}

    def fake_run(caps, cmd, **kw):
        captured["cmd"] = [str(c) for c in cmd]
        raise RuntimeError("stop after capture")

    monkeypatch.setattr(pm, "run_tool", fake_run)

    aligner = pm.ProgressiveMauveAligner()
    params = AlignParams(threads=1, extra={"seed_weight": 11})
    with pytest.raises(RuntimeError):
        aligner.align(_genomes(tmp_path, 3), None, tmp_path / "out", params, _LOG)

    cmd = captured["cmd"]
    assert "--seed-weight" in cmd
    assert cmd[cmd.index("--seed-weight") + 1] == "11"


def test_sibeliaz_forwards_kmer_and_filtermemory(monkeypatch, tmp_path: Path) -> None:
    captured: dict[str, list[str]] = {}

    def fake_run(caps, cmd, **kw):
        captured["cmd"] = [str(c) for c in cmd]
        raise RuntimeError("stop after capture")

    monkeypatch.setattr(sz, "run_tool", fake_run)
    # Avoid the macOS wrapper-patching path so the command stays plain "sibeliaz".
    monkeypatch.setattr(sz, "_sibeliaz_invocation", lambda out_dir, logger: ["sibeliaz"])

    aligner = sz.SibeliazAligner()
    params = AlignParams(threads=4, extra={"kmer": 15, "filtermemory": 12})
    with pytest.raises(RuntimeError):
        aligner.align(_genomes(tmp_path, 3), None, tmp_path / "out", params, _LOG)

    cmd = captured["cmd"]
    assert cmd[cmd.index("-k") + 1] == "15"
    assert cmd[cmd.index("-f") + 1] == "12"
    assert cmd[cmd.index("-t") + 1] == "4"
