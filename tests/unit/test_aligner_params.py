"""Aligner adapters forward tuning options into the tool command line."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from recomfi.aligners import progressivemauve as pm
from recomfi.aligners import sibeliaz as sz
from recomfi.aligners.base import AlignParams
from recomfi.core.errors import UserInputError

_LOG = logging.getLogger("recomfi.test")


def _genomes(tmp_path: Path, n: int) -> list[Path]:
    paths = []
    for i in range(n):
        p = tmp_path / f"g{i}.fasta"
        p.write_text(f">g{i}\nACGT\n")
        paths.append(p)
    return paths


def test_alignparams_flag_parses_boolean_strings() -> None:
    assert AlignParams(extra={"single": "true"}).flag("single") is True
    assert AlignParams(extra={"single": "1"}).flag("single") is True
    # The footgun: a non-empty "false" string must not read as truthy.
    assert AlignParams(extra={"single": "false"}).flag("single") is False
    assert AlignParams(extra={"single": "0"}).flag("single") is False
    assert AlignParams(extra={}).flag("single") is False
    assert AlignParams(extra={}).flag("single", default=True) is True


def test_alignparams_flag_rejects_unrecognised_value() -> None:
    with pytest.raises(UserInputError):
        AlignParams(extra={"single": "ture"}).flag("single")


def test_progressivemauve_single_false_does_not_serialise(monkeypatch, tmp_path: Path) -> None:
    captured: dict[str, int] = {}

    def fake_parallel_map(fn, items, workers, **kw):
        captured["workers"] = workers
        return []  # skip the real alignment work

    monkeypatch.setattr(pm, "parallel_map", fake_parallel_map)
    monkeypatch.setattr(pm, "_concatenate", lambda *a, **k: None)

    aligner = pm.ProgressiveMauveAligner()
    params = AlignParams(threads=4, extra={"single": "false"})
    aligner.align(_genomes(tmp_path, 3), None, tmp_path / "out", params, _LOG)
    assert captured["workers"] == 4  # single=false keeps parallelism

    params = AlignParams(threads=4, extra={"single": "true"})
    aligner.align(_genomes(tmp_path, 3), None, tmp_path / "out", params, _LOG)
    assert captured["workers"] == 1  # single=true serialises


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
