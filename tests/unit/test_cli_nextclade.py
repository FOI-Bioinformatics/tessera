"""CLI: the nextclade seed-source options reach FillParams."""

from __future__ import annotations

from typer.testing import CliRunner

from recomfi.cli import cmd_build_panel, cmd_detect, cmd_fill_references  # noqa: F401
from recomfi.cli.main import app

runner = CliRunner()


def _capture(monkeypatch):
    captured = {}
    import recomfi.discover.iterate as iterate

    def fake_fill(params, logger):
        captured["params"] = params
        return []

    # all three commands import fill_references from discover.iterate at call time
    monkeypatch.setattr(iterate, "fill_references", fake_fill)
    return captured


def test_detect_nextclade_flag_sets_seed_source(monkeypatch, tmp_path):
    captured = _capture(monkeypatch)
    q = tmp_path / "q.fasta"
    q.write_text(">q\nACGT\n")
    result = runner.invoke(app, ["detect", "-q", str(q), "-o", str(tmp_path / "o"),
                                 "--nextclade"])
    assert result.exit_code == 0, result.output
    assert captured["params"].seed_source == "nextclade"


def test_detect_nextclade_dataset_override(monkeypatch, tmp_path):
    captured = _capture(monkeypatch)
    q = tmp_path / "q.fasta"
    q.write_text(">q\nACGT\n")
    result = runner.invoke(app, ["detect", "-q", str(q), "-o", str(tmp_path / "o"),
                                 "--nextclade-dataset", "nextstrain/sars-cov-2/XBB"])
    assert result.exit_code == 0, result.output
    assert captured["params"].seed_source == "nextclade"
    assert captured["params"].nextclade_dataset == "nextstrain/sars-cov-2/XBB"


def test_fill_references_accepts_nextclade_seed_source(monkeypatch, tmp_path):
    captured = _capture(monkeypatch)
    q = tmp_path / "q.fasta"
    q.write_text(">q\nACGT\n")
    result = runner.invoke(app, ["fill-references", "-q", str(q), "-o", str(tmp_path / "o"),
                                 "--seed-source", "nextclade"])
    assert result.exit_code == 0, result.output
    assert captured["params"].seed_source == "nextclade"
