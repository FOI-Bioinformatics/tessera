"""Inputs the CLI must reject before any work starts.

Each case used to be accepted: the option was silently ignored, or the run failed
later under "Unexpected error". None of these tests may reach the network -- the
panel-building entry points are replaced with a function that fails the test if the
command gets that far.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

from tessera.cli.main import app

runner = CliRunner()

EXAMPLE = "example_data/divergent.msa.fasta"


def _must_not_run(*args, **kwargs):
    raise AssertionError("the command reached the pipeline; it should have been rejected")


@pytest.fixture
def no_pipeline(monkeypatch):
    """Fail the test if a command gets past validation into panel building or typing."""
    monkeypatch.setattr("tessera.discover.iterate.fill_references", _must_not_run)
    monkeypatch.setattr("tessera.discover.lineage_assign.assign_lineages", _must_not_run)
    monkeypatch.setattr("tessera.discover.run.find_references", _must_not_run)


def _query(tmp_path: Path) -> Path:
    path = tmp_path / "query.fasta"
    path.write_text(">query\nACGTACGTACGT\n")
    return path


def _collection(tmp_path: Path) -> Path:
    directory = tmp_path / "collection"
    directory.mkdir()
    (directory / "ref.fasta").write_text(">ref\nACGTACGTACGT\n")
    return directory


def _rejected(result, needle: str) -> None:
    assert result.exit_code == 1, result.output
    assert "Unexpected error" not in result.output
    assert needle in " ".join(result.output.split())


# --- --lineage-map must exist ---------------------------------------------

def test_recomb_rejects_a_missing_lineage_map(tmp_path: Path) -> None:
    """A mistyped path used to be ignored: the run exited 0 with an untyped report, and
    the barcode caller and donor re-attribution silently did nothing."""
    result = runner.invoke(app, [
        "recomb", "--msa", EXAMPLE, "--query", "query",
        "--output", str(tmp_path / "out"), "--lineage-map", str(tmp_path / "nope.tsv"),
    ])
    _rejected(result, "--lineage-map file not found")
    assert not (tmp_path / "out" / "report.html").exists()


def test_recomb_rejects_a_lineage_map_that_is_a_directory(tmp_path: Path) -> None:
    result = runner.invoke(app, [
        "recomb", "--msa", EXAMPLE, "--query", "query",
        "--output", str(tmp_path / "out"), "--lineage-map", str(tmp_path),
    ])
    _rejected(result, "is a directory, not a file")


def test_type_lineages_rejects_a_missing_lineage_map(tmp_path: Path, no_pipeline) -> None:
    result = runner.invoke(app, [
        "type-lineages", "--collection", str(_collection(tmp_path)),
        "--output", str(tmp_path / "out"), "--lineage-map", str(tmp_path / "nope.tsv"),
    ])
    _rejected(result, "--lineage-map file not found")


@pytest.mark.parametrize("command", ["detect", "fill-references", "build-panel"])
def test_panel_commands_reject_a_missing_lineage_map(
    tmp_path: Path, no_pipeline, command: str
) -> None:
    result = runner.invoke(app, [
        command, "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        "--lineage-map", str(tmp_path / "nope.tsv"),
    ])
    _rejected(result, "--lineage-map file not found")
