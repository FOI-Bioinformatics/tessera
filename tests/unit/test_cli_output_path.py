"""Every command that writes into an output directory rejects an existing file there.

``recomb`` checked this (audit B11); five other commands crashed with a raw traceback or
an "Unexpected error" from wherever the directory was first created. None of these tests
may reach the pipeline: the entry points are replaced with a function that fails the test.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

from tessera.cli.main import app

runner = CliRunner()


def _must_not_run(*args, **kwargs):
    raise AssertionError("the command reached the pipeline; it should have been rejected")


@pytest.fixture
def no_pipeline(monkeypatch):
    monkeypatch.setattr("tessera.discover.iterate.fill_references", _must_not_run)
    monkeypatch.setattr("tessera.discover.lineage_assign.assign_lineages", _must_not_run)
    monkeypatch.setattr("tessera.discover.run.find_references", _must_not_run)
    monkeypatch.setattr("tessera.discover.panel.curate_collection", _must_not_run)
    monkeypatch.setattr("tessera.cli.cmd_reassort.assign_segments", _must_not_run)


def _args(command: str, tmp_path: Path, occupied: Path) -> list[str]:
    query = tmp_path / "query.fasta"
    query.write_text(">query\nACGTACGTACGT\n>NA\nACGTACGTACGT\n")
    collection = tmp_path / "collection"
    collection.mkdir()
    (collection / "ref.fasta").write_text(">ref\nACGTACGTACGT\n")
    common = ["--output", str(occupied)]
    if command == "type-lineages":
        return [command, "--collection", str(collection), *common]
    if command == "curate-panel":
        return [command, "--query", str(query), "--collection", str(collection), *common]
    if command == "find-references":
        return [command, "--msa", "example_data/divergent.msa.fasta", "--query", "query",
                *common]
    return [command, "--query", str(query), *common]


@pytest.mark.parametrize("command", [
    "detect", "build-panel", "fill-references", "find-references", "curate-panel",
    "type-lineages", "reassort",
])
def test_command_rejects_an_output_path_that_is_a_file(
    tmp_path: Path, no_pipeline, command: str
) -> None:
    occupied = tmp_path / "out"
    occupied.write_text("not a directory\n")
    result = runner.invoke(app, _args(command, tmp_path, occupied))
    assert result.exit_code == 1, result.output
    assert result.exception is None or isinstance(result.exception, SystemExit), (
        f"raw exception: {result.exception!r}"
    )
    assert "Unexpected error" not in result.output
    assert "Output path exists and is not a directory" in " ".join(result.output.split())
    assert occupied.read_text() == "not a directory\n"  # left untouched


@pytest.mark.parametrize("command", ["detect", "build-panel", "type-lineages", "reassort"])
def test_a_new_nested_output_directory_is_still_created(
    tmp_path: Path, no_pipeline, command: str
) -> None:
    """Review focus 5, guard (passes before and after): the check must not stop a command
    from creating an output directory that does not exist yet (the run log goes there)."""
    out = tmp_path / "new" / "nested" / "out"
    runner.invoke(app, _args(command, tmp_path, out))
    assert out.is_dir()
    assert (out / "tessera.log").exists()
