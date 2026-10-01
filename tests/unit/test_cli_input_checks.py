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


# --- paths ------------------------------------------------------------------

def test_find_references_rejects_a_missing_msa(tmp_path: Path, no_pipeline) -> None:
    result = runner.invoke(app, [
        "find-references", "--msa", str(tmp_path / "nope.fasta"), "--query", "query",
        "--output", str(tmp_path / "out"),
    ])
    _rejected(result, "MSA file not found")


def test_recomb_rejects_an_output_path_that_is_a_file(tmp_path: Path) -> None:
    occupied = tmp_path / "out"
    occupied.write_text("not a directory\n")
    result = runner.invoke(app, [
        "recomb", "--msa", EXAMPLE, "--query", "query", "--output", str(occupied),
    ])
    _rejected(result, "is not a directory")
    assert occupied.read_text() == "not a directory\n"  # left untouched


# --- numeric options on the panel-building commands ------------------------

@pytest.mark.parametrize("command", ["detect", "fill-references", "build-panel"])
@pytest.mark.parametrize(
    ("option", "value"),
    [("--max-rounds", "0"), ("--window-size", "0"), ("--window-step", "0")],
)
def test_panel_commands_reject_out_of_range_options(
    tmp_path: Path, no_pipeline, command: str, option: str, value: str
) -> None:
    """`fill-references --max-rounds 0` exited 0 having built no alignment and no
    report; a bad window reached the scan only after seeding and alignment."""
    result = runner.invoke(app, [
        command, "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        option, value,
    ])
    _rejected(result, f"Invalid {option}")


@pytest.mark.parametrize(("option", "value"), [("--window-size", "0"), ("--window-step", "0")])
def test_find_references_rejects_out_of_range_windows(
    tmp_path: Path, no_pipeline, option: str, value: str
) -> None:
    result = runner.invoke(app, [
        "find-references", "--msa", EXAMPLE, "--query", "query",
        "--output", str(tmp_path / "out"), option, value,
    ])
    _rejected(result, f"Invalid {option}")


# --- reassort ---------------------------------------------------------------

def _segments(tmp_path: Path) -> Path:
    path = tmp_path / "segments.fasta"
    path.write_text(">HA\nACGTACGT\n>NA\nACGTACGT\n")
    return path


@pytest.fixture
def no_reassort(monkeypatch):
    monkeypatch.setattr("tessera.cli.cmd_reassort.assign_segments", _must_not_run)


@pytest.mark.parametrize(
    ("option", "value", "needle"),
    [("--ani-floor", "500", "Invalid --ani-floor"),
     ("--ani-floor", "-1", "Invalid --ani-floor"),
     ("--margin", "-3", "Invalid --margin")],
)
def test_reassort_rejects_out_of_range_options(
    tmp_path: Path, no_reassort, option: str, value: str, needle: str
) -> None:
    """A negative margin makes every near-best set empty (a clonal pair then reads as
    undetermined); an ANI floor above 100 leaves every segment unassigned. Both exited 0."""
    result = runner.invoke(app, [
        "reassort", "--query", str(_segments(tmp_path)), "--output", str(tmp_path / "out"),
        option, value,
    ])
    _rejected(result, needle)


def test_reassort_rejects_a_dataset_override_for_an_unknown_segment(
    tmp_path: Path, no_reassort
) -> None:
    """A typo in the segment name silently fell back to dataset auto-detection."""
    result = runner.invoke(app, [
        "reassort", "--query", str(_segments(tmp_path)), "--output", str(tmp_path / "out"),
        "--dataset", "HAA=nextstrain/flu/h3n2/ha",
    ])
    _rejected(result, "HAA")
    assert "HA, NA" in " ".join(result.output.split())  # names the segments that do exist


def test_segment_scan_tsv_keeps_one_row_per_segment(tmp_path: Path, monkeypatch) -> None:
    """An aligner error message spans several lines and may hold tabs; written verbatim
    into the note column it broke the table."""
    from tessera.reassort.assign import ReassortmentResult, SegmentAssignment
    from tessera.reassort.scan import SegmentScan

    def fake_assign(query, **kwargs):
        return ReassortmentResult(
            segments=[SegmentAssignment("HA", "dataset", None, None, 0.0, "unassigned")],
            verdict="undetermined", groups=[], pair_notes=[],
            scans=[SegmentScan("HA", False, False, 0,
                               "scan failed: mafft failed:\nline1\tx\nline2")],
        )

    monkeypatch.setattr("tessera.cli.cmd_reassort.assign_segments", fake_assign)
    out = tmp_path / "out"
    result = runner.invoke(app, [
        "reassort", "--query", str(_segments(tmp_path)), "--output", str(out),
    ])
    assert result.exit_code == 0, result.output
    lines = (out / "segment_scan.tsv").read_text().splitlines()
    assert len(lines) == 2  # header + one segment
    assert all(len(line.split("\t")) == 4 for line in lines)
    assert lines[1].split("\t")[3] == "scan failed: mafft failed: line1 x line2"


# --- type-lineages reads the same collection the other commands do ---------

def test_type_lineages_accepts_what_a_collection_may_hold(tmp_path: Path, monkeypatch) -> None:
    """Every other command treats each file in the collection as a genome, whatever its
    extension and gzip-compressed or not; type-lineages filtered on three suffixes and
    reported "No FASTA genomes found"."""
    import gzip

    collection = tmp_path / "collection"
    collection.mkdir()
    with gzip.open(collection / "a.fasta.gz", "wt") as fo:
        fo.write(">a\nACGTACGT\n")
    (collection / "b.fas").write_text(">b\nACGTACGT\n")
    (collection / ".DS_Store").write_text("not a genome")
    seen: list[str] = []

    def fake_assign(genomes, **kwargs):
        seen.extend(p.name for p in genomes)
        return [("a", "L1", "denovo"), ("b.fas", "L1", "denovo")]

    monkeypatch.setattr("tessera.discover.lineage_assign.assign_lineages", fake_assign)
    out = tmp_path / "out"
    result = runner.invoke(app, [
        "type-lineages", "--collection", str(collection), "--output", str(out),
    ])
    assert result.exit_code == 0, result.output
    assert seen == ["a.fasta.gz", "b.fas"]  # hidden files are not genomes
    assert (out / "lineages.tsv").exists()


def test_type_lineages_rejects_a_file_that_is_not_fasta(tmp_path: Path, no_pipeline) -> None:
    collection = _collection(tmp_path)
    (collection / "notes.txt").write_text("these are my notes\n")
    result = runner.invoke(app, [
        "type-lineages", "--collection", str(collection), "--output", str(tmp_path / "out"),
    ])
    _rejected(result, "does not look like a FASTA file")
