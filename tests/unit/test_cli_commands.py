"""Every CLI command: it is registered, its help renders, and it fails cleanly.

The CLI layer is what users actually touch and was the least-tested part of the
package. These do not re-test the science -- the callers have their own suites --
they check the layer around it: that a command exists, that its options parse, that
a bad argument produces a clean message and a non-zero exit rather than a traceback,
and that the one command needing no external tool runs end to end.

`example_data/` is used deliberately: those alignments need no aligner, no network
and no downloaded data, so the whole file runs in CI.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

from tessera.cli.main import app

runner = CliRunner()

COMMANDS = [
    "build-panel", "curate-panel", "detect", "fill-references",
    "find-references", "msa", "reassort", "recomb", "type-lineages",
]

EXAMPLE = Path("example_data/divergent.msa.fasta")


# --- registration and help ------------------------------------------------

def test_every_command_is_registered() -> None:
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    listed = " ".join(result.output.split())
    for name in COMMANDS:
        assert name in listed, f"{name} missing from `tessera --help`"


@pytest.mark.parametrize("command", COMMANDS)
def test_command_help_renders(command: str) -> None:
    result = runner.invoke(app, [command, "--help"])
    assert result.exit_code == 0
    assert result.exception is None


@pytest.mark.parametrize("command", COMMANDS)
def test_missing_required_options_exit_as_usage_errors(command: str) -> None:
    # Click's own usage error (exit 2), not a traceback from inside the command.
    result = runner.invoke(app, [command])
    assert result.exit_code != 0
    assert "Traceback" not in result.output


def test_version_flag() -> None:
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert "tessera" in result.output


# --- clean failures rather than tracebacks --------------------------------

def _clean_failure(result) -> None:
    """Exit non-zero without leaking an unexpected exception to the user."""
    assert result.exit_code != 0
    assert "Traceback" not in result.output
    assert "Unexpected error" not in result.output


def test_recomb_missing_msa(tmp_path: Path) -> None:
    result = runner.invoke(app, [
        "recomb", "--msa", str(tmp_path / "nope.fasta"), "--query", "q",
        "--output", str(tmp_path / "out"),
    ])
    _clean_failure(result)


def test_recomb_unknown_query_label_lists_what_is_available(tmp_path: Path) -> None:
    result = runner.invoke(app, [
        "recomb", "--msa", str(EXAMPLE), "--query", "not_in_the_alignment",
        "--output", str(tmp_path / "out"),
    ])
    _clean_failure(result)
    assert "parent_A" in result.output  # the message names the labels that do exist


def test_recomb_unknown_method(tmp_path: Path) -> None:
    result = runner.invoke(app, [
        "recomb", "--msa", str(EXAMPLE), "--query", "query",
        "--output", str(tmp_path / "out"), "--method", "nosuchcaller",
    ])
    _clean_failure(result)


def test_recomb_unknown_metric(tmp_path: Path) -> None:
    result = runner.invoke(app, [
        "recomb", "--msa", str(EXAMPLE), "--query", "query",
        "--output", str(tmp_path / "out"), "--metric", "jc69",
    ])
    _clean_failure(result)


def test_recomb_unknown_plot_format(tmp_path: Path) -> None:
    result = runner.invoke(app, [
        "recomb", "--msa", str(EXAMPLE), "--query", "query",
        "--output", str(tmp_path / "out"), "--plot-format", "bmp",
    ])
    _clean_failure(result)


def test_recomb_window_larger_than_the_alignment(tmp_path: Path) -> None:
    result = runner.invoke(app, [
        "recomb", "--msa", str(EXAMPLE), "--query", "query",
        "--output", str(tmp_path / "out"), "--window-size", "10000000",
    ])
    _clean_failure(result)


def test_msa_missing_query(tmp_path: Path) -> None:
    coll = tmp_path / "coll"
    coll.mkdir()
    (coll / "a.fasta").write_text(">a\nACGT\n")
    result = runner.invoke(app, [
        "msa", "--query", str(tmp_path / "absent.fasta"), "--collection", str(coll),
        "--output", str(tmp_path / "out.fasta"),
    ])
    _clean_failure(result)


def test_msa_empty_collection(tmp_path: Path) -> None:
    coll = tmp_path / "empty"
    coll.mkdir()
    query = tmp_path / "q.fasta"
    query.write_text(">q\nACGT\n")
    result = runner.invoke(app, [
        "msa", "--query", str(query), "--collection", str(coll),
        "--output", str(tmp_path / "out.fasta"),
    ])
    _clean_failure(result)


def test_msa_unknown_aligner(tmp_path: Path) -> None:
    coll = tmp_path / "coll"
    coll.mkdir()
    (coll / "a.fasta").write_text(">a\nACGT\n")
    query = tmp_path / "q.fasta"
    query.write_text(">q\nACGT\n")
    result = runner.invoke(app, [
        "msa", "--query", str(query), "--collection", str(coll),
        "--output", str(tmp_path / "out.fasta"), "--aligner", "nosuchaligner",
    ])
    _clean_failure(result)


def test_curate_panel_missing_collection(tmp_path: Path) -> None:
    query = tmp_path / "q.fasta"
    query.write_text(">q\nACGT\n")
    result = runner.invoke(app, [
        "curate-panel", "--query", str(query), "--collection", str(tmp_path / "absent"),
        "--output", str(tmp_path / "out"),
    ])
    _clean_failure(result)


def test_reassort_missing_query(tmp_path: Path) -> None:
    result = runner.invoke(app, [
        "reassort", "--query", str(tmp_path / "absent.fasta"),
        "--output", str(tmp_path / "out"),
    ])
    _clean_failure(result)


def test_type_lineages_missing_collection(tmp_path: Path) -> None:
    result = runner.invoke(app, [
        "type-lineages", "--collection", str(tmp_path / "absent"),
        "--output", str(tmp_path / "out"),
    ])
    _clean_failure(result)


# --- the one command that runs with no external tool ----------------------

def test_recomb_runs_end_to_end_and_writes_its_outputs(tmp_path: Path) -> None:
    out = tmp_path / "out"
    result = runner.invoke(app, [
        "recomb", "--msa", str(EXAMPLE), "--query", "query", "--output", str(out),
        "--window-size", "300", "--window-step", "30", "--plot-format", "png",
    ])
    assert result.exit_code == 0, result.output

    for name in (
        "recombination_regions.tsv", "recombination_profile.tsv",
        "similarity_windows.tsv", "similarity_stats.tsv", "window_winners.tsv",
        "coverage_gaps.tsv", "report.html", "run_provenance.json",
    ):
        assert (out / name).is_file(), f"{name} was not written"

    # and it found the region the example exists to demonstrate
    rows = (out / "recombination_regions.tsv").read_text().splitlines()
    assert len(rows) >= 2
    assert "parent_B" in rows[1]


def test_recomb_single_method_omits_the_ensemble_breakdown(tmp_path: Path) -> None:
    out = tmp_path / "out"
    result = runner.invoke(app, [
        "recomb", "--msa", str(EXAMPLE), "--query", "query", "--output", str(out),
        "--window-size", "300", "--window-step", "30", "--method", "hmm",
        "--plot-format", "png",
    ])
    assert result.exit_code == 0, result.output
    # recombination_methods.tsv is an ensemble artefact; one caller means no comparison
    assert not (out / "recombination_methods.tsv").exists()


def test_recomb_no_phi_skips_the_parent_free_profile(tmp_path: Path) -> None:
    out = tmp_path / "out"
    result = runner.invoke(app, [
        "recomb", "--msa", str(EXAMPLE), "--query", "query", "--output", str(out),
        "--window-size", "300", "--window-step", "30", "--no-phi", "--plot-format", "png",
    ])
    assert result.exit_code == 0, result.output
    assert not (out / "recombination_profile.tsv").exists()
