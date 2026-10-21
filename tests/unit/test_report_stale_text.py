"""Report and help text that describes what the run did, not an earlier version of it."""

from __future__ import annotations

import re
from pathlib import Path

from typer.testing import CliRunner

from tessera.cli.main import app
from tessera.recomb.report_assets import _REFERENCES
from tessera.recomb.run import RecombParams, run_recomb


def _footer(out: Path) -> str:
    report = (out / "report.html").read_text()
    return re.search(r"<footer>.*?</footer>", report, flags=re.S).group(0)


def _run(example_data: Path, out: Path, logger, **kwargs) -> Path:
    run_recomb(
        RecombParams(msa=example_data / "divergent.msa.fasta", output=out, query="query",
                     window_size=300, window_step=30, **kwargs),
        logger,
    )
    return out


def test_footer_lists_the_plots_in_the_format_they_were_written(
    example_data, tmp_path, logger
) -> None:
    out = tmp_path / "out"
    out.mkdir()
    (out / "similarity_pair.pdf").write_text("left over from an earlier run")
    footer = _footer(_run(example_data, out, logger, plot_format="png"))
    assert "similarity_pair.png" in footer
    assert "similarity_top3.png" in footer
    # Was hard-coded; and a file this run did not write is not listed even if present.
    assert ".pdf" not in footer


def test_footer_lists_only_files_that_exist(example_data, tmp_path, logger) -> None:
    out = _run(example_data, tmp_path / "out", logger, plot_format="png",
               methods=("hmm",), phi=False)
    footer = _footer(out)
    listed = re.findall(r"<code>(.*?)</code>", footer)
    assert listed, "the footer should list the companion files"
    for name in listed:
        assert (out / name).exists(), f"{name} is listed but was not written"
    # A single-caller run writes no method comparison; --no-phi writes no profile.
    assert "recombination_methods.tsv" not in listed
    assert "recombination_profile.tsv" not in listed
    assert "run_provenance.json" in listed


def test_methods_paragraph_describes_the_default_ensemble(
    example_data, tmp_path, logger
) -> None:
    report = (_run(example_data, tmp_path / "out", logger, plot_format="png")
              / "report.html").read_text()
    methods = re.search(r'<details class="methods">.*?</details>', report, flags=re.S).group(0)
    assert "(HMM and the 3SEQ triplet test)" not in methods  # the two-caller description
    for caller in ("3SEQ", "MaxChi", "Bootscan"):
        assert caller in methods


def test_references_cite_every_default_caller() -> None:
    titles = " ".join(title for title, _ in _REFERENCES)
    assert "MaxChi" in titles
    assert "Bootscan" in titles


def test_method_help_lists_every_caller_and_the_real_default() -> None:
    result = CliRunner().invoke(app, ["recomb", "--help"], env={"COLUMNS": "200"})
    assert result.exit_code == 0
    text = " ".join(result.output.replace("│", " ").split())
    assert "geneconv" in text  # run by `all`, and was missing from the list
    assert "all but the legacy heuristic" not in text  # geneconv and barcode are opt-in too
