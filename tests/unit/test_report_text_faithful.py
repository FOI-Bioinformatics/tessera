"""Report, glossary and provenance text that describes what the run did.

Audit 4, item B5 and the smaller reporting items: the provenance gave base-pair window
geometry under informative-site windowing, pointed at a coordinate convention it did not
record, and the methods text and glossary described the barcode caller, the breakpoint
interval and the major parent as something they are not.
"""

from __future__ import annotations

import json
from pathlib import Path

from tessera.recomb.report_assets import _GLOSSARY
from tessera.recomb.run import RecombParams, run_recomb


def _run(example_data: Path, tmp_path: Path, logger, **kwargs) -> tuple[dict, str]:
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=example_data / "divergent.msa.fasta", output=out, query="query",
                     window_size=300, window_step=30, plot_format="png", **kwargs),
        logger,
    )
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    return run, (out / "report.html").read_text()


def test_window_line_names_the_base_pair_scan_under_informative_windowing(
    example_data, tmp_path, logger
) -> None:
    run, _ = _run(example_data, tmp_path, logger, informative_sites=True)
    assert run["windowing"].startswith("informative-site")
    assert run["window / step"].startswith("300 / 30 (base-pair scan")


def test_window_line_is_unchanged_under_base_pair_windowing(
    example_data, tmp_path, logger
) -> None:
    """Guard (passes before and after)."""
    run, _ = _run(example_data, tmp_path, logger, informative_sites=False)
    assert run["window / step"] == "300 / 30"


def test_provenance_records_the_coordinate_convention(example_data, tmp_path, logger) -> None:
    run, _ = _run(example_data, tmp_path, logger)
    assert run["coordinates"] == "0-based, end-exclusive; MSA columns and query bases"


def test_methods_text_does_not_say_every_caller_applies_bh(
    example_data, tmp_path, logger
) -> None:
    _, report = _run(example_data, tmp_path, logger)
    assert "Each caller applies Benjamini-Hochberg" not in report
    assert "barcode caller has no significance test" in report


def test_glossary_describes_the_breakpoint_and_major_parent_faithfully() -> None:
    glossary = dict(_GLOSSARY)
    assert "site caller" in glossary["Breakpoint"]
    major = glossary["Major parent (backbone)"]
    assert "most windows overall" not in major
    assert "siblings" in major and "pooled" in major
