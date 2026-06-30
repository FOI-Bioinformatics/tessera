"""Reporting entry point: orchestrate the tables, plots and HTML for a scan.

The reporting code is split by responsibility:

- ``report_text``    console tables and the TSV companion files
- ``report_plots``   the shared colour palette and the static / interactive plots
- ``report_assets``  static presentation constants (stylesheet, glossary, references)
- ``report_html``    the self-contained HTML report

This module wires them together (``write_reports``) and re-exports the
console/HTML entry points that the rest of the package imports. All user-facing
values are *similarity* (1 = identical). Outputs:

- ``similarity_windows.tsv``     raw per-window matrix (MSA + query coords, winner)
- ``similarity_stats.tsv``       per-dataset similarity statistics
- ``window_winners.tsv``         per-dataset window-win counts (ties included)
- ``recombination_regions.tsv``  called recombinant regions
- ``similarity_top{N}.{fmt}``    static top-N similarity plot (regions shaded)
- ``similarity_pair.{fmt}``      static major-vs-minor pairwise plot
- ``report.html``                self-contained summary (tables + interactive plot)
"""

from __future__ import annotations

import logging
from pathlib import Path

from .analyze import AnalysisResult, rank_datasets
from .coverage import CoverageGap
from .diagnostics import RecombinationSignal
from .regions import Region
from .report_html import write_html_report
from .report_plots import build_interactive_figure, plot_pairwise, plot_top_n
from .report_text import (
    print_coverage,
    print_regions,
    print_summary,
    write_coverage_tsv,
    write_methods_tsv,
    write_profile_tsv,
    write_regions_tsv,
    write_stats_tsv,
    write_windows_tsv,
    write_winners_tsv,
)
from .similarity import WindowSimilarity
from .typing import LineageMap

__all__ = [
    "print_summary", "print_regions", "print_coverage",
    "build_interactive_figure", "plot_top_n", "plot_pairwise",
    "write_html_report", "write_reports",
]


def write_reports(
    result: WindowSimilarity,
    analysis: AnalysisResult,
    regions: list[Region],
    per_window_winners: list[list[str]],
    provenance: dict[str, str],
    output_dir: Path,
    top_n: int,
    plot_format: str,
    logger: logging.Logger,
    coverage_gaps: list[CoverageGap] | None = None,
    coverage_threshold: float = 0.0,
    extra_sections: list[tuple[str, str]] | None = None,
    lineage_map: LineageMap | None = None,
    query_lineage: str | None = None,
    signal: RecombinationSignal | None = None,
    organism: str | None = None,
    methods_run: tuple[str, ...] = (),
    method_breakdown: list[dict] | None = None,
    per_major: dict[str, str] | None = None,
) -> None:
    """Write every table, plot and the HTML report for a completed scan."""
    output_dir.mkdir(parents=True, exist_ok=True)
    gaps = coverage_gaps or []

    write_windows_tsv(result, per_window_winners, output_dir, logger)
    write_stats_tsv(analysis, output_dir, logger)
    write_winners_tsv(analysis, output_dir, logger)
    write_regions_tsv(regions, output_dir, logger)
    write_coverage_tsv(gaps, coverage_threshold, output_dir, logger)
    if signal is not None:
        write_profile_tsv(signal, output_dir, logger)
    if len(methods_run) > 1 and method_breakdown is not None:
        write_methods_tsv(method_breakdown, methods_run, output_dir, logger)

    top_datasets = rank_datasets(analysis, top_n)
    logger.info("Top %d nearest datasets: %s", len(top_datasets), ", ".join(top_datasets))
    plot_top_n(result, top_datasets, regions, output_dir, plot_format, logger)

    pair = rank_datasets(analysis, 2)
    plot_pairwise(result, pair, regions, output_dir, plot_format, logger)

    write_html_report(
        result, analysis, regions, top_datasets, provenance, output_dir, logger,
        coverage_gaps=gaps, coverage_threshold=coverage_threshold,
        extra_sections=extra_sections, lineage_map=lineage_map,
        query_lineage=query_lineage, signal=signal, organism=organism,
        methods_run=methods_run, method_breakdown=method_breakdown, per_major=per_major,
    )
