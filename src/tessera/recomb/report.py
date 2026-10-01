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
- ``similarity_pair.{fmt}``      static pairwise plot: the major parent against the
                                 donor of the longest called region (the two leading
                                 window winners when no donor was called)
- ``report.html``                self-contained summary (tables + interactive plot)

Under informative-site windowing (near-identical panels) the windows the HMM
segmented hold identity at polymorphic columns only. Those are written separately,
never as "similarity":

- ``informative_site_windows.tsv``      per-window identity at informative sites
- ``informative_sites_top{N}.{fmt}``    static plot of the same track
"""

from __future__ import annotations

import logging
from pathlib import Path

from .analyze import AnalysisResult, rank_datasets, winners_per_window
from .regions import Region
from .report_context import ReportContext
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
    write_site_windows_tsv,
    write_stats_tsv,
    write_windows_tsv,
    write_winners_tsv,
)
from .similarity import WindowSimilarity

__all__ = [
    "ReportContext",
    "print_summary", "print_regions", "print_coverage",
    "build_interactive_figure", "plot_top_n", "plot_pairwise",
    "write_html_report", "write_reports", "pair_datasets",
]


def pair_datasets(regions: list[Region], ranked: list[str]) -> list[str]:
    """The two datasets of the pairwise plot: major parent, then the leading donor.

    The leading donor is the donor of the longest called, donor-present region. Without
    one there is no minor parent to show and the plot falls back to ``ranked`` (the two
    leading window winners) -- which is not "major versus minor" when the panel holds a
    near-duplicate of the backbone, hence the region-based choice whenever possible.
    """
    present = [r for r in regions if not r.donor_absent]
    if not present:
        return ranked
    longest = max(present, key=lambda r: r.length_bp)
    return [longest.major_parent, longest.minor_parent]


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
    ctx: ReportContext,
) -> None:
    """Write every table, plot and the HTML report for a completed scan."""
    output_dir.mkdir(parents=True, exist_ok=True)
    gaps = ctx.gaps

    # The files written beside the report, in the order the footer lists them. Built as
    # they are written, so the footer never names a file that does not exist or a plot
    # in a format that was not asked for.
    companions: list[str] = []

    write_regions_tsv(regions, output_dir, logger)
    companions.append("recombination_regions.tsv")
    if len(ctx.methods_run) > 1 and ctx.method_breakdown is not None:
        write_methods_tsv(
            ctx.method_breakdown, ctx.methods_run, output_dir, logger, ctx.methods_not_run
        )
        companions.append("recombination_methods.tsv")
    write_coverage_tsv(gaps, ctx.coverage_threshold, output_dir, logger)
    companions.append("coverage_gaps.tsv")
    if ctx.signal is not None:
        write_profile_tsv(ctx.signal, output_dir, logger)
        companions.append("recombination_profile.tsv")
    write_winners_tsv(analysis, output_dir, logger)
    write_stats_tsv(analysis, output_dir, logger)
    write_windows_tsv(result, per_window_winners, output_dir, logger)
    companions += ["window_winners.tsv", "similarity_stats.tsv", "similarity_windows.tsv"]

    # Rank on the windows the caller segmented: on a near-identical panel base-pair
    # windows tie almost everywhere, so they order the references poorly.
    ranking = ctx.site_analysis or analysis
    top_datasets = rank_datasets(ranking, top_n)
    logger.info("Top %d nearest datasets: %s", len(top_datasets), ", ".join(top_datasets))
    plots = [plot_top_n(result, top_datasets, regions, output_dir, plot_format, logger)]
    if ctx.site_result is not None:
        write_site_windows_tsv(
            ctx.site_result, winners_per_window(ctx.site_result), output_dir, logger
        )
        companions.append("informative_site_windows.tsv")
        plots.append(plot_top_n(
            ctx.site_result, top_datasets, regions, output_dir, plot_format, logger,
            ylabel="Identity to query at informative sites",
            title=f"Identity at informative sites, query {result.query}",
            stem="informative_sites_top",
        ))

    pair = pair_datasets(regions, rank_datasets(ranking, 2))
    plots.append(plot_pairwise(result, pair, regions, output_dir, plot_format, logger))
    # A plot function returns None when it had nothing to draw and wrote no file.
    companions += [path.name for path in plots if path is not None]
    if (output_dir / "run_provenance.json").exists():
        companions.append("run_provenance.json")

    write_html_report(
        result, analysis, regions, top_datasets, provenance, output_dir, logger, ctx,
        companion_files=companions,
    )
