"""Console tables and TSV outputs for a completed scan.

The terminal-facing tables (``print_*``) and the tab-separated companion files
(``write_*_tsv``); the HTML report lives in ``report_html`` and the plots in
``report_plots``.
"""

from __future__ import annotations

import logging
from math import isnan
from pathlib import Path

from .analyze import AnalysisResult, stats_sort_key, winner_label
from .coverage import CoverageGap
from .diagnostics import RecombinationSignal
from .regions import Region
from .similarity import WindowSimilarity


# ---------------------------------------------------------------------------
# Terminal tables
# ---------------------------------------------------------------------------
def print_formatted_table(data, header: list[str], echo=print) -> None:
    """Print ``data`` (list of (key, value-or-list)) as an aligned table."""
    num_columns = len(header) - 1

    expanded_data = []
    for key, value in data:
        if isinstance(value, list):
            expanded_values = value + [None] * (num_columns - len(value))
            expanded_data.append([key, *expanded_values[:num_columns]])
        else:
            expanded_data.append([key, value] + [None] * (num_columns - 1))

    column_widths = [
        max((len(str(row[i])) for row in expanded_data if row[i] is not None), default=0)
        for i in range(num_columns + 1)
    ]
    header_widths = [len(str(header[i])) for i in range(num_columns + 1)]
    column_widths = [max(column_widths[i], header_widths[i]) + 2 for i in range(num_columns + 1)]

    offset = "  "
    echo(offset + "".join(f"{str(header[i]):<{column_widths[i]}}" for i in range(num_columns + 1)))
    echo(offset + "-" * sum(column_widths))
    for row in expanded_data:
        echo(
            offset
            + "".join(
                f"{str(row[i]) if row[i] is not None else '':<{column_widths[i]}}"
                for i in range(num_columns + 1)
            )
        )


def print_summary(analysis: AnalysisResult, echo=print) -> None:
    """Print the winners and per-dataset similarity-statistics tables."""
    echo("Window winners (no ties):")
    print_formatted_table(
        sorted(analysis.winners_no_ties.items(), key=lambda x: x[1], reverse=True),
        header=["Dataset", "Windows"],
        echo=echo,
    )
    echo("")
    echo("Window winners (ties allowed):")
    print_formatted_table(
        sorted(analysis.winners_with_ties.items(), key=lambda x: x[1], reverse=True),
        header=["Dataset", "Windows"],
        echo=echo,
    )
    echo("")
    echo("Per-dataset similarity stats:")
    print_formatted_table(
        sorted(analysis.stats.items(), key=stats_sort_key, reverse=True),
        header=analysis.stats_header,
        echo=echo,
    )
    echo("")


REGION_HEADER = [
    "Minor parent", "Major parent", "Query start", "Query end", "Length(bp)",
    "Sim minor", "Sim major", "Support", "q-value", "Breakpoint", "Method(s)",
]


def _methods_text(region: Region) -> str:
    """The caller(s) that found a region, for the console/TSV ('-' for donor-absent)."""
    return ", ".join(region.methods) if region.methods else "-"


def _breakpoint_str(region: Region) -> str:
    if region.breakpoint_lo is None:
        return "-"
    if region.breakpoint_lo == region.breakpoint_hi:
        return str(region.breakpoint_lo)
    return f"{region.breakpoint_lo}-{region.breakpoint_hi}"


def _region_row(region: Region) -> tuple[str, list]:
    support = "-" if region.support is None else f"{region.support:.2f}"
    qval = "-" if region.qvalue is None else f"{region.qvalue:.2g}"
    return (
        region.minor_parent,
        [
            region.major_parent, region.query_start, region.query_end,
            region.length_bp, region.mean_sim_minor, region.mean_sim_major,
            support, qval, _breakpoint_str(region), _methods_text(region),
        ],
    )


def print_regions(regions: list[Region], major_parent: str | None, echo=print) -> None:
    """Print the called recombinant regions (or a 'none found' note)."""
    echo(f"Recombination regions (major parent: {major_parent or 'n/a'}):")
    if not regions:
        echo("  none detected")
        echo("")
        return
    print_formatted_table([_region_row(r) for r in regions], header=REGION_HEADER, echo=echo)
    echo("")


COVERAGE_HEADER = [
    "Best ref", "Query start", "Query end", "Length(bp)", "Windows", "Mean best", "Kind",
]


def print_coverage(gaps: list[CoverageGap], threshold: float, echo=print) -> None:
    """Print under-covered regions (possible missing references), or a clear note."""
    echo(f"Reference coverage (best-similarity threshold {threshold:.3f}):")
    if not gaps:
        echo("  no coverage gaps -- every region has a close reference")
        echo("")
        return
    rows = [
        (g.best_label, [g.query_start, g.query_end, g.length_bp, g.n_windows,
                        g.mean_best, g.kind])
        for g in gaps
    ]
    print_formatted_table(rows, header=COVERAGE_HEADER, echo=echo)
    echo("  ^ the closest reference here is poor; the true source may be missing. "
         "Run 'tessera find-references' to search NCBI.")
    echo("")


# ---------------------------------------------------------------------------
# TSV outputs
# ---------------------------------------------------------------------------
def write_windows_tsv(
    result: WindowSimilarity,
    per_window_winners: list[list[str]],
    output_dir: Path,
    logger: logging.Logger,
) -> None:
    """Write the full per-window similarity matrix with MSA + query coordinates."""
    datasets = list(result.similarities)
    path = output_dir / "similarity_windows.tsv"
    logger.info("Writing per-window similarity matrix: %s", path)
    with open(path, "w") as fo:
        fo.write("\t".join(["msa_position", "query_position", "winner", *datasets]) + "\n")
        for i, (msa_pos, q_pos) in enumerate(
            zip(result.positions, result.query_positions, strict=True)
        ):
            row = [msa_pos, q_pos, winner_label(per_window_winners[i])]
            row += [
                "NA" if isnan(result.similarities[ds][i]) else f"{result.similarities[ds][i]:.4f}"
                for ds in datasets
            ]
            fo.write("\t".join(map(str, row)) + "\n")


def write_stats_tsv(analysis: AnalysisResult, output_dir: Path, logger: logging.Logger) -> None:
    path = output_dir / "similarity_stats.tsv"
    logger.info("Writing per-dataset similarity stats: %s", path)
    with open(path, "w") as fo:
        fo.write("\t".join(analysis.stats_header) + "\n")
        for dataset, values in sorted(analysis.stats.items(), key=stats_sort_key, reverse=True):
            fo.write("\t".join(map(str, [dataset, *values])) + "\n")


def write_winners_tsv(analysis: AnalysisResult, output_dir: Path, logger: logging.Logger) -> None:
    path = output_dir / "window_winners.tsv"
    logger.info("Writing window winners: %s", path)
    with open(path, "w") as fo:
        fo.write("Dataset\tNumber of windows (ties included)\n")
        for dataset, num in sorted(
            analysis.winners_with_ties.items(), key=lambda x: x[1], reverse=True
        ):
            fo.write(f"{dataset}\t{num}\n")


def write_regions_tsv(regions: list[Region], output_dir: Path, logger: logging.Logger) -> None:
    path = output_dir / "recombination_regions.tsv"
    logger.info("Writing recombination regions: %s", path)
    header = [
        "minor_parent", "major_parent", "msa_start", "msa_end",
        "query_start", "query_end", "length_bp", "n_windows",
        "mean_sim_minor", "mean_sim_major", "margin",
        "support", "pvalue", "qvalue", "posterior", "breakpoint_lo", "breakpoint_hi",
        "donor_undercovered", "donor_absent", "methods", "parent_free_support",
    ]
    with open(path, "w") as fo:
        fo.write("\t".join(header) + "\n")
        for r in regions:
            fo.write(
                "\t".join(map(str, [
                    r.minor_parent, r.major_parent, r.msa_start, r.msa_end,
                    r.query_start, r.query_end, r.length_bp, r.n_windows,
                    r.mean_sim_minor, r.mean_sim_major, r.margin,
                    "NA" if r.support is None else r.support,
                    "NA" if r.pvalue is None else r.pvalue,
                    "NA" if r.qvalue is None else r.qvalue,
                    r.posterior_support,
                    "NA" if r.breakpoint_lo is None else r.breakpoint_lo,
                    "NA" if r.breakpoint_hi is None else r.breakpoint_hi,
                    "yes" if r.donor_undercovered else "no",
                    "yes" if r.donor_absent else "no",
                    ",".join(r.methods) if r.methods else "NA",
                    "yes" if r.parent_free_support else "no",
                ])) + "\n"
            )


def write_methods_tsv(
    breakdown: list[dict], methods_run: tuple[str, ...], output_dir: Path,
    logger: logging.Logger,
) -> None:
    """Write the per-region x per-method agreement matrix (the ensemble breakdown)."""
    path = output_dir / "recombination_methods.tsv"
    logger.info("Writing method comparison: %s", path)
    with open(path, "w") as fo:
        fo.write("\t".join(["minor_parent", "query_start", "query_end",
                            *methods_run, "parent_free_support"]) + "\n")
        for b in breakdown:
            called = b["per_method_support"]
            cells = [("yes" if m in called else "no") for m in methods_run]
            fo.write("\t".join(map(str, [
                b["minor_parent"], b["query_start"], b["query_end"], *cells,
                "yes" if b["parent_free_support"] else "no",
            ])) + "\n")


def write_coverage_tsv(
    gaps: list[CoverageGap], threshold: float, output_dir: Path, logger: logging.Logger
) -> None:
    """Write the under-covered regions (possible missing references)."""
    path = output_dir / "coverage_gaps.tsv"
    logger.info("Writing reference-coverage gaps: %s", path)
    header = [
        "msa_start", "msa_end", "query_start", "query_end", "length_bp",
        "n_windows", "best_reference", "mean_best_similarity", "kind",
    ]
    with open(path, "w") as fo:
        fo.write(f"# best-similarity threshold\t{threshold:.4f}\n")
        fo.write("\t".join(header) + "\n")
        for g in gaps:
            fo.write("\t".join(map(str, [
                g.msa_start, g.msa_end, g.query_start, g.query_end, g.length_bp,
                g.n_windows, g.best_label, g.mean_best, g.kind,
            ])) + "\n")


def write_profile_tsv(
    signal: RecombinationSignal, output_dir: Path, logger: logging.Logger
) -> None:
    """Write the per-informative-site PHI profile (the parent-free signal track)."""
    path = output_dir / "recombination_profile.tsv"
    logger.info("Writing recombination signal profile: %s", path)
    with open(path, "w") as fo:
        fo.write(
            f"# PHI p-value\t{signal.phi_p:.4g}\t(window {signal.phi_window} "
            f"informative sites, {signal.n_informative} sites, Rmin {signal.rmin})\n"
        )
        fo.write("msa_pos\tquery_pos\tphi\n")
        for msa_pos, query_pos, value in signal.profile:
            fo.write(f"{msa_pos}\t{query_pos}\t{value:.4f}\n")
