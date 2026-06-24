"""Console tables, TSV outputs, plots, and the self-contained HTML report.

All user-facing values are *similarity* (1 = identical). Outputs:

- ``similarity_windows.tsv``     raw per-window matrix (MSA + query coords, winner)
- ``similarity_stats.tsv``       per-dataset similarity statistics
- ``window_winners.tsv``         per-dataset window-win counts (ties included)
- ``recombination_regions.tsv``  called recombinant regions
- ``similarity_top{N}.{fmt}``    static top-N similarity plot (regions shaded)
- ``similarity_pair.{fmt}``      static major-vs-minor pairwise plot
- ``report.html``                self-contained summary (tables + interactive plot)
"""

from __future__ import annotations

import html
import logging
from math import isnan
from pathlib import Path

from .analyze import AnalysisResult, rank_datasets, stats_sort_key, winner_label
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
    "Minor parent", "Major parent", "MSA start", "MSA end",
    "Query start", "Query end", "Length(bp)", "Windows", "Sim minor", "Sim major",
]


def _region_row(region: Region) -> tuple[str, list]:
    return (
        region.minor_parent,
        [
            region.major_parent, region.msa_start, region.msa_end,
            region.query_start, region.query_end, region.length_bp,
            region.n_windows, region.mean_sim_minor, region.mean_sim_major,
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
    ]
    with open(path, "w") as fo:
        fo.write("\t".join(header) + "\n")
        for r in regions:
            fo.write(
                "\t".join(map(str, [
                    r.minor_parent, r.major_parent, r.msa_start, r.msa_end,
                    r.query_start, r.query_end, r.length_bp, r.n_windows,
                    r.mean_sim_minor, r.mean_sim_major, r.margin,
                ])) + "\n"
            )


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
def _palette(datasets: list[str]):
    """Deterministic label -> colour mapping using a qualitative colormap."""
    import matplotlib as mpl

    cmap = mpl.colormaps["tab10"]
    return {label: cmap(i % cmap.N) for i, label in enumerate(datasets)}


def _ylim(values) -> tuple[float, float]:
    """Data-driven y-limits; ignores nan windows, with a sensible fallback."""
    import numpy as np

    finite_min = np.nanmin(values) if np.isfinite(values).any() else 0.8
    lower = max(0.0, float(finite_min) - 0.02)
    return lower, 1.005


def _shade_regions(ax, regions: list[Region], colors: dict) -> None:
    seen: set[str] = set()
    for region in regions:
        color = colors.get(region.minor_parent, "grey")
        label = f"recombinant: {region.minor_parent}" if region.minor_parent not in seen else None
        seen.add(region.minor_parent)
        ax.axvspan(region.msa_start, region.msa_end, color=color, alpha=0.12, label=label)


def plot_top_n(
    result: WindowSimilarity,
    datasets: list[str],
    regions: list[Region],
    output_dir: Path,
    plot_format: str,
    logger: logging.Logger,
) -> Path | None:
    """Static top-N similarity plot with called regions shaded."""
    from matplotlib import pyplot as plt

    df = result.to_dataframe()
    available = [d for d in datasets if d in df.index]
    if not available:
        logger.warning("No datasets available for the top-N plot; skipping.")
        return None
    subset = df.loc[available]
    colors = _palette(available)

    fig, ax = plt.subplots(figsize=(14, 7))
    _shade_regions(ax, regions, colors)
    for label in subset.index:
        ax.plot(subset.columns, subset.loc[label], lw=2, alpha=0.85,
                color=colors[label], label=label)

    ax.set_ylim(*_ylim(subset.to_numpy()))
    ax.set_xlabel("MSA position (bp)")
    ax.set_ylabel("Similarity to query")
    ax.set_title(f"Similarity to query {result.query}")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), title="Sequences", fontsize="small")
    fig.tight_layout()

    out = output_dir / f"similarity_top{len(available)}.{plot_format}"
    fig.savefig(out)
    plt.close(fig)
    logger.info("Plot saved: %s", out)
    return out


def plot_pairwise(
    result: WindowSimilarity,
    datasets: list[str],
    regions: list[Region],
    output_dir: Path,
    plot_format: str,
    logger: logging.Logger,
) -> Path | None:
    """Static pairwise plot of the two leading datasets, regions shaded."""
    from matplotlib import pyplot as plt

    if len(datasets) < 2:
        logger.warning("Need two datasets for the pairwise plot; skipping.")
        return None
    df = result.to_dataframe()
    seq1, seq2 = datasets[0], datasets[1]
    if seq1 not in df.index or seq2 not in df.index:
        logger.warning("Pairwise datasets not both present; skipping pairwise plot.")
        return None
    colors = _palette([seq1, seq2])

    fig, ax = plt.subplots(figsize=(16, 6))
    _shade_regions(ax, regions, colors)
    ax.plot(df.columns, df.loc[seq1], lw=4, alpha=0.8, color=colors[seq1], label=seq1)
    ax.plot(df.columns, df.loc[seq2], lw=4, alpha=0.8, color=colors[seq2], label=seq2)

    ax.set_ylim(*_ylim(df.loc[[seq1, seq2]].to_numpy()))
    ax.set_xlabel("MSA position (bp)")
    ax.set_ylabel("Similarity to query")
    ax.set_title(f"Similarity to query {result.query}: {seq1} vs {seq2}")
    ax.legend(loc="best")
    fig.tight_layout()

    out = output_dir / f"similarity_pair.{plot_format}"
    fig.savefig(out)
    plt.close(fig)
    logger.info("Plot saved: %s", out)
    return out


def build_interactive_figure(
    result: WindowSimilarity, datasets: list[str], regions: list[Region]
):
    """Plotly figure of the top-N similarities with called regions shaded."""
    import plotly.graph_objects as go

    df = result.to_dataframe()
    fig = go.Figure()
    for region in regions:
        fig.add_vrect(
            x0=region.msa_start, x1=region.msa_end,
            fillcolor="LightSalmon", opacity=0.25, line_width=0,
            annotation_text=region.minor_parent, annotation_position="top left",
        )
    for dataset in datasets:
        if dataset not in df.index:
            continue
        fig.add_trace(go.Scatter(
            x=df.columns, y=df.loc[dataset], mode="lines", name=dataset,
            customdata=result.query_positions,
            hovertemplate="MSA %{x}<br>query %{customdata}<br>similarity %{y:.3f}<extra>"
            + dataset + "</extra>",
        ))
    fig.update_layout(
        title=f"Similarity to query {result.query}",
        xaxis_title="MSA position (bp)", yaxis_title="Similarity",
        hovermode="closest", dragmode="zoom", showlegend=True,
    )
    return fig


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------
def _provenance_html(provenance: dict[str, str]) -> str:
    items = "".join(
        f"<tr><th>{html.escape(str(k))}</th><td>{html.escape(str(v))}</td></tr>"
        for k, v in provenance.items()
    )
    return f"<table class='kv'>{items}</table>"


def _regions_html(regions: list[Region]) -> str:
    if not regions:
        return "<p><em>No recombinant regions detected.</em></p>"
    head = (
        "<tr><th>Minor parent</th><th>Major parent</th><th>MSA start</th><th>MSA end</th>"
        "<th>Query start</th><th>Query end</th><th>Length (bp)</th><th>Windows</th>"
        "<th>Sim minor</th><th>Sim major</th></tr>"
    )
    rows = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(v))}</td>"
            for v in [
                r.minor_parent, r.major_parent, r.msa_start, r.msa_end,
                r.query_start, r.query_end, r.length_bp, r.n_windows,
                r.mean_sim_minor, r.mean_sim_major,
            ]
        )
        + "</tr>"
        for r in regions
    )
    return f"<table class='data'>{head}{rows}</table>"


def write_html_report(
    result: WindowSimilarity,
    analysis: AnalysisResult,
    regions: list[Region],
    datasets: list[str],
    provenance: dict[str, str],
    output_dir: Path,
    logger: logging.Logger,
) -> Path:
    """Write a single self-contained ``report.html``."""
    fig = build_interactive_figure(result, datasets, regions)
    plot_div = fig.to_html(full_html=False, include_plotlyjs="inline")

    stats_rows = "".join(
        "<tr><td>" + html.escape(dataset) + "</td>"
        + "".join(f"<td>{html.escape(str(v))}</td>" for v in values)
        + "</tr>"
        for dataset, values in sorted(analysis.stats.items(), key=stats_sort_key, reverse=True)
    )
    stats_head = "".join(f"<th>{html.escape(h)}</th>" for h in analysis.stats_header)

    doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>RecomFi report: {html.escape(result.query)}</title>
<style>
 body {{ font-family: system-ui, sans-serif; margin: 2rem; color: #222; }}
 h1 {{ font-size: 1.5rem; }} h2 {{ margin-top: 2rem; font-size: 1.2rem; }}
 table {{ border-collapse: collapse; margin: 0.5rem 0; }}
 table.data th, table.data td {{ border: 1px solid #ccc; padding: 4px 8px; text-align: right; }}
 table.data th:first-child, table.data td:first-child {{ text-align: left; }}
 table.kv th {{ text-align: left; padding: 2px 12px 2px 0; color: #555; }}
 table.kv td {{ padding: 2px 0; }}
 .note {{ color: #666; font-size: 0.9rem; }}
</style></head><body>
<h1>RecomFi recombination report</h1>
<h2>Run</h2>
{_provenance_html(provenance)}
<h2>Recombination regions</h2>
<p class="note">Heuristic screen: regions where the query is closer to a minor
parent than to the major parent. Coordinates are given in both MSA columns and
query bases. This is indicative, not a statistical significance test.</p>
{_regions_html(regions)}
<h2>Similarity across the alignment</h2>
{plot_div}
<h2>Per-dataset similarity statistics</h2>
<table class="data"><tr>{stats_head}</tr>{stats_rows}</table>
</body></html>"""

    out = output_dir / "report.html"
    out.write_text(doc)
    logger.info("Report written: %s", out)
    return out


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
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
) -> None:
    """Write every table, plot and the HTML report for a completed scan."""
    output_dir.mkdir(parents=True, exist_ok=True)

    write_windows_tsv(result, per_window_winners, output_dir, logger)
    write_stats_tsv(analysis, output_dir, logger)
    write_winners_tsv(analysis, output_dir, logger)
    write_regions_tsv(regions, output_dir, logger)

    top_datasets = rank_datasets(analysis, top_n)
    logger.info("Top %d nearest datasets: %s", len(top_datasets), ", ".join(top_datasets))
    plot_top_n(result, top_datasets, regions, output_dir, plot_format, logger)

    pair = rank_datasets(analysis, 2)
    plot_pairwise(result, pair, regions, output_dir, plot_format, logger)

    write_html_report(result, analysis, regions, top_datasets, provenance, output_dir, logger)
