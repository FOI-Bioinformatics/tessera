"""Plotting for the report: the shared colour palette, static matplotlib plots,
and the interactive plotly figure embedded in the HTML report.

The palette and ``_color_map`` live here because both the plots and the HTML
builders (``report_html``) draw a dataset in the same colour everywhere.
"""

from __future__ import annotations

import logging
from pathlib import Path

from .coverage import CoverageGap
from .regions import Region
from .similarity import WindowSimilarity

# One categorical palette, reused across the mosaic track, the interactive plot,
# the static plots and the table swatches, so a dataset reads as the same colour
# everywhere in the report.
PALETTE = (
    "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#6b46c1",
    "#2c7a7b", "#b7791f", "#3182ce", "#97266d", "#4a5568",
)
GREY = "#94a3b8"


def _color_map(datasets: list[str]) -> dict[str, str]:
    """Deterministic label -> hex colour, stable across the whole report."""
    return {label: PALETTE[i % len(PALETTE)] for i, label in enumerate(datasets)}


def _palette(datasets: list[str]):
    """Matplotlib RGBA map drawn from the shared categorical palette."""
    from matplotlib.colors import to_rgba

    hexmap = _color_map(datasets)
    return {label: to_rgba(hexmap[label]) for label in datasets}


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
    result: WindowSimilarity, datasets: list[str], regions: list[Region],
    coverage_gaps: list[CoverageGap] | None = None,
):
    """Plotly figure of the top-N similarities with called regions shaded."""
    import plotly.graph_objects as go

    df = result.to_dataframe()
    colors = _color_map(datasets)
    fig = go.Figure()
    for gap in coverage_gaps or []:
        fig.add_vrect(
            x0=gap.msa_start, x1=gap.msa_end,
            fillcolor="#94a3b8", opacity=0.18, line_width=0, layer="below",
            annotation_text="low coverage", annotation_position="bottom left",
            annotation_font_size=10, annotation_font_color="#5b6573",
        )
    seen: set[str] = set()
    for region in regions:
        color = colors.get(region.minor_parent, GREY)
        fig.add_vrect(
            x0=region.msa_start, x1=region.msa_end,
            fillcolor=color, opacity=0.12, line_width=0, layer="below",
            annotation_text=("" if region.minor_parent in seen else region.minor_parent),
            annotation_position="top left",
            annotation_font_size=11, annotation_font_color=color,
        )
        seen.add(region.minor_parent)
    for dataset in datasets:
        if dataset not in df.index:
            continue
        fig.add_trace(go.Scatter(
            x=df.columns, y=df.loc[dataset], mode="lines", name=dataset,
            line={"width": 2, "color": colors.get(dataset, GREY)},
            customdata=result.query_positions,
            hovertemplate="MSA %{x:,} bp<br>query %{customdata:,} bp<br>"
            "similarity %{y:.3f}<extra>" + dataset + "</extra>",
        ))
    fig.update_layout(
        template="plotly_white",
        margin={"l": 64, "r": 24, "t": 36, "b": 52},
        height=460,
        font={"family": "system-ui, -apple-system, sans-serif", "size": 12, "color": "#11161f"},
        xaxis={"title": "MSA position (bp)", "gridcolor": "#eef0f3", "zeroline": False,
               "tickformat": ","},
        yaxis={"title": "Similarity to query (1.0 = identical)", "gridcolor": "#eef0f3",
               "zeroline": False},
        hovermode="x unified", dragmode="zoom",
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "x": 0, "title": ""},
    )
    return fig
