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
from .coverage import CoverageGap
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
    "Sim minor", "Sim major", "Support", "q-value", "Breakpoint",
]


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
            support, qval, _breakpoint_str(region),
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
         "Run 'recomfi find-references' to search NCBI.")
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
        "donor_undercovered", "donor_absent",
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
                ])) + "\n"
            )


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


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------
_CSS = """
:root{
  --ink:#11161f;--muted:#5b6573;--faint:#8b94a3;--line:#e6e9ee;
  --panel:#f6f7f9;--paper:#fff;--bg:#fbfcfd;--accent:#2b6cb0;--hl:#eef4fb;
  --mono:ui-monospace,"SF Mono","JetBrains Mono",Menlo,Consolas,monospace;
  --sans:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font-family:var(--sans);
  font-size:15px;line-height:1.55;-webkit-font-smoothing:antialiased}
.wrap{max-width:1120px;margin:0 auto;padding:48px 28px 72px}
.mono{font-family:var(--mono);font-variant-numeric:tabular-nums}
.eyebrow{font-size:11px;font-weight:600;letter-spacing:.16em;text-transform:uppercase;color:var(--faint)}
.sw{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:6px;flex:none}
.sw.hatch{background:repeating-linear-gradient(45deg,#64748b 0 3px,#cbd5e1 3px 6px)}
header{border-bottom:1px solid var(--line);padding-bottom:26px;margin-bottom:30px}
header h1{font-size:30px;font-weight:600;margin:.35em 0 .55em;letter-spacing:-.01em;word-break:break-word}
.verdict{font-size:19px;line-height:1.5;margin:0;max-width:76ch}
.verdict strong{font-weight:650}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:14px;margin-bottom:40px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:15px 18px}
.card .k{font-size:11px;font-weight:600;letter-spacing:.07em;text-transform:uppercase;color:var(--faint);margin-bottom:9px}
.card .v{font-size:17px;font-weight:550;display:flex;align-items:center;flex-wrap:wrap;gap:3px}
.big{font-size:26px;font-weight:650;font-family:var(--mono)}
.card .sub{font-size:12px;color:var(--muted);font-weight:400;font-family:var(--mono)}
.section{margin:42px 0}
.section>.eyebrow{margin-bottom:14px}
.cap{color:var(--muted);font-size:13px;margin:.2em 0 1em;max-width:82ch}
.mosaic .track{position:relative;height:50px;border-radius:8px;overflow:hidden;border:1px solid var(--line);
  background:repeating-linear-gradient(90deg,transparent 0 3px,rgba(0,0,0,.015) 3px 6px),var(--bb);
  box-shadow:inset 0 1px 2px rgba(0,0,0,.05)}
.mosaic .seg{position:absolute;top:0;bottom:0;min-width:1px;box-shadow:0 0 0 1px rgba(255,255,255,.25) inset}
.mosaic .axis{display:flex;justify-content:space-between;margin-top:6px;font-family:var(--mono);font-size:11px;color:var(--faint)}
.mosaic .gap{position:absolute;top:0;bottom:0;min-width:1px;border-left:1px solid #64748b;border-right:1px solid #64748b;
  background:repeating-linear-gradient(45deg,rgba(100,116,139,.42) 0 4px,rgba(100,116,139,.14) 4px 8px)}
.mosaic .legend{display:flex;flex-wrap:wrap;gap:16px;margin-top:13px;font-size:13px}
.mosaic .leg{display:inline-flex;align-items:center;color:var(--muted)}
.mosaic .leg .hatch{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:6px;
  background:repeating-linear-gradient(45deg,#64748b 0 3px,#cbd5e1 3px 6px)}
.caveat{display:flex;gap:11px;align-items:flex-start;background:#fff7ed;border:1px solid #fed7aa;
  border-left:3px solid #dd6b20;border-radius:10px;padding:13px 16px;margin:20px 0 0;font-size:14px;color:#7c4a13;max-width:76ch}
.caveat .ic{font-weight:700;color:#dd6b20;font-size:16px;line-height:1.3}
.flag{display:inline-block;margin-left:7px;font-size:10.5px;color:#92600a;border:1px solid #e6c98a;
  background:#fdf6e7;border-radius:4px;padding:0 5px;font-weight:600;letter-spacing:.02em;vertical-align:middle}
.scroll{overflow-x:auto}
table.table{border-collapse:collapse;width:100%;font-size:14px}
table.table th{font-size:11px;font-weight:600;letter-spacing:.05em;text-transform:uppercase;color:var(--faint);
  text-align:right;padding:8px 14px;border-bottom:2px solid var(--line);white-space:nowrap}
table.table td{padding:9px 14px;border-bottom:1px solid var(--line);text-align:right}
table.table th:first-child,table.table td:first-child{text-align:left}
table.table .num{font-family:var(--mono);font-variant-numeric:tabular-nums;white-space:nowrap}
table.table .lbl{text-align:left;font-weight:500;white-space:nowrap}
table.table .strong{font-weight:650;color:var(--ink)}
table.table tr:hover td{background:var(--panel)}
table.table tr.hl td{background:var(--hl);font-weight:600}
.empty{color:var(--muted);font-style:italic}
.bars{display:flex;flex-direction:column;gap:7px;max-width:780px}
.barrow{display:grid;grid-template-columns:190px 1fr 60px;align-items:center;gap:12px}
.blabel{display:inline-flex;align-items:center;font-size:13px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.btrack{height:14px;background:var(--panel);border:1px solid var(--line);border-radius:7px;overflow:hidden}
.bfill{display:block;height:100%}
.bnum{font-family:var(--mono);font-size:13px;color:var(--muted);text-align:right}
details.methods{border:1px solid var(--line);border-radius:12px;padding:2px 18px;background:var(--paper)}
details.methods summary{cursor:pointer;font-weight:600;padding:13px 0;list-style:none;font-size:14px}
details.methods summary::-webkit-details-marker{display:none}
details.methods summary::before{content:"+";display:inline-block;width:1.3em;color:var(--faint);font-family:var(--mono)}
details.methods[open] summary::before{content:"-"}
details.methods[open]{padding-bottom:18px}
.glossary{display:grid;grid-template-columns:max-content 1fr;gap:7px 20px;margin:4px 0 20px}
.glossary dt{font-weight:600;font-size:13px}
.glossary dd{margin:0;color:var(--muted);font-size:13px}
details.methods h3{font-size:11px;text-transform:uppercase;letter-spacing:.08em;color:var(--faint);margin:10px 0 6px}
table.kv{border-collapse:collapse}
table.kv th{text-align:left;font-weight:500;color:var(--muted);padding:3px 18px 3px 0;vertical-align:top;font-size:13px;white-space:nowrap}
table.kv td{padding:3px 0;font-size:13px}
code{font-family:var(--mono);font-size:12px;background:var(--panel);padding:1px 5px;border-radius:4px}
footer{margin-top:50px;padding-top:20px;border-top:1px solid var(--line);color:var(--faint);font-size:12.5px;
  display:flex;flex-direction:column;gap:6px}
footer code{font-size:11px}
summary:focus-visible{outline:2px solid var(--accent);outline-offset:3px}
@media (max-width:640px){.wrap{padding:32px 16px 56px}header h1{font-size:23px}.verdict{font-size:16px}
  .barrow{grid-template-columns:120px 1fr 50px}}
"""


def _fmt_int(n) -> str:
    return f"{int(round(float(n))):,}"


def _fmt_kb(n) -> str:
    n = float(n)
    return f"{n / 1000:.1f} kb" if abs(n) >= 1000 else f"{int(round(n))} bp"


def _swatch(color: str) -> str:
    return f'<span class="sw" style="background:{color}"></span>'


def _summary(
    result: WindowSimilarity, regions: list[Region], datasets: list[str]
) -> dict:
    """Headline numbers: query length, backbone, donor list, recombinant fraction."""
    qcum = result.query_cumulative
    query_len = int(qcum[-1]) if len(qcum) else 0
    if regions:
        major = regions[0].major_parent
    elif datasets:
        major = datasets[0]
    else:
        major = "n/a"
    present = [r for r in regions if not r.donor_absent]
    absent = [r for r in regions if r.donor_absent]
    recomb_bp = sum(max(0, r.query_end - r.query_start) for r in present)
    minors: list[str] = []
    for r in present:
        if r.minor_parent not in minors:
            minors.append(r.minor_parent)
    return {
        "query_len": query_len, "major": major,
        "n_regions": len(present), "n_absent": len(absent),
        "recomb_bp": recomb_bp,
        "pct": (100.0 * recomb_bp / query_len) if query_len else 0.0,
        "minors": minors,
    }


def _absent_clause(s: dict) -> str:
    if not s["n_absent"]:
        return ""
    word = "region" if s["n_absent"] == 1 else "regions"
    return (f' {s["n_absent"]} further {word} diverge from every reference &mdash; a '
            f'<strong>donor may be missing</strong> (see Reference coverage).')


def _verdict_html(s: dict, query: str, colors: dict[str, str]) -> str:
    major = html.escape(s["major"])
    if s["n_regions"] == 0:
        lead = (f'No recombination among the present references &mdash; the query is most '
                f'similar to <strong>{major}</strong> throughout.')
        return f'<p class="verdict">{lead}{_absent_clause(s)}</p>'
    donors = ", ".join(
        f'{_swatch(colors.get(m, GREY))}<strong>{html.escape(m)}</strong>' for m in s["minors"]
    )
    word = "region" if s["n_regions"] == 1 else "regions"
    return (
        f'<p class="verdict">The query is a <strong>recombinant</strong>: a '
        f'<strong>{major}</strong> backbone carrying {s["n_regions"]} donor {word} from '
        f'{donors}, covering <span class="mono">{_fmt_kb(s["recomb_bp"])}</span> '
        f'(<span class="mono">{s["pct"]:.1f}%</span>) of the query.{_absent_clause(s)}</p>'
    )


def _cards_html(s: dict, colors: dict[str, str]) -> str:
    donors = "".join(
        f'{_swatch(colors.get(m, GREY))}{html.escape(m)} ' for m in s["minors"]
    ) or "&mdash;"
    regions_val = f'<span class="big">{s["n_regions"]}</span>'
    if s["n_absent"]:
        regions_val += f'<span class="sub">&nbsp;+{s["n_absent"]} donor-absent</span>'
    cards = [
        ("Recombinant regions", regions_val),
        ("Backbone &middot; major parent",
         f'{_swatch(colors.get(s["major"], GREY))}{html.escape(s["major"])}'),
        ("Donor &middot; minor parent", donors),
        ("Query recombinant",
         f'<span class="big">{s["pct"]:.1f}%</span>'
         f'<span class="sub">&nbsp;{_fmt_kb(s["recomb_bp"])} of {_fmt_kb(s["query_len"])}</span>'),
    ]
    items = "".join(
        f'<div class="card"><div class="k">{k}</div><div class="v">{v}</div></div>'
        for k, v in cards
    )
    return f'<div class="cards">{items}</div>'


def _mosaic_html(
    regions: list[Region], colors: dict[str, str], s: dict, gaps: list[CoverageGap]
) -> str:
    query_len = s["query_len"] or 1
    backbone = colors.get(s["major"], GREY)
    segs = "".join(
        f'<div class="seg" style="left:{100.0 * r.query_start / query_len:.3f}%;'
        f'width:{100.0 * max(0, r.query_end - r.query_start) / query_len:.3f}%;'
        f'background:{colors.get(r.minor_parent, GREY)}" '
        f'title="{html.escape(r.minor_parent)}: '
        f'{_fmt_int(r.query_start)}-{_fmt_int(r.query_end)} bp"></div>'
        for r in regions if not r.donor_absent  # absent-donor regions show as gaps
    )
    segs += "".join(
        f'<div class="gap" style="left:{100.0 * g.query_start / query_len:.3f}%;'
        f'width:{100.0 * max(0, g.query_end - g.query_start) / query_len:.3f}%" '
        f'title="low coverage ({g.kind}): closest reference {html.escape(g.best_label)} '
        f'~{g.mean_best:.2f}"></div>'
        for g in gaps
    )
    axis = "".join(
        f"<span>{_fmt_kb(frac * query_len)}</span>" for frac in (0, 0.25, 0.5, 0.75, 1.0)
    )
    legend = f'<span class="leg">{_swatch(backbone)}{html.escape(s["major"])} (backbone)</span>'
    legend += "".join(
        f'<span class="leg">{_swatch(colors.get(m, GREY))}{html.escape(m)} (donor)</span>'
        for m in s["minors"]
    )
    if gaps:
        legend += '<span class="leg"><span class="hatch"></span>low coverage</span>'
    return (
        f'<div class="mosaic"><div class="track" style="--bb:{backbone}">{segs}</div>'
        f'<div class="axis">{axis}</div><div class="legend">{legend}</div>'
        f'<p class="cap">The query genome painted by its closest reference per window: solid '
        f'backbone is the major parent; coloured segments are donor regions; hatched bands are '
        f'low-coverage stretches where even the closest reference is a poor match.</p></div>'
    )


def _regions_html(regions: list[Region], colors: dict[str, str], query_len: int) -> str:
    if not regions:
        return '<p class="empty">No recombinant regions were called.</p>'
    head = (
        "<tr><th>Donor (minor)</th><th>Backbone (major)</th><th>Query span (bp)</th>"
        "<th>Length</th><th>% query</th><th>Sim donor</th><th>Sim backbone</th>"
        "<th>Support</th><th>q-value</th><th>Breakpoint</th></tr>"
    )
    rows = ""
    for r in regions:
        qlen = max(0, r.query_end - r.query_start)
        pct = (100.0 * qlen / query_len) if query_len else 0.0
        if r.donor_absent:
            swatch = '<span class="sw hatch"></span>'
            donor = '<em>missing</em>'
            flag = ('<span class="flag" title="query is far from every reference here; '
                    'the true donor is likely absent">donor absent</span>')
        else:
            swatch = _swatch(colors.get(r.minor_parent, GREY))
            donor = html.escape(r.minor_parent)
            flag = ('<span class="flag" title="donor is itself a poor match">low conf</span>'
                    if r.donor_undercovered else "")
        support = "&ndash;" if r.support is None else f"{r.support:.2f}"
        qval = "&ndash;" if r.qvalue is None else f"{r.qvalue:.2g}"
        if r.breakpoint_lo is None:
            bp = "&ndash;"
        elif r.breakpoint_lo == r.breakpoint_hi:
            bp = _fmt_int(r.breakpoint_lo)
        else:
            bp = f"{_fmt_int(r.breakpoint_lo)}&ndash;{_fmt_int(r.breakpoint_hi)}"
        rows += (
            "<tr>"
            f'<td class="lbl">{swatch}{donor}{flag}</td>'
            f'<td class="lbl">{html.escape(r.major_parent)}</td>'
            f'<td class="num">{_fmt_int(r.query_start)}&ndash;{_fmt_int(r.query_end)}</td>'
            f'<td class="num">{_fmt_kb(qlen)}</td>'
            f'<td class="num">{pct:.1f}%</td>'
            f'<td class="num">{r.mean_sim_minor:.3f}</td>'
            f'<td class="num">{r.mean_sim_major:.3f}</td>'
            f'<td class="num strong">{support}</td>'
            f'<td class="num">{qval}</td>'
            f'<td class="num">{bp}</td>'
            "</tr>"
        )
    return f'<div class="scroll"><table class="table">{head}{rows}</table></div>'


def _winners_html(analysis: AnalysisResult, colors: dict[str, str]) -> str:
    items = sorted(analysis.winners_with_ties.items(), key=lambda x: x[1], reverse=True)
    if not items:
        return '<p class="empty">No window winners recorded.</p>'
    top = items[0][1] or 1
    bars = "".join(
        f'<div class="barrow"><span class="blabel">{_swatch(colors.get(label, GREY))}'
        f'{html.escape(label)}</span><span class="btrack">'
        f'<span class="bfill" style="width:{100.0 * count / top:.1f}%;'
        f'background:{colors.get(label, GREY)}"></span></span>'
        f'<span class="bnum">{_fmt_int(count)}</span></div>'
        for label, count in items
    )
    return f'<div class="bars">{bars}</div>'


def _stats_html(analysis: AnalysisResult, major: str) -> str:
    head = "".join(f"<th>{html.escape(h)}</th>" for h in analysis.stats_header)
    rows = ""
    for dataset, values in sorted(analysis.stats.items(), key=stats_sort_key, reverse=True):
        cls = ' class="hl"' if dataset == major else ""
        cells = f'<td class="lbl">{html.escape(dataset)}</td>' + "".join(
            f'<td class="num">{html.escape(str(v))}</td>' for v in values
        )
        rows += f"<tr{cls}>{cells}</tr>"
    return f'<div class="scroll"><table class="table"><tr>{head}</tr>{rows}</table></div>'


_GLOSSARY = [
    ("Major parent (backbone)",
     "The reference the query matches in the most windows overall."),
    ("Minor parent (donor)",
     "A reference the query matches better than the backbone over a stretch of the "
     "alignment -- a candidate recombination donor."),
    ("Similarity",
     "Per-window fraction of identical canonical bases (1.0 = identical). Windows with "
     "no comparable position are ignored."),
    ("Window / step",
     "The scan slides a fixed-width window along the alignment in fixed steps; each "
     "window is scored independently."),
    ("Support",
     "The share of distinguishing (discordant) sites -- where the query matches one "
     "candidate parent but not the other -- that favour the donor. 0.5 = no "
     "preference, 1.0 = every distinguishing site favours the donor."),
    ("q-value",
     "The sign-test p-value after Benjamini-Hochberg correction across all candidate "
     "segments (false-discovery-rate control). A region is reported when q <= alpha."),
    ("Breakpoint",
     "The query position where the source switches, with a posterior-derived "
     "uncertainty interval from the HMM."),
    ("Region calling (HMM)",
     "An HMM segments the query against the reference panel (a jump rate penalises "
     "switching reference); a segment is reported as recombinant only when its donor "
     "beats the major parent on the discordant sites by a sign test at level alpha."),
]


def _methods_html(provenance: dict[str, str]) -> str:
    glossary = "".join(
        f"<dt>{html.escape(t)}</dt><dd>{html.escape(d)}</dd>" for t, d in _GLOSSARY
    )
    params = "".join(
        f'<tr><th>{html.escape(str(k))}</th><td class="mono">{html.escape(str(v))}</td></tr>'
        for k, v in provenance.items()
    )
    return (
        '<details class="methods"><summary>Methods &amp; glossary</summary>'
        '<p class="cap">RecomFi segments the query against the reference panel with an HMM '
        '(jpHMM-style) and reports a region only when its donor beats the major parent on the '
        'sites that distinguish them (a sign test on discordant sites, immune to window '
        'overlap; Benjamini-Hochberg FDR across segments), with a posterior breakpoint '
        'interval. It remains an indicative screen, not a full phylogenetic test (e.g. 3SEQ, '
        'GARD) -- confirm strong candidates.</p>'
        f'<dl class="glossary">{glossary}</dl>'
        f'<h3>Run parameters</h3><table class="kv">{params}</table></details>'
    )


def _footer_html(provenance: dict[str, str]) -> str:
    files = [
        "recombination_regions.tsv", "coverage_gaps.tsv", "window_winners.tsv",
        "similarity_stats.tsv", "similarity_windows.tsv",
        "similarity_top*.pdf", "similarity_pair.pdf",
    ]
    flist = ", ".join(f"<code>{f}</code>" for f in files)
    ver = html.escape(provenance.get("recomfi version", ""))
    date = html.escape(provenance.get("date (UTC)", ""))
    return (
        f'<footer><div>Generated by RecomFi <span class="mono">{ver}</span> &middot; '
        f'<span class="mono">{date}</span> UTC</div>'
        f'<div>Companion files in this folder: {flist}.</div></footer>'
    )


def _caveat_html(gaps: list[CoverageGap], threshold: float) -> str:
    if not gaps:
        return ""
    n = len(gaps)
    word = "region" if n == 1 else "regions"
    verb = "is" if n == 1 else "are"
    total = sum(max(0, g.query_end - g.query_start) for g in gaps)
    worst = min(gaps, key=lambda g: g.mean_best)
    return (
        f'<div class="caveat"><span class="ic">&#9888;</span><div>'
        f'<strong>Possible missing reference.</strong> {n} {word} '
        f'(<span class="mono">{_fmt_kb(total)}</span> of the query) {verb} poorly covered: '
        f'even the closest reference stays below <span class="mono">{threshold:.2f}</span>. '
        f'The weakest is query <span class="mono">{_fmt_int(worst.query_start)}&ndash;'
        f'{_fmt_int(worst.query_end)}</span>, where <strong>{html.escape(worst.best_label)}'
        f'</strong> reaches only <span class="mono">{worst.mean_best:.2f}</span> &mdash; its '
        f'true source may not be in the collection. See <em>Reference coverage</em> below, or '
        f'run <span class="mono">recomfi find-references</span> to search NCBI.</div></div>'
    )


def _coverage_html(gaps: list[CoverageGap], threshold: float) -> str:
    intro = (
        f'<p class="cap">Stretches where even the closest reference is below the '
        f'<span class="mono">{threshold:.3f}</span> best-similarity threshold. '
        f'<strong>divergent</strong> = the query is genuinely far from every reference '
        f'(a likely missing reference); <strong>low information</strong> = too few comparable '
        f'bases to judge.</p>'
    )
    if not gaps:
        return (
            '<p class="cap">Every region has a close reference &mdash; no coverage gaps at the '
            f'<span class="mono">{threshold:.3f}</span> threshold.</p>'
        )
    head = (
        "<tr><th>Closest reference</th><th>Query span (bp)</th><th>Length</th>"
        "<th>Windows</th><th>Mean best sim</th><th>Kind</th></tr>"
    )
    rows = "".join(
        "<tr>"
        f'<td class="lbl">{html.escape(g.best_label)}</td>'
        f'<td class="num">{_fmt_int(g.query_start)}&ndash;{_fmt_int(g.query_end)}</td>'
        f'<td class="num">{_fmt_kb(g.length_bp)}</td>'
        f'<td class="num">{_fmt_int(g.n_windows)}</td>'
        f'<td class="num strong">{g.mean_best:.3f}</td>'
        f'<td class="lbl">{html.escape(g.kind.replace("_", " "))}</td>'
        "</tr>"
        for g in gaps
    )
    return f'{intro}<div class="scroll"><table class="table">{head}{rows}</table></div>'


def write_html_report(
    result: WindowSimilarity,
    analysis: AnalysisResult,
    regions: list[Region],
    datasets: list[str],
    provenance: dict[str, str],
    output_dir: Path,
    logger: logging.Logger,
    coverage_gaps: list[CoverageGap] | None = None,
    coverage_threshold: float = 0.0,
    extra_sections: list[tuple[str, str]] | None = None,
) -> Path:
    """Write a single self-contained ``report.html``."""
    gaps = coverage_gaps or []
    fig = build_interactive_figure(result, datasets, regions, gaps)
    plot_div = fig.to_html(full_html=False, include_plotlyjs="inline")

    colors = _color_map(datasets)
    s = _summary(result, regions, datasets)
    extras = "".join(
        f'<section class="section"><div class="eyebrow">{html.escape(title)}</div>{body}</section>'
        for title, body in (extra_sections or [])
    )

    doc = (
        '<!DOCTYPE html>\n<html lang="en"><head><meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '<link rel="icon" href="data:,">\n'
        f"<title>RecomFi report: {html.escape(result.query)}</title>\n"
        f"<style>{_CSS}</style></head><body><div class=\"wrap\">"
        '<header><div class="eyebrow">RecomFi &middot; recombination report</div>'
        f'<h1 class="mono">{html.escape(result.query)}</h1>'
        f"{_verdict_html(s, result.query, colors)}"
        f"{_caveat_html(gaps, coverage_threshold)}</header>"
        f"{_cards_html(s, colors)}"
        '<section class="section"><div class="eyebrow">Query mosaic</div>'
        f"{_mosaic_html(regions, colors, s, gaps)}</section>"
        '<section class="section"><div class="eyebrow">Recombinant regions</div>'
        f'{_regions_html(regions, colors, s["query_len"])}</section>'
        '<section class="section"><div class="eyebrow">Reference coverage</div>'
        f"{_coverage_html(gaps, coverage_threshold)}</section>"
        f"{extras}"
        '<section class="section"><div class="eyebrow">Similarity across the alignment</div>'
        '<p class="cap">Each line is one reference\'s similarity to the query along the '
        'alignment; coloured bands are called donor regions, hatched bands are low-coverage '
        'stretches. Drag to zoom, hover for values.</p>'
        f"{plot_div}</section>"
        '<section class="section"><div class="eyebrow">Window winners</div>'
        '<p class="cap">Windows in which each reference is the query\'s closest match '
        '(ties included).</p>'
        f"{_winners_html(analysis, colors)}</section>"
        '<section class="section"><div class="eyebrow">Per-dataset similarity statistics</div>'
        f'{_stats_html(analysis, s["major"])}</section>'
        f'<section class="section">{_methods_html(provenance)}</section>'
        f"{_footer_html(provenance)}"
        "</div></body></html>"
    )

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
    coverage_gaps: list[CoverageGap] | None = None,
    coverage_threshold: float = 0.0,
    extra_sections: list[tuple[str, str]] | None = None,
) -> None:
    """Write every table, plot and the HTML report for a completed scan."""
    output_dir.mkdir(parents=True, exist_ok=True)
    gaps = coverage_gaps or []

    write_windows_tsv(result, per_window_winners, output_dir, logger)
    write_stats_tsv(analysis, output_dir, logger)
    write_winners_tsv(analysis, output_dir, logger)
    write_regions_tsv(regions, output_dir, logger)
    write_coverage_tsv(gaps, coverage_threshold, output_dir, logger)

    top_datasets = rank_datasets(analysis, top_n)
    logger.info("Top %d nearest datasets: %s", len(top_datasets), ", ".join(top_datasets))
    plot_top_n(result, top_datasets, regions, output_dir, plot_format, logger)

    pair = rank_datasets(analysis, 2)
    plot_pairwise(result, pair, regions, output_dir, plot_format, logger)

    write_html_report(
        result, analysis, regions, top_datasets, provenance, output_dir, logger,
        coverage_gaps=gaps, coverage_threshold=coverage_threshold,
        extra_sections=extra_sections,
    )
