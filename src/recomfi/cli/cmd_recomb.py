"""The ``recomfi recomb`` command: scan an MSA for recombination events."""

from __future__ import annotations

from pathlib import Path

import typer

from .main import _require_choice, app, get_logger, stage_errors


@app.command()
def recomb(
    msa: Path = typer.Option(
        ..., "-i", "--input", "--msa", help="Path to the multiple sequence alignment file."
    ),
    output: Path = typer.Option(..., "-o", "--output", help="Path to the output directory."),
    query: str = typer.Option(
        ..., "-q", "--query",
        help="Query label as it appears in the MSA (common sequence extensions are stripped).",
    ),
    window_size: int = typer.Option(1000, "--window-size", help="Sliding window width (columns)."),
    window_step: int = typer.Option(100, "--window-step", help="Sliding window step (columns)."),
    metric: str = typer.Option("pdist", "--metric", help="Similarity metric: pdist."),
    method: str = typer.Option(
        "hmm", "--method",
        help="Region caller: hmm (default, HMM segmentation + a discordant-site "
        "significance test) or heuristic (legacy margin/merge).",
    ),
    jump_rate: float = typer.Option(
        1e-3, "--jump-rate",
        help="HMM prior probability of switching reference per window (lower = "
        "fewer, longer segments).",
    ),
    alpha: float = typer.Option(
        0.05, "--alpha",
        help="Significance level for the donor-vs-major discordant-site test (hmm).",
    ),
    top_n: int = typer.Option(5, "--top-n", help="Number of nearest datasets to plot."),
    plot_format: str = typer.Option(
        "pdf", "--plot-format", help="Static plot format: pdf, png, or svg."
    ),
    min_region: int | None = typer.Option(
        None, "--min-region",
        help="Minimum recombinant region length in bp (default: the window size).",
    ),
    margin: float = typer.Option(
        0.0, "--margin",
        help="Minimum similarity a minor parent must exceed the major parent by to "
        "mark a window recombinant.",
    ),
    merge_gap: int | None = typer.Option(
        None, "--merge-gap",
        help="Merge same-parent regions separated by at most this many bp "
        "(default: the window size).",
    ),
    coverage_floor: float | None = typer.Option(
        None, "--coverage-floor",
        help="Flag windows whose closest reference is below this similarity as a "
        "possible missing reference (default: adaptive, relative to the query).",
    ),
    coverage_rel_drop: float = typer.Option(
        0.05, "--coverage-rel-drop",
        help="Adaptive coverage: flag windows more than this far below the query's "
        "typical best similarity (ignored when --coverage-floor is set).",
    ),
    exclude_siblings: bool = typer.Option(
        True, "--exclude-siblings/--no-exclude-siblings",
        help="Set aside the query's whole-genome siblings (its own lineage) before "
        "calling, so a sibling cannot win every window and mask recombination (hmm).",
    ),
    cluster_lineages: bool = typer.Option(
        True, "--cluster-lineages/--no-cluster-lineages",
        help="Pool near-duplicate references into lineages and compete the lineages, so "
        "duplicates do not tie every window and fragment the call (hmm).",
    ),
) -> None:
    """Identify recombination events in a multiple sequence alignment."""
    from ..recomb.run import RecombParams, run_recomb

    logger = get_logger()
    with stage_errors(logger):
        _require_choice(plot_format, {"pdf", "png", "svg"}, "--plot-format")
        _require_choice(method, {"hmm", "heuristic"}, "--method")
        params = RecombParams(
            msa=msa,
            output=output,
            query=query,
            window_size=window_size,
            window_step=window_step,
            metric=metric,
            top_n=top_n,
            plot_format=plot_format,
            method=method,
            jump_rate=jump_rate,
            alpha=alpha,
            min_region=min_region,
            margin=margin,
            merge_gap=merge_gap,
            coverage_floor=coverage_floor,
            coverage_rel_drop=coverage_rel_drop,
            exclude_siblings=exclude_siblings,
            cluster_lineages=cluster_lineages,
        )
        run_recomb(params, logger)
