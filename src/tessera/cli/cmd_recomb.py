"""The ``tessera recomb`` command: scan an MSA for recombination events."""

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
        "hmm,3seq", "--method",
        help="Region caller(s), comma-separated, or 'all'. Several run as an ensemble "
        "and their regions are merged into a consensus (agreement raises confidence); "
        "the default is hmm,3seq. Callers: hmm (HMM segmentation + a discordant-site "
        "significance test), 3seq (scan-aware triplet max-drawdown test; strong at low "
        "divergence), maxchi (chi-square triplet test, complementary to 3seq), bootscan "
        "(distance + bootstrap support for the closest parent), heuristic (legacy "
        "margin/merge). Pass a single name (e.g. --method hmm) for one caller.",
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
    lineage_map: Path | None = typer.Option(
        None, "--lineage-map",
        help="TSV of reference genotypes (accession<TAB>genotype) to name parents by "
        "lineage in the report. Defaults to a 'lineages.tsv' beside the output or MSA.",
    ),
    informative_sites: bool | None = typer.Option(
        None, "--informative-sites/--no-informative-sites",
        help="Window over polymorphic (informative) sites instead of fixed base-pair "
        "windows -- recovers signal on near-identical panels (intra-species sets, DNA "
        "viruses: mpox, VZV, ebola). Default: auto, on when a base-pair window would "
        "hold too few discriminating sites (hmm).",
    ),
    informative_window: int = typer.Option(
        40, "--informative-window",
        help="Informative-site mode: polymorphic sites per window.",
    ),
    informative_step: int = typer.Option(
        5, "--informative-step",
        help="Informative-site mode: step in informative-site index space.",
    ),
    phi: bool = typer.Option(
        True, "--phi/--no-phi",
        help="Report the parent-free recombination signal (PHI test + Hudson-Kaplan "
        "Rmin) alongside the called regions. Independent of --method.",
    ),
    phi_window: int = typer.Option(
        100, "--phi-window",
        help="PHI test window width, in informative-site ranks.",
    ),
) -> None:
    """Identify recombination events in a multiple sequence alignment."""
    from ..recomb.regions import parse_methods
    from ..recomb.run import RecombParams, run_recomb

    logger = get_logger()
    with stage_errors(logger):
        _require_choice(plot_format, {"pdf", "png", "svg"}, "--plot-format")
        methods = parse_methods(method)
        params = RecombParams(
            msa=msa,
            output=output,
            query=query,
            window_size=window_size,
            window_step=window_step,
            metric=metric,
            top_n=top_n,
            plot_format=plot_format,
            methods=methods,
            jump_rate=jump_rate,
            alpha=alpha,
            min_region=min_region,
            margin=margin,
            merge_gap=merge_gap,
            coverage_floor=coverage_floor,
            coverage_rel_drop=coverage_rel_drop,
            exclude_siblings=exclude_siblings,
            cluster_lineages=cluster_lineages,
            informative_sites=informative_sites,
            informative_window=informative_window,
            informative_step=informative_step,
            phi=phi,
            phi_window=phi_window,
            lineage_map_path=lineage_map,
        )
        run_recomb(params, logger)
