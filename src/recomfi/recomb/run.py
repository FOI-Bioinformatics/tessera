"""Orchestrate the recombination scan: similarity -> analysis -> regions -> reports."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from .. import __version__
from ..core.io import strip_sequence_extension
from .analyze import analyze, winners_per_window
from .regions import RegionParams, call_regions
from .report import print_regions, print_summary, write_reports
from .similarity import compute_similarity


@dataclass
class RecombParams:
    msa: Path
    output: Path
    query: str
    window_size: int = 1000
    window_step: int = 100
    metric: str = "pdist"
    top_n: int = 5
    plot_format: str = "pdf"
    # Region calling (None -> derive from the window size).
    min_region: int | None = None
    margin: float = 0.0
    merge_gap: int | None = None


def run_recomb(params: RecombParams, logger: logging.Logger) -> None:
    """Run the full recombination scan described by ``params``."""
    query_label = strip_sequence_extension(params.query)

    logger.info(
        "Scanning MSA %s for query '%s' (window=%d, step=%d, metric=%s)",
        params.msa, query_label, params.window_size, params.window_step, params.metric,
    )
    result = compute_similarity(
        str(params.msa),
        query_label,
        window_size=params.window_size,
        window_step=params.window_step,
        metric=params.metric,
    )
    analysis = analyze(result)
    per_window_winners = winners_per_window(result)

    region_params = RegionParams.with_defaults(
        params.window_size,
        min_region=params.min_region,
        margin=params.margin,
        merge_gap=params.merge_gap,
    )
    regions, major_parent = call_regions(result, analysis, params.window_size, region_params)
    logger.info(
        "Major parent: %s; %d recombinant region(s) called.",
        major_parent or "n/a", len(regions),
    )

    # Summary + regions go to the logger (so they reach the run log) and stdout.
    print_summary(analysis, echo=logger.info)
    print_regions(regions, major_parent, echo=logger.info)

    provenance = {
        "recomfi version": __version__,
        "date (UTC)": datetime.now(UTC).strftime("%Y-%m-%d %H:%M:%S"),
        "query": query_label,
        "MSA": str(params.msa),
        "datasets": str(len(result.similarities)),
        "window / step": f"{params.window_size} / {params.window_step}",
        "metric": params.metric,
        "region min length / margin / merge gap":
            f"{region_params.min_region} / {region_params.margin} / {region_params.merge_gap}",
        "major parent": major_parent or "n/a",
    }

    output_dir = Path(params.output)
    logger.info("Writing outputs to %s", output_dir)
    write_reports(
        result, analysis, regions, per_window_winners, provenance, output_dir,
        top_n=params.top_n, plot_format=params.plot_format, logger=logger,
    )
    logger.info("All done.")
