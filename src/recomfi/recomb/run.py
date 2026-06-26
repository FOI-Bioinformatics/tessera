"""Orchestrate the recombination scan: similarity -> analysis -> regions -> reports."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from .. import __version__
from ..core.io import strip_sequence_extension
from .analyze import analyze, winners_per_window
from .coverage import CoverageParams, call_coverage_gaps, flag_undercovered_regions
from .regions import RegionParams, call_regions
from .report import print_coverage, print_regions, print_summary, write_reports
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
    # Reference-coverage diagnostic.
    coverage_floor: float | None = None  # None -> adaptive baseline
    coverage_rel_drop: float = 0.05


def run_recomb(
    params: RecombParams,
    logger: logging.Logger,
    extra_sections: list[tuple[str, str]] | None = None,
) -> None:
    """Run the full recombination scan described by ``params``.

    ``extra_sections`` are ``(title, html)`` blocks inserted into the HTML report
    after the coverage section -- used by ``fill-references`` to show its rounds.
    """
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

    coverage_params = CoverageParams.with_defaults(
        params.window_size,
        floor=params.coverage_floor,
        rel_drop=params.coverage_rel_drop,
    )
    coverage_gaps, coverage_threshold = call_coverage_gaps(
        result, params.window_size, coverage_params
    )
    flag_undercovered_regions(regions, coverage_threshold)
    if coverage_gaps:
        logger.info(
            "Reference coverage: %d region(s) where the closest reference is below "
            "%.3f -- a better reference may be missing.",
            len(coverage_gaps), coverage_threshold,
        )

    # Summary + regions go to the logger (so they reach the run log) and stdout.
    print_summary(analysis, echo=logger.info)
    print_regions(regions, major_parent, echo=logger.info)
    print_coverage(coverage_gaps, coverage_threshold, echo=logger.info)

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
        "coverage threshold / gaps": f"{coverage_threshold:.3f} / {len(coverage_gaps)}",
    }

    output_dir = Path(params.output)
    logger.info("Writing outputs to %s", output_dir)
    write_reports(
        result, analysis, regions, per_window_winners, provenance, output_dir,
        top_n=params.top_n, plot_format=params.plot_format, logger=logger,
        coverage_gaps=coverage_gaps, coverage_threshold=coverage_threshold,
        extra_sections=extra_sections,
    )
    logger.info("All done.")
