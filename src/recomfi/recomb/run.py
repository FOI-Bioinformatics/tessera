"""Orchestrate the recombination scan: similarity -> analysis -> regions -> reports."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from .. import __version__
from ..core.io import strip_sequence_extension
from .analyze import analyze, winners_per_window
from .coverage import (
    CoverageParams,
    call_coverage_gaps,
    flag_undercovered_regions,
    gaps_as_regions,
)
from .hmm import DEFAULT_JUMP_RATE
from .regions import RegionParams, call_regions
from .report import print_coverage, print_regions, print_summary, write_reports
from .similarity import compute_similarity
from .typing import LINEAGES_TSV, LineageMap, lineage_of, load_lineage_map


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
    # Region calling. method="hmm" (default, statistical) or "heuristic" (legacy).
    method: str = "hmm"
    jump_rate: float = DEFAULT_JUMP_RATE  # HMM prior switch probability per window
    alpha: float = 0.05  # significance level for the donor-vs-major site test
    exclude_siblings: bool = True  # set aside the query's own-lineage siblings first
    cluster_lineages: bool = True  # pool near-duplicate references into lineages first
    # Heuristic-only thresholds (None -> derive from the window size).
    min_region: int | None = None
    margin: float = 0.0
    merge_gap: int | None = None
    # Reference-coverage diagnostic.
    coverage_floor: float | None = None  # None -> adaptive baseline
    coverage_rel_drop: float = 0.05
    # Typed genotype/lineage names for references. A map keyed by reference label takes
    # precedence; otherwise a lineages.tsv beside the output or the MSA is read.
    lineage_map: LineageMap | None = None
    lineage_map_path: Path | None = None
    # The query's own typed lineage (shown in the verdict). Defaults to the query's
    # entry in the lineage map when present.
    query_lineage: str | None = None


def _discover_lineage_tsv(output: Path, msa: Path) -> Path | None:
    """Find a ``lineages.tsv`` written by panel building, beside the output or the MSA."""
    candidates = [
        Path(output) / LINEAGES_TSV,
        Path(msa).parent / LINEAGES_TSV,
        Path(msa).with_suffix(".lineages.tsv"),
    ]
    return next((p for p in candidates if p.exists()), None)


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

    lineage_map = params.lineage_map
    if lineage_map is None:
        tsv = params.lineage_map_path or _discover_lineage_tsv(params.output, params.msa)
        lineage_map = load_lineage_map(tsv) or None
        if lineage_map:
            logger.info("Typed names for %d reference(s) from %s.", len(lineage_map), tsv)
    query_lineage = params.query_lineage or lineage_of(query_label, lineage_map)

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
        method=params.method,
        jump_rate=params.jump_rate,
        alpha=params.alpha,
        exclude_siblings=params.exclude_siblings,
        cluster_lineages=params.cluster_lineages,
    )
    regions, major_parent, excluded_siblings = call_regions(
        result, analysis, params.window_size, region_params
    )
    if excluded_siblings:
        logger.info(
            "Excluded %d whole-genome sibling(s) of the query (its own lineage) from the "
            "competition: %s.",
            len(excluded_siblings),
            ", ".join(f"{ev.label} (leads {ev.lead_frac:.0%})" for ev in excluded_siblings),
        )
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

    # Bridge: a divergent coverage gap (query far from every reference) is a
    # putative recombination whose donor is absent from the collection. Report
    # those that a confident donor-present region does not already cover.
    absent = [
        r for r in gaps_as_regions(coverage_gaps, result, major_parent)
        if not any(r.msa_start < p.msa_end and p.msa_start < r.msa_end for p in regions)
    ]
    regions = sorted(regions + absent, key=lambda r: r.msa_start)
    if absent:
        logger.info("Added %d donor-absent region(s) (likely a missing reference).",
                    len(absent))

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
        "caller": (
            f"hmm (jump-rate {params.jump_rate:g}, alpha {params.alpha:g})"
            if params.method == "hmm"
            else f"heuristic (min {region_params.min_region} / margin "
                 f"{region_params.margin} / merge {region_params.merge_gap})"
        ),
        "major parent": major_parent or "n/a",
        "coverage threshold / gaps": f"{coverage_threshold:.3f} / {len(coverage_gaps)}",
    }
    if excluded_siblings:
        provenance["excluded siblings (query's own lineage)"] = ", ".join(
            ev.label for ev in excluded_siblings
        )

    output_dir = Path(params.output)
    logger.info("Writing outputs to %s", output_dir)
    write_reports(
        result, analysis, regions, per_window_winners, provenance, output_dir,
        top_n=params.top_n, plot_format=params.plot_format, logger=logger,
        coverage_gaps=coverage_gaps, coverage_threshold=coverage_threshold,
        extra_sections=extra_sections, lineage_map=lineage_map,
        query_lineage=query_lineage,
    )
    logger.info("All done.")
