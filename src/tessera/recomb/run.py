"""Orchestrate the recombination scan: similarity -> analysis -> regions -> reports."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from .. import __version__
from ..core.errors import UserInputError
from ..core.io import strip_sequence_extension
from .analyze import analyze, winners_per_window
from .coverage import (
    CoverageParams,
    call_coverage_gaps,
    flag_undercovered_regions,
    gaps_as_regions,
)
from .diagnostics import recombination_signal
from .ensemble import consensus_regions, reconcile_major
from .hmm import DEFAULT_JUMP_RATE
from .regions import DEFAULT_METHODS, RegionParams, call_regions
from .report import print_coverage, print_regions, print_summary, write_reports
from .similarity import (
    compute_similarity,
    compute_similarity_informative,
    informative_site_count,
)
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
    # Region calling. The default ensemble runs the hmm and 3seq callers and merges
    # their regions into a consensus; a single-method tuple reproduces one caller.
    methods: tuple[str, ...] = DEFAULT_METHODS
    jump_rate: float = DEFAULT_JUMP_RATE  # HMM prior switch probability per window
    alpha: float = 0.05  # significance level for the donor-vs-major site test
    exclude_siblings: bool = True  # set aside the query's own-lineage siblings first
    cluster_lineages: bool = True  # pool near-duplicate references into lineages first
    # Informative-site windowing for near-identical panels (intra-species sets, DNA
    # viruses: mpox, VZV, ebola). None -> auto: switch when a base-pair window would
    # hold fewer than min_info_per_window polymorphic sites. True/False force it.
    informative_sites: bool | None = None
    informative_window: int = 40  # polymorphic sites per window
    informative_step: int = 5  # step in informative-site index space
    # Parent-free recombination signal (PHI test + Hudson-Kaplan Rmin), reported
    # alongside the parent-attributed regions. Independent of --method; cheap.
    phi: bool = True
    phi_window: int = 100  # PHI window width, in informative-site ranks
    # Auto-trigger: switch to informative-site windowing when the references differ at
    # fewer than this fraction of columns (low inter-reference divergence, so a
    # base-pair window dilutes the few discriminating sites). 8% sits below the
    # genotype-level panels that work in base-pair mode (measles ~10%, HIV ~37%).
    informative_max_fraction: float = 0.08
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


def _select_windowing(bp_result, params: RecombParams, query_label: str, logger):
    """Choose base-pair or informative-site windowing; return ``(result, label)``.

    Informative-site windowing is for near-identical panels where a base-pair window
    holds too few polymorphic columns. It is auto-enabled (``informative_sites=None``)
    when the expected informative sites per base-pair window fall below
    ``min_info_per_window``, and only for the HMM caller (the heuristic path assumes
    base-pair window geometry). Falls back to base-pair windowing if the panel has
    fewer than ``informative_window`` informative sites.
    """
    forced = params.informative_sites
    # Informative-site windowing is the HMM caller's; it only applies when hmm runs.
    if "hmm" not in params.methods:
        if forced:
            logger.warning(
                "--informative-sites applies to the hmm caller; using base-pair windowing."
            )
        return bp_result, "base-pair"

    n_info = informative_site_count(bp_result.rows, query_label)
    fraction = n_info / bp_result.width if bp_result.width else 0.0
    use_informative = forced
    if use_informative is None:
        use_informative = (
            n_info >= params.informative_window
            and fraction < params.informative_max_fraction
        )
    if not use_informative:
        return bp_result, "base-pair"

    try:
        result = compute_similarity_informative(
            str(params.msa), query_label,
            info_window=params.informative_window,
            info_step=params.informative_step,
            metric=params.metric,
        )
    except UserInputError as exc:
        logger.warning(
            "Informative-site windowing unavailable (%s); using base-pair windowing.", exc
        )
        return bp_result, "base-pair"

    spans = result.window_spans
    mean_bp = (spans[-1][1] - spans[0][0]) / max(len(result.positions), 1) if spans else 0
    logger.info(
        "Near-identical panel (%d informative site(s), %.2f%% of columns): using "
        "informative-site windowing -- %d window(s), ~%.0f bp/window.",
        n_info, 100.0 * fraction, len(result.positions), mean_bp,
    )
    return result, (
        f"informative-site ({params.informative_window} sites/window, "
        f"step {params.informative_step})"
    )


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
) -> str:
    """Run the full recombination scan described by ``params``; return the windowing
    label used ("base-pair" or "informative-site (...)").

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
    bp_result = compute_similarity(
        str(params.msa),
        query_label,
        window_size=params.window_size,
        window_step=params.window_step,
        metric=params.metric,
    )
    # On a near-identical panel a base-pair window holds too few discriminating
    # sites for the HMM emission to separate the references; switch to windows that
    # span a fixed number of informative (polymorphic) columns instead. Auto unless
    # forced; the base-pair result is kept for the coverage diagnostic either way.
    result, windowing = _select_windowing(bp_result, params, query_label, logger)

    # Parent-free recombination signal: does the alignment carry recombination at all,
    # and where -- computed once from every sequence, independent of the caller. Most
    # informative where divergence is low or the true donor is absent from the panel.
    signal = None
    if params.phi:
        signal = recombination_signal(
            bp_result.rows, query_label, bp_result.column_to_query,
            window=params.phi_window,
        )
        if signal is not None:
            logger.info(
                "Recombination signal (parent-free): PHI p=%.4g, Rmin=%d (%d informative "
                "sites).", signal.phi_p, signal.rmin, signal.n_informative,
            )

    analysis = analyze(result)
    per_window_winners = winners_per_window(result)
    # 3seq / heuristic assume base-pair window geometry; the HMM uses the selected
    # (possibly informative-site) result. Reuse the same analysis when they coincide.
    analysis_bp = analysis if result is bp_result else analyze(bp_result)

    def _region_params(method: str) -> RegionParams:
        return RegionParams.with_defaults(
            params.window_size,
            min_region=params.min_region,
            margin=params.margin,
            merge_gap=params.merge_gap,
            method=method,
            jump_rate=params.jump_rate,
            alpha=params.alpha,
            exclude_siblings=params.exclude_siblings,
            cluster_lineages=params.cluster_lineages,
        )

    # Run each selected caller (sharing the one similarity scan + analysis), then merge
    # their regions into a single consensus that records per-region method agreement.
    per_method: dict[str, list] = {}
    majors: dict[str, str | None] = {}
    excluded_siblings: list = []
    for method in params.methods:
        res, ana = (result, analysis) if method == "hmm" else (bp_result, analysis_bp)
        regs, major_m, sibs = call_regions(res, ana, params.window_size, _region_params(method))
        per_method[method] = regs
        majors[method] = major_m
        if method == "hmm":
            excluded_siblings = sibs

    major_parent, per_major = reconcile_major(majors)
    regions, method_breakdown = consensus_regions(
        per_method, major=major_parent,
        rmin_intervals=signal.rmin_intervals if signal else None,
        lineage_map=lineage_map,
    )
    if excluded_siblings:
        logger.info(
            "Excluded %d whole-genome sibling(s) of the query (its own lineage) from the "
            "competition: %s.",
            len(excluded_siblings),
            ", ".join(f"{ev.label} (leads {ev.lead_frac:.0%})" for ev in excluded_siblings),
        )
    n_agree = sum(1 for r in regions if len(r.methods) >= 2)
    logger.info(
        "Caller(s): %s -- major parent %s; %d recombinant region(s)%s.",
        ", ".join(params.methods), major_parent or "n/a", len(regions),
        f", {n_agree} called by >1 method" if len(params.methods) > 1 else "",
    )

    coverage_params = CoverageParams.with_defaults(
        params.window_size,
        floor=params.coverage_floor,
        rel_drop=params.coverage_rel_drop,
    )
    coverage_gaps, coverage_threshold = call_coverage_gaps(
        bp_result, params.window_size, coverage_params
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

    def _caller_desc(method: str) -> str:
        if method == "hmm":
            return f"hmm (jump-rate {params.jump_rate:g}, alpha {params.alpha:g})"
        if method == "3seq":
            return f"3seq (triplet max-descent test, alpha {params.alpha:g})"
        min_region = params.min_region if params.min_region is not None else params.window_size
        merge_gap = params.merge_gap if params.merge_gap is not None else params.window_size
        return f"heuristic (min {min_region} / margin {params.margin} / merge {merge_gap})"

    caller_desc = " + ".join(_caller_desc(m) for m in params.methods)
    provenance = {
        "tessera version": __version__,
        "date (UTC)": datetime.now(UTC).strftime("%Y-%m-%d %H:%M:%S"),
        "query": query_label,
        "MSA": str(params.msa),
        "datasets": str(len(result.similarities)),
        "window / step": f"{params.window_size} / {params.window_step}",
        "metric": params.metric,
        "caller": f"ensemble: {caller_desc}" if len(params.methods) > 1 else caller_desc,
        "windowing": windowing,
        "major parent": major_parent or "n/a",
        "coverage threshold / gaps": f"{coverage_threshold:.3f} / {len(coverage_gaps)}",
    }
    if excluded_siblings:
        provenance["excluded siblings (query's own lineage)"] = ", ".join(
            ev.label for ev in excluded_siblings
        )
    if signal is not None:
        provenance["recombination signal (PHI)"] = (
            f"p={signal.phi_p:.4g} ({signal.n_informative} informative sites, "
            f"window {signal.phi_window})"
        )
        provenance["min recombination events (Rmin)"] = str(signal.rmin)

    output_dir = Path(params.output)
    logger.info("Writing outputs to %s", output_dir)
    write_reports(
        result, analysis, regions, per_window_winners, provenance, output_dir,
        top_n=params.top_n, plot_format=params.plot_format, logger=logger,
        coverage_gaps=coverage_gaps, coverage_threshold=coverage_threshold,
        extra_sections=extra_sections, lineage_map=lineage_map,
        query_lineage=query_lineage, signal=signal,
        methods_run=params.methods, method_breakdown=method_breakdown, per_major=per_major,
    )
    logger.info("All done.")
    return windowing
