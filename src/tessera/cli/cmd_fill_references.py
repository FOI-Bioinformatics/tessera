"""The ``tessera fill-references`` command: iteratively add missing references."""

from __future__ import annotations

from pathlib import Path

import typer

from .main import _require_choice, app, get_logger, stage_errors


@app.command(name="fill-references")
def fill_references(
    query: Path = typer.Option(..., "-q", "--query", help="Path to the query FASTA."),
    collection: Path | None = typer.Option(
        None, "-c", "--collection",
        help="Starting reference collection directory (left untouched; a copy is grown). "
        "Omit to start fresh: the collection is seeded from a whole-query NCBI search.",
    ),
    output: Path = typer.Option(
        ..., "-o", "--output", help="Output directory (expanded collection, reports, summary)."
    ),
    aligner: str = typer.Option("mafft", "--aligner", help="Aligner backend for each round."),
    reference: str | None = typer.Option(
        None, "--reference", help="Backbone genome (label or filename)."
    ),
    max_rounds: int = typer.Option(3, "--max-rounds", help="Maximum search/download rounds."),
    seed_source: str = typer.Option(
        "blast", "--seed-source",
        help="Fresh-start source (no --collection): 'blast' (NCBI nt), 'local' (a genome pool "
        "directory, --candidate-pool), or 'ncbi-virus' (a taxon-scoped set, --taxon).",
    ),
    candidate_pool: Path | None = typer.Option(
        None, "--candidate-pool",
        help="Genome pool directory for --seed-source local (regional selection picks parents).",
    ),
    taxon: str | None = typer.Option(
        None, "--taxon",
        help="Taxon for --seed-source ncbi-virus (e.g. 'HIV-1'). Auto-detected if omitted.",
    ),
    nextclade_dataset: str | None = typer.Option(
        None, "--nextclade-dataset",
        help="Nextclade dataset path for --seed-source nextclade "
        "(e.g. nextstrain/sars-cov-2/XBB). Auto-detected from the query if omitted.",
    ),
    source_complete: bool = typer.Option(
        False, "--source-complete",
        help="NCBI Virus: fetch all complete genomes (then dereplicate) instead of RefSeq only.",
    ),
    seed_keep_siblings: bool = typer.Option(
        False, "--seed-keep-siblings",
        help="Keep the query's siblings when selecting from a pool (do not drop near-identical "
        "genomes). Use for closely related parents with no masking sibling.",
    ),
    auto_diversify: bool = typer.Option(
        True, "--auto-diversify/--no-auto-diversify",
        help="When BLAST seeding finds only siblings (a saturated lineage), switch to "
        "recruiting parental lineages by NCBI Virus taxonomy diversity (needs the datasets CLI).",
    ),
    negative_lineage: bool = typer.Option(
        True, "--negative-lineage/--no-negative-lineage",
        help="On saturation, re-BLAST each region excluding the query's own lineage (from "
        "the hit titles) to surface the divergent parents -- works for mega-taxa, no download.",
    ),
    fetch_limit: int = typer.Option(
        2000, "--fetch-limit",
        help="Notify (and dereplicate locally) when a broad NCBI Virus fetch returns "
        "more than this many complete genomes; the full set is downloaded either way.",
    ),
    cache_dir: Path | None = typer.Option(
        None, "--cache-dir",
        help="Where to cache fetched NCBI Virus panels (default ~/.cache/tessera or "
        "$TESSERA_CACHE), so repeat runs skip the network.",
    ),
    seed_mode: str = typer.Option(
        "windowed", "--seed-mode",
        help="BLAST seeding (--seed-source blast): 'whole' (closest whole-genome relatives), "
        "'windowed' (per-window best hits), or 'parents' (suppress siblings to recruit the "
        "parental lineages of a recombinant query).",
    ),
    seed_hits: int = typer.Option(
        10, "--seed-hits", help="BLAST hits to keep per seed search (per window, or whole-query)."
    ),
    seed_window: int = typer.Option(
        1500, "--seed-window", help="Window width (bp) for windowed/parents seeding."
    ),
    window_size: int = typer.Option(1000, "--window-size", help="Sliding window width (columns)."),
    window_step: int = typer.Option(100, "--window-step", help="Sliding window step (columns)."),
    coverage_floor: float | None = typer.Option(
        None, "--coverage-floor", help="Absolute best-similarity floor (default: adaptive)."
    ),
    coverage_rel_drop: float = typer.Option(
        0.05, "--coverage-rel-drop", help="Adaptive coverage drop below typical best similarity."
    ),
    min_improvement: float = typer.Option(
        0.01, "--min-improvement",
        help="Stop once a round improves the worst gap's best similarity by less than this.",
    ),
    max_hits: int = typer.Option(5, "--max-hits", help="BLAST hits to keep per gap."),
    top_gaps: int = typer.Option(3, "--top-gaps", help="Search the N longest gaps each round."),
    subtile: int = typer.Option(
        400, "--subtile",
        help="Search gaps longer than this (bp) in overlapping sub-intervals, so a short "
        "divergent tract surfaces its own donor instead of being diluted (0 disables).",
    ),
    email: str | None = typer.Option(
        None, "--email", help="Contact email for NCBI (or set NCBI_EMAIL)."
    ),
    exclude: list[str] = typer.Option(
        [], "--exclude", help="Accession(s) to never add (repeatable; version-insensitive)."
    ),
    keep_self_hits: bool = typer.Option(
        False, "--keep-self-hits", help="Keep near-identical hits (the query's own record)."
    ),
    curate: bool = typer.Option(
        False, "--curate",
        help="Drop the query's siblings and dereplicate each round (needs skani/skDER).",
    ),
    sibling_margin: float = typer.Option(
        3.0, "--sibling-margin",
        help="A candidate is a sibling if its query-ANI beats the backbone's by this many %.",
    ),
    af_min: float = typer.Option(
        80.0, "--af-min", help="... over at least this %% of the query (whole-genome match)."
    ),
    derep_ani: float = typer.Option(
        99.0, "--derep-ani", help="skDER: collapse references at or above this ANI %%."
    ),
    report: bool = typer.Option(
        True, "--report/--no-report",
        help="Run recombination detection and write the report after building the panel. "
        "Use --no-report to stop at the panel and run 'tessera recomb' separately.",
    ),
    method: str = typer.Option(
        "hmm,3seq,maxchi,bootscan", "--method",
        help="Region caller(s) for the detection step: a comma-separated list of "
        "hmm/3seq/maxchi/bootscan/heuristic, or 'all'. Several run as an ensemble "
        "(default hmm,3seq,maxchi,bootscan).",
    ),
    pool_consensus: bool = typer.Option(
        False, "--pool-consensus/--no-pool-consensus",
        help="With a Nextclade pool (--seed-source nextclade), use one denoised consensus "
        "genome per clade (a stable per-lineage reference) instead of every tree tip.",
    ),
    organism: str | None = typer.Option(
        None, "--organism",
        help="Organism / species name for the report header; defaults to the taxon.",
    ),
    lineage_map: Path | None = typer.Option(
        None, "--lineage-map",
        help="TSV of reference genotypes (accession<TAB>genotype) to override the typed "
        "names mined from genome headers; written into the panel's lineages.tsv.",
    ),
    threads: int = typer.Option(4, "-t", "--threads", help="Aligner worker threads."),
) -> None:
    """Iteratively search NCBI and add references until coverage stops improving.

    Each round rebuilds the MSA, scans for coverage gaps, and downloads the best
    new reference per gap. Needs an aligner and Entrez Direct; contacts NCBI.
    """
    import os

    from ..aligners.base import registry as aligner_registry
    from ..discover.iterate import FillParams, fill_references
    from ..recomb.regions import parse_methods

    logger = get_logger()
    with stage_errors(logger):
        _require_choice(aligner, set(aligner_registry.names()), "--aligner")
        _require_choice(seed_mode, {"whole", "windowed", "parents"}, "--seed-mode")
        _require_choice(seed_source, {"blast", "local", "ncbi-virus", "nextclade"}, "--seed-source")
        params = FillParams(
            query=query, collection=collection, output=output,
            aligner=aligner, reference=reference, max_rounds=max_rounds,
            seed_source=seed_source, candidate_pool=candidate_pool,
            taxon=taxon, nextclade_dataset=nextclade_dataset,
            source_refseq=not source_complete, seed_keep_siblings=seed_keep_siblings,
            auto_diversify=auto_diversify, negative_lineage=negative_lineage,
            fetch_limit=fetch_limit, cache_dir=cache_dir,
            seed_mode=seed_mode, seed_hits=seed_hits, seed_window=seed_window,
            window_size=window_size, window_step=window_step,
            coverage_floor=coverage_floor, coverage_rel_drop=coverage_rel_drop,
            min_improvement=min_improvement, max_hits=max_hits, top_gaps=top_gaps, subtile=subtile,
            email=email or os.environ.get("NCBI_EMAIL"),
            exclude=tuple(exclude), keep_self_hits=keep_self_hits, threads=threads,
            curate=curate, sibling_margin=sibling_margin, af_min=af_min, derep_ani=derep_ani,
            report=report, methods=parse_methods(method), pool_consensus=pool_consensus,
            organism=organism, lineage_map=lineage_map,
        )
        fill_references(params, logger)
