"""The ``tessera detect`` command: one-shot recombination detection from a query alone.

Detects the query's taxon, recruits a diverse reference panel from NCBI (no genomes
supplied by the user), aligns, and calls recombination with the sibling- and
lineage-aware caller -- a thin detection-tuned preset over ``fill-references``.
"""

from __future__ import annotations

from pathlib import Path

import typer

from .main import _require_choice, app, get_logger, stage_errors


def _seed_source(candidate_pool, nextclade: bool, nextclade_dataset: str | None) -> str:
    """Local pool wins, then Nextclade, then the default BLAST recruitment."""
    if candidate_pool:
        return "local"
    if nextclade or nextclade_dataset:
        return "nextclade"
    return "blast"


@app.command(name="detect")
def detect(
    query: Path = typer.Option(
        ..., "-q", "--query", help="Path to the query FASTA (the only input needed)."
    ),
    output: Path = typer.Option(
        ..., "-o", "--output", help="Output directory (report, recruited panel, run log)."
    ),
    taxon: str | None = typer.Option(
        None, "--taxon", help="Taxon for NCBI Virus recruitment (auto-detected if omitted)."
    ),
    aligner: str = typer.Option("mafft", "--aligner", help="Aligner backend."),
    window_size: int = typer.Option(1000, "--window-size", help="Sliding window width (columns)."),
    window_step: int = typer.Option(100, "--window-step", help="Sliding window step (columns)."),
    max_rounds: int = typer.Option(2, "--max-rounds", help="Maximum recruit/align rounds."),
    email: str | None = typer.Option(
        None, "--email", help="Contact email for NCBI (or set NCBI_EMAIL)."
    ),
    candidate_pool: Path | None = typer.Option(
        None, "--candidate-pool",
        help="Use a local genome pool instead of fetching from NCBI (for mega-taxa).",
    ),
    nextclade: bool = typer.Option(
        False, "--nextclade/--no-nextclade",
        help="Recruit the panel from a Nextclade dataset auto-detected from the query.",
    ),
    nextclade_dataset: str | None = typer.Option(
        None, "--nextclade-dataset",
        help="Nextclade dataset path (implies --nextclade; e.g. nextstrain/sars-cov-2/XBB).",
    ),
    cache_dir: Path | None = typer.Option(
        None, "--cache-dir", help="Where to cache fetched panels (so repeat runs are fast)."
    ),
    lineage_map: Path | None = typer.Option(
        None, "--lineage-map",
        help="TSV of reference genotypes (accession<TAB>genotype) to override the typed "
        "names mined from genome headers; the report names parents by lineage.",
    ),
    method: str = typer.Option(
        "hmm,3seq", "--method",
        help="Region caller(s): a comma-separated list of hmm/3seq/maxchi/bootscan/"
        "heuristic, or 'all'. Several run as an ensemble and their regions are merged "
        "(default hmm,3seq).",
    ),
    pool_consensus: bool = typer.Option(
        False, "--pool-consensus/--no-pool-consensus",
        help="With a Nextclade pool, use one denoised consensus genome per clade (a "
        "stable per-lineage reference) instead of every tree tip.",
    ),
    threads: int = typer.Option(4, "-t", "--threads", help="Aligner worker threads."),
) -> None:
    """Detect recombination in a query with no reference genomes supplied.

    Recruits the parental lineages from NCBI (negative-lineage BLAST + taxonomy
    diversity), drops the query's siblings, and competes lineages, so a recombinant
    query whose own lineage is common in NCBI is not masked. Needs an aligner and
    Entrez Direct; for a heavily sequenced taxon (e.g. SARS-CoV-2) supply
    ``--candidate-pool``.
    """
    import os

    from ..aligners.base import registry as aligner_registry
    from ..discover.iterate import FillParams, fill_references
    from ..recomb.regions import parse_methods

    logger = get_logger(output)
    with stage_errors(logger):
        _require_choice(aligner, set(aligner_registry.names()), "--aligner")
        params = FillParams(
            query=query, collection=None, output=output,
            aligner=aligner, max_rounds=max_rounds,
            window_size=window_size, window_step=window_step,
            email=email or os.environ.get("NCBI_EMAIL"),
            threads=threads, cache_dir=cache_dir,
            # Detection-tuned preset: recruit parents, drop siblings, curate the panel.
            seed_source=_seed_source(candidate_pool, nextclade, nextclade_dataset),
            nextclade_dataset=nextclade_dataset,
            candidate_pool=candidate_pool, taxon=taxon,
            seed_mode="parents", curate=True,
            auto_diversify=True, negative_lineage=True,
            methods=parse_methods(method), pool_consensus=pool_consensus,
            lineage_map=lineage_map,
        )
        fill_references(params, logger)
