"""The ``recomfi detect`` command: one-shot recombination detection from a query alone.

Detects the query's taxon, recruits a diverse reference panel from NCBI (no genomes
supplied by the user), aligns, and calls recombination with the sibling- and
lineage-aware caller -- a thin detection-tuned preset over ``fill-references``.
"""

from __future__ import annotations

from pathlib import Path

import typer

from .main import _require_choice, app, get_logger, stage_errors


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
    cache_dir: Path | None = typer.Option(
        None, "--cache-dir", help="Where to cache fetched panels (so repeat runs are fast)."
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
            seed_source="local" if candidate_pool else "blast",
            candidate_pool=candidate_pool, taxon=taxon,
            seed_mode="parents", curate=True,
            auto_diversify=True, negative_lineage=True,
        )
        fill_references(params, logger)
