"""The ``recomfi build-panel`` command: recruit a donor panel without detecting.

Builds (or grows) the query's donor/reference panel and stops before recombination
detection -- a panel-only counterpart to ``detect``. The donors are recruited with the
same parent-recruiting, sibling-dropping preset as ``detect``; the result is a curated
``collection/``, a stable ``panel.msa.fasta`` alignment, and ``panel_lineages.tsv``.
Run detection separately with ``recomfi recomb -i <output>/panel.msa.fasta``.
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


@app.command(name="build-panel")
def build_panel(
    query: Path = typer.Option(
        ..., "-q", "--query", help="Path to the query FASTA."
    ),
    output: Path = typer.Option(
        ..., "-o", "--output", help="Output directory (panel.msa.fasta, collection, run log)."
    ),
    collection: Path | None = typer.Option(
        None, "-c", "--collection",
        help="Starting reference collection to grow (left untouched; a copy is grown). "
        "Omit to recruit a panel from NCBI from the query alone.",
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
        "names mined from genome headers; written into the panel's lineages.tsv.",
    ),
    threads: int = typer.Option(4, "-t", "--threads", help="Aligner worker threads."),
) -> None:
    """Recruit a donor panel for the query and stop before detection.

    Uses the same recruitment as ``detect`` (recruit parental lineages, drop the
    query's siblings, curate the panel), but writes only the panel -- no report.
    The follow-up command to call recombination is logged on completion.
    """
    import os

    from ..aligners.base import registry as aligner_registry
    from ..discover.iterate import FillParams, fill_references

    logger = get_logger(output)
    with stage_errors(logger):
        _require_choice(aligner, set(aligner_registry.names()), "--aligner")
        params = FillParams(
            query=query, collection=collection, output=output,
            aligner=aligner, max_rounds=max_rounds,
            window_size=window_size, window_step=window_step,
            email=email or os.environ.get("NCBI_EMAIL"),
            threads=threads, cache_dir=cache_dir,
            # Same parent-recruiting preset as `detect`, but stop at the panel.
            seed_source=_seed_source(candidate_pool, nextclade, nextclade_dataset),
            nextclade_dataset=nextclade_dataset,
            candidate_pool=candidate_pool, taxon=taxon,
            seed_mode="parents", curate=True,
            auto_diversify=True, negative_lineage=True,
            report=False, lineage_map=lineage_map,
        )
        fill_references(params, logger)
