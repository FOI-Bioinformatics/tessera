"""The ``recomfi msa`` command: build a reference-anchored MSA."""

from __future__ import annotations

from pathlib import Path

import typer

from .main import DEFAULT_THREADS, _parse_key_values, _require_choice, app, get_logger, stage_errors


@app.command()
def msa(
    query: Path = typer.Option(..., "-q", "--query", help="Path to the query file."),
    collection: Path = typer.Option(
        ..., "-c", "--collection", help="Path to the reference sequence collection directory."
    ),
    output: Path = typer.Option(
        ..., "-o", "--output", help="Path to the output MSA file (e.g. msa.fasta)."
    ),
    aligner: str = typer.Option(
        "sibeliaz", "--aligner",
        help="Aligner backend: sibeliaz (default), mafft, minimap2, progressivemauve, cactus.",
    ),
    reference: str | None = typer.Option(
        None, "--reference",
        help="Backbone genome (label or filename). Default: first in collection.",
    ),
    query_as_backbone: bool = typer.Option(
        False, "--query-as-backbone",
        help="Use the query as the MSA backbone. Only for a single-contig query.",
    ),
    aligner_arg: list[str] = typer.Option(
        [], "--aligner-arg",
        help="Aligner tuning as key=value (repeatable), e.g. seed_weight=11 "
        "(progressivemauve) or kmer=15 (sibeliaz).",
    ),
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads", help="Worker threads."),
) -> None:
    """Generate a multiple sequence alignment from a query and a sequence collection."""
    from ..aligners.base import registry as aligner_registry
    from ..msa.build import MsaParams, build_msa

    logger = get_logger()
    with stage_errors(logger):
        _require_choice(aligner, set(aligner_registry.names()), "--aligner")
        params = MsaParams(
            query=query,
            collection=collection,
            output=output,
            aligner=aligner,
            reference=reference,
            query_as_backbone=query_as_backbone,
            threads=threads,
            extra=_parse_key_values(aligner_arg, "--aligner-arg"),
        )
        build_msa(params, logger)
