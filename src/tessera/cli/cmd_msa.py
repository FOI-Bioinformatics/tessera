"""The ``tessera msa`` command: build a reference-anchored MSA."""

from __future__ import annotations

from pathlib import Path

import typer

from ..core.errors import UserInputError
from .main import (
    DEFAULT_THREADS,
    _parse_key_values,
    _require_choice,
    _require_range,
    app,
    get_logger,
    stage_errors,
)


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
        _require_range(threads, "--threads", lo=1)
        # Both name the backbone. --reference used to win silently, so a user who asked
        # for the query as backbone got a different coordinate system without a word.
        if reference is not None and query_as_backbone:
            raise UserInputError(
                "--reference and --query-as-backbone both choose the backbone; give one."
            )
        if output.exists():
            logger.warning("Overwriting the existing alignment %s.", output)
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
