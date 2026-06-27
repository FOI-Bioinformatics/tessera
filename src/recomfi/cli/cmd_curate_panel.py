"""The ``recomfi curate-panel`` command: build a diverse, sibling-free panel."""

from __future__ import annotations

from pathlib import Path

import typer

from .main import app, get_logger, stage_errors


@app.command(name="curate-panel")
def curate_panel(
    query: Path = typer.Option(..., "-q", "--query", help="Path to the query FASTA."),
    collection: Path = typer.Option(
        ..., "-c", "--collection",
        help="Reference collection directory (left untouched; a curated copy is written).",
    ),
    output: Path = typer.Option(
        ..., "-o", "--output", help="Output directory (curated collection + panel_lineages.tsv)."
    ),
    reference: str | None = typer.Option(
        None, "--reference",
        help="Backbone genome (label or filename). Default: the query's closest match.",
    ),
    sibling_margin: float = typer.Option(
        3.0, "--sibling-margin",
        help="A reference is a sibling if its query-ANI beats the backbone's by this many %.",
    ),
    af_min: float = typer.Option(
        80.0, "--af-min", help="... over at least this %% of the query (whole-genome match)."
    ),
    derep_ani: float = typer.Option(
        99.0, "--derep-ani", help="skDER: collapse references at or above this ANI %%."
    ),
) -> None:
    """Curate a collection into a diverse panel for recombination detection.

    Drops the query's own-lineage relatives (siblings), which would mask
    recombination, and dereplicates near-duplicates, keeping distinct parental
    lineages. Needs skani (and skDER for dereplication).
    """
    from ..discover.panel import curate_collection

    logger = get_logger()
    with stage_errors(logger):
        curate_collection(
            query, collection, output,
            reference=reference, ani_margin=sibling_margin,
            af_min=af_min, derep_ani=derep_ani, logger=logger,
        )
