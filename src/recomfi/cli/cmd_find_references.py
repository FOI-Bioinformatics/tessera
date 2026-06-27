"""The ``recomfi find-references`` command: find references missing from a collection."""

from __future__ import annotations

from pathlib import Path

import typer

from .main import app, get_logger, stage_errors


@app.command(name="find-references")
def find_references(
    msa: Path = typer.Option(
        ..., "-i", "--input", "--msa", help="Path to the multiple sequence alignment file."
    ),
    query: Path = typer.Option(
        ..., "-q", "--query", help="Query label as it appears in the MSA."
    ),
    output: Path = typer.Option(
        ..., "-o", "--output", help="Output directory for the candidate report."
    ),
    collection: Path | None = typer.Option(
        None, "-c", "--collection",
        help="Existing collection directory (so already-present hits are marked).",
    ),
    window_size: int = typer.Option(1000, "--window-size", help="Sliding window width (columns)."),
    window_step: int = typer.Option(100, "--window-step", help="Sliding window step (columns)."),
    coverage_floor: float | None = typer.Option(
        None, "--coverage-floor",
        help="Treat windows whose closest reference is below this similarity as gaps "
        "(default: adaptive, relative to the query).",
    ),
    coverage_rel_drop: float = typer.Option(
        0.05, "--coverage-rel-drop",
        help="Adaptive coverage drop below the query's typical best similarity.",
    ),
    max_hits: int = typer.Option(5, "--max-hits", help="BLAST hits to keep per gap."),
    top_gaps: int = typer.Option(3, "--top-gaps", help="Search the N longest gaps."),
    email: str | None = typer.Option(
        None, "--email", help="Contact email for NCBI (etiquette; or set NCBI_EMAIL)."
    ),
    exclude: list[str] = typer.Option(
        [], "--exclude",
        help="Accession(s) to drop from the candidates, e.g. the query's own record "
        "(repeatable; version-insensitive).",
    ),
    keep_self_hits: bool = typer.Option(
        False, "--keep-self-hits",
        help="Keep near-identical hits; by default a hit matching the query almost "
        "exactly (the query's own record) is auto-skipped.",
    ),
    download: Path | None = typer.Option(
        None, "--download",
        help="Download the best new reference per gap into this collection directory.",
    ),
    curate: bool = typer.Option(
        False, "--curate",
        help="After download, drop the query's siblings and dereplicate (needs skani/skDER, "
        "--collection as the backbone source).",
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
    """Find (and optionally download) reference genomes missing from the collection.

    Flags query regions the collection covers poorly, BLASTs them against NCBI nt,
    and reports candidate references. This contacts NCBI over the network.
    """
    import os

    from ..discover.run import FindRefParams, find_references

    logger = get_logger()
    with stage_errors(logger):
        params = FindRefParams(
            msa=msa, query=str(query), output=output, collection=collection,
            window_size=window_size, window_step=window_step,
            coverage_floor=coverage_floor, coverage_rel_drop=coverage_rel_drop,
            max_hits=max_hits, top_gaps=top_gaps,
            email=email or os.environ.get("NCBI_EMAIL"),
            exclude=tuple(exclude), keep_self_hits=keep_self_hits,
            download=download, curate=curate, sibling_margin=sibling_margin,
            af_min=af_min, derep_ani=derep_ani,
        )
        find_references(params, logger)
