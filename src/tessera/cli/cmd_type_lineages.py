"""The ``tessera type-lineages`` command: assign lineages to a genome collection.

Types an NCBI/local set of genomes (header/metadata, then nearest-neighbour to a
Nextclade dataset, then de-novo ANI clustering) and writes ``lineages.tsv``, so the
collection is lineage-ready for detection and panel building. Needs skani; the
nearest-neighbour step contacts the Nextclade dataset server on first use.
"""

from __future__ import annotations

import os
from pathlib import Path

import typer

from ..core.errors import UserInputError
from ..core.io import _require_fasta, collection_genomes
from .main import _require_directory, _require_lineage_map, app, get_logger, stage_errors


@app.command(name="type-lineages")
def type_lineages(
    collection: Path = typer.Option(
        ..., "-c", "--collection", help="Directory of genome FASTAs to type."
    ),
    output: Path = typer.Option(
        ..., "-o", "--output", help="Output directory (writes lineages.tsv)."
    ),
    taxon: str | None = typer.Option(
        None, "--taxon", help="Taxon for the Nextclade reference step (auto-detected if omitted)."
    ),
    nextclade_dataset: str | None = typer.Option(
        None, "--nextclade-dataset", help="Nextclade dataset path (else auto-detected)."
    ),
    lineage_map: Path | None = typer.Option(
        None, "--lineage-map", help="User TSV (accession<TAB>genotype) taking precedence."
    ),
    ref_ani_floor: float = typer.Option(
        90.0, "--ref-ani-floor", help="Min ANI to inherit a nearest reference tip's clade."
    ),
    cluster_ani: float = typer.Option(
        95.0, "--cluster-ani", help="De-novo: genomes at/above this ANI share a lineage."
    ),
    email: str | None = typer.Option(
        None, "--email", help="Contact email for NCBI (or set NCBI_EMAIL)."
    ),
) -> None:
    """Assign a lineage to every genome in a collection and write lineages.tsv."""
    from ..discover.lineage_assign import assign_lineages
    from ..recomb.typing import LINEAGES_TSV, write_lineage_map

    logger = get_logger(output)
    with stage_errors(logger):
        _require_directory(collection, "Collection directory")
        _require_lineage_map(lineage_map)
        # The same reading of a collection as every other command: each non-hidden
        # file is a genome, whatever its extension and gzip-compressed or not, and a
        # file that is not FASTA is an error rather than something to skip.
        genomes = collection_genomes(collection)
        if not genomes:
            raise UserInputError(f"No genome files found in {collection}")
        for genome in genomes:
            _require_fasta(genome)
        rows = assign_lineages(
            genomes, user_lineage_map=lineage_map, taxon=taxon,
            nextclade_dataset=nextclade_dataset, ref_ani_floor=ref_ani_floor,
            cluster_ani=cluster_ani, email=email or os.environ.get("NCBI_EMAIL"),
            logger=logger,
        )
        output.mkdir(parents=True, exist_ok=True)
        write_lineage_map(output / LINEAGES_TSV, rows)
        logger.info("Wrote %d lineage row(s): %s", len(rows), output / LINEAGES_TSV)
