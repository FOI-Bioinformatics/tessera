"""The ``tessera reassort`` command: per-segment reassortment detection.

Types each segment of a multi-FASTA query against its Nextclade dataset (nearest
reference tip by skani ANI) and calls reassortment when the segments trace to different
parent strains/lineages. Needs skani; the typing step contacts the Nextclade dataset
server on first use.
"""

from __future__ import annotations

import os
from pathlib import Path

import typer

from ..reassort import assign_segments
from ..reassort.assign import DEFAULT_ANI_FLOOR
from .main import app, get_logger, stage_errors


@app.command(name="reassort")
def reassort(
    query: Path = typer.Option(
        ..., "-q", "--query", help="Segmented query: a multi-FASTA, one record per segment."
    ),
    output: Path = typer.Option(
        ..., "-o", "--output", help="Output directory (writes reassortment.tsv)."
    ),
    email: str | None = typer.Option(
        None, "--email", help="Contact email for NCBI dataset auto-detection (or NCBI_EMAIL)."
    ),
    dataset: list[str] = typer.Option(
        None, "--dataset",
        help="Override a segment's Nextclade dataset, as SEGMENT=path (repeatable).",
    ),
    ani_floor: float = typer.Option(
        DEFAULT_ANI_FLOOR, "--ani-floor",
        help="A segment below this ANI to every reference tip is left unassigned.",
    ),
) -> None:
    """Detect reassortment: assign each segment to its nearest reference lineage and
    report the per-segment genotype plus a clonal/reassortant verdict."""
    logger = get_logger(output)
    with stage_errors(logger):
        overrides: dict[str, str] = {}
        for item in dataset or []:
            if "=" not in item:
                raise typer.BadParameter(f"--dataset must be SEGMENT=path, got {item!r}")
            seg, path = item.split("=", 1)
            overrides[seg.strip()] = path.strip()

        result = assign_segments(
            query, dataset_overrides=overrides,
            email=email or os.environ.get("NCBI_EMAIL"),
            ani_floor=ani_floor, logger=logger,
        )

        output.mkdir(parents=True, exist_ok=True)
        tsv = output / "reassortment.tsv"
        with open(tsv, "w") as fo:
            fo.write("segment\tdataset\tnearest_strain\tclade\tani\tstatus\n")
            for s in result.segments:
                fo.write(f"{s.segment}\t{s.dataset}\t{s.strain or ''}\t{s.clade or ''}\t"
                         f"{s.ani:.1f}\t{s.status}\n")

        mosaic = " | ".join(
            f"{s.segment}:{s.clade} ({s.strain})" if s.status == "assigned"
            else f"{s.segment}:unassigned"
            for s in result.segments
        )
        logger.info("Reassortment verdict: %s", result.verdict.upper())
        logger.info("Per-segment genotype: %s", mosaic)
        logger.info("Wrote %s", tsv)
