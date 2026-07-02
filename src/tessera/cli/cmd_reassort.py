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
from ..reassort.constellation import DEFAULT_MARGIN
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
    dataset: list[str] | None = typer.Option(
        None, "--dataset",
        help="Override a segment's Nextclade dataset, as SEGMENT=path (repeatable).",
    ),
    ani_floor: float = typer.Option(
        DEFAULT_ANI_FLOOR, "--ani-floor",
        help="A segment below this ANI to every reference tip is left unassigned.",
    ),
    margin: float = typer.Option(
        DEFAULT_MARGIN, "--margin",
        help="ANI window (percentage points) for a segment's near-best parents.",
    ),
    scan_segments: bool = typer.Option(
        False, "--scan-segments",
        help="Also scan each assigned segment for intragenic recombination (needs an aligner).",
    ),
    aligner: str = typer.Option(
        "mafft", "--aligner", help="Aligner backend for --scan-segments."),
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
            ani_floor=ani_floor, margin=margin, output=output,
            scan_segments=scan_segments, aligner=aligner, logger=logger,
        )

        output.mkdir(parents=True, exist_ok=True)
        group_of = {seg: i for i, g in enumerate(result.groups) for seg in g.segments}

        tsv = output / "reassortment.tsv"
        with open(tsv, "w") as fo:
            fo.write("segment\tdataset\tnearest_strain\tclade\tani\tstatus\tparent_group\n")
            for s in result.segments:
                grp = group_of.get(s.segment, "")
                fo.write(f"{s.segment}\t{s.dataset}\t{s.strain or ''}\t{s.clade or ''}\t"
                         f"{s.ani:.1f}\t{s.status}\t{grp}\n")

        ctsv = output / "constellation.tsv"
        with open(ctsv, "w") as fo:
            fo.write("group_index\tsegments\tparent_strains\n")
            for i, g in enumerate(result.groups):
                fo.write(f"{i}\t{','.join(g.segments)}\t{','.join(g.parent_strains)}\n")

        mosaic = " | ".join(
            f"group {i} {{{','.join(g.segments)}}} = {','.join(g.parent_strains) or '?'}"
            for i, g in enumerate(result.groups)
        )
        logger.info("Reassortment verdict: %s", result.verdict.upper())
        logger.info("Parent constellation: %s", mosaic or "(no segments assigned)")
        if result.verdict == "undetermined" and result.pair_notes:
            logger.info("Why undetermined: %s", "; ".join(result.pair_notes))
        logger.info("Wrote %s and %s", tsv, ctsv)

        if result.scans:
            stsv = output / "segment_scan.tsv"
            with open(stsv, "w") as fo:
                fo.write("segment\tintragenic_recombination\tn_regions\tnote\n")
                for sc in result.scans:
                    flag = "yes" if sc.recombinant else ("no" if sc.scanned else "n/a")
                    fo.write(f"{sc.segment}\t{flag}\t{sc.n_regions}\t{sc.note}\n")
            rollup = " | ".join(f"{sc.segment}: {sc.note}" for sc in result.scans)
            logger.info("Intragenic scan: %s", rollup)
            logger.info("Wrote %s", stsv)
