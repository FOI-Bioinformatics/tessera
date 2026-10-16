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

from ..core.errors import UserInputError
from ..core.io import read_fasta
from ..reassort import assign_segments
from ..reassort.assign import DEFAULT_ANI_FLOOR
from ..reassort.constellation import DEFAULT_MARGIN, MAX_MARGIN
from .main import _require_file, _require_range, app, get_logger, stage_errors


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
        help="A segment below this ANI (a percentage, 0-100) to every reference tip is "
        "left unassigned.",
    ),
    margin: float = typer.Option(
        DEFAULT_MARGIN, "--margin",
        help="ANI window in percentage points (above 0, at most 5) for a segment's "
        "near-best parents. Not the 0-1 similarity fraction 'recomb --margin' takes. A "
        "wide window lets every pair share a strain, so every query reads clonal.",
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
        _require_file(query, "Query file")
        # ANI is a percentage. Above 100 nothing can be assigned; a value in (0, 1] is
        # almost certainly a fraction meant as a percentage (0.8 for 80), and would
        # accept every tip skani reports.
        _require_range(ani_floor, "--ani-floor", lo=0.0, hi=100.0)
        if 0.0 < ani_floor <= 1.0:
            raise UserInputError(
                f"Invalid --ani-floor {ani_floor:g}: ANI is a percentage (0-100), not a "
                f"fraction. Did you mean --ani-floor {100 * ani_floor:g}?"
            )
        # A margin of zero makes a 0.01-point tie-break between two tips read as
        # different parents; a wide one makes every pair share a strain, so the verdict
        # is clonal whatever the data (0.5 gave reassortant where 40 gave clonal).
        _require_range(margin, "--margin", lo=0.0, hi=MAX_MARGIN, lo_open=True)
        overrides: dict[str, str] = {}
        for item in dataset or []:
            if "=" not in item:
                raise UserInputError(f"--dataset must be SEGMENT=path, got {item!r}")
            seg, path = item.split("=", 1)
            overrides[seg.strip()] = path.strip()
        if overrides:
            # An override is looked up by segment name; one that matches no record would
            # be ignored and that segment's dataset auto-detected instead.
            segments = [name for name, _seq in read_fasta(query) if name]
            unknown = sorted(set(overrides) - set(segments))
            if unknown:
                raise UserInputError(
                    f"--dataset names segment(s) not in the query: {', '.join(unknown)}. "
                    f"Segments in {query.name}: {', '.join(segments) or '(none)'}."
                )

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
        # The pairwise evidence behind the verdict, for every verdict: which pair
        # disagreed, and on which strains, is what a reader needs to check a call.
        ptsv = output / "segment_pairs.tsv"
        with open(ptsv, "w") as fo:
            fo.write("segment_a\tsegment_b\trelation\tnote\n")
            for pair in result.pairs:
                fo.write(f"{pair.segment_a}\t{pair.segment_b}\t{pair.relation}\t"
                         f"{' '.join(pair.note.split())}\n")

        logger.info("Reassortment verdict: %s", result.verdict.upper())
        logger.info("Parent constellation: %s", mosaic or "(no segments assigned)")
        if result.pair_notes:
            logger.info("Pairwise evidence: %s", "; ".join(result.pair_notes))
        logger.info("Wrote %s, %s and %s", tsv, ctsv, ptsv)

        if result.scans:
            stsv = output / "segment_scan.tsv"
            with open(stsv, "w") as fo:
                fo.write("segment\tintragenic_recombination\tn_regions\tnote\n")
                for sc in result.scans:
                    flag = "yes" if sc.recombinant else ("no" if sc.scanned else "n/a")
                    # A failure note quotes the aligner's own message, which can span
                    # lines and hold tabs; keep it to one cell.
                    note = " ".join(sc.note.split())
                    fo.write(f"{sc.segment}\t{flag}\t{sc.n_regions}\t{note}\n")
            rollup = " | ".join(f"{sc.segment}: {sc.note}" for sc in result.scans)
            logger.info("Intragenic scan: %s", rollup)
            logger.info("Wrote %s", stsv)
