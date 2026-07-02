"""Intragenic per-segment recombination scan (the ``reassort --scan-segments`` flag).

For each assigned segment, build a per-clade-consensus panel from its Nextclade dataset,
align the segment query to it, and run the ordinary single-backbone recombination scan
(:func:`tessera.recomb.run.run_recomb`) inside that one segment. This is orthogonal to the
whole-segment reassortment call: reassortment asks which parent each segment came from; this
asks whether a single segment is itself a within-segment mosaic of two lineages.
"""

from __future__ import annotations

import logging
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

from ..core.cache import nextclade_cache
from ..core.errors import UserInputError
from ..core.io import strip_sequence_extension, write_fasta_record
from ..discover.nextclade import build_pool
from ..msa.build import MsaParams, build_msa
from ..recomb.regions import DEFAULT_METHODS
from ..recomb.run import RecombParams, run_recomb
from ..recomb.typing import (
    LINEAGES_TSV,
    first_header,
    lineage_map_from_rows,
    write_lineage_map,
)


@dataclass
class SegmentScan:
    segment: str
    scanned: bool          # False if it could not be scanned (single-clade panel or error)
    recombinant: bool      # True if run_recomb called at least one present region
    n_regions: int
    note: str


def window_params(genome_len: int) -> tuple[int, int, int]:
    """Adaptive (recomb window, step, selection window) scaled to a short segment length,
    so a ~1-2.5 kb flu segment is not given a window wider than its alignment."""
    window = max(120, min(500, genome_len // 12))
    step = max(20, window // 10)
    select_window = max(window, min(1500, genome_len // 4))
    return window, step, select_window


def require_aligner(aligner: str) -> None:
    """Fail up front if ``aligner`` is unknown or its binary is not on PATH."""
    from ..aligners.base import registry
    names = set(registry.names())
    if aligner not in names:
        raise UserInputError(
            f"Unknown aligner '{aligner}'. Available: {', '.join(sorted(names))}."
        )
    registry.create(aligner).preflight()  # raises MissingBinaryError if the binary is absent


def _clade_of_header(path: Path) -> str:
    """The clade token from a consensus genome's ``>{label} {clade}`` header."""
    parts = first_header(path).split(None, 1)
    return parts[1].strip() if len(parts) > 1 and parts[1].strip() else "?"


def _summarize_regions(path: Path) -> tuple[int, bool]:
    """Count present (non ``donor_absent``) regions in a recomb TSV -> ``(n, n > 0)``."""
    if not path.exists():
        return 0, False
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        return 0, False
    header = lines[0].split("\t")
    absent_idx = header.index("donor_absent") if "donor_absent" in header else None
    n = 0
    for ln in lines[1:]:
        fields = ln.split("\t")
        if absent_idx is not None and absent_idx < len(fields) and fields[absent_idx] == "yes":
            continue
        n += 1
    return n, n > 0


def scan_segment(
    segment: str, seq: str, dataset, out_dir: Path, *,
    aligner: str, cache_dir: Path | None, logger: logging.Logger,
) -> SegmentScan:
    """Scan one segment for intragenic recombination. Never raises: a failure is recorded as
    ``scanned=False`` so the caller can continue with the other segments."""
    seg_dir = out_dir / re.sub(r"[^\w.-]+", "_", segment)
    try:
        panel = build_pool(
            dataset,
            cache_dir=nextclade_cache(dataset.path, dataset.tag, override=cache_dir),
            logger=logger, per_clade_consensus=True,
        )
        if len(panel) < 2:
            return SegmentScan(segment, False, False, 0, "single-clade panel")

        seg_dir.mkdir(parents=True, exist_ok=True)
        collection = seg_dir / "collection"
        if collection.exists():
            shutil.rmtree(collection)
        collection.mkdir(parents=True)
        rows = []
        for p in panel:
            shutil.copy(p, collection / p.name)
            rows.append((strip_sequence_extension(p.name), _clade_of_header(p), "consensus"))
        query = seg_dir / "query.fasta"
        with open(query, "w") as fo:
            write_fasta_record(fo, segment, seq)
        rows.append((segment, "query", "query"))
        write_lineage_map(seg_dir / LINEAGES_TSV, rows)
        lineage_map = lineage_map_from_rows(rows)

        window, step, _sel = window_params(len(seq))
        msa = seg_dir / "panel.msa.fasta"
        build_msa(
            MsaParams(query=query, collection=collection, output=msa, aligner=aligner), logger)
        run_recomb(RecombParams(msa=msa, output=seg_dir, query=segment,
                                window_size=window, window_step=step, organism=segment,
                                methods=DEFAULT_METHODS, lineage_map=lineage_map), logger)
    except Exception as exc:  # noqa: BLE001 - a per-segment scan failure is non-fatal
        logger.info("[%s] intragenic scan failed (%s); not scanned.", segment, exc)
        return SegmentScan(segment, False, False, 0, f"scan failed: {exc}")

    n_regions, recombinant = _summarize_regions(seg_dir / "recombination_regions.tsv")
    return SegmentScan(segment, True, recombinant, n_regions,
                       f"{n_regions} region(s)" if recombinant else "none")
