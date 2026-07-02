"""Assign each segment of a multi-FASTA query to its nearest Nextclade reference and
call reassortment by nearest-reference-strain consistency.

Reassortment is a per-segment phenomenon: for each segment, which parent did the query
get it from? A reassortant is a query whose segments trace to different parents. Each
segment is typed independently against its own Nextclade dataset, so this sidesteps the
single-backbone model that the intragenic recombination scan uses.
"""

from __future__ import annotations

import logging
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from ..core.cache import nextclade_cache
from ..core.errors import UserInputError
from ..core.io import read_fasta, strip_sequence_extension, write_fasta_record
from ..discover.nextclade import build_pool, resolve_dataset
from ..discover.panel import skani_available, skani_query_ani
from ..recomb.typing import first_header

DEFAULT_ANI_FLOOR = 80.0  # a segment below this ANI to every tip is left unassigned
TOP_K = 5                 # nearest strains kept per segment for the consistency test


@dataclass
class SegmentAssignment:
    segment: str
    dataset: str
    strain: str | None
    clade: str | None
    ani: float
    status: str  # "assigned" | "unassigned"


@dataclass
class ReassortmentResult:
    segments: list[SegmentAssignment] = field(default_factory=list)
    verdict: str = "undetermined"


def _call_reassortment(per_segment_topk: dict[str, set[str]], assigned: list[str]) -> str:
    """clonal iff one reference strain is in the top-k nearest of EVERY assigned segment
    (one parent explains all segments); reassortant iff no strain is shared across them;
    undetermined for fewer than two assigned segments."""
    considered = [s for s in assigned if s in per_segment_topk]
    if len(considered) < 2:
        return "undetermined"
    shared = set.intersection(*(per_segment_topk[s] for s in considered))
    return "clonal" if shared else "reassortant"


def _clade_of_tip(tip: Path) -> str | None:
    """The clade from a reconstructed tip's header ``>{strain} {clade}`` (2nd token)."""
    parts = first_header(tip).split(None, 1)
    clade = parts[1].strip() if len(parts) > 1 else ""
    return clade or None


def assign_segments(
    query: Path,
    *,
    dataset_overrides: dict[str, str] | None = None,
    ani_floor: float = DEFAULT_ANI_FLOOR,
    email: str | None = None,
    cache_dir: Path | None = None,
    logger: logging.Logger,
) -> ReassortmentResult:
    """Type each segment record of ``query`` against its Nextclade dataset and call
    reassortment. Each record's first-token name is the segment name (e.g. ``HA``)."""
    if not skani_available():
        raise UserInputError(
            "reassort needs skani. Install with: conda install -c bioconda skani skder"
        )
    overrides = dataset_overrides or {}
    records = read_fasta(query)
    if not records:
        raise UserInputError(f"Query FASTA {query} has no sequence.")

    result = ReassortmentResult()
    per_segment_topk: dict[str, set[str]] = {}
    with tempfile.TemporaryDirectory() as tmp:
        for name, seq in records:
            seg = name or f"segment_{len(result.segments) + 1}"
            seg_fasta = Path(tmp) / f"{strip_sequence_extension(seg)}.fasta"
            with open(seg_fasta, "w") as fo:
                write_fasta_record(fo, seg, seq)
            try:
                dataset = resolve_dataset(seg_fasta, overrides.get(seg), email=email,
                                          logger=logger)
            except Exception as exc:  # noqa: BLE001 - no dataset is a per-segment skip
                logger.info("[%s] no Nextclade dataset resolved (%s); unassigned.", seg, exc)
                result.segments.append(SegmentAssignment(seg, "?", None, None, 0.0, "unassigned"))
                continue
            tips = build_pool(
                dataset, cache_dir=nextclade_cache(dataset.path, dataset.tag, override=cache_dir),
                logger=logger,
            )
            ani = skani_query_ani(seg_fasta, tips, logger)
            ranked = sorted((t for t in tips if ani.get(t, (0.0, 0.0))[0] >= ani_floor),
                            key=lambda t: ani[t][0], reverse=True)
            if not ranked:
                result.segments.append(
                    SegmentAssignment(seg, dataset.path, None, None, 0.0, "unassigned"))
                continue
            per_segment_topk[seg] = {strip_sequence_extension(t.name) for t in ranked[:TOP_K]}
            best = ranked[0]
            result.segments.append(SegmentAssignment(
                seg, dataset.path, strip_sequence_extension(best.name),
                _clade_of_tip(best), ani[best][0], "assigned"))
            logger.info("[%s] %s -> %s (%s), ANI %.1f%%", seg, dataset.path,
                        strip_sequence_extension(best.name), _clade_of_tip(best), ani[best][0])

    assigned = [s.segment for s in result.segments if s.status == "assigned"]
    result.verdict = _call_reassortment(per_segment_topk, assigned)
    return result
