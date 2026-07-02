"""Assign each segment of a multi-FASTA query to its nearest Nextclade references and
call reassortment with the coverage-aware constellation model.

Reassortment is a per-segment phenomenon: for each segment, which parent did the query get
it from? A reassortant is a query whose segments trace to different parents. Each segment is
typed independently against its own Nextclade dataset, so this sidesteps the single-backbone
model that the intragenic recombination scan uses.
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
from .constellation import DEFAULT_MARGIN, ParentGroup, call_constellation

DEFAULT_ANI_FLOOR = 80.0  # a segment below this ANI to every tip is left unassigned
MIN_AF = 0.5              # a tip aligning over less than this fraction of the segment is ignored
TOP_K = 25                # internal cap on candidate strains kept per segment


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
    groups: list[ParentGroup] = field(default_factory=list)
    pair_notes: list[str] = field(default_factory=list)


def _clade_of_tip(tip: Path) -> str | None:
    """The clade from a reconstructed tip's header ``>{strain} {clade}`` (2nd token)."""
    parts = first_header(tip).split(None, 1)
    clade = parts[1].strip() if len(parts) > 1 else ""
    return clade or None


def _type_segment(seg, seq, overrides, ani_floor, email, cache_dir, tmp, logger):
    """Type one segment. Returns ``(SegmentAssignment, candidates, universe)`` where
    ``candidates`` is ``[(strain, ani)]`` best-first (empty if unassigned) and ``universe``
    is every strain in the dataset (empty if unassigned). Never raises for this segment."""
    seg_fasta = Path(tmp) / f"{strip_sequence_extension(seg)}.fasta"
    with open(seg_fasta, "w") as fo:
        write_fasta_record(fo, seg, seq)
    try:
        dataset = resolve_dataset(seg_fasta, overrides.get(seg), email=email, logger=logger)
        tips = build_pool(
            dataset, cache_dir=nextclade_cache(dataset.path, dataset.tag, override=cache_dir),
            logger=logger,
        )
        ani = skani_query_ani(seg_fasta, tips, logger)
    except Exception as exc:  # noqa: BLE001 - a per-segment failure is non-fatal
        logger.info("[%s] could not type segment (%s); unassigned.", seg, exc)
        return SegmentAssignment(seg, "?", None, None, 0.0, "unassigned"), [], set()

    eligible = [t for t in tips
                if ani.get(t, (0.0, 0.0))[0] >= ani_floor and ani.get(t, (0.0, 0.0))[1] >= MIN_AF]
    ranked = sorted(eligible, key=lambda t: ani[t][0], reverse=True)
    if not ranked:
        return SegmentAssignment(seg, dataset.path, None, None, 0.0, "unassigned"), [], set()

    candidates = [(strip_sequence_extension(t.name), ani[t][0]) for t in ranked[:TOP_K]]
    universe = {strip_sequence_extension(t.name) for t in tips}
    best = ranked[0]
    logger.info("[%s] %s -> %s (%s), ANI %.1f%%", seg, dataset.path,
                strip_sequence_extension(best.name), _clade_of_tip(best), ani[best][0])
    assignment = SegmentAssignment(seg, dataset.path, strip_sequence_extension(best.name),
                                   _clade_of_tip(best), ani[best][0], "assigned")
    return assignment, candidates, universe


def assign_segments(
    query: Path,
    *,
    dataset_overrides: dict[str, str] | None = None,
    ani_floor: float = DEFAULT_ANI_FLOOR,
    margin: float = DEFAULT_MARGIN,
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
    candidates: dict[str, list[tuple[str, float]]] = {}
    universes: dict[str, set[str]] = {}
    with tempfile.TemporaryDirectory() as tmp:
        for name, seq in records:
            seg = name or f"segment_{len(result.segments) + 1}"
            assignment, cand, universe = _type_segment(
                seg, seq, overrides, ani_floor, email, cache_dir, tmp, logger)
            result.segments.append(assignment)
            if assignment.status == "assigned":
                candidates[seg] = cand
                universes[seg] = universe

    call = call_constellation(candidates, universes, margin)
    result.verdict = call.verdict
    result.groups = call.groups
    result.pair_notes = call.pair_notes
    return result
