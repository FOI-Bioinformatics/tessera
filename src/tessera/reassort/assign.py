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
from ..core.errors import ToolExecutionError, UserInputError
from ..core.io import read_fasta, strip_sequence_extension, write_fasta_record
from ..discover.nextclade import NON_CLADE_MARKERS, build_pool, resolve_dataset
from ..discover.panel import skani_available, skani_query_ani
from ..recomb.typing import first_header
from .constellation import DEFAULT_MARGIN, ParentGroup, call_constellation
from .scan import SegmentScan, require_aligner, scan_segment

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
    scans: list[SegmentScan] = field(default_factory=list)


def cap_candidates(
    ranked: list[tuple[str, float]], *, margin: float, top_k: int = TOP_K
) -> list[tuple[str, float]]:
    """Bound the candidate list without ever cutting inside the near-best ANI window.

    ``top_k`` exists only to keep the list from growing without limit. The constellation
    call, though, decides concordance from the strains within ``margin`` ANI of the best
    (``constellation._near_best``), so a rank cap applied first can drop a genuinely
    shared parent before the margin is consulted -- and two segments of one clonal strain
    then look like they have different parents. Densely sampled Nextclade trees routinely
    put more than ``top_k`` strains inside a half-point window, so the cap is widened to
    cover the window whenever that happens.
    """
    if not ranked:
        return []
    floor = ranked[0][1] - margin
    within = sum(1 for _strain, ani in ranked if ani >= floor)
    return ranked[: max(top_k, within)]


def _clade_of_tip(tip: Path) -> str | None:
    """The clade from a reconstructed tip's header ``>{strain} {clade}`` (2nd token), or None
    when the tip is not clade-typed (an example genome or an untyped tip)."""
    parts = first_header(tip).split(None, 1)
    clade = parts[1].strip() if len(parts) > 1 else ""
    return None if clade in NON_CLADE_MARKERS else clade


def _type_segment(seg, seq, overrides, ani_floor, margin, email, cache_dir, tmp, logger):
    """Type one segment. Returns ``(SegmentAssignment, candidates, universe, dataset)`` where
    ``candidates`` is ``[(strain, ani)]`` best-first (empty if unassigned), ``universe`` is every
    strain in the dataset (empty if unassigned), and ``dataset`` is the resolved Nextclade dataset
    (``None`` if no dataset maps). Returns an unassigned row for genuine per-segment skips (no
    dataset, or a skani rejection of a short/odd segment); a transient failure (a download/network
    error, or any unexpected error) propagates so it surfaces rather than reading as unassigned."""
    seg_fasta = Path(tmp) / f"{strip_sequence_extension(seg)}.fasta"
    with open(seg_fasta, "w") as fo:
        write_fasta_record(fo, seg, seq)
    try:
        dataset = resolve_dataset(seg_fasta, overrides.get(seg), email=email, logger=logger)
    except UserInputError as exc:  # no Nextclade dataset maps to this segment (a genuine skip)
        logger.info("[%s] no Nextclade dataset resolved (%s); unassigned.", seg, exc)
        return SegmentAssignment(seg, "?", None, None, 0.0, "unassigned"), [], set(), None
    try:
        tips = build_pool(
            dataset, cache_dir=nextclade_cache(dataset.path, dataset.tag, override=cache_dir),
            logger=logger,
        )
        ani = skani_query_ani(seg_fasta, tips, logger)
    except ToolExecutionError as exc:  # skani rejects a very short / odd segment (a genuine skip)
        logger.info("[%s] could not type against %s (%s); unassigned.", seg, dataset.path, exc)
        skipped = SegmentAssignment(seg, dataset.path, None, None, 0.0, "unassigned")
        return skipped, [], set(), dataset
    # A download / network failure (raised as UserInputError) or any unexpected error is NOT
    # swallowed here: it propagates so a transient failure surfaces loudly instead of being
    # reported as a biological "unassigned".

    eligible = [t for t in tips
                if ani.get(t, (0.0, 0.0))[0] >= ani_floor and ani.get(t, (0.0, 0.0))[1] >= MIN_AF]
    ranked = sorted(eligible, key=lambda t: ani[t][0], reverse=True)
    if not ranked:
        unassigned = SegmentAssignment(seg, dataset.path, None, None, 0.0, "unassigned")
        return unassigned, [], set(), dataset

    scored = [(strip_sequence_extension(t.name), ani[t][0]) for t in ranked]
    candidates = cap_candidates(scored, margin=margin)
    universe = {strip_sequence_extension(t.name) for t in tips}
    best = ranked[0]
    logger.info("[%s] %s -> %s (%s), ANI %.1f%%", seg, dataset.path,
                strip_sequence_extension(best.name), _clade_of_tip(best), ani[best][0])
    assignment = SegmentAssignment(seg, dataset.path, strip_sequence_extension(best.name),
                                   _clade_of_tip(best), ani[best][0], "assigned")
    return assignment, candidates, universe, dataset


def assign_segments(
    query: Path,
    *,
    dataset_overrides: dict[str, str] | None = None,
    ani_floor: float = DEFAULT_ANI_FLOOR,
    margin: float = DEFAULT_MARGIN,
    email: str | None = None,
    cache_dir: Path | None = None,
    output: Path | None = None,
    scan_segments: bool = False,
    aligner: str = "mafft",
    logger: logging.Logger,
) -> ReassortmentResult:
    """Type each segment record of ``query`` against its Nextclade dataset and call
    reassortment. Each record's first-token name is the segment name (e.g. ``HA``). With
    ``scan_segments`` each assigned segment is additionally scanned for intragenic
    recombination (needs an ``aligner`` and an ``output`` directory)."""
    if not skani_available():
        raise UserInputError(
            "reassort needs skani. Install with: conda install -c bioconda skani skder"
        )
    if scan_segments:
        if output is None:
            raise UserInputError("--scan-segments needs an output directory (-o/--output).")
        require_aligner(aligner)
    overrides = dataset_overrides or {}
    records = read_fasta(query)
    if not records:
        raise UserInputError(f"Query FASTA {query} has no sequence.")
    names = [name for name, _ in records if name]
    dupes = sorted({n for n in names if names.count(n) > 1})
    if dupes:
        raise UserInputError(
            f"Query has duplicate segment names ({', '.join(dupes)}); each record must be a "
            f"distinct segment with a unique name."
        )

    result = ReassortmentResult()
    candidates: dict[str, list[tuple[str, float]]] = {}
    universes: dict[str, set[str]] = {}
    to_scan: dict[str, tuple[str, object]] = {}  # segment -> (seq, dataset) for assigned segments
    with tempfile.TemporaryDirectory() as tmp:
        for name, seq in records:
            seg = name or f"segment_{len(result.segments) + 1}"
            assignment, cand, universe, dataset = _type_segment(
                seg, seq, overrides, ani_floor, margin, email, cache_dir, tmp, logger)
            result.segments.append(assignment)
            if assignment.status == "assigned":
                candidates[seg] = cand
                universes[seg] = universe
                to_scan[seg] = (seq, dataset)

    call = call_constellation(candidates, universes, margin)
    result.verdict = call.verdict
    result.groups = call.groups
    result.pair_notes = call.pair_notes

    if scan_segments:
        for s in result.segments:
            if s.status == "assigned":
                seq, dataset = to_scan[s.segment]
                assert output is not None  # guarded above when scan_segments is set
                result.scans.append(scan_segment(
                    s.segment, seq, dataset, output,
                    aligner=aligner, cache_dir=cache_dir, logger=logger))
            else:
                result.scans.append(SegmentScan(s.segment, False, False, 0, "unassigned"))
    return result
