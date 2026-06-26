"""`recomfi find-references`: find and download references missing from a collection.

Flags under-covered stretches of the query (where even the closest reference is a
poor match), BLASTs each against NCBI to find better-matching published
sequences, reports the candidates, and -- on request -- downloads them into a
collection directory ready to re-run with.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

from ..core.errors import UserInputError
from ..core.io import read_fasta, strip_sequence_extension
from ..recomb.coverage import CoverageGap, CoverageParams, call_coverage_gaps
from ..recomb.similarity import compute_similarity
from .blast import BlastError, Hit, blast_subsequence
from .fetch import efetch_available, efetch_fasta

MIN_SUBSEQ = 50  # too short to BLAST meaningfully


@dataclass
class FindRefParams:
    msa: Path
    query: str
    output: Path
    collection: Path | None = None
    window_size: int = 1000
    window_step: int = 100
    coverage_floor: float | None = None
    coverage_rel_drop: float = 0.05
    max_hits: int = 5
    top_gaps: int = 3
    email: str | None = None
    download: Path | None = None


@dataclass
class Candidate:
    gap: CoverageGap
    hit: Hit
    in_collection: bool


def find_references(params: FindRefParams, logger: logging.Logger) -> list[Candidate]:
    """Run the full find-references workflow; returns the candidate list."""
    if params.download is not None and not efetch_available():
        raise UserInputError(
            "Downloading needs Entrez Direct. Install with: "
            "conda install -c bioconda entrez-direct"
        )

    query_label = strip_sequence_extension(params.query)
    result = compute_similarity(
        str(params.msa), query_label,
        window_size=params.window_size, window_step=params.window_step,
    )
    cov = CoverageParams.with_defaults(
        params.window_size, floor=params.coverage_floor, rel_drop=params.coverage_rel_drop,
    )
    gaps, threshold = call_coverage_gaps(result, params.window_size, cov)
    logger.info("Coverage threshold %.3f; %d under-covered region(s).", threshold, len(gaps))

    params.output.mkdir(parents=True, exist_ok=True)
    if not gaps:
        logger.info(
            "No coverage gaps -- the collection already covers the query well. "
            "Nothing to search for."
        )
        _write_candidates(params.output, [], logger)
        return []

    records = dict(read_fasta(str(params.msa)))
    query_row = records.get(query_label)
    if query_row is None:
        raise UserInputError(
            f"Query '{query_label}' not found in the MSA. Available: "
            f"{', '.join(sorted(records))}"
        )
    existing = _existing_labels(result, params.collection, query_label)

    targets = sorted(gaps, key=lambda g: g.length_bp, reverse=True)[: params.top_gaps]
    candidates: list[Candidate] = []
    for gap in targets:
        subseq = query_row[gap.msa_start : gap.msa_end].replace("-", "")
        if len(subseq) < MIN_SUBSEQ:
            logger.info(
                "Skipping gap query %d-%d: only %d bp after de-gapping.",
                gap.query_start, gap.query_end, len(subseq),
            )
            continue
        logger.info(
            "Searching gap query %d-%d (%d bp, closest current ref %s ~%.2f)...",
            gap.query_start, gap.query_end, len(subseq), gap.best_label, gap.mean_best,
        )
        try:
            hits = blast_subsequence(
                subseq, max_hits=params.max_hits, logger=logger, email=params.email,
            )
        except BlastError as exc:
            logger.warning("  BLAST failed for this gap, skipping: %s", exc)
            continue
        for hit in hits:
            candidates.append(Candidate(gap, hit, hit.accession in existing))

    _write_candidates(params.output, candidates, logger)
    _print_candidates(candidates, logger)

    if params.download is not None:
        _download(candidates, params.download, logger)
    elif candidates:
        logger.info(
            "Re-run with --download <collection_dir> to add the new references, "
            "then rebuild the MSA with 'recomfi msa'."
        )
    return candidates


def _existing_labels(result, collection: Path | None, query_label: str) -> set[str]:
    labels = set(result.similarities) | {query_label}
    if collection and collection.is_dir():
        labels |= {strip_sequence_extension(p.name) for p in collection.iterdir() if p.is_file()}
    return labels


def _write_candidates(output: Path, candidates: list[Candidate], logger: logging.Logger) -> None:
    path = output / "candidate_references.tsv"
    logger.info("Writing candidate references: %s", path)
    header = [
        "gap_query_start", "gap_query_end", "gap_kind", "accession",
        "pct_identity", "query_coverage", "evalue", "in_collection", "title",
    ]
    with open(path, "w") as fo:
        fo.write("\t".join(header) + "\n")
        for c in candidates:
            fo.write("\t".join(map(str, [
                c.gap.query_start, c.gap.query_end, c.gap.kind, c.hit.accession,
                c.hit.pct_identity, c.hit.query_coverage, c.hit.evalue,
                "yes" if c.in_collection else "no", c.hit.title,
            ])) + "\n")


def _print_candidates(candidates: list[Candidate], logger: logging.Logger) -> None:
    if not candidates:
        logger.info("No candidate references found.")
        return
    logger.info("Candidate references (new = not already in the collection):")
    logger.info("  %-12s %-6s %-5s %-4s %s", "Accession", "Ident", "Cov", "New", "Title")
    logger.info("  %s", "-" * 76)
    for c in candidates:
        logger.info(
            "  %-12s %5.1f%% %4.0f%% %-4s %s",
            c.hit.accession, c.hit.pct_identity, c.hit.query_coverage,
            "no" if c.in_collection else "yes", c.hit.title[:48],
        )


def _download(candidates: list[Candidate], dest: Path, logger: logging.Logger) -> None:
    # The best new hit per gap (dedup by accession), so the collection grows by the
    # most relevant sequences rather than every hit.
    chosen: list[str] = []
    for c in candidates:
        if not c.in_collection and c.hit.accession not in chosen:
            chosen.append(c.hit.accession)
    if not chosen:
        logger.info("No new references to download (all hits already in the collection).")
        return
    logger.info("Downloading %d new reference(s) to %s", len(chosen), dest)
    for accession in chosen:
        try:
            out = efetch_fasta(accession, dest, logger)
            logger.info("  + %s -> %s", accession, out.name)
        except Exception as exc:  # noqa: BLE001 - report and continue the batch
            logger.warning("  ! failed to download %s: %s", accession, exc)
