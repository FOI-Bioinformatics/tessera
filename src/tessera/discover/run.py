"""`tessera find-references`: find and download references missing from a collection.

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
# A hit this close over near-full coverage is almost certainly the query's own
# record (the MSA labels the query by lineage, so BLAST can't recognise it).
SELF_HIT_IDENTITY = 99.5
SELF_HIT_COVERAGE = 95.0


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
    exclude: tuple[str, ...] = ()  # accessions to drop from candidates
    keep_self_hits: bool = False  # keep near-identical hits (the query's own record)
    curate: bool = False  # after download, drop the query's siblings and dereplicate
    sibling_margin: float = 3.0  # query-ANI must beat the backbone's by this many %
    af_min: float = 80.0  # ... over at least this % of the query (whole-genome match)
    derep_ani: float = 99.0  # skDER: collapse references >= this ANI to one representative
    subtile: int = 400  # search gaps longer than this in overlapping sub-intervals (0 = off)


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
    candidates = collect_candidates(
        targets, query_row, existing,
        max_hits=params.max_hits, email=params.email,
        exclude={_base_accession(e) for e in params.exclude},
        keep_self_hits=params.keep_self_hits, logger=logger, subtile=params.subtile,
    )

    _write_candidates(params.output, candidates, logger)
    _print_candidates(candidates, logger)

    if params.download is not None:
        downloaded = _download(candidates, params.download, logger)
        _write_downloaded(params.output, downloaded, logger)
        if params.curate and downloaded:
            _curate_download(params, query_label, query_row, logger)
    elif candidates:
        logger.info(
            "Re-run with --download <collection_dir> to add the new references, "
            "then rebuild the MSA with 'tessera msa'."
        )
    return candidates


def _tiles(length: int, tile: int) -> list[tuple[int, int]]:
    """Overlapping sub-intervals (~``tile`` long, half-overlap) covering ``length``.

    A long gap is searched in pieces so a short divergent tract is not diluted by
    high-matching flanks (which would otherwise dominate one whole-gap BLAST and hide
    the tract's donor). Overlap keeps a tract that straddles a boundary intact in at
    least one tile.
    """
    if tile <= 0 or length <= tile:
        return [(0, length)]
    step = max(MIN_SUBSEQ, tile // 2)
    starts = list(range(0, length - tile + 1, step))
    if starts[-1] != length - tile:
        starts.append(length - tile)
    return [(s, s + tile) for s in starts]


def collect_candidates(
    targets: list[CoverageGap],
    query_row: str,
    existing: set[str],
    *,
    max_hits: int,
    email: str | None,
    exclude: set[str],
    keep_self_hits: bool,
    logger: logging.Logger,
    subtile: int = 0,
) -> list[Candidate]:
    """BLAST each gap's query subsequence and return the kept candidate hits.

    When ``subtile`` > 0, a gap longer than that is searched in overlapping
    sub-intervals so a short divergent tract surfaces its own donor instead of being
    diluted by the gap's high-matching flanks. Hits in ``exclude`` (version-insensitive)
    or recognised as the query's own record are dropped and logged.
    """
    candidates: list[Candidate] = []
    skipped_self: list[str] = []
    skipped_excluded: list[str] = []
    for gap in targets:
        subseq = query_row[gap.msa_start : gap.msa_end].replace("-", "")
        if len(subseq) < MIN_SUBSEQ:
            logger.info(
                "Skipping gap query %d-%d: only %d bp after de-gapping.",
                gap.query_start, gap.query_end, len(subseq),
            )
            continue
        tiles = _tiles(len(subseq), subtile)
        logger.info(
            "Searching gap query %d-%d (%d bp, closest current ref %s ~%.2f)%s...",
            gap.query_start, gap.query_end, len(subseq), gap.best_label, gap.mean_best,
            f" in {len(tiles)} sub-intervals" if len(tiles) > 1 else "",
        )
        seen_here: set[str] = set()  # dedup an accession found in several overlapping tiles
        for tstart, tend in tiles:
            tile_seq = subseq[tstart:tend]
            if len(tile_seq) < MIN_SUBSEQ:
                continue
            try:
                hits = blast_subsequence(tile_seq, max_hits=max_hits, logger=logger, email=email)
            except BlastError as exc:
                logger.warning("  BLAST failed for a sub-interval, skipping: %s", exc)
                continue
            for hit in hits:
                if _base_accession(hit.accession) in exclude:
                    skipped_excluded.append(hit.accession)
                elif _is_self_hit(hit, keep_self_hits):
                    skipped_self.append(hit.accession)
                elif hit.accession not in seen_here:
                    seen_here.add(hit.accession)
                    candidates.append(Candidate(gap, hit, hit.accession in existing))

    if skipped_self:
        logger.info(
            "Skipped %d near-identical hit(s) as the query's own record (>= %.1f%% "
            "identity): %s. Use --keep-self-hits to keep them.",
            len(skipped_self), SELF_HIT_IDENTITY, ", ".join(skipped_self),
        )
    if skipped_excluded:
        logger.info("Excluded %d hit(s) by --exclude: %s",
                    len(skipped_excluded), ", ".join(skipped_excluded))
    return candidates


def _base_accession(accession: str) -> str:
    """Drop the version suffix so 'U54771' and 'U54771.1' compare equal."""
    return accession.split(".")[0].strip()


def _is_self_hit(hit: Hit, keep_self_hits: bool) -> bool:
    """A near-identical, near-full-length hit -- almost certainly the query itself."""
    return (
        not keep_self_hits
        and hit.pct_identity >= SELF_HIT_IDENTITY
        and hit.query_coverage >= SELF_HIT_COVERAGE
    )


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


def _download(
    candidates: list[Candidate], dest: Path, logger: logging.Logger
) -> list[Candidate]:
    """Fetch the best new hit per gap (dedup by accession). Returns what was saved."""
    chosen: dict[str, Candidate] = {}
    for c in candidates:
        if not c.in_collection and c.hit.accession not in chosen:
            chosen[c.hit.accession] = c
    if not chosen:
        logger.info("No new references to download (all hits already in the collection).")
        return []
    logger.info("Downloading %d new reference(s) to %s", len(chosen), dest)
    downloaded: list[Candidate] = []
    for accession, candidate in chosen.items():
        try:
            out = efetch_fasta(accession, dest, logger)
            logger.info("  + %s -> %s", accession, out.name)
            downloaded.append(candidate)
        except Exception as exc:  # noqa: BLE001 - report and continue the batch
            logger.warning("  ! failed to download %s: %s", accession, exc)
    return downloaded


def _curate_download(
    params: FindRefParams, query_label: str, query_row: str, logger: logging.Logger
) -> None:
    """Drop the query's siblings and dereplicate the download directory in place.

    The backbone (the query's whole-genome anchor) is chosen from the existing
    ``--collection``, so a freshly-downloaded sibling cannot be mistaken for it. The
    query is reconstructed from its (de-gapped) MSA row, as skani needs a FASTA.
    """
    from .panel import (
        curate_collection_dir,
        pick_backbone,
        skani_available,
        write_panel_tsv,
    )

    if not skani_available():
        raise UserInputError(
            "--curate needs skani. Install with: conda install -c bioconda skani skder"
        )
    if params.collection is None or not params.collection.is_dir():
        logger.warning("--curate needs --collection (the backbone source); skipping curation.")
        return
    qfasta = params.output / "query.degapped.fasta"
    qfasta.write_text(f">{query_label}\n{query_row.replace('-', '')}\n")
    backbone = pick_backbone(
        qfasta, [p for p in params.collection.iterdir() if p.is_file()],
        af_min=params.af_min, logger=logger,
    )
    if backbone is None:
        logger.warning("Could not determine a backbone from --collection; skipping curation.")
        return
    curation = curate_collection_dir(
        qfasta, params.download, backbone,
        ani_margin=params.sibling_margin, af_min=params.af_min,
        derep_ani=params.derep_ani, logger=logger,
    )
    write_panel_tsv(params.output / "panel_lineages.tsv", curation.table, logger)


def _write_downloaded(
    output: Path, downloaded: list[Candidate], logger: logging.Logger
) -> None:
    """Record exactly which references were added (a precise 'what was downloaded')."""
    path = output / "downloaded_references.tsv"
    logger.info("Writing downloaded references: %s", path)
    header = ["accession", "filename", "pct_identity", "gap_query_start", "gap_query_end", "title"]
    with open(path, "w") as fo:
        fo.write("\t".join(header) + "\n")
        for c in downloaded:
            fo.write("\t".join(map(str, [
                c.hit.accession, f"{c.hit.accession}.fasta", c.hit.pct_identity,
                c.gap.query_start, c.gap.query_end, c.hit.title,
            ])) + "\n")
