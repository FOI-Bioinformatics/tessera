"""Remote BLAST of an under-covered query subsequence against NCBI nt.

Uses Biopython's ``NCBIWWW.qblast`` (no local BLAST install, no extra dependency)
to find published sequences that match a stretch of the query the current
collection covers poorly -- candidate references to add. This contacts NCBI over
the network and is slow; it is only invoked by ``tessera find-references``.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass
from pathlib import Path

from ..core.cache import blast_cache_path, read_blast_cache, write_blast_cache
from ..core.errors import TesseraError
from ..core.ncbi import blast_cache_days, blast_throttle, resolve_email, with_retries


class BlastError(TesseraError):
    """A remote BLAST request failed (network, NCBI, or parse error)."""


@dataclass
class Hit:
    accession: str
    title: str
    pct_identity: float
    query_coverage: float
    evalue: float


def _organism(title: str) -> str:
    """Best-effort organism from an NCBI defline (``... [Organism]`` if present)."""
    if title.endswith("]") and "[" in title:
        return title[title.rfind("[") + 1 : -1]
    return title.split(",")[0]


def blast_subsequence(
    seq: str,
    *,
    max_hits: int,
    logger: logging.Logger,
    email: str | None = None,
    program: str = "blastn",
    database: str = "nt",
    entrez_query: str | None = None,
    cache_dir: str | Path | None = None,
) -> list[Hit]:
    """BLAST ``seq`` against NCBI ``database`` and return up to ``max_hits`` hits.

    ``entrez_query`` restricts the searched set with an Entrez expression (e.g.
    ``'Norovirus[Organism] NOT GII.P16-GII.1'``), used to exclude the query's own
    over-represented lineage so the divergent parents surface.

    Results are cached on disk, keyed by the sequence and every parameter that changes
    what comes back. This is the slowest step in recruitment -- minutes per search, and
    paced at one submission per 10 seconds -- so a re-run, a resumed run, or a second
    round asking about the same window would otherwise pay for it again. Entries expire
    (see :data:`BLAST_CACHE_DAYS`) because ``nt`` grows.

    Raises :class:`BlastError` on any network/NCBI/parse failure so the caller can
    skip a gap and carry on.
    """
    cache_path = blast_cache_path(
        (program, database, max_hits, entrez_query or "", seq), override=cache_dir
    )
    cached = read_blast_cache(cache_path, max_age_days=blast_cache_days())
    if cached is not None:
        payload, age_days = cached
        logger.info(
            "Reusing a cached BLAST result (%d hit(s), %.1f day(s) old). Delete %s to "
            "force a fresh search.", len(payload), age_days, cache_path,
        )
        return [Hit(**row) for row in payload]

    try:
        from Bio.Blast import NCBIWWW, NCBIXML
    except ImportError as exc:  # pragma: no cover - biopython is a hard dependency
        raise BlastError(f"Biopython BLAST is unavailable: {exc}") from exc

    resolved_email = resolve_email(email)
    if resolved_email:
        from Bio import Entrez

        # Sets the contact address for Entrez calls. Note that Biopython's `qblast`
        # takes neither an email nor an api_key parameter, so the BLAST URL API below
        # is unauthenticated whatever is passed here; pacing is the only lever we have.
        Entrez.email = resolved_email

    logger.info("Submitting %d bp to NCBI %s/%s (this can take minutes)...",
                len(seq), program, database)
    qblast_kwargs = {"hitlist_size": max_hits}
    if entrez_query:
        qblast_kwargs["entrez_query"] = entrez_query

    def _submit():
        # Space submissions: this runs on every machine Tessera is installed on, and a
        # gap-by-subtile sweep would otherwise fire them back to back.
        blast_throttle.wait(logger)
        handle = NCBIWWW.qblast(program, database, seq, **qblast_kwargs)
        return NCBIXML.read(handle)

    try:
        record = with_retries(_submit, what="NCBI BLAST request", logger=logger)
    except Exception as exc:  # noqa: BLE001 - any failure becomes a clean skip
        raise BlastError(f"NCBI BLAST request failed: {exc}") from exc

    qlen = record.query_length or len(seq)
    hits: list[Hit] = []
    for alignment in record.alignments[:max_hits]:
        hsp = alignment.hsps[0]
        hits.append(Hit(
            accession=alignment.accession,
            title=alignment.hit_def,
            pct_identity=round(100.0 * hsp.identities / hsp.align_length, 1),
            query_coverage=round(100.0 * hsp.align_length / qlen, 1),
            evalue=hsp.expect,
        ))
    # Cached even when empty: "this search returned nothing" is a real answer and cost
    # the same minutes to obtain as any other.
    write_blast_cache(cache_path, [asdict(hit) for hit in hits])
    return hits
