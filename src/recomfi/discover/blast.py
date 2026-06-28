"""Remote BLAST of an under-covered query subsequence against NCBI nt.

Uses Biopython's ``NCBIWWW.qblast`` (no local BLAST install, no extra dependency)
to find published sequences that match a stretch of the query the current
collection covers poorly -- candidate references to add. This contacts NCBI over
the network and is slow; it is only invoked by ``recomfi find-references``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

from ..core.errors import RecomFiError


class BlastError(RecomFiError):
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
) -> list[Hit]:
    """BLAST ``seq`` against NCBI ``database`` and return up to ``max_hits`` hits.

    ``entrez_query`` restricts the searched set with an Entrez expression (e.g.
    ``'Norovirus[Organism] NOT GII.P16-GII.1'``), used to exclude the query's own
    over-represented lineage so the divergent parents surface.

    Raises :class:`BlastError` on any network/NCBI/parse failure so the caller can
    skip a gap and carry on.
    """
    try:
        from Bio.Blast import NCBIWWW, NCBIXML
    except ImportError as exc:  # pragma: no cover - biopython is a hard dependency
        raise BlastError(f"Biopython BLAST is unavailable: {exc}") from exc

    if email:
        from Bio import Entrez

        Entrez.email = email

    logger.info("Submitting %d bp to NCBI %s/%s (this can take minutes)...",
                len(seq), program, database)
    qblast_kwargs = {"hitlist_size": max_hits}
    if entrez_query:
        qblast_kwargs["entrez_query"] = entrez_query
    try:
        handle = NCBIWWW.qblast(program, database, seq, **qblast_kwargs)
        record = NCBIXML.read(handle)
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
    return hits
