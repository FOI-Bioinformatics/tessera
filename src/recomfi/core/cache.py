"""A small on-disk cache for expensive, reusable recruitment artifacts.

Fresh-start recruitment fetches a taxon-scoped genome set from NCBI Virus every run
-- the dominant cost. Caching the fetched panel per taxon makes a second ``detect``
or ``fill-references`` run, and every iterative round, skip the network entirely.

The cache is a plain directory tree keyed by a sanitized taxon name plus a hash of
the cache key; there is no database and no expiry (genome sets are append-only at the
source, and a stale panel is harmless -- delete the directory to refresh).
"""

from __future__ import annotations

import hashlib
import os
import re
from pathlib import Path

_FASTA_SUFFIXES = {".fasta", ".fa", ".fna", ".gz"}


def cache_root(override: str | Path | None = None) -> Path:
    """The cache base directory: ``override``, ``$RECOMFI_CACHE``, or ``~/.cache/recomfi``."""
    if override:
        return Path(override)
    env = os.environ.get("RECOMFI_CACHE")
    return Path(env) if env else Path.home() / ".cache" / "recomfi"


def _slug(text: str) -> str:
    return re.sub(r"\W+", "_", text).strip("_")[:48] or "x"


def ncbi_virus_cache(taxon: str, *, override: str | Path | None = None) -> Path:
    """The cache directory for a taxon's recruited NCBI Virus panel (created on demand)."""
    key = hashlib.sha1(taxon.encode()).hexdigest()[:12]  # noqa: S324 - non-cryptographic
    path = cache_root(override) / "ncbi_virus" / f"{_slug(taxon)}_{key}"
    return path


def cached_genomes(directory: Path) -> list[Path]:
    """Genome FASTA files already present in a cache directory (empty if none/missing)."""
    if not directory.is_dir():
        return []
    return sorted(
        p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in _FASTA_SUFFIXES
    )
