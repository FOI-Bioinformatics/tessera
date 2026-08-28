"""A small on-disk cache for expensive, reusable recruitment artifacts.

Fresh-start recruitment fetches a taxon-scoped genome set from NCBI Virus every run
-- the dominant cost. Caching the fetched panel per taxon makes a second ``detect``
or ``fill-references`` run, and every iterative round, skip the network entirely.

The cache is a plain directory tree keyed by a sanitized taxon name plus a hash of
the cache key; there is no database and no expiry (genome sets are append-only at the
source, and a stale panel is harmless -- delete the directory to refresh).

A *partial* panel is not harmless, though: panel composition decides which donors can
be found at all, so a fetch interrupted half way through would otherwise leave a
truncated genome set that the next run reuses in silence. Two mechanisms prevent that.
:func:`atomic_cache_dir` stages a build in a sibling directory and moves it into place
only once it is complete, and :func:`write_cache_manifest` records what was fetched;
``cached_genomes(..., manifest_required=True)`` then ignores any directory without a
manifest, which also rejects partial caches left by earlier versions.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import tempfile
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path

_FASTA_SUFFIXES = {".fasta", ".fa", ".fna", ".gz"}

# Marks a cache directory as completely written. Its absence means "do not trust".
CACHE_MANIFEST = "tessera_cache.json"


def cache_root(override: str | Path | None = None) -> Path:
    """The cache base directory: ``override``, ``$TESSERA_CACHE``, or ``~/.cache/tessera``."""
    if override:
        return Path(override)
    env = os.environ.get("TESSERA_CACHE")
    return Path(env) if env else Path.home() / ".cache" / "tessera"


def _slug(text: str) -> str:
    return re.sub(r"\W+", "_", text).strip("_")[:48] or "x"


def _scope_key(scope: Mapping[str, object] | None) -> str:
    """Canonical string for the parameters that decide *what* a fetch returns."""
    if not scope:
        return ""
    return ";".join(f"{k}={scope[k]!r}" for k in sorted(scope))


def ncbi_virus_cache(
    taxon: str,
    *,
    scope: Mapping[str, object] | None = None,
    override: str | Path | None = None,
) -> Path:
    """The cache directory for a taxon's recruited NCBI Virus panel.

    ``scope`` carries the request parameters that change which genomes come back
    (RefSeq versus all complete genomes, any fetch cap). They belong in the key:
    keyed on the taxon alone, a ``--source-refseq`` run and a default run share one
    slot, and whichever ran second would silently inherit the other's panel.
    """
    payload = f"{taxon}|{_scope_key(scope)}" if scope else taxon
    key = hashlib.sha1(payload.encode()).hexdigest()[:12]  # noqa: S324 - non-cryptographic
    path = cache_root(override) / "ncbi_virus" / f"{_slug(taxon)}_{key}"
    return path


def nextclade_cache(path: str, tag: str, *, override: str | Path | None = None) -> Path:
    """Cache directory for a reconstructed Nextclade pool, keyed by ``path@tag``."""
    key = hashlib.sha1(f"{path}@{tag}".encode()).hexdigest()[:12]  # noqa: S324 - non-cryptographic
    return cache_root(override) / "nextclade" / f"{_slug(path)}_{key}"


def pango_alias_path(override: str | Path | None = None) -> Path:
    """The cached Pango ``alias_key.json`` path (fetched once, reused across runs)."""
    return cache_root(override) / "pango" / "alias_key.json"


@contextmanager
def atomic_cache_dir(cache_dir: Path) -> Iterator[Path]:
    """Stage a cache build next to its destination; install it only once complete.

    Yields a private staging directory. On a clean exit it is moved onto
    ``cache_dir`` with :func:`os.replace`, which is atomic on the same filesystem,
    so a reader either sees no cache or sees the finished one -- never a fetch
    caught half way. Any exception (``KeyboardInterrupt`` included, hence
    ``BaseException``) removes the staging tree instead.

    Paths handed out during the build point into the staging directory and stop
    existing once it is installed; re-list from ``cache_dir`` afterwards.
    """
    cache_dir.parent.mkdir(parents=True, exist_ok=True)
    build_dir = Path(tempfile.mkdtemp(dir=cache_dir.parent))
    try:
        yield build_dir
    except BaseException:
        shutil.rmtree(build_dir, ignore_errors=True)
        raise
    try:
        os.replace(build_dir, cache_dir)
    except OSError:
        # A concurrent run finished the same fetch first. Its copy answers the same
        # request as ours, so keep it rather than fail; only re-raise if the
        # destination is not in fact a populated cache.
        if not cache_dir.is_dir():
            shutil.rmtree(build_dir, ignore_errors=True)
            raise
        shutil.rmtree(build_dir, ignore_errors=True)


def write_cache_manifest(directory: Path, **fields: object) -> Path:
    """Mark ``directory`` as a complete cache entry and record how it was built."""
    manifest = {
        "written (UTC)": datetime.now(UTC).strftime("%Y-%m-%d %H:%M:%S"),
        **fields,
    }
    path = directory / CACHE_MANIFEST
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return path


def read_cache_manifest(directory: Path) -> dict | None:
    """The cache manifest for ``directory``, or ``None`` if absent or unreadable."""
    path = directory / CACHE_MANIFEST
    try:
        loaded = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    return loaded if isinstance(loaded, dict) else None


def cached_genomes(directory: Path, *, manifest_required: bool = False) -> list[Path]:
    """Genome FASTA files already present in a cache directory (empty if none/missing).

    With ``manifest_required``, a directory carrying no manifest is treated as
    absent rather than as a usable panel -- the interrupted-fetch case, and any
    cache written before manifests existed.
    """
    if not directory.is_dir():
        return []
    if manifest_required and read_cache_manifest(directory) is None:
        return []
    return sorted(
        p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in _FASTA_SUFFIXES
    )
