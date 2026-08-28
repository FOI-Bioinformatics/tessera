"""On-disk cache key helpers, atomic installation, and completeness manifests."""

from __future__ import annotations

from pathlib import Path

import pytest

from tessera.core.cache import (
    CACHE_MANIFEST,
    atomic_cache_dir,
    cached_genomes,
    ncbi_virus_cache,
    nextclade_cache,
    read_cache_manifest,
    write_cache_manifest,
)


def test_nextclade_cache_path_is_under_override_and_stable(tmp_path: Path):
    a = nextclade_cache("community/neherlab/hiv-1/hxb2", "2026-04-14--11-55-23Z", override=tmp_path)
    b = nextclade_cache("community/neherlab/hiv-1/hxb2", "2026-04-14--11-55-23Z", override=tmp_path)
    assert a == b  # deterministic
    assert tmp_path in a.parents
    assert a.parent.name == "nextclade"


def test_nextclade_cache_tag_changes_key(tmp_path: Path):
    old = nextclade_cache("nextstrain/sars-cov-2/XBB", "2024-01-01--00-00-00Z", override=tmp_path)
    new = nextclade_cache("nextstrain/sars-cov-2/XBB", "2026-01-01--00-00-00Z", override=tmp_path)
    assert old != new  # a newer tag is a distinct cache dir


# --- NCBI Virus: scope belongs in the key ---------------------------------

def test_ncbi_virus_scope_changes_key(tmp_path: Path):
    # A --source-refseq run and a default run ask for different genome sets; sharing
    # one slot would hand the second run the first one's panel.
    refseq = ncbi_virus_cache("Orthopoxvirus", scope={"source_refseq": True}, override=tmp_path)
    complete = ncbi_virus_cache(
        "Orthopoxvirus", scope={"source_refseq": False}, override=tmp_path
    )
    assert refseq != complete


def test_ncbi_virus_scope_key_is_order_insensitive_and_stable(tmp_path: Path):
    a = ncbi_virus_cache("Norovirus", scope={"source_refseq": True, "fetch_limit": 2000},
                         override=tmp_path)
    b = ncbi_virus_cache("Norovirus", scope={"fetch_limit": 2000, "source_refseq": True},
                         override=tmp_path)
    assert a == b


# --- completeness: a partial cache must not be reused ---------------------

def test_partial_cache_without_manifest_is_refused(tmp_path: Path):
    """The interrupted-fetch case: genomes on disk, but the fetch never finished.

    Reusing these silently narrows the donor panel, so detection quietly loses
    events it would otherwise find.
    """
    cache = tmp_path / "ncbi_virus" / "partial"
    cache.mkdir(parents=True)
    (cache / "NC_000001.fasta").write_text(">NC_000001\nACGT\n")
    (cache / "NC_000002.fasta").write_text(">NC_000002\nACGT\n")

    assert len(cached_genomes(cache)) == 2  # a plain listing still sees them
    assert cached_genomes(cache, manifest_required=True) == []  # but it is not trusted

    write_cache_manifest(cache, taxon="test", genomes=2)
    assert len(cached_genomes(cache, manifest_required=True)) == 2


def test_manifest_roundtrip_and_unreadable_manifest(tmp_path: Path):
    cache = tmp_path / "c"
    cache.mkdir()
    write_cache_manifest(cache, taxon="Norovirus", source="refseq", genomes=7)
    loaded = read_cache_manifest(cache)
    assert loaded["taxon"] == "Norovirus"
    assert loaded["genomes"] == 7
    assert "written (UTC)" in loaded

    (cache / CACHE_MANIFEST).write_text("{ truncated")
    assert read_cache_manifest(cache) is None
    assert cached_genomes(cache, manifest_required=True) == []


# --- atomic installation --------------------------------------------------

def test_atomic_cache_dir_installs_only_on_success(tmp_path: Path):
    cache = tmp_path / "ncbi_virus" / "panel"
    with atomic_cache_dir(cache) as staging:
        (staging / "a.fasta").write_text(">a\nACGT\n")
        assert not cache.exists()  # nothing visible to a reader mid-build
        write_cache_manifest(staging, genomes=1)
    assert len(cached_genomes(cache, manifest_required=True)) == 1


def test_atomic_cache_dir_leaves_nothing_behind_on_failure(tmp_path: Path):
    cache = tmp_path / "ncbi_virus" / "panel"
    cache.parent.mkdir(parents=True)
    with pytest.raises(RuntimeError):
        with atomic_cache_dir(cache) as staging:
            (staging / "a.fasta").write_text(">a\nACGT\n")
            raise RuntimeError("network dropped mid-fetch")
    assert not cache.exists()
    # and no staging tree is orphaned beside it
    assert list(cache.parent.iterdir()) == []


def test_atomic_cache_dir_keeps_a_concurrent_winner(tmp_path: Path):
    """Two runs fetching the same panel: the first to finish wins, the second yields.

    Both answer the same request, so failing the loser would be gratuitous.
    """
    cache = tmp_path / "ncbi_virus" / "panel"
    with atomic_cache_dir(cache) as staging:
        (staging / "mine.fasta").write_text(">mine\nACGT\n")
        write_cache_manifest(staging, genomes=1)
        # simulate the race: another process installs a complete cache first
        cache.mkdir(parents=True)
        (cache / "theirs.fasta").write_text(">theirs\nACGT\n")
        write_cache_manifest(cache, genomes=1)

    kept = [p.name for p in cached_genomes(cache, manifest_required=True)]
    assert kept == ["theirs.fasta"]
