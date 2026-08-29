"""On-disk cache key helpers, atomic installation, and completeness manifests."""

from __future__ import annotations

from pathlib import Path

import pytest

from tessera.core.cache import (
    CACHE_MANIFEST,
    atomic_cache_dir,
    blast_cache_path,
    cached_genomes,
    ncbi_virus_cache,
    nextclade_cache,
    read_blast_cache,
    read_cache_manifest,
    write_blast_cache,
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


# --- BLAST result cache ---------------------------------------------------

def test_blast_key_covers_everything_that_changes_the_answer(tmp_path: Path):
    """A cached hit list must never be served for a different question."""
    base = ("blastn", "nt", 5, "", "ACGTACGT")
    same = blast_cache_path(base, override=tmp_path)
    assert blast_cache_path(base, override=tmp_path) == same  # deterministic

    for changed in (
        ("blastp", "nt", 5, "", "ACGTACGT"),          # program
        ("blastn", "refseq", 5, "", "ACGTACGT"),      # database
        ("blastn", "nt", 10, "", "ACGTACGT"),         # hit limit
        ("blastn", "nt", 5, "Norovirus[Organism]", "ACGTACGT"),  # entrez filter
        ("blastn", "nt", 5, "", "TTTTTTTT"),          # the sequence itself
    ):
        assert blast_cache_path(changed, override=tmp_path) != same


def test_blast_cache_roundtrip(tmp_path: Path):
    path = blast_cache_path(("blastn", "nt", 5, "", "ACGT"), override=tmp_path)
    hits = [{"accession": "NC_1", "title": "a hit", "pct_identity": 98.0,
             "query_coverage": 95.0, "evalue": 1e-30}]
    write_blast_cache(path, hits)

    cached = read_blast_cache(path, max_age_days=14.0)
    assert cached is not None
    payload, age_days = cached
    assert payload == hits
    assert age_days < 1.0


def test_blast_cache_stores_an_empty_result(tmp_path: Path):
    # "this search found nothing" is a real answer that cost the same minutes.
    path = blast_cache_path(("blastn", "nt", 5, "", "ACGT"), override=tmp_path)
    write_blast_cache(path, [])
    cached = read_blast_cache(path, max_age_days=14.0)
    assert cached is not None and cached[0] == []


def test_blast_cache_expires(tmp_path: Path):
    """nt grows, so an indefinitely reused hit list would quietly stop reflecting it."""
    import json
    from datetime import UTC, datetime, timedelta

    path = blast_cache_path(("blastn", "nt", 5, "", "ACGT"), override=tmp_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    old = datetime.now(UTC) - timedelta(days=30)
    path.write_text(json.dumps({"fetched (UTC)": old.isoformat(), "hits": []}))

    assert read_blast_cache(path, max_age_days=14.0) is None  # too old
    assert read_blast_cache(path, max_age_days=60.0) is not None  # within a longer window


def test_blast_cache_zero_days_disables_reuse(tmp_path: Path):
    path = blast_cache_path(("blastn", "nt", 5, "", "ACGT"), override=tmp_path)
    write_blast_cache(path, [{"accession": "NC_1", "title": "x", "pct_identity": 1.0,
                              "query_coverage": 1.0, "evalue": 0.0}])
    assert read_blast_cache(path, max_age_days=0.0) is None


def test_blast_cache_survives_a_corrupt_entry(tmp_path: Path):
    path = blast_cache_path(("blastn", "nt", 5, "", "ACGT"), override=tmp_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    for junk in ("{ truncated", "[]", '{"hits": "not a list"}', '{"fetched (UTC)": "nope"}'):
        path.write_text(junk)
        assert read_blast_cache(path, max_age_days=14.0) is None  # falls back to searching


def test_blast_cache_write_failure_is_not_fatal(tmp_path: Path):
    # A cache that cannot be written must not fail the search it was serving.
    blocked = tmp_path / "blocked"
    blocked.write_text("i am a file, not a directory")
    write_blast_cache(blocked / "sub" / "x.json", [])  # must not raise
