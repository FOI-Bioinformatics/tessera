"""On-disk cache key helpers."""

from __future__ import annotations

from pathlib import Path

from tessera.core.cache import nextclade_cache


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
