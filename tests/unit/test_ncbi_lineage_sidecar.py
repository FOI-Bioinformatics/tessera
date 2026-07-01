"""The NCBI-datasets lineage sidecar written during fetch."""

from __future__ import annotations

from pathlib import Path

from tessera.discover.pool import NCBI_LINEAGES_TSV, write_ncbi_lineages


def test_write_ncbi_lineages_round_trip(tmp_path: Path):
    dest = tmp_path / "cache"
    dest.mkdir()
    out = write_ncbi_lineages(dest, {"ACC1": "DENV1", "ACC2": "DENV4"})
    assert out == dest / NCBI_LINEAGES_TSV
    lines = sorted(out.read_text().splitlines())
    assert lines == ["ACC1\tDENV1", "ACC2\tDENV4"]


def test_write_ncbi_lineages_empty_writes_nothing(tmp_path: Path):
    dest = tmp_path / "cache"
    dest.mkdir()
    out = write_ncbi_lineages(dest, {})
    assert out is None
    assert not (dest / NCBI_LINEAGES_TSV).exists()
