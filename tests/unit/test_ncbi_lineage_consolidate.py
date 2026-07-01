"""Consolidating and reading the NCBI-datasets lineage sidecar in recruitment."""

from __future__ import annotations

from pathlib import Path

from tessera.discover.iterate import _read_ncbi_lineages, _record_ncbi_lineages
from tessera.discover.pool import NCBI_LINEAGES_TSV


def _sidecar(d: Path, rows: dict[str, str]) -> None:
    d.mkdir(parents=True, exist_ok=True)
    (d / NCBI_LINEAGES_TSV).write_text("".join(f"{a}\t{lin}\n" for a, lin in rows.items()))


def test_read_absent_returns_none(tmp_path: Path):
    assert _read_ncbi_lineages(tmp_path) is None


def test_record_then_read_merges_and_dedups(tmp_path: Path):
    out = tmp_path / "out"
    out.mkdir()
    src1 = tmp_path / "cache1"
    src2 = tmp_path / "cache2"
    _sidecar(src1, {"ACC1": "DENV1", "ACC2": "DENV4"})
    _sidecar(src2, {"ACC2": "DENV4-updated", "ACC3": "DENV2"})  # ACC2 re-recorded
    _record_ncbi_lineages(out, src1)
    _record_ncbi_lineages(out, src2)
    rows = dict(_read_ncbi_lineages(out))
    # The later record wins for the duplicate accession.
    assert rows == {"ACC1": "DENV1", "ACC2": "DENV4-updated", "ACC3": "DENV2"}


def test_record_missing_source_is_noop(tmp_path: Path):
    out = tmp_path / "out"
    out.mkdir()
    _record_ncbi_lineages(out, tmp_path / "does_not_exist")
    assert _read_ncbi_lineages(out) is None
