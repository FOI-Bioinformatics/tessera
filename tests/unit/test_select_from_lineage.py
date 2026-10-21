"""_select_from should type the pool from its headers and forward lineage
selection into select_regional."""

from __future__ import annotations

import logging
from pathlib import Path

from tessera.discover import iterate as it
from tessera.discover.iterate import FillParams
from tessera.discover.pool import PoolSelection

_LOG = logging.getLogger("tessera.test")


def _write(dir_: Path, acc: str, note: str) -> Path:
    p = dir_ / f"{acc}.fasta"
    p.write_text(f">{acc} {note}\nACGT\n")
    return p


def test_select_from_passes_lineage_map(monkeypatch, tmp_path):
    pool = tmp_path / "pool"
    pool.mkdir()
    g1 = _write(pool, "ACC1", "HIV-1 subtype A1bc complete genome")
    g2 = _write(pool, "ACC2", "HIV-1 subtype B2ab complete genome")
    captured = {}

    def fake_select_regional(query, genomes, **kwargs):
        captured.update(kwargs)
        return PoolSelection(selected=list(genomes), table=[])

    monkeypatch.setattr("tessera.discover.pool.select_regional", fake_select_regional)

    params = FillParams(
        query=tmp_path / "q.fasta", collection=None, output=tmp_path / "out",
        taxon="Human immunodeficiency virus 1", keep_recombinant=False,
    )
    it._select_from(params, [g1, g2], _LOG)

    assert captured["keep_recombinant"] is False
    assert captured["lineage_of"]  # a non-empty map was built and forwarded
    assert set(captured["lineage_of"]) == {"ACC1", "ACC2"}



def _captured_lineage_map(monkeypatch, tmp_path, *, pool_headers: bool):
    pool = tmp_path / "pool"
    pool.mkdir()
    genomes = [_write(pool, "ACC1", "B"), _write(pool, "ACC2", "XBB"),
               _write(pool, "ACC3", "NA")]
    captured = {}

    def fake_select_regional(query, genomes, **kwargs):
        captured.update(kwargs)
        return PoolSelection(selected=list(genomes), table=[])

    monkeypatch.setattr("tessera.discover.pool.select_regional", fake_select_regional)
    params = FillParams(query=tmp_path / "q.fasta", collection=None, output=tmp_path / "out")
    it._select_from(params, genomes, _LOG, pool_headers=pool_headers)
    return captured["lineage_of"]


def test_select_from_reads_structured_pool_labels(monkeypatch, tmp_path):
    # Short clade labels as a Nextclade pool writes them. Mined as free text they are all
    # dropped, so lineage selection and the recombinant exclusion never see them.
    lineage_of = _captured_lineage_map(monkeypatch, tmp_path, pool_headers=True)
    assert lineage_of == {"ACC1": "B", "ACC2": "XBB"}  # "NA" is not a clade


def test_select_from_still_mines_a_user_pool(monkeypatch, tmp_path):
    # A --candidate-pool holds the user's own files; their headers are free text and a
    # trailing word must not be taken for a clade.
    assert _captured_lineage_map(monkeypatch, tmp_path, pool_headers=False) is None
