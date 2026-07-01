"""_type_panel routing: tier-1 by default, the full ladder under --deep-typing."""

from __future__ import annotations

import logging
from pathlib import Path

from tessera.discover import iterate as it
from tessera.discover.iterate import FillParams, _type_panel
from tessera.recomb.typing import LINEAGES_TSV

_LOG = logging.getLogger("tessera.test")


def _panel(tmp_path: Path) -> tuple[FillParams, Path]:
    out = tmp_path / "out"
    coll = out / "collection"
    coll.mkdir(parents=True)
    (coll / "ACC1.fasta").write_text(">ACC1 genotype 1\nACGT\n")
    q = tmp_path / "q.fasta"
    q.write_text(">QRY genotype 2\nACGT\n")
    return FillParams(query=q, collection=None, output=out, taxon=None), out


def test_default_uses_tier1_not_ladder(tmp_path, monkeypatch):
    params, out = _panel(tmp_path)
    called = {"ladder": False}
    monkeypatch.setattr(it, "assign_lineages",
                        lambda *a, **k: called.__setitem__("ladder", True) or [])
    _type_panel(params, out / "collection", "QRY", _LOG)
    assert called["ladder"] is False
    assert (out / LINEAGES_TSV).exists()


def test_deep_typing_calls_ladder_with_datasets_rows(tmp_path, monkeypatch):
    params, out = _panel(tmp_path)
    params.deep_typing = True
    (out / "ncbi_lineages.tsv").write_text("ACC1\tDENV1\n")
    seen = {}

    def fake_ladder(genomes, **kwargs):
        seen.update(kwargs)
        return [("ACC1", "DENV1", "ncbi-datasets")]

    monkeypatch.setattr(it, "assign_lineages", fake_ladder)
    lineage_map, _ = _type_panel(params, out / "collection", "QRY", _LOG)
    assert seen["datasets_rows"] == [("ACC1", "DENV1")]
    assert lineage_map["ACC1"] == "DENV1"
