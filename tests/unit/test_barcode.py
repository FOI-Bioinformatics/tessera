"""Clade-barcode caller: per-clade markers and lineage-attributed region calling."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from tessera.recomb.analyze import analyze
from tessera.recomb.barcode import clade_markers
from tessera.recomb.regions import RegionParams, call_regions
from tessera.recomb.similarity import compute_similarity

from ..conftest import write_fasta

L = 3000
LMAP = {"a1": "A", "a2": "A", "b1": "B", "b2": "B", "c1": "C", "c2": "C"}


def _clade_seqs() -> dict[str, list[str]]:
    """Three clades over a monomorphic 'T' background, each fixed for a distinctive base at
    its own marker columns (A at 0,30,..; C at 10,40,..; G at 20,50,..)."""
    a = ["T"] * L
    b = ["T"] * L
    c = ["T"] * L
    for p in range(0, L, 30):
        a[p] = "A"
    for p in range(10, L, 30):
        b[p] = "C"
    for p in range(20, L, 30):
        c[p] = "G"
    return {"A": a, "B": b, "C": c}


def _rows(seqs: dict[str, str]) -> dict[str, np.ndarray]:
    return {k: np.frombuffer(s.encode(), dtype=np.uint8) for k, s in seqs.items()}


def test_clade_markers_picks_characteristic_columns() -> None:
    cl = _clade_seqs()
    rows = _rows({"query": "".join(cl["A"]), "a1": "".join(cl["A"]), "a2": "".join(cl["A"]),
                  "b1": "".join(cl["B"]), "b2": "".join(cl["B"]),
                  "c1": "".join(cl["C"]), "c2": "".join(cl["C"])})
    cols, alleles, rep = clade_markers(rows, "query", LMAP)
    assert set(cols) == {"A", "B", "C"}
    # A's markers are exactly the ≡0 (mod 30) columns, each carrying 'A'
    assert cols["A"].tolist() == list(range(0, L, 30))
    assert set(alleles["A"].tolist()) == {ord("A")}
    assert rep["A"] in ("a1", "a2")


def test_clade_markers_empty_when_untyped_or_single_clade() -> None:
    cl = _clade_seqs()
    rows = _rows({"query": "".join(cl["A"]), "a1": "".join(cl["A"]), "a2": "".join(cl["A"])})
    assert clade_markers(rows, "query", LMAP) == ({}, {}, {})        # one clade
    assert clade_markers(rows, "query", None) == ({}, {}, {})        # untyped


def _panel(tmp_path: Path) -> Path:
    cl = _clade_seqs()
    a, b = "".join(cl["A"]), "".join(cl["B"])
    query = a[:1000] + b[1000:2000] + a[2000:]  # A backbone, B insert in the middle
    records = {
        "query": query,
        "a1": a, "a2": a, "b1": b, "b2": b, "c1": "".join(cl["C"]), "c2": "".join(cl["C"]),
    }
    return write_fasta(tmp_path / "panel.fasta", records)


def test_barcode_attributes_the_insert_to_its_clade(tmp_path: Path) -> None:
    result = compute_similarity(str(_panel(tmp_path)), "query", window_size=300, window_step=30)
    params = RegionParams.with_defaults(300, method="barcode", lineage_map=LMAP)
    regions, major, _ = call_regions(result, analyze(result), 300, params)
    assert major in ("a1", "a2")  # backbone clade A
    assert len(regions) == 1
    r = regions[0]
    assert r.minor_parent in ("b1", "b2")  # donor clade B
    assert r.query_start < 1200 and r.query_end > 1800  # spans the true insert
    assert r.methods == ("barcode",)


def test_barcode_silent_on_untyped_panel(tmp_path: Path) -> None:
    result = compute_similarity(str(_panel(tmp_path)), "query", window_size=300, window_step=30)
    params = RegionParams.with_defaults(300, method="barcode")  # no lineage map
    regions, _, _ = call_regions(result, analyze(result), 300, params)
    assert regions == []
