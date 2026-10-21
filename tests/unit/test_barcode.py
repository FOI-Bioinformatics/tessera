"""Clade-barcode caller: per-clade markers and lineage-attributed region calling."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path

import numpy as np
import pytest

from tessera.core.errors import UserInputError
from tessera.recomb.analyze import analyze
from tessera.recomb.barcode import clade_markers
from tessera.recomb.regions import RegionParams, call_regions
from tessera.recomb.run import RecombParams, run_recomb
from tessera.recomb.similarity import compute_similarity

from ..conftest import recombinant_msa, write_fasta

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


# --- "could not run" is not "found nothing" --------------------------------

def test_barcode_names_no_major_parent_when_it_cannot_run(tmp_path: Path) -> None:
    """It used to return the first record of the alignment as the major parent, which a
    barcode-only run then reported as the backbone."""
    result = compute_similarity(str(_panel(tmp_path)), "query", window_size=300, window_step=30)
    params = RegionParams.with_defaults(300, method="barcode")  # no lineage map
    regions, major, _ = call_regions(result, analyze(result), 300, params)
    assert regions == []
    assert major is None


def _run(tmp_path: Path, logger, **kwargs) -> Path:
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=recombinant_msa(tmp_path, recombinant=True), output=out,
                     query="query", plot_format="png", **kwargs),
        logger,
    )
    return out


def test_barcode_only_run_on_an_untyped_panel_is_refused(tmp_path: Path, logger) -> None:
    """Nothing was tested, so there is no result to report -- not a clean negative."""
    with pytest.raises(UserInputError, match="typed references"):
        _run(tmp_path, logger, methods=("barcode",))
    assert not (tmp_path / "out" / "report.html").exists()


def _collecting_logger() -> tuple[logging.Logger, list[str]]:
    """A logger outside the ``tessera`` hierarchy that keeps its warnings in a list.

    Not ``caplog``: once a CLI test has configured the ``tessera`` logger it stops
    propagating to the root logger, and whether that has happened depends on test order.
    """
    messages: list[str] = []

    class _Collect(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            messages.append(record.getMessage())

    log = logging.getLogger("test_barcode.collect")
    log.handlers = [_Collect(level=logging.WARNING)]
    log.propagate = False
    return log, messages


def test_ensemble_reports_barcode_as_not_run_on_an_untyped_panel(tmp_path: Path) -> None:
    log, warnings = _collecting_logger()
    out = _run(tmp_path, log, methods=("3seq", "barcode"))
    assert any("barcode" in message and "could not run" in message for message in warnings)

    rows = list(csv.DictReader((out / "recombination_methods.tsv").open(), delimiter="\t"))
    assert rows, "3seq should still call the insert"
    assert {row["3seq"] for row in rows} == {"yes"}
    assert {row["barcode"] for row in rows} == {"not run"}  # was "no"

    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert "barcode" in run["callers not run"]
    assert "not run" in (out / "report.html").read_text()


def test_agreement_gate_counts_only_the_callers_that_ran(tmp_path: Path, logger) -> None:
    """--min-methods is clamped to the callers actually run. A caller that could not run
    cannot agree, so it must not make the gate unreachable."""
    out = _run(tmp_path, logger, methods=("3seq", "barcode"), min_methods=2)
    rows = list(csv.DictReader((out / "recombination_regions.tsv").open(), delimiter="\t"))
    assert [row["methods"] for row in rows] == ["3seq"]


def test_a_lowered_agreement_gate_is_logged_and_recorded(tmp_path: Path) -> None:
    """The user asked for two agreeing callers and only one could run. The regions are
    then single-caller regions; the log and the provenance must say the gate was 1."""
    log, warnings = _collecting_logger()
    out = _run(tmp_path, log, methods=("3seq", "barcode"), min_methods=2)
    assert any("--min-methods 2" in m and "1 caller" in m for m in warnings)
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["min methods (agreement gate)"] == "1 (requested 2; 1 caller ran)"


def test_barcode_not_run_on_a_typed_panel_without_marked_clades(tmp_path: Path) -> None:
    """Typed, but every reference is in one clade: there are no clade markers to
    compete, so the caller cannot run and the warning says why."""
    log, warnings = _collecting_logger()
    one_clade = {"A": "L1", "B": "L1", "other": "L1"}
    out = _run(tmp_path, log, methods=("3seq", "barcode"), lineage_map=one_clade)
    assert any("fewer than two typed clades" in message for message in warnings)
    rows = list(csv.DictReader((out / "recombination_methods.tsv").open(), delimiter="\t"))
    assert {row["barcode"] for row in rows} == {"not run"}
    # The panel is typed, so "needs typed references" would be the wrong reason.
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert "fewer than two typed clades" in run["callers not run"]
    assert "needs typed references" not in run["callers not run"]
    report = (out / "report.html").read_text()
    assert "fewer than two typed clades" in report
    assert "none were available" not in report
