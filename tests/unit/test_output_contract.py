"""The TSV output contract: column schemas and what each column means.

These are characterization tests. They pin the shape of the files downstream analysis
consumes, so a change to the contract has to be made deliberately rather than by
accident -- and so the columns that used to mean two different things cannot silently
go back to doing so.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from tessera.recomb.ensemble import consensus_regions
from tessera.recomb.regions import Region
from tessera.recomb.report_text import write_regions_tsv

REGION_COLUMNS = [
    "minor_parent", "major_parent", "msa_start", "msa_end",
    "query_start", "query_end", "length_bp", "length_msa", "n_windows",
    "mean_sim_minor", "mean_sim_major", "margin",
    "support", "pvalue", "qvalue", "test", "statistic",
    "posterior", "breakpoint_lo", "breakpoint_hi",
    "donor_undercovered", "donor_absent", "methods", "parent_free_support",
]


def mk(minor="B", major="A", *, msa=(0, 300), query=(0, 100), method="hmm",
       pvalue=None, qvalue=None, support=None, test="", statistic="") -> Region:
    return Region(
        minor_parent=minor, major_parent=major,
        msa_start=msa[0], msa_end=msa[1],
        query_start=query[0], query_end=query[1],
        n_windows=1, mean_sim_minor=0.99, mean_sim_major=0.90, margin=0.09,
        support=support, pvalue=pvalue, qvalue=qvalue,
        test=test, statistic=statistic, methods=(method,),
    )


# --- lengths have exactly one definition each -----------------------------

def test_length_bp_is_query_bases_and_length_msa_is_columns() -> None:
    """The two were one column that meant the MSA span for caller rows and the query
    span for merged ensemble rows."""
    r = mk(msa=(1000, 4000), query=(900, 2000))
    assert r.length_bp == 1100  # query bases
    assert r.length_msa == 3000  # alignment columns, gaps included


def test_lengths_cannot_drift_from_the_coordinates() -> None:
    # Derived, not stored, so no caller can set them inconsistently.
    r = mk(msa=(0, 500), query=(10, 60))
    assert r.length_bp == r.query_end - r.query_start
    assert r.length_msa == r.msa_end - r.msa_start


# --- each row names the test behind its own p-value -----------------------

def test_merged_row_reports_the_test_its_pvalue_came_from() -> None:
    """A merged region lists every agreeing caller but keeps one member's p-value.

    Without the test column travelling with it, the row shows four methods beside a
    p-value from one of them and a reader cannot tell which.
    """
    weak = mk(method="hmm", qvalue=0.04, support=0.6,
              test="sign test on discordant sites (one-sided)",
              statistic="fraction of discordant sites favouring the donor")
    strong = mk(method="3seq", qvalue=1e-9, support=0.9,
                test="3SEQ max-descent (exact)",
                statistic="fraction of tract sites matching the minor parent")

    merged, _breakdown = consensus_regions(
        {"hmm": [weak], "3seq": [strong]}, major="A"
    )
    assert len(merged) == 1
    row = merged[0]
    assert set(row.methods) == {"hmm", "3seq"}
    assert row.qvalue == 1e-9  # the most significant member's
    assert row.test == "3SEQ max-descent (exact)"  # ... and its test, not the other one
    assert "tract sites" in row.statistic


# --- the file itself ------------------------------------------------------

def _write(tmp_path: Path, regions: list[Region]) -> Path:
    write_regions_tsv(regions, tmp_path, logging.getLogger("tessera"))
    return tmp_path / "recombination_regions.tsv"


def test_regions_tsv_column_schema(tmp_path: Path) -> None:
    path = _write(tmp_path, [mk()])
    assert path.read_text().splitlines()[0].split("\t") == REGION_COLUMNS


def test_regions_tsv_parses_with_a_plain_reader(tmp_path: Path) -> None:
    """No comment preamble: this is the most widely parsed output, and a leading `#`
    would break every reader that does not know to skip it."""
    path = _write(tmp_path, [mk(), mk(minor="C")])
    assert not path.read_text().startswith("#")
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    assert len(rows) == 2
    assert rows[0]["minor_parent"] == "B"
    assert rows[1]["minor_parent"] == "C"


def test_regions_tsv_keeps_full_pvalue_precision(tmp_path: Path) -> None:
    # Rounding to two significant figures threw away precision a reader cannot recover
    # and made an underflow indistinguishable from an exact zero.
    tiny = 1.422084672242363e-12
    path = _write(tmp_path, [mk(pvalue=tiny, qvalue=tiny * 2)])
    row = next(csv.DictReader(path.open(), delimiter="\t"))
    assert float(row["pvalue"]) == tiny
    assert float(row["qvalue"]) == tiny * 2


def test_regions_tsv_marks_absent_values_na(tmp_path: Path) -> None:
    path = _write(tmp_path, [mk()])  # no p-value, no test: the barcode/heuristic case
    row = next(csv.DictReader(path.open(), delimiter="\t"))
    assert row["pvalue"] == "NA"
    assert row["qvalue"] == "NA"
    assert row["test"] == "NA"
    assert row["statistic"] == "NA"
