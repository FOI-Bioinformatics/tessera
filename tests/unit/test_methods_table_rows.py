"""recombination_methods.tsv keeps one row per row of recombination_regions.tsv.

Audit 4, smaller reporting item: donor-absent regions were appended after the method
breakdown was final, so the methods file had fewer rows than the regions file and a
positional join paired the wrong rows.
"""

from __future__ import annotations

import csv
import random
from pathlib import Path

from tessera.recomb.run import RecombParams, run_recomb


def _panel_with_an_absent_donor(tmp_path: Path) -> Path:
    """A parent_B tract at 2000-4000 and a 2.5 kb tract from a donor not in the panel."""
    rng = random.Random(31)
    anc = [rng.choice("ACGT") for _ in range(12000)]

    def mut(seq: list[str], rate: float) -> list[str]:
        return [rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c
                for c in seq]

    pa, pb, og = mut(anc, 0.04), mut(anc, 0.04), mut(anc, 0.10)
    ghost = mut(anc, 0.30)
    query = list(pa)
    query[2000:4000] = pb[2000:4000]
    query[7000:9500] = ghost[7000:9500]
    msa = tmp_path / "absent.fasta"
    msa.write_text("".join(
        f">{n}\n{''.join(s)}\n"
        for n, s in (("query", query), ("parent_A", pa), ("parent_B", pb), ("outgroup", og))
    ))
    return msa


def _rows(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(encoding="utf-8"), delimiter="\t"))


def test_methods_rows_match_regions_rows_with_a_donor_absent_region(
    tmp_path: Path, logger
) -> None:
    out = tmp_path / "out"
    run_recomb(RecombParams(msa=_panel_with_an_absent_donor(tmp_path), output=out,
                            query="query", plot_format="png"), logger)
    regions = _rows(out / "recombination_regions.tsv")
    methods = _rows(out / "recombination_methods.tsv")
    assert any(r["donor_absent"] == "yes" for r in regions), "fixture needs an absent row"
    assert len(methods) == len(regions)
    for region, row in zip(regions, methods, strict=True):
        assert (row["query_start"], row["query_end"]) == (
            region["query_start"], region["query_end"])
        if region["donor_absent"] == "yes":
            assert {row[m] for m in ("hmm", "3seq", "maxchi", "bootscan")} == {"no"}
    assert "donor absent" in (out / "report.html").read_text()


def test_methods_rows_match_regions_rows_under_the_agreement_gate(
    tmp_path: Path, logger
) -> None:
    """Review focus 4: --min-methods 2 removes rows before the donor-absent rows are
    added; the two files must still pair row for row."""
    out = tmp_path / "out"
    run_recomb(RecombParams(msa=_panel_with_an_absent_donor(tmp_path), output=out,
                            query="query", plot_format="png", min_methods=2), logger)
    regions = _rows(out / "recombination_regions.tsv")
    methods = _rows(out / "recombination_methods.tsv")
    assert any(r["donor_absent"] == "yes" for r in regions)
    assert [(r["query_start"], r["query_end"]) for r in regions] == [
        (m["query_start"], m["query_end"]) for m in methods]
