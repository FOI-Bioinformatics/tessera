"""A low-information coverage gap is not reported as a possible missing reference.

Audit 4, item B3: a stretch of the query that is mostly ``N`` has too few comparable
bases to judge, which ``coverage.py`` labels ``low_information``. The report headline
still called it a possible missing reference and advised searching NCBI.
"""

from __future__ import annotations

import csv
import random
from pathlib import Path

from tessera.recomb.run import RecombParams, run_recomb


def _query_with_an_unsequenced_stretch(tmp_path: Path) -> Path:
    """2.5 kb of very poor sequencing: 97 % N, the surviving bases random miscalls."""
    rng = random.Random(11)
    anc = [rng.choice("ACGT") for _ in range(8000)]

    def mut(seq: list[str], rate: float) -> list[str]:
        return [rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c
                for c in seq]

    pa, pb, og = mut(anc, 0.05), mut(anc, 0.05), mut(anc, 0.12)
    query = list(pa)
    for i in range(3000, 5500):
        query[i] = "N" if rng.random() < 0.97 else rng.choice("ACGT")
    msa = tmp_path / "lowinfo.fasta"
    msa.write_text("".join(
        f">{name}\n{''.join(seq)}\n"
        for name, seq in (("query", query), ("parent_A", pa), ("parent_B", pb),
                          ("outgroup", og))
    ))
    return msa


def test_low_information_gap_is_not_a_missing_reference(tmp_path: Path, logger) -> None:
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=_query_with_an_unsequenced_stretch(tmp_path), output=out,
                     query="query", plot_format="png"),
        logger,
    )
    lines = (out / "coverage_gaps.tsv").read_text().splitlines()[1:]  # skip the preamble
    kinds = {row["kind"] for row in csv.DictReader(lines, delimiter="\t")}
    assert kinds == {"low_information"}, kinds  # the fixture holds no divergent stretch
    report = (out / "report.html").read_text()
    assert "Possible missing reference" not in report
    assert "find-references" not in report.split("</header>")[0]
    assert "too few comparable bases" in report.split("</header>")[0]


def test_console_table_does_not_advise_a_search_for_a_low_information_gap() -> None:
    from tessera.recomb.coverage import CoverageGap
    from tessera.recomb.report_text import print_coverage

    lines: list[str] = []
    gap = CoverageGap(msa_start=0, msa_end=1000, query_start=0, query_end=1000,
                      length_bp=1000, n_windows=5, best_label="ref", mean_best=0.4,
                      kind="low_information")
    print_coverage([gap], 0.95, echo=lines.append)
    text = "\n".join(lines)
    assert "find-references" not in text
    assert "too few comparable bases" in text


def test_divergent_and_low_information_gaps_get_separate_notes() -> None:
    """Review focus 3: with both kinds present the missing-reference caveat counts only
    the divergent gap, and the low-information note still appears."""
    from tessera.recomb.coverage import CoverageGap
    from tessera.recomb.report_context import ReportContext
    from tessera.recomb.report_html import _caveat_html, _low_information_html

    divergent = CoverageGap(msa_start=0, msa_end=1500, query_start=0, query_end=1500,
                            length_bp=1500, n_windows=6, best_label="ref", mean_best=0.70,
                            kind="divergent")
    low = CoverageGap(msa_start=3000, msa_end=5500, query_start=3000, query_end=5500,
                      length_bp=2500, n_windows=20, best_label="ref", mean_best=0.40,
                      kind="low_information")
    ctx = ReportContext(coverage_gaps=[divergent, low], coverage_threshold=0.95)
    caveat = _caveat_html(ctx.missing_reference_gaps, 0.95)
    assert "1 region" in caveat and "1.5 kb" in caveat
    assert "Too little sequence to judge" in _low_information_html(ctx.low_information_gaps)
    assert ctx.caveat_gaps == [divergent, low]  # the mosaic and plots still hatch both


def test_console_note_for_a_local_dip_does_not_advise_a_search() -> None:
    """Issue #82: a local identity dip is not a missing reference.

    The closest reference there still matches the query closely, so the note must
    say what the stretch is without sending the reader to ``find-references``.
    """
    from tessera.recomb.coverage import LOCAL_DIP_KIND, CoverageGap
    from tessera.recomb.report_text import print_coverage

    lines: list[str] = []
    gap = CoverageGap(msa_start=0, msa_end=1500, query_start=0, query_end=1500,
                      length_bp=1500, n_windows=6, best_label="ref", mean_best=0.997,
                      kind=LOCAL_DIP_KIND)
    print_coverage([gap], 0.998, echo=lines.append)
    text = "\n".join(lines)

    assert "find-references" not in text
    assert "local_dip" in text
    assert "still matches" in text


def test_a_local_dip_gets_its_own_note_and_no_missing_reference_headline() -> None:
    """It is reported, but not as a gap in the panel and not as uncertainty."""
    from tessera.recomb.coverage import LOCAL_DIP_KIND, CoverageGap
    from tessera.recomb.report_context import ReportContext
    from tessera.recomb.report_html import (
        _caveat_html,
        _local_dip_html,
        _low_information_html,
    )

    dip = CoverageGap(msa_start=1000, msa_end=2500, query_start=1000, query_end=2500,
                      length_bp=1500, n_windows=6, best_label="ref", mean_best=0.997,
                      kind=LOCAL_DIP_KIND)
    ctx = ReportContext(coverage_gaps=[dip], coverage_threshold=0.998)

    assert ctx.local_dip_gaps == [dip]
    assert ctx.missing_reference_gaps == []
    assert ctx.low_information_gaps == []
    assert _caveat_html(ctx.missing_reference_gaps, 0.998) == ""
    assert _low_information_html(ctx.low_information_gaps) == ""
    note = _local_dip_html(ctx.local_dip_gaps)
    assert "1.5 kb" in note
    assert "missing" not in note.lower().replace("not a missing", "")


def test_the_plot_does_not_label_a_local_dip_as_low_coverage() -> None:
    """Coverage is high there; what is low is the local alignment quality."""
    from tessera.recomb.coverage import LOCAL_DIP_KIND, CoverageGap
    from tessera.recomb.report_plots import _gap_annotation

    dip = CoverageGap(msa_start=0, msa_end=1500, query_start=0, query_end=1500,
                      length_bp=1500, n_windows=6, best_label="ref", mean_best=0.997,
                      kind=LOCAL_DIP_KIND)
    divergent = CoverageGap(msa_start=0, msa_end=1500, query_start=0, query_end=1500,
                            length_bp=1500, n_windows=6, best_label="ref", mean_best=0.70,
                            kind="divergent")

    assert _gap_annotation(divergent) == "low coverage"
    assert _gap_annotation(dip) == "local dip"
