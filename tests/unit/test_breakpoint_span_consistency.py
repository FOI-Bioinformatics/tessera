"""A region's breakpoint interval must be reconcilable with its own span (issue #96).

The interval comes from the HMM's window grid and posterior; the reported span starts
at the first column that distinguishes the parents. On a panel where those sites are
sparse the two can part company, and on the mpox run behind the issue the interval sat
8.6 kb before the region it belonged to and was 5 kb wide. Both numbers were
individually defensible; printing them side by side was what misled.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from tessera.recomb.regions import breakpoint_meets_span
from tessera.recomb.run import RecombParams, run_recomb

# --- the rule itself -------------------------------------------------------

def test_an_interval_short_of_its_span_is_extended_to_reach_it() -> None:
    # The mpox numbers from the issue, in query coordinates.
    assert breakpoint_meets_span(50680, 55657, 64305) == (50680, 64305)


def test_an_interval_that_already_reaches_the_span_is_untouched() -> None:
    # The ordinary case: 59 of the 70 regions on disk when the issue was filed.
    assert breakpoint_meets_span(3999, 4004, 4000) == (3999, 4004)
    # Adjacent-but-outside, the benign case: extended by the few bp that separate them.
    assert breakpoint_meets_span(3950, 3998, 4000) == (3950, 4000)


def test_the_interval_is_never_narrowed() -> None:
    """Clipping into the span would discard real uncertainty before the first
    informative column, so the lower bound stays where the posterior put it."""
    lo, hi = breakpoint_meets_span(1000, 9000, 4000)
    assert (lo, hi) == (1000, 9000)


def test_a_first_segment_has_no_interval_to_extend() -> None:
    assert breakpoint_meets_span(None, None, 4000) == (None, None)


# --- the property, on real runs -------------------------------------------

def _regions(example_data: Path, tmp_path: Path, msa: str, **kwargs) -> list[dict[str, str]]:
    out = tmp_path / msa
    run_recomb(
        RecombParams(msa=example_data / msa, output=out, query="query",
                     window_size=300, window_step=30, plot_format="png", **kwargs),
        logging.getLogger("tessera.test"),
    )
    with (out / "recombination_regions.tsv").open(encoding="utf-8") as fo:
        return list(csv.DictReader(fo, delimiter="\t"))


def test_no_reported_interval_falls_outside_its_own_region(example_data, tmp_path) -> None:
    """Whatever the caller and whatever the panel, an interval that does not even
    touch the region it belongs to cannot be reconciled by a reader."""
    rows: list[dict[str, str]] = []
    for msa in ("divergent.msa.fasta", "cryptic_insert.msa.fasta"):
        rows += _regions(example_data, tmp_path, msa, methods=("hmm", "3seq", "maxchi", "bootscan"))
    assert rows, "no regions to check"

    checked = 0
    for row in rows:
        if row["breakpoint_lo"] == "NA":
            continue
        lo, hi = int(row["breakpoint_lo"]), int(row["breakpoint_hi"])
        start, end = int(row["query_start"]), int(row["query_end"])
        assert lo <= hi, row
        # The interval and the span must overlap or at least abut: hi below start, or
        # lo above end, is the inconsistency issue #96 reported.
        assert hi >= start, f"interval ends {start - hi} bp before its region: {row}"
        assert lo <= end, f"interval starts {lo - end} bp after its region: {row}"
        checked += 1
    assert checked, "every region reported NA, so the property was not exercised"
