"""A coverage gap produced by a window straddling a called breakpoint is not a
missing reference.

A window that spans a breakpoint between two divergent parents matches neither parent
well on its own, so its best similarity dips below the coverage threshold even though
both parents are in the panel. Such a gap is relabelled ``breakpoint``; it must not
caveat the region, become a donor-absent region, or reach the report headline.
"""

from __future__ import annotations

import csv
import re

import numpy as np

from tessera.recomb.coverage import (
    BREAKPOINT_KIND,
    CoverageGap,
    gaps_as_regions,
    mark_breakpoint_gaps,
)
from tessera.recomb.regions import Region
from tessera.recomb.run import RecombParams, run_recomb

WIDTH = 2000
WINDOW = 400


def _enc(seq: str) -> np.ndarray:
    return np.frombuffer(seq.encode("ascii"), dtype=np.uint8)


def _rows(query: str) -> dict[str, np.ndarray]:
    """Parents that differ at every tenth column; the query is given by the caller."""
    major = "A" * WIDTH
    minor = "".join("C" if i % 10 == 0 else "A" for i in range(WIDTH))
    return {"query": _enc(query), "major": _enc(major), "minor": _enc(minor)}


def _mosaic() -> str:
    """Major up to column 1000, minor after it: a clean breakpoint at 1000."""
    minor = "".join("C" if i % 10 == 0 else "A" for i in range(WIDTH))
    return "A" * 1000 + minor[1000:]


def _region() -> Region:
    return Region(
        minor_parent="minor", major_parent="major", msa_start=1000, msa_end=WIDTH,
        query_start=1000, query_end=WIDTH, n_windows=10,
        mean_sim_minor=1.0, mean_sim_major=0.9, margin=0.1,
    )


def _gap(start: int, end: int, kind: str = "divergent") -> CoverageGap:
    return CoverageGap(
        msa_start=start, msa_end=end, query_start=start, query_end=end,
        length_bp=end - start, n_windows=2, best_label="major", mean_best=0.94, kind=kind,
    )


def test_gap_straddling_a_called_breakpoint_is_relabelled() -> None:
    gaps = [_gap(800, 1200)]
    n = mark_breakpoint_gaps(gaps, [_region()], _rows(_mosaic()), "query", WINDOW, 0.95)
    assert n == 1
    assert gaps[0].kind == BREAKPOINT_KIND


def test_breakpoint_gap_is_not_bridged_to_a_donor_absent_region() -> None:
    gaps = [_gap(800, 1200)]
    mark_breakpoint_gaps(gaps, [_region()], _rows(_mosaic()), "query", WINDOW, 0.95)

    class _Result:
        positions: list[int] = []
        similarities: dict[str, list[float]] = {}

    assert gaps_as_regions(gaps, _Result(), "major") == []


def test_gap_the_two_parents_do_not_explain_stays_divergent() -> None:
    # Next to the breakpoint, but the query carries a base neither parent has at one
    # column in five: a third source, not a straddling window.
    query = list(_mosaic())
    for i in range(800, 1200, 5):
        query[i] = "G"
    gaps = [_gap(800, 1200)]
    n = mark_breakpoint_gaps(gaps, [_region()], _rows("".join(query)), "query", WINDOW, 0.95)
    assert n == 0
    assert gaps[0].kind == "divergent"


def test_gap_far_from_every_region_boundary_stays_divergent() -> None:
    gaps = [_gap(100, 500)]  # more than one window away from 1000 and from 2000
    n = mark_breakpoint_gaps(gaps, [_region()], _rows(_mosaic()), "query", WINDOW, 0.95)
    assert n == 0
    assert gaps[0].kind == "divergent"


def test_low_information_gap_is_left_alone() -> None:
    gaps = [_gap(800, 1200, kind="low_information")]
    mark_breakpoint_gaps(gaps, [_region()], _rows(_mosaic()), "query", WINDOW, 0.95)
    assert gaps[0].kind == "low_information"


def test_donor_absent_region_does_not_explain_a_gap() -> None:
    region = _region()
    region.donor_absent = True
    gaps = [_gap(800, 1200)]
    assert mark_breakpoint_gaps(gaps, [region], _rows(_mosaic()), "query", WINDOW, 0.95) == 0


def test_region_naming_a_label_outside_the_alignment_is_skipped() -> None:
    # A region can carry a placeholder parent ("n/a" when no backbone was resolved).
    region = _region()
    region.major_parent = "n/a"
    gaps = [_gap(800, 1200)]
    assert mark_breakpoint_gaps(gaps, [region], _rows(_mosaic()), "query", WINDOW, 0.95) == 0
    assert gaps[0].kind == "divergent"


# --- the shipped example, end to end --------------------------------------

def test_clean_recombinant_is_not_reported_as_a_missing_reference(
    example_data, tmp_path, logger
) -> None:
    """`example_data/divergent.msa.fasta`: the donor matches the query exactly over the
    insert and four callers agree. The only coverage gaps are the windows straddling
    the two breakpoints."""
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=example_data / "divergent.msa.fasta", output=out, query="query",
                     window_size=300, window_step=30, plot_format="png"),
        logger,
    )
    regions = list(csv.DictReader((out / "recombination_regions.tsv").open(), delimiter="\t"))
    assert len(regions) == 1
    assert regions[0]["donor_undercovered"] == "no"
    assert regions[0]["donor_absent"] == "no"

    lines = (out / "coverage_gaps.tsv").read_text().splitlines()
    gaps = list(csv.DictReader([ln for ln in lines if not ln.startswith("#")], delimiter="\t"))
    assert gaps and {g["kind"] for g in gaps} == {BREAKPOINT_KIND}

    report = (out / "report.html").read_text()
    verdict = re.sub(r"<[^>]+>", "", re.search(r'<p class="verdict">.*?</p>', report).group(0))
    assert "high confidence" in " ".join(verdict.split())
    assert "Possible missing reference" not in report
