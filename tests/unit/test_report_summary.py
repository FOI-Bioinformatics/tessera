"""The report's headline numbers agree with the regions they summarise."""

from __future__ import annotations

import numpy as np

from tessera.recomb.regions import Region
from tessera.recomb.report_html import _summary, _verdict_html


class _Result:
    """The one attribute ``_summary`` reads: a 6000-base, gap-free query."""

    query = "query"
    query_cumulative = np.arange(6001)


def _region(minor: str, start: int, end: int, *, donor_absent: bool = False) -> Region:
    return Region(
        minor_parent=minor, major_parent="backbone", msa_start=start, msa_end=end,
        query_start=start, query_end=end, n_windows=5,
        mean_sim_minor=0.99, mean_sim_major=0.92, margin=0.07,
        qvalue=1e-9, support=0.95, methods=("hmm", "3seq"), donor_absent=donor_absent,
    )


def test_overlapping_regions_are_counted_once() -> None:
    """The ensemble keeps overlapping regions that name different donors separate.
    Adding their lengths reported a 2.2 kb union as 3.3 kb (54.6 % of the query)."""
    regions = [_region("donorB", 1900, 3100), _region("donorA", 2025, 4100)]
    s = _summary(_Result(), regions, ["backbone", "donorA", "donorB"])
    assert s["recomb_bp"] == 2200  # the union 1900-4100, not 1200 + 2075
    assert round(s["pct"], 1) == 36.7
    assert s["n_regions"] == 2  # both regions are still listed


def test_recombinant_fraction_never_exceeds_the_query() -> None:
    regions = [_region(f"donor{i}", 0, 6000) for i in range(3)]
    s = _summary(_Result(), regions, ["backbone"])
    assert s["recomb_bp"] == 6000
    assert s["pct"] == 100.0


def test_touching_and_nested_regions() -> None:
    touching = [_region("donorA", 1000, 2000), _region("donorB", 2000, 3000)]
    assert _summary(_Result(), touching, ["backbone"])["recomb_bp"] == 2000
    nested = [_region("donorA", 1000, 4000), _region("donorB", 2000, 2500)]
    assert _summary(_Result(), nested, ["backbone"])["recomb_bp"] == 3000


def test_disjoint_regions_still_add_up() -> None:
    regions = [_region("donorA", 500, 1500), _region("donorB", 3000, 3500)]
    assert _summary(_Result(), regions, ["backbone"])["recomb_bp"] == 1500


def test_donor_absent_regions_stay_out_of_the_recombinant_span() -> None:
    regions = [_region("donorA", 500, 1500), _region("x", 3000, 5000, donor_absent=True)]
    s = _summary(_Result(), regions, ["backbone"])
    assert s["recomb_bp"] == 1000
    assert s["n_absent"] == 1


def test_verdict_states_the_union() -> None:
    regions = [_region("donorB", 1900, 3100), _region("donorA", 2025, 4100)]
    s = _summary(_Result(), regions, ["backbone", "donorA", "donorB"])
    verdict = _verdict_html(s, "query", {})
    assert "2.2 kb" in verdict
    assert "36.7%" in verdict
