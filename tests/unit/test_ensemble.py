"""Ensemble consensus merge and the --method parser."""

from __future__ import annotations

import pytest

from tessera.core.errors import UserInputError
from tessera.recomb.ensemble import (
    consensus_regions,
    filter_by_agreement,
    reconcile_major,
)
from tessera.recomb.regions import CALLERS, Region, parse_methods


def mk(minor: str, qs: int, qe: int, method: str, *,
       major: str = "A", support: float | None = None,
       qvalue: float | None = None, breakpoint_lo: int | None = None) -> Region:
    return Region(
        minor_parent=minor, major_parent=major, msa_start=qs, msa_end=qe,
        query_start=qs, query_end=qe, n_windows=1,
        mean_sim_minor=0.99, mean_sim_major=0.90, margin=0.0,
        support=support, qvalue=qvalue, breakpoint_lo=breakpoint_lo,
        breakpoint_hi=breakpoint_lo, methods=(method,),
    )


# --- parse_methods ---------------------------------------------------------

def test_parse_methods_all_expands_to_every_caller() -> None:
    assert parse_methods("all") == CALLERS


def test_parse_methods_orders_by_callers_and_dedups() -> None:
    # input order ignored; CALLERS order imposed; duplicates collapsed
    assert parse_methods("3seq,hmm") == ("hmm", "3seq")
    assert parse_methods("hmm,hmm") == ("hmm",)


def test_parse_methods_rejects_unknown_token() -> None:
    with pytest.raises(UserInputError):
        parse_methods("hmm,bogus")
    with pytest.raises(UserInputError):
        parse_methods("")


# --- consensus_regions -----------------------------------------------------

def test_single_method_round_trips_unchanged() -> None:
    regions = [mk("B", 100, 200, "hmm", support=0.9, qvalue=1e-9)]
    merged, breakdown = consensus_regions({"hmm": regions}, major="A")
    assert len(merged) == 1
    assert merged[0].methods == ("hmm",)
    assert merged[0].query_start == 100 and merged[0].query_end == 200
    assert breakdown[0]["methods"] == ("hmm",)


def test_overlapping_same_minor_merges_with_agreement_and_best_support() -> None:
    hmm = mk("B", 100, 220, "hmm", support=0.8, qvalue=1e-6, breakpoint_lo=110)
    seq = mk("B", 120, 200, "3seq", support=0.95, qvalue=1e-12)
    merged, breakdown = consensus_regions(
        {"hmm": [hmm], "3seq": [seq]}, major="A",
    )
    assert len(merged) == 1
    r = merged[0]
    assert r.methods == ("hmm", "3seq")          # both callers, CALLERS order
    assert (r.query_start, r.query_end) == (100, 220)  # union span
    assert r.support == 0.95 and r.qvalue == 1e-12     # most significant member
    assert r.breakpoint_lo == 110                # HMM breakpoint preferred
    assert breakdown[0]["methods"] == ("hmm", "3seq")
    assert set(breakdown[0]["per_method_support"]) == {"hmm", "3seq"}


def test_lineage_map_merges_same_lineage_different_genomes() -> None:
    # Two callers pick different genomes (G1, G2) of one lineage "B" over the same span.
    g1 = mk("G1", 100, 200, "hmm", qvalue=1e-6)
    g2 = mk("G2", 110, 190, "3seq", qvalue=1e-9)
    lmap = {"G1": "B", "G2": "B"}
    # Without typing, different genome labels -> two separate single-method regions.
    plain, _ = consensus_regions({"hmm": [g1], "3seq": [g2]}, major="A")
    assert len(plain) == 2
    # With the lineage map, they are one event called by both methods (agreement).
    merged, _ = consensus_regions({"hmm": [g1], "3seq": [g2]}, major="A", lineage_map=lmap)
    assert len(merged) == 1
    assert merged[0].methods == ("hmm", "3seq")


def test_overlapping_different_minor_stays_separate() -> None:
    a = mk("B", 100, 200, "hmm", qvalue=1e-9)
    b = mk("C", 150, 250, "3seq", qvalue=1e-9)
    merged, _ = consensus_regions({"hmm": [a], "3seq": [b]}, major="A")
    assert len(merged) == 2
    assert {r.minor_parent for r in merged} == {"B", "C"}
    assert all(len(r.methods) == 1 for r in merged)


def test_parent_free_support_set_when_overlapping_rmin_interval() -> None:
    regions = [mk("B", 100, 200, "hmm", qvalue=1e-9)]
    merged, _ = consensus_regions(
        {"hmm": regions}, major="A", rmin_intervals=[(150, 160)],
    )
    assert merged[0].parent_free_support is True
    merged2, _ = consensus_regions(
        {"hmm": regions}, major="A", rmin_intervals=[(900, 950)],
    )
    assert merged2[0].parent_free_support is False


def test_agreement_raises_confidence_to_high() -> None:
    from tessera.recomb.report_html import _confidence
    # a single marginal caller is not high...
    one = [mk("B", 100, 200, "3seq", support=0.55, qvalue=0.03)]
    assert _confidence(one, []) == "moderate"
    # ...but the same region called by two methods is high on agreement alone
    two = [mk("B", 100, 200, "hmm", support=0.55, qvalue=0.03)]
    two[0].methods = ("hmm", "3seq")
    assert _confidence(two, []) == "high"


def test_reconcile_major_prefers_hmm() -> None:
    canonical, per_major = reconcile_major({"hmm": "A1", "3seq": "A2"})
    assert canonical == "A1"
    assert per_major == {"hmm": "A1", "3seq": "A2"}


def test_reconcile_major_available_with_no_regions() -> None:
    # a caller reports its backbone even when it finds no recombination
    canonical, per_major = reconcile_major({"hmm": "backbone", "3seq": "backbone"})
    assert canonical == "backbone"


def test_reconcile_major_overrides_zero_win_hmm_major() -> None:
    # HMM dropped the true backbone A1 as a sibling and picked K (0 windows); the windowed
    # vote clearly favours A1, so the reconciled backbone is A1, not the zero-win HMM major.
    canonical, per_major = reconcile_major(
        {"hmm": "K", "3seq": "A1", "maxchi": "A1"}, window_wins={"A1": 298}
    )
    assert canonical == "A1"
    assert per_major["hmm"] == "K"   # per-method disagreement still recorded


def test_reconcile_major_keeps_hmm_major_that_wins_windows() -> None:
    # the masking-twin case: HMM's major wins its own segment's windows -> left untouched.
    canonical, _ = reconcile_major(
        {"hmm": "parent", "3seq": "twin"}, window_wins={"parent": 40, "twin": 260}
    )
    assert canonical == "parent"


def test_reconcile_major_without_window_wins_keeps_hmm() -> None:
    canonical, _ = reconcile_major({"hmm": "A1", "3seq": "A2"})
    assert canonical == "A1"


# --- agreement gate (filter_by_agreement) ----------------------------------

def test_agreement_gate_drops_regions_below_the_threshold() -> None:
    # B over 1000-2000 is found by two callers; C over 5000-6000 by one only.
    per_method = {
        "hmm": [mk("B", 1000, 2000, "hmm")],
        "maxchi": [mk("B", 1000, 2000, "maxchi")],
        "bootscan": [mk("C", 5000, 6000, "bootscan")],
    }
    regions, breakdown = consensus_regions(per_method, major="A")
    kept, kept_breakdown, suppressed = filter_by_agreement(regions, breakdown, 2)

    assert [r.minor_parent for r in kept] == ["B"]
    assert kept[0].methods == ("hmm", "maxchi")
    assert suppressed == 1
    # the breakdown must not advertise a region the region table no longer has
    assert [b["minor_parent"] for b in kept_breakdown] == ["B"]


def test_agreement_gate_of_one_keeps_every_region() -> None:
    per_method = {
        "hmm": [mk("B", 1000, 2000, "hmm")],
        "bootscan": [mk("C", 5000, 6000, "bootscan")],
    }
    regions, breakdown = consensus_regions(per_method, major="A")
    kept, kept_breakdown, suppressed = filter_by_agreement(regions, breakdown, 1)
    assert sorted(r.minor_parent for r in kept) == ["B", "C"]
    assert len(kept_breakdown) == 2
    assert suppressed == 0
