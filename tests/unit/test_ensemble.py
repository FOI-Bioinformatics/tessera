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


def test_region_is_not_relabelled_as_a_donor_into_itself() -> None:
    # The callers disagree on the backbone: the HMM reports H (it set the sibling S
    # aside), 3SEQ kept S and found an H tract in it. H is the canonical backbone, so
    # relabelling the 3SEQ region's major to H would read "H donated into H".
    hmm = mk("X", 100, 200, "hmm", major="H", qvalue=1e-6)
    seq = mk("H", 500, 600, "3seq", major="S", qvalue=1e-9)
    merged, _ = consensus_regions({"hmm": [hmm], "3seq": [seq]}, major="H")
    assert all(r.minor_parent != r.major_parent for r in merged)
    by_minor = {r.minor_parent: r.major_parent for r in merged}
    assert by_minor == {"X": "H", "H": "S"}  # the tested pair is kept for the odd one out


def test_same_lineage_regions_merge_transitively_across_labels() -> None:
    # A1 and A2 are one lineage. A1 has two separate regions; the A2 region overlaps
    # both. All three are one event; sorting by label must not leave two groups.
    lmap = {"A1": "L", "A2": "L"}
    first = mk("A1", 0, 10, "hmm", qvalue=1e-6)
    second = mk("A1", 20, 30, "hmm", qvalue=1e-6)
    bridge = mk("A2", 5, 25, "3seq", qvalue=1e-9)
    merged, breakdown = consensus_regions(
        {"hmm": [first, second], "3seq": [bridge]}, major="M", lineage_map=lmap
    )
    assert len(merged) == 1
    assert (merged[0].query_start, merged[0].query_end) == (0, 30)
    assert merged[0].methods == ("hmm", "3seq")
    assert len(breakdown) == 1


def _scan_result():
    """A 100-column alignment: B matches the query on [10, 40), A on [40, 100). Five
    20-column base-pair windows centred at 10, 30, 50, 70, 90."""
    import numpy as np

    from tessera.recomb.similarity import WindowSimilarity

    width = 100
    query = np.frombuffer(b"A" * width, dtype=np.uint8)
    minor = np.frombuffer(b"C" * 10 + b"A" * 30 + b"C" * 60, dtype=np.uint8)
    major = np.frombuffer(b"A" * 10 + b"C" * 30 + b"A" * 60, dtype=np.uint8)
    positions = [10, 30, 50, 70, 90]
    return WindowSimilarity(
        positions=positions, query_positions=positions,
        similarities={"A": [1.0] * 5, "B": [1.0] * 5}, query="q", width=width,
        query_cumulative=np.arange(width + 1),
        rows={"q": query, "A": major, "B": minor},
    )


def test_merged_identity_is_measured_over_the_merged_span() -> None:
    # The HMM member covers [10, 40) where B is identical to the query; the 3SEQ member
    # covers [30, 80). The merged row spans [10, 80) and must not report the HMM
    # member's 1.0 over a span of which it describes less than half.
    hmm = mk("B", 10, 40, "hmm", qvalue=1e-9)
    hmm.mean_sim_minor, hmm.mean_sim_major = 1.0, 0.0
    seq = mk("B", 30, 80, "3seq", qvalue=1e-3)
    merged, _ = consensus_regions({"hmm": [hmm], "3seq": [seq]}, major="A",
                                  scan=(_scan_result(), 20))
    assert len(merged) == 1
    region = merged[0]
    assert (region.msa_start, region.msa_end) == (10, 80)
    assert region.mean_sim_minor == round(30 / 70, 4)
    assert region.mean_sim_major == round(40 / 70, 4)
    assert region.margin == round(30 / 70 - 40 / 70, 4)


def test_merged_identity_is_column_weighted_when_parts_differ_tenfold() -> None:
    # A 100-column alignment, B identical to the query on [0, 10) and different on
    # [10, 100). A 10-column member and a 90-column member merge to [0, 100); the
    # identity is 10/100, not the mean of the parts' own figures (1.0 and 0.0 -> 0.5)
    # and not the longer part's.
    import numpy as np

    from tessera.recomb.similarity import WindowSimilarity

    width = 100
    query = np.frombuffer(b"A" * width, dtype=np.uint8)
    minor = np.frombuffer(b"A" * 10 + b"C" * 90, dtype=np.uint8)
    major = np.frombuffer(b"A" * width, dtype=np.uint8)
    positions = [10, 30, 50, 70, 90]
    scan = WindowSimilarity(
        positions=positions, query_positions=positions,
        similarities={"A": [1.0] * 5, "B": [1.0] * 5}, query="q", width=width,
        query_cumulative=np.arange(width + 1),
        rows={"q": query, "A": major, "B": minor},
    )
    short = mk("B", 0, 10, "hmm", qvalue=1e-9)
    short.mean_sim_minor, short.mean_sim_major = 1.0, 1.0
    long = mk("B", 5, 100, "3seq", qvalue=1e-3)
    long.mean_sim_minor, long.mean_sim_major = 0.0, 1.0
    merged, _ = consensus_regions({"hmm": [short], "3seq": [long]}, major="A",
                                  scan=(scan, 20))
    assert (merged[0].msa_start, merged[0].msa_end) == (0, 100)
    assert merged[0].mean_sim_minor == 0.1
    assert merged[0].mean_sim_major == 1.0


def test_merged_window_count_is_in_one_unit() -> None:
    # 3SEQ counts discriminating sites in n_windows, the HMM counts windows; the merged
    # row restates it as base-pair windows overlapping the merged span.
    hmm = mk("B", 10, 40, "hmm", qvalue=1e-9)
    seq = mk("B", 30, 80, "3seq", qvalue=1e-3)
    seq.n_windows = 37  # discriminating sites, not windows
    merged, _ = consensus_regions({"hmm": [hmm], "3seq": [seq]}, major="A",
                                  scan=(_scan_result(), 20))
    # windows [0,20) [20,40) [40,60) [60,80) [80,100): four overlap [10, 80)
    assert merged[0].n_windows == 4


def test_merge_without_a_scan_keeps_the_best_members_values() -> None:
    # Guard: callers that do not pass the scan (tests, the API) get the old behaviour.
    hmm = mk("B", 10, 40, "hmm", qvalue=1e-9)
    seq = mk("B", 30, 80, "3seq", qvalue=1e-3)
    merged, _ = consensus_regions({"hmm": [hmm], "3seq": [seq]}, major="A")
    assert merged[0].mean_sim_minor == hmm.mean_sim_minor
