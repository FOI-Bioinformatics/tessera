"""Unit tests for the reassortment call (per-segment inputs; no binaries)."""

from __future__ import annotations

from tessera.reassort.assign import _call_reassortment


def test_clonal_when_a_strain_is_near_across_all_segments():
    # HA and NA both have strain S1 in their top-k -> one parent explains both -> clonal.
    topk = {"HA": {"S1", "S2"}, "NA": {"S1", "S3"}}
    assert _call_reassortment(topk, assigned=["HA", "NA"]) == "clonal"


def test_reassortant_when_topk_strain_sets_are_disjoint():
    topk = {"HA": {"S1", "S2"}, "NA": {"S3", "S4"}}
    assert _call_reassortment(topk, assigned=["HA", "NA"]) == "reassortant"


def test_undetermined_with_fewer_than_two_assigned_segments():
    assert _call_reassortment({"HA": {"S1"}}, assigned=["HA"]) == "undetermined"


def test_unassigned_segments_excluded_from_the_decision():
    # NA unassigned -> decision on HA + PB2; they share S1 -> clonal.
    topk = {"HA": {"S1"}, "PB2": {"S1"}}
    assert _call_reassortment(topk, assigned=["HA", "PB2"]) == "clonal"
