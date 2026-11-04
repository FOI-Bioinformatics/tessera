"""Setting references aside from a scan, and expanding a dropped sibling to its cluster.

Issue #120: the HMM sets the query's whole-genome siblings aside, but 3SEQ, MaxChi and
Bootscan scored the unfiltered scan, where a near-copy of the query wins every window and
leaves them nothing to find. They now receive the scan without those references.
"""

from __future__ import annotations

import numpy as np

from tessera.recomb.siblings import SiblingEvidence, sibling_members, unbeaten_labels
from tessera.recomb.similarity import WindowSimilarity, without_references

NAN = float("nan")


def _scan() -> WindowSimilarity:
    sims = {
        "sib": [1.00, 1.00, 1.00],
        "refA": [0.95, 0.96, 0.97],
        "refB": [0.90, 0.99, NAN],
    }
    return WindowSimilarity(
        positions=[10, 20, 30],
        query_positions=[10, 20, 30],
        similarities=sims,
        query="q",
        width=40,
        query_cumulative=np.arange(41),
        best_sim=[1.0, 1.0, 1.0],
        best_label=["sib", "sib", "sib"],
        informative=[8, 8, 8],
        numerators={k: [int(100 * x) if x == x else 0 for x in v] for k, v in sims.items()},
        denominators={k: [100, 100, 100] for k in sims},
        rows={k: np.frombuffer(b"ACGT", dtype=np.uint8) for k in ("q", *sims)},
    )


def test_references_are_removed_from_every_per_label_mapping() -> None:
    view = without_references(_scan(), {"sib"})

    assert list(view.similarities) == ["refA", "refB"]
    assert list(view.numerators) == ["refA", "refB"]
    assert list(view.denominators) == ["refA", "refB"]


def test_the_best_reference_is_recomputed_over_the_survivors() -> None:
    view = without_references(_scan(), {"sib"})

    assert view.best_label == ["refA", "refB", "refA"]
    assert view.best_sim == [0.95, 0.99, 0.97]


def test_rows_and_the_scan_geometry_are_kept() -> None:
    """The discordant-site tests read ``rows`` by label, so every genome stays in it."""
    scan = _scan()
    view = without_references(scan, {"sib"})

    assert set(view.rows) == {"q", "sib", "refA", "refB"}
    assert view.positions == scan.positions
    assert view.informative == scan.informative
    assert view.query == scan.query


def test_the_input_scan_is_not_modified() -> None:
    scan = _scan()
    without_references(scan, {"sib"})

    assert list(scan.similarities) == ["sib", "refA", "refB"]
    assert scan.best_label == ["sib", "sib", "sib"]


def test_unknown_labels_are_ignored() -> None:
    view = without_references(_scan(), {"nope"})

    assert list(view.similarities) == ["sib", "refA", "refB"]
    assert view.best_label == ["sib", "sib", "sib"]


def test_removing_nothing_returns_an_equivalent_scan() -> None:
    scan = _scan()
    view = without_references(scan, set())

    assert view.similarities == scan.similarities
    assert view.best_label == scan.best_label


def _evidence(label: str) -> SiblingEvidence:
    return SiblingEvidence(label, coverage=1.0, lead_frac=1.0, beaten_frac=0.0, max_beaten_run=0)


def test_a_dropped_sibling_stands_for_its_whole_cluster() -> None:
    """The HMM competes pooled lineages, so the label it drops is a cluster's
    representative; the near-duplicates of a sibling mask the event just as it does."""
    clusters = [["sibA", "sibB", "sibC"], ["refA", "refA2"]]

    assert sibling_members([_evidence("sibA")], clusters) == {"sibA", "sibB", "sibC"}


def test_a_dropped_label_in_no_cluster_stands_for_itself() -> None:
    assert sibling_members([_evidence("sib")], [["refA", "refA2"]]) == {"sib"}
    assert sibling_members([_evidence("sib")], None) == {"sib"}


def test_no_siblings_means_nothing_is_removed() -> None:
    assert sibling_members([], [["a", "b"]]) == set()


def _counts(matches: dict[str, list[int]], sites: int = 1000) -> WindowSimilarity:
    """A scan whose per-window match counts are given; ``sites`` comparable sites each."""
    n_win = len(next(iter(matches.values())))
    return WindowSimilarity(
        positions=list(range(n_win)),
        query_positions=list(range(n_win)),
        similarities={k: [m / sites for m in v] for k, v in matches.items()},
        query="q",
        width=n_win,
        query_cumulative=np.arange(n_win + 1),
        numerators=matches,
        denominators={k: [sites] * n_win for k in matches},
    )


def test_a_true_twin_is_never_beaten() -> None:
    scan = _counts({"twin": [995, 995, 995], "other": [900, 900, 900]})
    assert unbeaten_labels(scan, 0.05) == {"twin"}


def test_a_backbone_beaten_inside_the_tract_is_not_unbeaten() -> None:
    """The backbone leads everywhere except one window, where the donor beats it."""
    scan = _counts({
        "backbone": [990, 990, 800, 990],
        "donor": [900, 900, 990, 900],
    })
    assert unbeaten_labels(scan, 0.05) == set()


def test_ties_do_not_count_as_beating() -> None:
    scan = _counts({"a": [990, 990], "b": [990, 990]})
    assert unbeaten_labels(scan, 0.05) == {"a", "b"}
