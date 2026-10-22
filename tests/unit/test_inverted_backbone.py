"""A region whose major parent is a backbone the run has rejected (issue #86).

``reconcile_major`` already guards one artefact: the HMM segments on a view that can
disagree with the genome-wide evidence, and when the major it returns wins no base-pair
window at all, the windowed vote is trusted instead. That guard corrects the run's
backbone but not the regions the caller had already built on the rejected one, so the
report could name one backbone in its headline and the opposite in its table -- with the
genome-wide closest reference cast as a donor into a reference that is closest nowhere.

Measured on seed 4 of ``make_gap_degraded_msa``, the only seed that crosses the
informative-site auto-trigger: the caller's major is ``outclade0`` with zero base-pair
window wins, the reconciled backbone is ``sister`` with 129, and both rows read
``minor=sister major=outclade0``.
"""

from __future__ import annotations

from tessera.recomb.ensemble import (
    is_inverted_backbone,
    withhold_inverted_backbone,
)
from tessera.recomb.regions import Region


def _region(minor: str, major: str, lo: int = 100, hi: int = 900) -> Region:
    return Region(
        minor_parent=minor, major_parent=major, msa_start=lo, msa_end=hi,
        query_start=lo, query_end=hi, n_windows=5,
        mean_sim_minor=0.99, mean_sim_major=0.98, margin=0.01,
        methods=("hmm",),
    )


def test_a_row_naming_the_backbone_as_its_donor_is_withheld() -> None:
    """The inversion: the reconciled backbone appears as the region's minor parent."""
    region = _region("sister", "outclade0")
    wins = {"sister": 129, "cousin": 67}
    assert is_inverted_backbone(region, backbone="sister", window_wins=wins)


def test_a_row_is_kept_when_its_major_wins_windows_of_its_own() -> None:
    """Callers may legitimately disagree about the backbone.

    The guard is tied to the artefact that prompted it -- a major that is the closest
    reference in not one window -- so a disagreement backed by real window evidence is
    left alone rather than being overruled by the headline.
    """
    region = _region("sister", "cousin")
    wins = {"sister": 129, "cousin": 67}
    assert not is_inverted_backbone(region, backbone="sister", window_wins=wins)


def test_an_ordinary_donor_region_is_untouched() -> None:
    """The usual case: the row's major is the backbone and its minor is a real donor."""
    region = _region("cousin", "sister")
    wins = {"sister": 129, "cousin": 67}
    assert not is_inverted_backbone(region, backbone="sister", window_wins=wins)


def test_nothing_is_withheld_without_a_reconciled_backbone() -> None:
    """With no backbone there is nothing for a row to contradict."""
    region = _region("sister", "outclade0")
    assert not is_inverted_backbone(region, backbone=None, window_wins={"sister": 129})


def test_nothing_is_withheld_without_window_evidence() -> None:
    """An empty vote cannot establish that a major won no window.

    Without this the guard would withhold every row on a panel whose windows are all
    uninformative, which is the opposite of what the evidence supports.
    """
    region = _region("sister", "outclade0")
    assert not is_inverted_backbone(region, backbone="sister", window_wins={})


def test_the_breakdown_rows_stay_aligned_with_the_regions_they_describe() -> None:
    """The parallel list is filtered with the regions, not left to slip out of step.

    A breakdown row describes a particular region; dropping a region without its row
    would shift every row after it onto the wrong region, which is worse than the
    inverted call this guard exists to remove.
    """
    inverted = _region("sister", "outclade0", 100, 900)
    ordinary = _region("cousin", "sister", 2000, 3000)
    rows = [{"tag": "inverted"}, {"tag": "ordinary"}]
    kept, kept_rows, withheld = withhold_inverted_backbone(
        [inverted, ordinary], rows,
        backbone="sister", window_wins={"sister": 129, "cousin": 67},
    )
    assert kept == [ordinary]
    assert kept_rows == [{"tag": "ordinary"}]
    assert withheld == [inverted]


def test_the_report_note_names_the_rejected_backbone_and_the_real_one() -> None:
    """The reader is told the rows existed and why, not left with a shorter table."""
    from tessera.recomb.report_html import _inverted_html

    note = _inverted_html([_region("sister", "outclade0")], "sister")

    assert "outclade0" in note
    assert "sister" in note
    assert "closest in no window" in note


def test_there_is_no_note_when_nothing_was_withheld() -> None:
    from tessera.recomb.report_html import _inverted_html

    assert _inverted_html([], "sister") == ""
    assert _inverted_html([_region("sister", "outclade0")], None) == ""
