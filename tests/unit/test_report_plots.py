"""The plots label what the regions table says, in the same colours."""

from __future__ import annotations

from pathlib import Path

from matplotlib.figure import Figure

from tessera.recomb.regions import Region
from tessera.recomb.report import pair_datasets
from tessera.recomb.report_plots import (
    GREY,
    _color_map,
    _palette,
    _shade_regions,
    build_interactive_figure,
    region_labels,
)
from tessera.recomb.similarity import compute_similarity

from ..conftest import recombinant_msa


def _region(minor: str, major: str = "A", *, start: int = 2000, end: int = 4000,
            donor_absent: bool = False) -> Region:
    return Region(
        minor_parent=minor, major_parent=major, msa_start=start, msa_end=end,
        query_start=start, query_end=end, n_windows=10,
        mean_sim_minor=0.99, mean_sim_major=0.95, margin=0.04, donor_absent=donor_absent,
    )


# --- colours cover every label a region names ------------------------------

def test_region_labels_appends_region_parents_to_the_plotted_datasets() -> None:
    labels = region_labels(["A"], [_region("B"), _region("other", donor_absent=True)])
    assert labels == ["A", "B"]  # the donor is added; a donor-absent stand-in is not
    assert region_labels(["A", "B"], [_region("B")]) == ["A", "B"]  # no duplicates


def test_donor_outside_the_top_n_keeps_a_colour(tmp_path: Path) -> None:
    """With --top-n 1 only the backbone was in the colour map, so the donor band and
    its label were drawn grey."""
    result = compute_similarity(str(recombinant_msa(tmp_path, recombinant=True)), "query")
    fig = build_interactive_figure(result, ["A"], [_region("B")])
    donor_colour = _color_map(["A", "B"])["B"]
    band = next(a for a in fig.layout.annotations if a.text == "B")
    assert band.font.color == donor_colour
    assert donor_colour != GREY


# --- a donor-absent region is not a recombinant from the backbone ----------

def test_static_plot_does_not_call_a_donor_absent_region_recombinant() -> None:
    """The stand-in label of a donor-absent region is the closest reference, often the
    backbone itself: the legend read "recombinant: <backbone>" in the backbone's colour."""
    ax = Figure().subplots()  # no pyplot state, no display backend needed
    _shade_regions(ax, [_region("A", donor_absent=True)], _palette(["A", "B"]))
    labels = ax.get_legend_handles_labels()[1]
    assert labels == ["donor absent (no close reference)"]


def test_static_plot_still_labels_a_called_donor() -> None:
    ax = Figure().subplots()
    _shade_regions(ax, [_region("B")], _palette(["A", "B"]))
    labels = ax.get_legend_handles_labels()[1]
    assert labels == ["recombinant: B"]


def test_interactive_plot_marks_a_donor_absent_region(tmp_path: Path) -> None:
    result = compute_similarity(str(recombinant_msa(tmp_path, recombinant=True)), "query")
    fig = build_interactive_figure(result, ["A", "B"], [_region("A", donor_absent=True)])
    texts = [a.text for a in fig.layout.annotations]
    assert texts == ["donor absent"]
    assert fig.layout.annotations[0].font.color == GREY


# --- the pair plot is major versus leading minor ---------------------------

def test_pair_plot_shows_the_backbone_and_the_leading_donor() -> None:
    """It showed the top two window winners, which are the backbone and a
    near-duplicate of it when the panel holds one."""
    regions = [_region("variola", "cowpox", start=2000, end=2600),
               _region("camelpox", "cowpox", start=4000, end=4100)]
    assert pair_datasets(regions, ["cowpox", "cowpox2"]) == ["cowpox", "variola"]


def test_pair_plot_falls_back_to_the_window_ranking_without_a_called_donor() -> None:
    assert pair_datasets([], ["cowpox", "cowpox2"]) == ["cowpox", "cowpox2"]
    absent = [_region("cowpox", "cowpox", donor_absent=True)]
    assert pair_datasets(absent, ["cowpox", "cowpox2"]) == ["cowpox", "cowpox2"]
