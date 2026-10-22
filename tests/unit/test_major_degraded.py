"""The report says when the backbone itself is locally degraded (issue #82).

``donor_undercovered`` reports on the minor parent. A major parent that is locally
gap-degraded was never flagged, so on the mpox run behind the issue all three regions
read ``donor_undercovered = no`` while the backbone was the degraded sequence, and a
reader had nothing in the output saying the call rested on a poor backbone.

Why the gap rate and not the backbone's identity. Over a true donor region the backbone
is *by definition* the poorer match -- that is what makes it a donor region -- so
flagging a low ``mean_sim_major`` fires on every genuine recombinant and tells the
reader nothing. Issue #81 established the discriminator that does separate the two: the
backbone's gap rate over the region against its own genome-wide rate. A backbone
degraded enough to lose the comparison measured 2.27x to 19.06x its genome-wide rate at
5.02 % to 31.21 % locally, while the backbone under a real recombinant region stayed at
or below 0.74x and 3.09 %.

#81 withholds a region above 2.0x and 5 %. This flag covers the band below that bar and
above anything a real region shows, so a mildly degraded backbone is caveated rather
than silently reported or silently removed.
"""

from __future__ import annotations

import numpy as np

from tessera.recomb.regions import Region, flag_degraded_backbone
from tessera.recomb.similarity import GAP


class _Result:
    """The one thing the flag reads: the backbone's alignment row."""

    def __init__(self, rows: dict[str, np.ndarray]) -> None:
        self.rows = rows


def _row(width: int, gap_spans: list[tuple[int, int]]) -> np.ndarray:
    row = np.full(width, ord("A"), dtype=np.uint8)
    for lo, hi in gap_spans:
        row[lo:hi] = GAP
    return row


def _region(lo: int = 1000, hi: int = 2000) -> Region:
    return Region(
        minor_parent="donor", major_parent="backbone",
        msa_start=lo, msa_end=hi, query_start=lo, query_end=hi,
        n_windows=5, mean_sim_minor=0.99, mean_sim_major=0.98, margin=0.01,
        methods=("hmm",),
    )


def test_a_mildly_degraded_backbone_is_flagged() -> None:
    """Locally gappier than it is elsewhere, but under #81's withholding bar."""
    region = _region()
    # 4 % over the region, 2 % genome-wide: ratio 2.0 is the withholding bar, so keep
    # the local rate under 5 % and the ratio at 2.0 -- #81 leaves this one reported.
    row = _row(10000, [(1000, 1040), (5000, 5160)])
    flag_degraded_backbone([region], _Result({"backbone": row}))

    assert region.major_degraded is True


def test_a_backbone_aligned_normally_over_the_region_is_not_flagged() -> None:
    """A true donor region: the backbone keeps its columns, it is just not the source."""
    region = _region()
    row = _row(10000, [(5000, 5100)])  # 1 % genome-wide, none over the region
    flag_degraded_backbone([region], _Result({"backbone": row}))

    assert region.major_degraded is False


def test_a_uniformly_gappy_backbone_is_not_flagged() -> None:
    """A draft assembly gappy throughout is not degraded *here*.

    The test is relative to the backbone itself, as in issue #81, so being gappy
    everywhere is not evidence about this region.
    """
    region = _region()
    row = _row(10000, [(i, i + 30) for i in range(0, 10000, 100)])  # 30 % everywhere
    flag_degraded_backbone([region], _Result({"backbone": row}))

    assert region.major_degraded is False


def test_a_donor_absent_row_is_exempt() -> None:
    """It says the query is far from every reference, not that the backbone failed.

    The same exemption issue #81's withholding makes, for the same reason.
    """
    region = _region()
    region.donor_absent = True
    row = _row(10000, [(1000, 1040), (5000, 5160)])
    flag_degraded_backbone([region], _Result({"backbone": row}))

    assert region.major_degraded is False


def test_a_backbone_absent_from_the_alignment_is_not_flagged() -> None:
    """No row means nothing was measured, which is not evidence of degradation."""
    region = _region()
    flag_degraded_backbone([region], _Result({}))

    assert region.major_degraded is False


def test_the_flag_is_independent_of_the_donor_caveat() -> None:
    """Separate columns: which parent is at fault changes the reader's conclusion."""
    region = _region()
    region.donor_undercovered = True
    row = _row(10000, [(5000, 5100)])
    flag_degraded_backbone([region], _Result({"backbone": row}))

    assert region.donor_undercovered is True
    assert region.major_degraded is False
