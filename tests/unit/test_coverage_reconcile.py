"""reconcile_gaps: a coverage gap caveats an overlapping called region; a
non-overlapping gap stays a donor-absent region."""

from __future__ import annotations

from tessera.recomb.coverage import reconcile_gaps
from tessera.recomb.regions import Region


def _region(msa_start: int, msa_end: int, *, donor_absent: bool = False) -> Region:
    return Region(
        minor_parent="m", major_parent="M", msa_start=msa_start, msa_end=msa_end,
        query_start=msa_start, query_end=msa_end, length_bp=msa_end - msa_start,
        n_windows=1, mean_sim_minor=0.98, mean_sim_major=0.93, margin=0.05,
        donor_absent=donor_absent,
    )


def test_overlapping_gap_caveats_region_and_is_not_emitted():
    called = _region(5000, 9000)
    gap = _region(5200, 5800, donor_absent=True)
    absent = reconcile_gaps([called], [gap])
    assert called.donor_undercovered is True   # the region is caveated
    assert absent == []                         # the gap is not double-emitted


def test_non_overlapping_gap_stays_donor_absent():
    called = _region(1000, 4000)
    gap = _region(5200, 5800, donor_absent=True)
    absent = reconcile_gaps([called], [gap])
    assert called.donor_undercovered is False
    assert absent == [gap]                       # genuine donor-absent region survives


def test_gap_overlapping_several_regions_caveats_all():
    a, b = _region(5000, 5500), _region(5400, 6000)
    absent = reconcile_gaps([a, b], [_region(5300, 5600, donor_absent=True)])
    assert a.donor_undercovered and b.donor_undercovered and absent == []
