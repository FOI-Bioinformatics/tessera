"""A region that rests on a locally misaligned backbone is withheld (issue #81).

The scan always crowns the closest reference genome-wide, then asks window by window
whether some other reference beats it. Where the backbone is locally gap-degraded it
loses that comparison to any reference that stayed intact, and the callers read the
switch as a minor parent. On the mpox panel behind the issue the backbone carried
31.2 % gaps over the called locus against 1.64 % genome-wide, and three method
families converged on the same locus naming three different donors, which points at
the locus rather than at any one donor.

The ratio of the backbone's gap rate over a region to its own genome-wide rate
separates the two populations. Measured over the five seeds of the fixture, the mpox
panel, and every region of the eight curated datasets in ``validation/``:

    false, backbone locally degraded   2.27x  3.20x  4.34x  4.61x  4.67x  19.06x
    true, real recombinants            0.00x ... 0.74x

One curated row sits at 2.00x and is a donor-absent coverage row rather than a
caller's region: it states that the true donor is missing, not that a donor beat the
backbone, so such rows are exempt.
"""

from __future__ import annotations

from pathlib import Path

from tessera.recomb.regions import Region, withhold_misaligned_backbone
from tessera.recomb.similarity import compute_similarity

from ..conftest import write_fasta
from ..factories import make_gap_degraded_msa, make_recombination_msa

WINDOW, STEP = 1000, 100
MAX_RATIO = 2.0


def _scan(tmp_path: Path, records: dict[str, str], query: str = "query"):
    path = write_fasta(tmp_path / "msa.fasta", records)
    return compute_similarity(str(path), query, window_size=WINDOW, window_step=STEP)


def _region(major: str, lo: int, hi: int, minor: str = "other", **kw) -> Region:
    return Region(
        minor_parent=minor, major_parent=major, msa_start=lo, msa_end=hi,
        query_start=lo, query_end=hi, n_windows=10,
        mean_sim_minor=0.999, mean_sim_major=0.995, margin=0.004, **kw
    )


def test_a_region_over_a_locally_degraded_backbone_is_withheld(tmp_path: Path) -> None:
    """The comparison that produced it was not between comparably aligned genomes."""
    msa = make_gap_degraded_msa(seed=1)
    result = _scan(tmp_path, msa.records)
    lo, hi = msa.locus
    region = _region(msa.nearest, lo, hi)

    kept, withheld = withhold_misaligned_backbone([region], result, MAX_RATIO)

    assert kept == []
    assert withheld == [region]


def test_a_real_recombinant_region_is_kept(tmp_path: Path) -> None:
    """A donor block on a normally aligned panel must survive untouched."""
    msa = make_recombination_msa(seed=0)
    result = _scan(tmp_path, msa.records)
    assert msa.donor_block is not None
    lo, hi = msa.donor_block
    region = _region(msa.backbone, lo, hi, minor=msa.donor)

    kept, withheld = withhold_misaligned_backbone([region], result, MAX_RATIO)

    assert kept == [region]
    assert withheld == []


def test_a_donor_absent_row_is_exempt(tmp_path: Path) -> None:
    """It states that the true donor is missing, not that a donor beat the backbone.

    One curated dataset carries such a row at twice its backbone's genome-wide gap
    rate; withholding it would drop a coverage statement on the strength of a
    comparison it never made.
    """
    msa = make_gap_degraded_msa(seed=1)
    result = _scan(tmp_path, msa.records)
    lo, hi = msa.locus
    region = _region(msa.nearest, lo, hi, donor_absent=True)

    kept, withheld = withhold_misaligned_backbone([region], result, MAX_RATIO)

    assert kept == [region]
    assert withheld == []


def test_a_uniformly_gappy_backbone_is_not_a_local_defect(tmp_path: Path) -> None:
    """The test is local: the ratio is against the backbone's own genome-wide rate.

    A backbone gappy everywhere -- a draft assembly, say -- is no worse over the
    region than it is anywhere else, so there is nothing local to withhold for.
    """
    msa = make_gap_degraded_msa(seed=1, locus_gap_rate=0.09, background_gap_rate=0.09)
    result = _scan(tmp_path, msa.records)
    lo, hi = msa.locus
    region = _region(msa.nearest, lo, hi)

    kept, withheld = withhold_misaligned_backbone([region], result, MAX_RATIO)

    assert kept == [region]
    assert withheld == []


def test_a_ratio_of_zero_disables_the_test(tmp_path: Path) -> None:
    """An escape hatch, so a user can see what the filter removed."""
    msa = make_gap_degraded_msa(seed=1)
    result = _scan(tmp_path, msa.records)
    lo, hi = msa.locus

    kept, withheld = withhold_misaligned_backbone(
        [_region(msa.nearest, lo, hi)], result, 0.0
    )

    assert len(kept) == 1
    assert withheld == []
