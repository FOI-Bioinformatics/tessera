"""A clonal query whose nearest reference is locally gap-degraded.

This is the failure mode behind issues #81 and #82, reduced from the mpox run that
prompted them. There the called regions rested on a locally poor alignment rather
than on recombination: at the locus the major parent carried 8.65 % gaps against
1.64 % genome-wide, and the query 5.85 % against 1.59 %, while the minor parent was
unchanged (3.36 % against 3.00 %).

The fixture reproduces that geometry with a query that is clonal by construction --
a tip inside one clade, carrying no donor block from anywhere -- so every region
reported on it is a false positive, and the coverage guard has something it ought to
flag. Issue #81 is the false region; issue #82 is that nothing flags it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from tessera.recomb.analyze import analyze
from tessera.recomb.coverage import (
    LOCAL_DIP_KIND,
    CoverageParams,
    call_coverage_gaps,
)
from tessera.recomb.regions import RegionParams, call_regions
from tessera.recomb.similarity import compute_similarity

from ..conftest import write_fasta
from ..factories import make_gap_degraded_msa

WINDOW, STEP = 1000, 100


def _scan(tmp_path: Path, msa):
    path = write_fasta(tmp_path / "msa.fasta", msa.records)
    result = compute_similarity(str(path), "query", window_size=WINDOW, window_step=STEP)
    analysis = analyze(result)
    regions, major, _ = call_regions(
        result, analysis, WINDOW, RegionParams.with_defaults(WINDOW)
    )
    return result, regions, major


def test_fixture_is_clonal_with_a_locally_degraded_nearest_reference() -> None:
    """The data carries the defect's two preconditions and no recombination."""
    msa = make_gap_degraded_msa(seed=1)
    lo, hi = msa.locus

    # The query is a clade tip: no block of it was copied from another reference.
    assert msa.donor_block is None
    # The nearest reference is several-fold gappier inside the locus than outside,
    # the contrast measured on the real locus.
    assert msa.locus_gap_rate > 4 * msa.background_gap_rate
    assert 0.05 < msa.locus_gap_rate < 0.20
    assert msa.background_gap_rate < 0.05
    # Every other reference is aligned normally across the same columns.
    for label, seq in msa.records.items():
        if label in ("query", msa.nearest):
            continue
        assert seq[lo:hi].count("-") / (hi - lo) < 0.05


def test_nearest_reference_matches_the_query_everywhere_but_the_locus(tmp_path: Path) -> None:
    """The degradation is local: a similarity dip at the locus, no dip outside it."""
    msa = make_gap_degraded_msa(seed=1)
    result, _, major = _scan(tmp_path, msa)
    lo, hi = msa.locus

    assert major == msa.nearest
    sims = result.similarities[msa.nearest]
    # Windows are centred, so only centres a half-window clear of the locus are
    # untouched by it; the ones between overlap it partially.
    half = WINDOW // 2
    inside = [s for pos, s in zip(result.positions, sims, strict=True) if lo <= pos < hi]
    outside = [
        s
        for pos, s in zip(result.positions, sims, strict=True)
        if pos + half < lo or pos - half >= hi
    ]
    assert min(outside) > 0.95
    assert max(inside) < min(outside)


@pytest.mark.xfail(
    strict=True,
    reason="issue #81: a locally gap-degraded major parent is called as a recombinant",
)
@pytest.mark.parametrize("seed", [1, 2, 3, 4, 5])
def test_clonal_query_yields_no_region(tmp_path: Path, seed: int) -> None:
    """A query that is a clade tip has no recombinant region to find.

    Fails today: the locally misaligned major parent loses the locus to an intact
    competitor, which reads as a mosaic. Parametrized because the call is systematic
    rather than one seed's luck, and strict, so fixing #81 turns these red until the
    marker comes off.
    """
    msa = make_gap_degraded_msa(seed=seed)
    _, regions, _ = _scan(tmp_path, msa)

    assert regions == []


def test_the_degraded_locus_is_flagged(tmp_path: Path) -> None:
    """Something marks the locus rather than letting it pass unremarked.

    Issue #82, first half: the best-similarity threshold was a fixed drop below the
    query's typical best match. On a panel this close that landed below every value
    the panel shows, so no window could reach it, the gaps table was empty, and the
    false region of issue #81 reached the report with nothing attached. The
    threshold now falls back to the panel's own scale when the fixed one cannot
    resolve anything, and the locus is reported.

    It is reported as a local dip, not as a possible missing reference: the closest
    reference over the locus still matches the query closely, since an intact
    competitor covers the columns the nearest reference lost. The second half of
    #82 -- a flag that describes the major parent rather than the donor -- is not
    addressed here, so ``donor_undercovered`` still speaks only of the donor.
    """
    msa = make_gap_degraded_msa(seed=1)
    result, _, _ = _scan(tmp_path, msa)
    lo, hi = msa.locus
    gaps, _ = call_coverage_gaps(result, WINDOW, CoverageParams.with_defaults(WINDOW))

    overlapping = [g for g in gaps if g.msa_start < hi and g.msa_end > lo]
    assert overlapping
    assert {g.kind for g in overlapping} == {LOCAL_DIP_KIND}
