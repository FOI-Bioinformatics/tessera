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

import csv
import json
import logging
from pathlib import Path

import pytest

from tessera.recomb.analyze import analyze
from tessera.recomb.coverage import (
    LOCAL_DIP_KIND,
    CoverageParams,
    call_coverage_gaps,
)
from tessera.recomb.regions import RegionParams, call_regions
from tessera.recomb.run import RecombParams, run_recomb
from tessera.recomb.similarity import compute_similarity

from ..conftest import write_fasta
from ..factories import make_gap_degraded_msa

WINDOW, STEP = 1000, 100


def _read_tsv(path: Path) -> list[dict[str, str]]:
    lines = [x for x in path.read_text().splitlines() if x and not x.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


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


def _present_regions(out: Path) -> list[dict[str, str]]:
    return [r for r in _read_tsv(out / "recombination_regions.tsv")
            if r["donor_absent"] == "no"]


def _run(tmp_path: Path, msa, name: str = "msa.fasta") -> Path:
    path = write_fasta(tmp_path / name, msa.records)
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=path, output=out, query="query", window_size=WINDOW,
                     window_step=STEP, plot_format="png"),
        logging.getLogger("tessera.test"),
    )
    return out


@pytest.mark.parametrize("seed", [1, 2, 3, 4, 5])
def test_no_region_is_reported_over_the_degraded_locus(tmp_path: Path, seed: int) -> None:
    """The locus itself is no longer called, which is the symptom issue #81 reports.

    End to end, because what the reader sees is the regions table and the scan sets
    such a region aside on the way there rather than inside a caller. Parametrized
    because the false call was systematic rather than one seed's luck: on the main
    branch every seed reports a region overlapping the locus.
    """
    msa = make_gap_degraded_msa(seed=seed)
    lo, hi = msa.locus
    out = _run(tmp_path, msa)

    overlapping = [
        r for r in _present_regions(out)
        if int(r["msa_start"]) < hi and int(r["msa_end"]) > lo
    ]
    assert overlapping == []


@pytest.mark.xfail(
    strict=True,
    reason="issue #86: the degraded locus also corrupts the genome-wide statistics, "
           "so the flanks are reported as a mosaic of the wrong backbone",
)
def test_a_clonal_query_yields_no_region_at_all(tmp_path: Path) -> None:
    """The locus is not the only thing its degradation breaks.

    On seed 4 the scan reports two regions that do not overlap the locus at all: it
    names ``outclade0``, of the other clade, as the backbone and the query's own
    sister as a donor across the flanks, which is the mosaic read the wrong way
    round. Ranking the references by their mean similarity lets the degraded locus
    drag the sister's genome-wide score below a competitor's, and the median picks
    the sister on all five seeds where the mean does not. That is a second mechanism
    with the same cause, it predates this fix (the same two regions are on main),
    and it needs the genome-wide statistics to be robust rather than one locus to be
    withheld, so it is filed separately as #86.
    """
    msa = make_gap_degraded_msa(seed=4)
    out = _run(tmp_path, msa)

    assert _present_regions(out) == []


def test_the_withheld_region_is_reported_not_silently_dropped(tmp_path: Path) -> None:
    """A region removed from the table is still accounted for.

    Reported numbers have to be faithful, so the run records how many regions it
    set aside and why, rather than leaving the reader with a table that is shorter
    for no stated reason. The assertions name a provenance key and a phrase of the
    report: the input path carries this test's own name, so matching a bare word
    anywhere in the output would pass on the path alone.
    """
    msa = make_gap_degraded_msa(seed=1)
    path = write_fasta(tmp_path / "msa.fasta", msa.records)
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=path, output=out, query="query", window_size=WINDOW,
                     window_step=STEP, plot_format="png"),
        logging.getLogger("tessera.test"),
    )

    run = json.loads((out / "run_provenance.json").read_text())["run"]
    key = "regions withheld (misaligned backbone)"
    assert key in run
    assert int(str(run[key]).split()[0]) >= 1
    assert "locally gap-degraded" in (out / "report.html").read_text()


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
