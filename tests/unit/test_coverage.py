"""Reference-coverage diagnostics: flagging possible missing references."""

from __future__ import annotations

import numpy as np

from tessera.recomb.coverage import (
    CoverageGap,
    CoverageParams,
    call_coverage_gaps,
    flag_subclade_identifiability,
    flag_undercovered_regions,
    gaps_as_regions,
)
from tessera.recomb.regions import Region
from tessera.recomb.similarity import WindowSimilarity


def _rows(**seqs):
    """Alignment rows as the callers see them: label -> uint8 array of bases."""
    return {label: np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
            for label, seq in seqs.items()}


def _result(best_sim, best_label, informative, *, window=100):
    """A minimal WindowSimilarity with non-overlapping windows over an ungapped query."""
    n = len(best_sim)
    width = n * window
    positions = [i * window + window // 2 for i in range(n)]
    return WindowSimilarity(
        positions=positions,
        query_positions=positions,
        similarities={},
        query="q",
        width=width,
        query_cumulative=np.arange(width + 1),  # ungapped: column == query base
        best_sim=best_sim,
        best_label=best_label,
        informative=informative,
    )


def test_divergent_gap_is_flagged():
    # high coverage everywhere except a contiguous block where the best ref is poor
    best = [0.98] * 10
    labels = ["refA"] * 10
    info = [100] * 10
    for i in (4, 5, 6):
        best[i], labels[i] = 0.80, "refB"
    result = _result(best, labels, info)

    gaps, threshold = call_coverage_gaps(result, 100, CoverageParams.with_defaults(100))

    assert 0.90 < threshold < 0.95  # adaptive: ~p90(best) - 0.05
    assert len(gaps) == 1
    gap = gaps[0]
    assert (gap.query_start, gap.query_end) == (400, 700)
    assert gap.best_label == "refB"
    assert gap.kind == "divergent"
    assert gap.mean_best == 0.80


def test_gap_length_bp_is_query_bases_not_alignment_columns():
    """`length_bp` was the MSA column span, gaps included, under a bp header."""
    best = [0.98] * 10
    labels = ["refA"] * 10
    info = [100] * 10
    for i in (4, 5, 6):
        best[i], labels[i] = 0.80, "refB"
    result = _result(best, labels, info)
    # 100 gap columns inside the gap's span: columns 500-599 carry no query base.
    present = np.ones(result.width, dtype=int)
    present[500:600] = 0
    result.query_cumulative = np.concatenate([[0], np.cumsum(present)])

    gaps, _ = call_coverage_gaps(result, 100, CoverageParams.with_defaults(100))

    gap = gaps[0]
    assert (gap.msa_start, gap.msa_end) == (400, 700)
    assert (gap.query_start, gap.query_end) == (400, 600)
    assert gap.length_bp == 200 == gap.query_end - gap.query_start


def test_low_information_gap_is_labelled_distinctly():
    best = [0.98] * 10
    labels = ["refA"] * 10
    info = [100] * 10
    for i in (4, 5, 6):
        best[i], info[i] = 0.80, 3  # poor best, but almost no comparable bases
    result = _result(best, labels, info)

    gaps, _ = call_coverage_gaps(result, 100, CoverageParams.with_defaults(100))

    assert len(gaps) == 1
    assert gaps[0].kind == "low_information"


def test_full_coverage_yields_no_gaps():
    result = _result([0.99] * 8, ["refA"] * 8, [100] * 8)
    gaps, _ = call_coverage_gaps(result, 100, CoverageParams.with_defaults(100))
    assert gaps == []


def test_absolute_floor_overrides_adaptive():
    best = [0.90] * 6  # uniformly moderate: adaptive would flag nothing
    result = _result(best, ["refA"] * 6, [100] * 6)

    none_flagged, thr_adaptive = call_coverage_gaps(
        result, 100, CoverageParams.with_defaults(100)
    )
    assert none_flagged == []

    flagged, thr_floor = call_coverage_gaps(
        result, 100, CoverageParams.with_defaults(100, floor=0.95)
    )
    assert thr_floor == 0.95
    assert len(flagged) == 1  # all six windows now under the 0.95 floor -> one gap


def test_gaps_as_regions_bridges_divergent_only():
    # a divergent gap becomes a donor-absent region; a low-information gap does not
    n, window = 10, 100
    width = n * window
    positions = [i * window + window // 2 for i in range(n)]
    result = WindowSimilarity(
        positions=positions, query_positions=positions,
        similarities={"major": [0.95] * n}, query="q", width=width,
        query_cumulative=np.arange(width + 1),
    )
    divergent = CoverageGap(400, 700, 400, 700, 300, 3, "refX", 0.82, "divergent")
    low_info = CoverageGap(800, 900, 800, 900, 100, 1, "refY", 0.80, "low_information")

    regions = gaps_as_regions([divergent, low_info], result, "major")
    assert len(regions) == 1
    r = regions[0]
    assert r.donor_absent is True
    assert r.donor_undercovered is True
    assert r.minor_parent == "refX"  # the closest (still poor) present reference
    assert r.major_parent == "major"
    assert (r.query_start, r.query_end) == (400, 700)


def test_flag_undercovered_regions():
    regions = [
        Region("variola", "cowpox", 0, 100, 0, 100, 5, 0.998, 0.97, 0.028),
        Region("env", "A1", 200, 300, 200, 300, 5, 0.80, 0.79, 0.01),
    ]
    flag_undercovered_regions(regions, threshold=0.90)
    assert regions[0].donor_undercovered is False  # donor matches well
    assert regions[1].donor_undercovered is True  # donor itself is a poor match


# The query over the region: 12 informative sites. The donor D clearly wins; the near
# lineage N shares most of them, the far lineage F shares almost none.
_IDENT_QUERY = "ACGTACGTACGT"
_IDENT_ROWS = _rows(
    q=_IDENT_QUERY,
    backbone="TTTTTTTTTTTT",
    donor=_IDENT_QUERY,                 # matches the query everywhere
    near="ACGTACGTACGA",                # 11/12 -- the query does not clearly prefer donor
    far="TGCATGCATGCA",                 # 0/12 -- the query clearly prefers donor over it
)
_LINEAGE = {"backbone": "A", "donor": "D", "near": "D.near", "far": "F"}


def _ident_region(minor="donor"):
    return Region(minor, "backbone", 0, 12, 0, 12, 1, 0.99, 0.10, 0.89, methods=("hmm",))


def test_subclade_identifiable_when_donor_clearly_beats_other_lineages():
    # F is the only competing lineage (near nests under the donor's own lineage D), and
    # the query prefers the donor over F at every site, so the sub-clade is resolved.
    lineage = {"backbone": "A", "donor": "D", "near": "D.near", "far": "F"}
    regions = [_ident_region()]
    flag_subclade_identifiability(regions, _IDENT_ROWS, "q", lineage.get, alpha=0.05)
    assert regions[0].subclade_identifiable is True
    assert regions[0].tied_lineages == ()


def test_subclade_unresolved_when_another_lineage_ties_the_donor():
    # Make F a near-equal match too, so no present lineage is clearly beaten: the panel
    # cannot resolve which sub-clade the tract came from (issue #90).
    rows = _rows(
        q=_IDENT_QUERY, backbone="TTTTTTTTTTTT",
        donor="ACGTACGTACGA",   # 11/12
        far="ACGTACGTACTA",     # 10/12 -- within a site of the donor
    )
    lineage = {"backbone": "A", "donor": "D", "far": "F"}
    regions = [_ident_region()]
    flag_subclade_identifiability(regions, rows, "q", lineage.get, alpha=0.05)
    assert regions[0].subclade_identifiable is False
    assert regions[0].tied_lineages == ("D", "F")


def test_subclade_identifiability_skips_donor_absent_and_unlabelled():
    absent = _ident_region()
    absent.donor_absent = True
    unlabelled = _ident_region()
    flag_subclade_identifiability([absent], _IDENT_ROWS, "q", _LINEAGE.get, alpha=0.05)
    flag_subclade_identifiability([unlabelled], _IDENT_ROWS, "q", lambda _l: None, alpha=0.05)
    assert absent.subclade_identifiable is None  # already flagged absent
    assert unlabelled.subclade_identifiable is None  # no lineage labels


def test_subclade_identifiability_not_evaluated_without_a_competing_lineage():
    # Only the donor's own lineage (and the backbone) are present: nothing to be
    # ambiguous with, so the region is left unevaluated rather than called resolved.
    lineage = {"backbone": "A", "donor": "D", "near": "D.near"}
    regions = [_ident_region()]
    flag_subclade_identifiability(regions, _IDENT_ROWS, "q", lineage.get, alpha=0.05)
    assert regions[0].subclade_identifiable is None
