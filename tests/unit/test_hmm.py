"""HMM segmentation of the query against the reference panel."""

from __future__ import annotations

import numpy as np

from recomfi.recomb.hmm import segment_query
from recomfi.recomb.similarity import WindowSimilarity


def _two_state_result(n=20, window=50):
    """refA matches the query in the first half, refB in the second -- one breakpoint."""
    width = n * window
    positions = [i * window + window // 2 for i in range(n)]
    half = n // 2

    def counts(hi_first, hi_second):
        return [hi_first if i < half else hi_second for i in range(n)]

    num = {"refA": counts(99, 80), "refB": counts(80, 99)}
    den = {"refA": [100] * n, "refB": [100] * n}
    sim = {s: [num[s][i] / den[s][i] for i in range(n)] for s in num}
    best = [max(sim["refA"][i], sim["refB"][i]) for i in range(n)]
    return WindowSimilarity(
        positions=positions, query_positions=positions, similarities=sim,
        query="q", width=width, query_cumulative=np.arange(width + 1),
        best_sim=best, best_label=["refA"] * n, informative=[100] * n,
        numerators=num, denominators=den,
    )


def test_hmm_finds_one_clean_breakpoint():
    result = _two_state_result()
    segments, major = segment_query(result, jump_rate=1e-3)

    assert [s.state for s in segments] == ["refA", "refB"]
    assert major in {"refA", "refB"}
    # the switch is near the middle of the 1000 bp query
    boundary = segments[1].query_start
    assert 400 <= boundary <= 600
    # the second segment carries a breakpoint uncertainty interval
    assert segments[1].breakpoint_lo is not None
    assert segments[1].breakpoint_lo <= segments[1].breakpoint_hi


def test_hmm_does_not_split_a_single_source():
    # query matches refA everywhere; refB is uniformly worse -> one segment, no flips
    result = _two_state_result()
    for i in range(len(result.positions)):
        result.numerators["refA"][i] = 99
        result.numerators["refB"][i] = 80
        result.similarities["refA"][i] = 0.99
        result.similarities["refB"][i] = 0.80
    segments, major = segment_query(result, jump_rate=1e-3)
    assert len(segments) == 1
    assert major == "refA"
