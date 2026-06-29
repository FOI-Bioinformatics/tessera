"""Bootscan / SimPlot caller (Salminen 1995; SimPlot++ 2022).

The published method of choice for identifying *which parent* a recombinant region came
from, and the only caller here that yields a **bootstrap support** for the call. For each
window it measures the query's identity to every candidate parent, then resamples the
window's alignment columns with replacement many times; the bootstrap support of a parent
is the fraction of resamples in which it is the query's closest match. A run of windows
where a non-major parent's support clears a threshold -- and beats the backbone -- is a
recombinant region, the support travelling onto the region as a confidence the HMM / 3SEQ
calls express only as a p-value. Dependency-free numpy.
"""

from __future__ import annotations

import numpy as np

from .clusters import _window_bounds
from .similarity import WindowSimilarity, _canonical_mask

_BOOTSTRAPS = 100  # column resamples per window
_SUPPORT = 0.70  # bootstrap support a donor must reach to claim a window
_MAX_CANDIDATES = 8  # cap candidate donors (by window wins) to bound cost
_MIN_WINDOWS = 2  # a region must span at least this many windows


def bootstrap_support(
    match: np.ndarray, comparable: np.ndarray, *, n: int = _BOOTSTRAPS, seed: int = 0
) -> np.ndarray:
    """Per-reference bootstrap support over one window's columns.

    ``match`` / ``comparable`` are ``(refs, columns)`` boolean arrays (query matches the
    reference / both carry a canonical base). Returns the fraction of ``n`` column
    resamples in which each reference has the highest identity (matches / comparable).
    """
    n_refs, n_cols = match.shape
    if n_cols == 0 or n_refs == 0:
        return np.zeros(n_refs)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n_cols, size=(n, n_cols))  # (n, columns) with replacement
    m = match[:, idx].sum(axis=2).astype(np.float64)  # (refs, n)
    c = comparable[:, idx].sum(axis=2).astype(np.float64)
    identity = np.divide(m, c, out=np.zeros_like(m), where=c > 0)
    winners = np.argmax(identity, axis=0)  # (n,) best reference per resample
    return np.bincount(winners, minlength=n_refs) / n


def call_regions_bootscan(result: WindowSimilarity, analysis, window_size: int, params):
    """Call recombinant regions by bootscanning. The major parent wins the most windows;
    a run of windows where a non-major candidate's bootstrap support reaches ``_SUPPORT``
    and exceeds the major's is a region for that donor. Returns ``(regions, major, [])``.
    """
    from .analyze import rank_by_wins
    from .regions import Region

    labels = list(result.similarities)
    if len(labels) < 2:
        return [], (labels[0] if labels else None), []
    ranked = rank_by_wins(analysis.winners_with_ties, len(labels)) or labels
    major = ranked[0]
    candidates = [major] + [c for c in ranked if c != major][: _MAX_CANDIDATES]
    major_idx = 0

    query = result.rows[result.query]
    q_canon = _canonical_mask(query)
    match_full = {}
    comp_full = {}
    for ref in candidates:
        r = result.rows[ref]
        comp = q_canon & _canonical_mask(r)
        comp_full[ref] = comp
        match_full[ref] = comp & (query == r)
    match_stack = np.array([match_full[c] for c in candidates])  # (refs, width)
    comp_stack = np.array([comp_full[c] for c in candidates])

    bounds = _window_bounds(result, window_size)  # window column spans (per position)
    width = match_stack.shape[1]
    # winner[w] = the candidate index with top bootstrap support at window w (or major)
    winner = np.full(len(bounds), major_idx)
    support = np.zeros(len(bounds))
    for w, (start, end) in enumerate(bounds):
        s, e = max(0, start), min(width, end)
        if e - s < 1:
            continue
        sup = bootstrap_support(match_stack[:, s:e], comp_stack[:, s:e], seed=w)
        best = int(np.argmax(sup))
        if best != major_idx and sup[best] >= _SUPPORT and sup[best] > sup[major_idx]:
            winner[w] = best
            support[w] = sup[best]

    # Contiguous runs of the same non-major winner -> regions.
    regions: list[Region] = []
    w = 0
    while w < len(winner):
        if winner[w] == major_idx:
            w += 1
            continue
        donor = winner[w]
        j = w
        while j < len(winner) and winner[j] == donor:
            j += 1
        if j - w >= _MIN_WINDOWS:
            msa_start = int(result.positions[w])
            msa_end = int(result.positions[j - 1]) + 1
            q_start = result.column_to_query(msa_start)
            q_end = result.column_to_query(msa_end)
            minor = candidates[donor]
            sim_minor = float(np.nanmean([result.similarities[minor][k] for k in range(w, j)]))
            sim_major = float(np.nanmean([result.similarities[major][k] for k in range(w, j)]))
            regions.append(Region(
                minor_parent=minor, major_parent=major,
                msa_start=msa_start, msa_end=msa_end,
                query_start=q_start, query_end=q_end,
                length_bp=msa_end - msa_start, n_windows=j - w,
                mean_sim_minor=round(sim_minor, 4), mean_sim_major=round(sim_major, 4),
                margin=round(sim_minor - sim_major, 4),
                support=round(float(np.mean(support[w:j])), 3),
                breakpoint_lo=q_start, breakpoint_hi=q_end,
            ))
        w = j
    return regions, major, []
