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
from .stats import benjamini_hochberg

_BOOTSTRAPS = 100  # column resamples per window
_SUPPORT = 0.70  # bootstrap support a donor must reach to claim a window
_MAX_CANDIDATES = 8  # cap candidate donors (by window wins) to bound cost
_MIN_WINDOWS = 2  # a region must span at least this many windows
_PERMUTATIONS = 2_000  # block permutations for the run-length null
_MAX_PERM_CELLS = 20_000_000  # cap permutations * windows
_MIN_PERMUTATIONS = 500


def _overlap_block(bounds: list[tuple[int, int]]) -> int:
    """Number of consecutive windows spanned by one window width.

    Sliding windows overlap (``window_step`` < ``window_size``), so neighbouring
    windows share most of their columns and their winners are strongly
    autocorrelated *even with no recombination*. That correlation length is the
    block size the null below must preserve; medians keep it right under
    informative-site windowing, where windows are not uniformly spaced.
    """
    if len(bounds) < 2:
        return 1
    widths = np.array([end - start for start, end in bounds], dtype=float)
    starts = np.array([start for start, _ in bounds], dtype=float)
    spacing = float(np.median(np.diff(starts)))
    if spacing <= 0:
        return 1
    return max(1, min(len(bounds), int(round(float(np.median(widths)) / spacing))))


def run_pvalue(wins: np.ndarray, run_len: int, *, block: int = 1, seed: int = 0) -> float:
    """Permutation ``P(longest run of donor wins >= run_len)``, win count held fixed.

    The null asks whether this donor's winning windows are *contiguous* beyond what
    their number alone explains -- a recombinant tract is one run, chance flips are
    scattered. Windows are permuted in blocks of ``block`` so the autocorrelation
    that overlapping windows induce survives into the null; a plain element-wise
    shuffle would destroy it and make every short run look significant. The block
    grid's phase is drawn at random per permutation, so the p-value depends on the
    run and not on where the run sits relative to a fixed grid.
    """
    n_win = int(wins.size)
    k = int(np.count_nonzero(wins))
    if run_len <= 0 or k == 0 or run_len > k:
        return 1.0
    block = max(1, min(block, n_win))
    perms = _PERMUTATIONS
    if n_win * perms > _MAX_PERM_CELLS:
        perms = max(_MIN_PERMUTATIONS, _MAX_PERM_CELLS // max(1, n_win))

    shuffled = _block_permutations(wins, block, perms, np.random.default_rng(seed))

    longest = np.zeros(perms, dtype=np.int32)
    cur = np.zeros(perms, dtype=np.int32)
    for j in range(n_win):  # vectorised over permutations; reset-cumulative run length
        cur = np.where(shuffled[:, j], cur + 1, 0)
        longest = np.maximum(longest, cur)
    return float((np.count_nonzero(longest >= run_len) + 1) / (perms + 1))


def _block_permutations(wins: np.ndarray, block: int, perms: int, rng) -> np.ndarray:
    """``(perms, n_win)`` block permutations of ``wins``, every row holding its wins.

    The block grid is laid down at a random phase for each permutation (offset drawn
    from ``0 .. block-1``), the blocks are shuffled, and the partial blocks the phase
    leaves at either end travel as blocks of their own. With a fixed grid the null
    depended on where a run happened to fall: a run inside one block was intact in
    every permutation and scored p = 1 whatever its length, while the same run
    shifted by a few windows was cut in two and scored p near 1/n_blocks. Randomising
    the phase makes the test a function of the run, not of the grid.

    The vector is padded to whole blocks at both ends and the padding cells are
    removed wherever they land, so every row keeps exactly ``n_win`` real cells and
    all of its wins (cutting the shuffled row back to ``n_win`` instead dropped the
    wins a shuffle moved past the end).
    """
    n_win = int(wins.size)
    block = max(1, min(block, n_win))
    offsets = rng.integers(0, block, size=perms) if block > 1 else np.zeros(perms, dtype=int)
    out = np.empty((perms, n_win), dtype=bool)
    for offset in np.unique(offsets):
        rows = np.flatnonzero(offsets == offset)
        front = int((block - offset) % block)  # cells before the first grid line
        n_blocks = int(np.ceil((n_win + front) / block))
        padded = np.zeros(n_blocks * block, dtype=bool)
        padded[front:front + n_win] = wins
        real = np.zeros(n_blocks * block, dtype=bool)
        real[front:front + n_win] = True
        order = np.argsort(rng.random((len(rows), n_blocks)), axis=1)
        values = padded.reshape(n_blocks, block)[order].reshape(len(rows), -1)
        keep = real.reshape(n_blocks, block)[order].reshape(len(rows), -1)
        # Every row has exactly n_win real cells, so the boolean selection reshapes cleanly.
        out[rows] = values[keep].reshape(len(rows), n_win)
    return out


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


def _candidate_groups(candidates: list[str], clusters: list[list[str]] | None) -> np.ndarray:
    """Group index per candidate: candidates in one lineage cluster share an index."""
    of = {label: i for i, members in enumerate(clusters or []) for label in members}
    keys: dict[object, int] = {}
    out = np.empty(len(candidates), dtype=np.int64)
    for k, label in enumerate(candidates):
        key: object = ("cluster", of[label]) if label in of else ("genome", label)
        out[k] = keys.setdefault(key, len(keys))
    return out


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
    # Near-duplicates of one lineage split the bootstrap support between them: two
    # genomes 0.4 % apart each win about half the resamples, and neither reaches the
    # 70 % threshold the lineage would. Support is therefore counted per lineage cluster
    # (the run's shared clusters; a genome outside every cluster is its own), and the
    # cluster's best-supported member names the donor.
    group = _candidate_groups(candidates, getattr(params, "reference_clusters", None))
    n_groups = int(group.max()) + 1
    major_group = int(group[major_idx])
    # winner[w] = the candidate index with top bootstrap support at window w (or major)
    winner = np.full(len(bounds), major_idx)
    support = np.zeros(len(bounds))
    for w, (start, end) in enumerate(bounds):
        s, e = max(0, start), min(width, end)
        if e - s < 1:
            continue
        sup = bootstrap_support(match_stack[:, s:e], comp_stack[:, s:e], seed=w)
        by_group = np.bincount(group, weights=sup, minlength=n_groups)
        best_group = int(np.argmax(by_group))
        if (best_group != major_group and by_group[best_group] >= _SUPPORT
                and by_group[best_group] > by_group[major_group]):
            members = np.flatnonzero(group == best_group)
            winner[w] = int(members[np.argmax(sup[members])])
            support[w] = by_group[best_group]

    # Contiguous runs of the same non-major winner -> candidate regions, each scored
    # against the block-permutation null so a run of chance flips cannot pass as an
    # event. Bootstrap support says which reference is closest here; it does not say
    # that the query recombined, which is what the run length is tested for.
    block = _overlap_block(bounds)
    # Runs are of one lineage, not one genome: the member with the top support can change
    # from window to window within a cluster.
    won = np.where(winner == major_idx, -1, group[winner])
    runs: list[tuple[int, int, int, float]] = []  # (donor, start_w, end_w, p)
    w = 0
    while w < len(won):
        if won[w] < 0:
            w += 1
            continue
        lineage = won[w]
        j = w
        while j < len(won) and won[j] == lineage:
            j += 1
        if j - w >= _MIN_WINDOWS:
            # The run's donor is the member that won most of its windows.
            donor = int(np.bincount(winner[w:j]).argmax())
            p = run_pvalue(won == lineage, j - w, block=block, seed=donor)
            runs.append((donor, w, j, p))
        w = j
    if not runs:
        return [], major, []

    qvalues = benjamini_hochberg([r[3] for r in runs])
    regions: list[Region] = []
    for (donor, w, j, p), q in zip(runs, qvalues, strict=True):
        if q > params.alpha:  # FDR-controlled gate (BH q, not raw p)
            continue
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
            n_windows=j - w,
            mean_sim_minor=round(sim_minor, 4), mean_sim_major=round(sim_major, 4),
            margin=round(sim_minor - sim_major, 4),
            support=round(float(np.mean(support[w:j])), 3),
            pvalue=p, qvalue=q,
            test="Bootscan run length (block permutation)",
            statistic="mean bootstrap support for the donor",
            breakpoint_lo=q_start, breakpoint_hi=q_end,
        ))
    return regions, major, []
