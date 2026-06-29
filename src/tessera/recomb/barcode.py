"""Clade-barcode recombination caller (RecombinHunt / rebar-style lineage attribution).

The other callers compete *individual genomes* on pairwise discriminating sites, which is
fragile when genomes are near-identical -- the source of Tessera's one real weak spot,
parental attribution on close lineages. This caller competes **lineages** instead, via
their *characteristic markers*: columns where one clade is near-fixed for a base that is
rare in every other clade. A marker is denoised across the clade's members and the
comparison is multi-way, so the call is robust to a single near-identical adjacent-clade
genome winning by chance. It is the shared primitive of the recent lineage-attribution
methods (RecombinHunt, rebar, VirusRecom), and dependency-free here because a Nextclade
panel already types every reference by clade.

Active only when the panel is typed (a lineage map); otherwise it returns no regions.
"""

from __future__ import annotations

import numpy as np

from .similarity import CANONICAL_BASES, WindowSimilarity
from .typing import lineage_of

_F_IN = 0.80  # a marker allele is carried by >= this fraction of its clade's members
_F_OUT = 0.20  # ... and by <= this fraction of all other clades' members
_MIN_MARKERS = 5  # a clade needs >= this many markers genome-wide to compete
_MIN_IN_WINDOW = 2  # ... and >= this many in a window to be scored there
_FLOOR = 0.60  # a non-major clade must match >= this fraction of its window markers
_MARGIN = 0.15  # ... and beat the major clade's local match-rate by this
_MIN_RUN = 2  # a region must span at least this many windows


def clade_markers(
    rows: dict[str, np.ndarray], query: str, lineage_map,
    *, f_in: float = _F_IN, f_out: float = _F_OUT, min_markers: int = _MIN_MARKERS,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, str]]:
    """Per-clade characteristic markers from a typed panel alignment.

    Returns ``(cols, alleles, rep)``: for each clade with enough markers, ``cols[clade]``
    are its marker MSA columns (sorted) and ``alleles[clade]`` the clade's near-fixed
    canonical base at each; ``rep[clade]`` is a representative member label. Empty when
    fewer than two clades clear ``min_markers``.
    """
    by_clade: dict[str, list[str]] = {}
    for label in rows:
        if label == query:
            continue
        clade = lineage_of(label, lineage_map)
        if clade:
            by_clade.setdefault(clade, []).append(label)
    if len(by_clade) < 2:
        return {}, {}, {}

    width = rows[query].size
    bases = np.array(sorted(CANONICAL_BASES), dtype=np.uint8)  # A C G T
    counts = {c: np.stack([(np.array([rows[m] for m in members]) == b).sum(axis=0)
                           for b in bases])  # (4, width)
              for c, members in by_clade.items()}
    total = sum(counts.values())  # (4, width) over all clades
    cols: dict[str, np.ndarray] = {}
    alleles: dict[str, np.ndarray] = {}
    idx = np.arange(width)
    for c, cnt in counts.items():
        in_canon = cnt.sum(axis=0)  # comparable members of c per column
        other = total - cnt
        out_canon = other.sum(axis=0)
        maj = np.argmax(cnt, axis=0)  # clade's majority canonical allele index per column
        within = np.divide(cnt[maj, idx], in_canon, out=np.zeros(width), where=in_canon > 0)
        out = np.divide(other[maj, idx], out_canon, out=np.zeros(width), where=out_canon > 0)
        is_marker = (in_canon > 0) & (out_canon > 0) & (within >= f_in) & (out <= f_out)
        marked = np.nonzero(is_marker)[0]
        if marked.size >= min_markers:
            cols[c] = marked
            alleles[c] = bases[maj[marked]]
    if len(cols) < 2:
        return {}, {}, {}
    rep = {c: by_clade[c][0] for c in cols}
    return cols, alleles, rep


def call_regions_barcode(result: WindowSimilarity, analysis, window_size: int, params):
    """Call recombinant regions by clade-marker scanning (typed panels only).

    The query's backbone clade matches the most markers genome-wide; a run of windows
    where another clade matches its local markers better (by ``_MARGIN``) is a region for
    that donor clade. Returns ``(regions, major, [])`` to match ``call_regions``.
    """
    from .clusters import _window_bounds
    from .regions import Region

    labels = list(result.similarities)
    default_major = labels[0] if labels else None
    lineage_map = getattr(params, "lineage_map", None)
    if not lineage_map:
        return [], default_major, []
    cols, alleles, rep = clade_markers(result.rows, result.query, lineage_map)
    if not cols:
        return [], default_major, []

    query = result.rows[result.query]
    clades = list(cols)
    # hit[c]: 1 where the query carries clade c's marker allele, with a cumulative sum for
    # fast per-window match counts; global_rate picks the backbone clade.
    hit = {c: (query[cols[c]] == alleles[c]).astype(np.int64) for c in clades}
    cumhit = {c: np.concatenate([[0], np.cumsum(hit[c])]) for c in clades}
    global_rate = {c: hit[c].sum() / hit[c].size for c in clades}
    major = max(clades, key=lambda c: global_rate[c])

    def window_rate(c: str, start: int, end: int):
        lo = int(np.searchsorted(cols[c], start))
        hi = int(np.searchsorted(cols[c], end))
        n = hi - lo
        return ((cumhit[c][hi] - cumhit[c][lo]) / n, n) if n else (0.0, 0)

    bounds = _window_bounds(result, window_size)
    winners: list[str] = []
    supports: list[float] = []
    for start, end in bounds:
        maj_rate, _ = window_rate(major, start, end)
        best, best_rate = major, -1.0
        for c in clades:
            if c == major:
                continue
            rate, n = window_rate(c, start, end)
            if (n >= _MIN_IN_WINDOW and rate >= _FLOOR
                    and rate >= maj_rate + _MARGIN and rate > best_rate):
                best, best_rate = c, rate
        winners.append(best)
        supports.append(best_rate if best != major else 0.0)

    regions: list[Region] = []
    w = 0
    while w < len(winners):
        if winners[w] == major:
            w += 1
            continue
        donor = winners[w]
        j = w
        while j < len(winners) and winners[j] == donor:
            j += 1
        if j - w >= _MIN_RUN:
            msa_start = int(result.positions[w])
            msa_end = int(result.positions[j - 1]) + 1
            q_start = result.column_to_query(msa_start)
            q_end = result.column_to_query(msa_end)
            minor_lbl, major_lbl = rep[donor], rep[major]
            sim_minor = float(np.nanmean([result.similarities[minor_lbl][k] for k in range(w, j)]))
            sim_major = float(np.nanmean([result.similarities[major_lbl][k] for k in range(w, j)]))
            regions.append(Region(
                minor_parent=minor_lbl, major_parent=major_lbl,
                msa_start=msa_start, msa_end=msa_end,
                query_start=q_start, query_end=q_end,
                length_bp=msa_end - msa_start, n_windows=j - w,
                mean_sim_minor=round(sim_minor, 4), mean_sim_major=round(sim_major, 4),
                margin=round(sim_minor - sim_major, 4),
                support=round(float(np.mean(supports[w:j])), 3),
                breakpoint_lo=q_start, breakpoint_hi=q_end,
            ))
        w = j
    return regions, rep[major], []
