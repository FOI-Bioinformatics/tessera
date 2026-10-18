"""Group panel references into lineages, for cluster-aware recombination calling.

Cold-start recruits many near-duplicate references per lineage. The HMM competes
individual genomes, so duplicates tie every window -- defeating the sibling test and
fragmenting the call. This module clusters the references into lineages and builds a
pooled "view" of the similarity result whose labels are one representative per
lineage, with each cluster's per-window evidence averaged over its members (so a pooled
lineage weighs as one genome). The existing HMM and region caller consume that view
verbatim, so the change is a front-end transform: with one genome per cluster the pooled
view equals the input and behaviour is unchanged.

**Clustering is conservative: it merges only clearly near-identical duplicates.**
Two references join the same lineage iff their per-window agreement (from the aligned
``rows``) stays at or above ``cluster_merge_identity`` with **no region-sized
contiguous run below it**. The absolute floor keeps distinct parents separate even
when they diverge only diffusely (a few % apart, no structural tract -- the case no
relative or structural rule can otherwise catch, e.g. distinct orthopox species);
the contiguous-run test additionally separates *structurally* distinct close parents
(genome-wide above the floor but differing over one tract, e.g. SARS-CoV-2
sublineages) -- the direct analogue of the sibling ``max_beaten_run`` test applied
reference-to-reference. Only genomes that are near-identical *everywhere* merge,
which is exactly the set that ties every window and fragments the call. Setting the
floor conservatively means same-lineage genomes that have drifted a few % apart stay
separate, but those are distinguishable by the caller and do not tie.
"""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np

from .similarity import (
    WindowSimilarity,
    _best_per_window,
    _canonical_mask,
    _informative_column_mask,
)


def _window_bounds(result: WindowSimilarity, window_size: int) -> list[tuple[int, int]]:
    """The ``[start, end)`` column span of each window.

    Informative-site windowing carries each window's exact column span in
    ``result.window_spans``; otherwise the base-pair window is centred on each
    position (matches ``compute_similarity``).
    """
    if result.window_spans:
        return list(result.window_spans)
    half = window_size // 2
    return [(p - half, p - half + window_size) for p in result.positions]


def median_window_step(result) -> int:
    """Typical column distance between consecutive windows, as a positive int.

    Converting a base-pair threshold like ``--min-region`` into a window count needs
    the spacing between windows. Under base-pair windowing the gaps equal
    ``window_step`` except the last, which is shorter when a final window is placed
    flush with the alignment end; the median is unmoved by that one gap. Under
    **informative-site windowing** -- automatic below ~8 % divergence, which is
    exactly where the low-divergence pathogens (mpox, VZV, ebola, within-Omicron) are
    analysed -- ``positions`` are the midpoints of informative-site windows, and their
    spacing tracks how densely polymorphic sites happen to fall. Reading one sampled gap
    and applying it genome-wide then over- or under-states the threshold wherever site
    density differs from that sample; the median is the honest central estimate.
    """
    positions = result.positions
    if len(positions) < 2:
        return 1
    gaps = np.diff(np.asarray(positions, dtype=np.int64))
    return max(1, int(np.median(gaps)))


def _max_low_run(
    comp_cumsum: np.ndarray, match_cumsum: np.ndarray,
    bounds: list[tuple[int, int]], t_floor: float,
) -> int:
    """Longest contiguous run of windows where pairwise agreement is below ``t_floor``.

    Windows with no comparable site (the pair does not overlap there) break the run --
    they carry no evidence either way.
    """
    run = max_run = 0
    for start, end in bounds:
        denom = int(comp_cumsum[end] - comp_cumsum[start])
        if denom == 0:
            run = 0
            continue
        if (match_cumsum[end] - match_cumsum[start]) / denom < t_floor:
            run += 1
            max_run = max(max_run, run)
        else:
            run = 0
    return max_run


class _UnionFind:
    def __init__(self, items: list[str]) -> None:
        self._parent = {x: x for x in items}

    def find(self, x: str) -> str:
        root = x
        while self._parent[root] != root:
            root = self._parent[root]
        while self._parent[x] != root:  # path compression
            self._parent[x], x = root, self._parent[x]
        return root

    def union(self, a: str, b: str) -> None:
        self._parent[self.find(a)] = self.find(b)


def cluster_references(
    result: WindowSimilarity, window_size: int, params
) -> list[list[str]]:
    """Group the reference labels of ``result`` into lineage clusters (list of lists)."""
    labels = list(result.similarities)
    if len(labels) < 2:
        return [[label] for label in labels]
    bounds = _window_bounds(result, window_size)
    step = median_window_step(result) if len(result.positions) > 1 else window_size
    min_run = max(1, math.ceil(params.min_region / max(step, 1)))
    t_floor = params.cluster_merge_identity

    canon = {label: _canonical_mask(result.rows[label]) for label in labels}
    # Informative-site windowing is what a near-identical panel is scanned under, and
    # on such a panel identity over all columns is above the floor for every pair: the
    # whole panel would pool into one lineage and leave the HMM a single state. Compare
    # the pair where the panel varies instead. There the floor again separates true
    # duplicates, which agree, from distinct lineages, which differ at a large share of
    # the informative columns.
    informative = (
        _informative_column_mask(result.rows, labels) if result.window_spans else None
    )
    uf = _UnionFind(labels)
    for i, a in enumerate(labels):
        a_row, a_canon = result.rows[a], canon[a]
        for b in labels[i + 1:]:
            comp = a_canon & canon[b]
            if informative is not None:
                comp = comp & informative
            match = comp & (a_row == result.rows[b])
            comp_cs = np.concatenate(([0], np.cumsum(comp)))
            match_cs = np.concatenate(([0], np.cumsum(match)))
            if _max_low_run(comp_cs, match_cs, bounds, t_floor) < min_run:
                uf.union(a, b)

    components: dict[str, list[str]] = {}
    for label in labels:
        components.setdefault(uf.find(label), []).append(label)
    return list(components.values())


def all_singletons(clusters: list[list[str]]) -> bool:
    return all(len(c) == 1 for c in clusters)


def lineage_clusters(result: WindowSimilarity, window_size: int, params) -> list[list[str]]:
    """The panel's lineage clusters with more than one member; ``[]`` when none apply.

    Applies the same gate as the HMM caller (``cluster_lineages`` on, and between
    ``cluster_min_refs`` and ``cluster_max_refs`` references). Computed once per run and
    shared, so that every caller and the ensemble merge see the same lineages: the HMM
    pools them, Bootscan sums their bootstrap support, and the merge treats two members
    of one cluster as the same donor.
    """
    n_refs = len(result.similarities)
    in_range = params.cluster_min_refs <= n_refs <= params.cluster_max_refs
    if not (params.cluster_lineages and in_range):
        return []
    return [c for c in cluster_references(result, window_size, params) if len(c) > 1]


def _representative(result: WindowSimilarity, members: list[str]) -> str:
    """The cluster's most-covering, best-matching member (deterministic tie-break)."""
    return min(
        members,
        key=lambda m: (-sum(result.denominators[m]), -sum(result.numerators[m]), m),
    )


def clustered_view(
    result: WindowSimilarity, clusters: list[list[str]]
) -> tuple[WindowSimilarity, dict[str, int]]:
    """A representative-labelled view of ``result`` with per-cluster pooled evidence.

    Each cluster becomes one label (its representative) whose per-window
    ``numerators``/``denominators`` are the element-wise **mean** over the cluster's
    members, rounded to whole sites, so a pooled lineage carries one genome's worth of
    evidence. The similarity is the ratio of the unrounded means, which equals the
    pooled ratio of the sums. ``rows`` is left untouched (all member
    genomes remain), so the downstream discordant-site test stays at genome level.
    Returns ``(view, cluster_size)`` mapping each representative to its member count.

    The mean, not the sum: the HMM emission is a log-likelihood over these counts and
    pads every state to the largest count in the window. Summed counts gave a cluster of
    k genomes k times the evidence of a single genome, and charged every singleton state
    a padding penalty for the difference; three near-duplicates of the backbone were
    enough to stop the HMM calling a 1 kb tract from a singleton donor.
    """
    n_win = len(result.positions)
    pooled_num: dict[str, list[int]] = {}
    pooled_den: dict[str, list[int]] = {}
    pooled_sim: dict[str, list[float]] = {}
    cluster_size: dict[str, int] = {}
    for members in clusters:
        rep = _representative(result, members)
        cluster_size[rep] = len(members)
        num = np.mean([result.numerators[m] for m in members], axis=0)
        den = np.mean([result.denominators[m] for m in members], axis=0)
        # Rounding moves a count by at most half a site per window, well inside the
        # binomial noise of the window, and keeps the counts integers for every consumer.
        pooled_num[rep] = [int(round(float(x))) for x in num]
        pooled_den[rep] = [int(round(float(x))) for x in den]
        pooled_sim[rep] = [
            float(num[i] / den[i]) if den[i] > 0 else float("nan") for i in range(n_win)
        ]
    best_sim, best_label = _best_per_window(pooled_sim, n_win)
    view = replace(
        result,
        similarities=pooled_sim,
        numerators=pooled_num,
        denominators=pooled_den,
        best_sim=best_sim,
        best_label=best_label,
    )
    return view, cluster_size
