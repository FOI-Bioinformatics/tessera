"""Post-hoc donor re-attribution.

The other callers attribute a recombinant region's donor by competing individual genomes;
at low divergence a near-identical adjacent-clade genome can win by chance. This step
refines a *called* region's donor -- after detection and the backbone are fixed -- to the
clade whose denoised consensus best matches the query over that region, guarded by a
similarity margin. It never changes detection or the backbone (attribution fixes are
post-hoc), and is a no-op on an untyped panel.
"""

from __future__ import annotations

import logging
from dataclasses import replace

import numpy as np

from .regions import Region
from .similarity import CANONICAL_BASES, WindowSimilarity
from .typing import LineageMap, lineage_of

_BASES = np.array(sorted(CANONICAL_BASES), dtype=np.uint8)  # A C G T byte codes


def _clade_consensus(
    rows: dict[str, np.ndarray], query: str, lineage_map: LineageMap
) -> tuple[dict[str, np.ndarray], dict[str, str]]:
    """Per-clade majority canonical base per column (0 where no canonical base), and a
    representative member label per clade."""
    by_clade: dict[str, list[str]] = {}
    for label in rows:
        if label == query:
            continue
        clade = lineage_of(label, lineage_map)
        if clade:
            by_clade.setdefault(clade, []).append(label)
    cons: dict[str, np.ndarray] = {}
    rep: dict[str, str] = {}
    for clade, members in by_clade.items():
        stack = np.array([rows[m] for m in members])  # (n_members, width)
        counts = np.stack([(stack == b).sum(axis=0) for b in _BASES])  # (4, width)
        majority = _BASES[np.argmax(counts, axis=0)]
        cons[clade] = np.where(counts.sum(axis=0) > 0, majority, 0).astype(np.uint8)
        rep[clade] = members[0]
    return cons, rep


def reattribute_donors(
    regions: list[Region],
    result: WindowSimilarity,
    lineage_map: LineageMap | None,
    major_clade: str | None,
    *,
    margin: float = 0.03,
    min_sites: int = 20,
    logger: logging.Logger | None = None,
) -> list[Region]:
    """Refine each region's donor to the best-matching non-backbone clade consensus.

    For each region, score every clade's consensus against the query over the region's MSA
    columns (fraction of matching canonical positions). Re-label the donor to the best
    non-``major_clade`` clade when it beats the region's current donor clade by ``margin``;
    otherwise leave the region unchanged. A no-op without a lineage map or with fewer than
    two scorable clades. Returns a new list; input regions are not mutated.
    """
    if not lineage_map:
        return regions
    cons, rep = _clade_consensus(result.rows, result.query, lineage_map)
    if len(cons) < 2:
        return regions
    query = result.rows[result.query]
    canon = np.isin(query, _BASES)

    def score(clade: str, lo: int, hi: int) -> float | None:
        c = cons[clade]
        comparable = canon[lo:hi] & (c[lo:hi] != 0)
        n = int(comparable.sum())
        if n < min_sites:
            return None
        match = int((query[lo:hi][comparable] == c[lo:hi][comparable]).sum())
        return match / n

    out: list[Region] = []
    for region in regions:
        lo, hi = region.msa_start, region.msa_end
        current = lineage_of(region.minor_parent, lineage_map)
        cur = score(current, lo, hi) if current in cons else None
        best, best_score = None, -1.0
        for clade in cons:
            if clade == major_clade:
                continue
            s = score(clade, lo, hi)
            if s is not None and s > best_score:
                best, best_score = clade, s
        if best is not None and best != current and best_score - (cur or 0.0) >= margin:
            if logger is not None:
                logger.info(
                    "Re-attributed donor %s -> %s over query %d-%d (consensus sim %.3f "
                    "vs %.3f).", region.minor_parent, rep[best], region.query_start,
                    region.query_end, best_score, cur or 0.0,
                )
            out.append(replace(region, minor_parent=rep[best],
                               mean_sim_minor=round(best_score, 4)))
        else:
            out.append(region)
    return out
