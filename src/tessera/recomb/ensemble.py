"""Combine several callers' regions into one transparent consensus.

The HMM and 3SEQ callers are complementary -- 3SEQ recovers short low-divergence tracts
the HMM segmentation dilutes, the HMM localizes breakpoints better on divergent parents.
Running them together raises recall (the union of what each finds) and, where they agree,
confidence. This module merges the per-method region sets into one annotated list, the
ensemble's output.

The merge is deliberately transparent: there is no fabricated combined score. Two regions
are *the same event* when they overlap in query coordinates and name the same minor
parent; a merged region records exactly **which** callers found it (``Region.methods``)
and whether the parent-free Hudson-Kaplan Rmin signal corroborates it
(``Region.parent_free_support``). Overlapping regions that name different minor parents
are a genuine disagreement and stay separate.
"""

from __future__ import annotations

from .regions import CALLERS, Region
from .typing import LineageMap, lineage_of


def _rank(method: str) -> int:
    return CALLERS.index(method) if method in CALLERS else len(CALLERS)


def _overlap(a: Region, b: Region) -> bool:
    """Do two regions overlap in query coordinates?"""
    return a.query_start < b.query_end and b.query_start < a.query_end


def _same_donor(a: Region, b: Region, lineage_map: LineageMap | None) -> bool:
    """Do two regions name the same donor? Equal genome labels, or -- when references
    are typed -- the same lineage, so two callers that pick different representative
    genomes of one lineage still count as the same event (and as agreement)."""
    if a.minor_parent == b.minor_parent:
        return True
    if lineage_map:
        la = lineage_of(a.minor_parent, lineage_map)
        return la is not None and la == lineage_of(b.minor_parent, lineage_map)
    return False


def reconcile_major(
    majors: dict[str, str | None],
) -> tuple[str | None, dict[str, str]]:
    """The canonical backbone (major parent) plus each method's own major.

    ``majors`` maps each method to the major parent its caller reported -- which exists
    even when the caller found no recombinant region (the no-recombination case), so the
    backbone is always available for the verdict. Callers usually agree on the backbone
    (the genome-wide closest reference); when they differ the HMM's major wins if HMM ran,
    else the most frequently reported one. The per-method mapping is returned so the
    report can surface a disagreement.
    """
    per = {m: major for m, major in majors.items() if major is not None}
    if not per:
        return None, per
    if "hmm" in per:
        return per["hmm"], per
    values = list(per.values())
    canonical = max(values, key=lambda m: (values.count(m), -values.index(m)))
    return canonical, per


def _best(group: list[Region]) -> Region:
    """The most significant member of a group (smallest q-value, then highest support)."""
    return min(group, key=lambda r: (r.qvalue if r.qvalue is not None else 1.0,
                                      -(r.support if r.support is not None else 0.0)))


def _merge_group(group: list[Region], major: str | None) -> Region:
    """One consensus region from a set of overlapping same-minor-parent regions."""
    best = _best(group)
    methods = tuple(sorted({m for r in group for m in r.methods}, key=_rank))
    q_start = min(r.query_start for r in group)
    q_end = max(r.query_end for r in group)
    # Prefer the HMM member's breakpoint (best localization) when present.
    hmm_member = next((r for r in group if "hmm" in r.methods), None)
    bp_lo = (hmm_member or best).breakpoint_lo
    bp_hi = (hmm_member or best).breakpoint_hi
    return Region(
        minor_parent=best.minor_parent,
        major_parent=major or best.major_parent,
        msa_start=min(r.msa_start for r in group),
        msa_end=max(r.msa_end for r in group),
        query_start=q_start,
        query_end=q_end,
        length_bp=q_end - q_start,
        n_windows=max(r.n_windows for r in group),
        mean_sim_minor=best.mean_sim_minor,
        mean_sim_major=best.mean_sim_major,
        margin=best.margin,
        posterior_support=best.posterior_support,
        breakpoint_lo=bp_lo,
        breakpoint_hi=bp_hi,
        support=best.support,
        pvalue=best.pvalue,
        qvalue=best.qvalue,
        donor_undercovered=any(r.donor_undercovered for r in group),
        donor_absent=any(r.donor_absent for r in group),
        minor_cluster_size=best.minor_cluster_size,
        major_cluster_size=best.major_cluster_size,
        methods=methods,
    )


def _group(regions: list[Region], lineage_map: LineageMap | None) -> list[list[Region]]:
    """Group regions that overlap (transitively) and name the same donor (genome or,
    when typed, lineage)."""
    groups: list[list[Region]] = []
    for region in sorted(regions, key=lambda r: (r.minor_parent, r.query_start)):
        placed = False
        for group in groups:
            if _same_donor(group[0], region, lineage_map) and any(
                _overlap(region, member) for member in group
            ):
                group.append(region)
                placed = True
                break
        if not placed:
            groups.append([region])
    return groups


def consensus_regions(
    per_method: dict[str, list[Region]],
    *,
    major: str | None,
    rmin_intervals: list[tuple[int, int]] | None = None,
    lineage_map: LineageMap | None = None,
) -> tuple[list[Region], list[dict]]:
    """Merge per-method regions into one consensus list plus a per-region breakdown.

    ``rmin_intervals`` are the Hudson-Kaplan recombination intervals in query
    coordinates; a consensus region overlapping one is flagged ``parent_free_support``.
    When ``lineage_map`` is given, two callers' regions over the same span merge if they
    name the same *lineage* (not only the same genome), so picking different
    representative genomes of one lineage still counts as agreement. The breakdown lists,
    per consensus region, which callers found it and their support, for the report's
    "method comparison" section.
    """
    all_regions = [r for regions in per_method.values() for r in regions]
    intervals = rmin_intervals or []
    consensus: list[Region] = []
    breakdown: list[dict] = []
    for group in _group(all_regions, lineage_map):
        region = _merge_group(group, major)
        region.parent_free_support = any(
            region.query_start < hi and lo < region.query_end for lo, hi in intervals
        )
        consensus.append(region)
        per: dict[str, float | None] = {}
        for member in group:
            for m in member.methods:
                per[m] = member.support
        breakdown.append({
            "minor_parent": region.minor_parent,
            "query_start": region.query_start,
            "query_end": region.query_end,
            "methods": region.methods,
            "per_method_support": per,
            "parent_free_support": region.parent_free_support,
        })
    consensus.sort(key=lambda r: r.query_start)
    breakdown.sort(key=lambda b: b["query_start"])
    return consensus, breakdown
