"""Assign lineage labels to a genome collection (the ``type-lineages`` command).

Three methods, most-authoritative first: (1) metadata/header mining, (2) nearest
neighbour to a Nextclade dataset's clade-labelled tips, (3) de-novo ANI clustering.
Produces the shared ``lineages.tsv`` so any collection becomes lineage-ready. The only
external tool is the already-optional ``skani``.
"""

from __future__ import annotations

import logging
from pathlib import Path

from ..core.io import strip_sequence_extension
from .panel import skani_query_ani

DEFAULT_REF_ANI_FLOOR = 90.0   # a genome joins a reference clade only above this ANI
DEFAULT_CLUSTER_ANI = 95.0     # de-novo: genomes at/above this ANI share a lineage


def cluster_by_ani(
    genomes: list[Path], edges: list[tuple[Path, Path]]
) -> dict[str, list[Path]]:
    """Connected components of the ANI graph -> de-novo lineages.

    ``edges`` are the genome pairs whose ANI cleared the clustering threshold. Every
    genome is a node; each connected component becomes one ``denovo_N`` lineage,
    numbered in order of first appearance. Returns ``{lineage_label: [genome, ...]}``.
    """
    parent: dict[Path, Path] = {g: g for g in genomes}

    def find(x: Path) -> Path:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        if a in parent and b in parent:
            parent[find(a)] = find(b)

    components: dict[Path, list[Path]] = {}
    for g in genomes:
        components.setdefault(find(g), []).append(g)
    return {f"denovo_{i}": members for i, members in enumerate(components.values(), start=1)}


def _ani_edges(
    genomes: list[Path], *, threshold: float, logger: logging.Logger
) -> list[tuple[Path, Path]]:
    """All-vs-all genome pairs with ANI >= ``threshold`` (pairwise ``skani_query_ani``;
    one skani call per genome vs the rest -- fine for a one-off typing pass)."""
    edges: list[tuple[Path, Path]] = []
    for i, g in enumerate(genomes):
        others = genomes[:i] + genomes[i + 1:]
        if not others:
            continue
        ani = skani_query_ani(g, others, logger)
        edges += [(g, other) for other, (a, _af) in ani.items() if a >= threshold]
    return edges


def assign_by_nearest_reference(
    genomes: list[Path], labeled_tips: dict[Path, str], *, ani_floor: float,
    logger: logging.Logger,
) -> dict[str, str]:
    """Assign each genome the clade of its nearest labelled reference tip.

    ``labeled_tips`` maps a reference-tip FASTA path to its clade. For each genome the
    tip with the highest whole-genome ANI wins, but only when it clears ``ani_floor``;
    genomes below the floor are left unassigned (absent from the result). Returns
    ``{genome_label: clade}``.
    """
    tips = list(labeled_tips)
    assigned: dict[str, str] = {}
    if not tips:
        return assigned
    for g in genomes:
        ani = skani_query_ani(g, tips, logger)
        best = max(ani, key=lambda t: ani[t][0], default=None)
        if best is not None and ani[best][0] >= ani_floor:
            assigned[strip_sequence_extension(g.name)] = labeled_tips[best]
    return assigned
