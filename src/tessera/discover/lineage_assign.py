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
from ..recomb.typing import (
    build_lineage_map,
    first_header,
    titles_from_collection,
)
from .nextclade import build_pool, resolve_dataset
from .panel import dereplicate, skani_query_ani

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


def _reference_tips(
    *, query: Path, nextclade_dataset: str | None, email: str | None,
    cache_dir: Path | None, logger: logging.Logger,
) -> dict[Path, str]:
    """Resolve a Nextclade dataset and return its clade-labelled tips, dereplicated per
    clade. Empty when no dataset resolves. ``resolve_dataset`` auto-detects the dataset
    from ``query`` (a representative collection genome) unless ``nextclade_dataset`` is
    given. Each tip FASTA header is ``>{accession} {clade}``, so the clade is the
    header's remainder."""
    from ..core.cache import nextclade_cache

    try:
        # resolve_dataset(query, override, *, email, logger): the override (an explicit
        # dataset path) wins; otherwise the dataset is auto-detected from the query genome.
        dataset = resolve_dataset(query, nextclade_dataset, email=email, logger=logger)
    except Exception as exc:  # noqa: BLE001 - no dataset is a normal, non-fatal outcome
        logger.info("No Nextclade dataset for reference typing (%s); skipping.", exc)
        return {}
    cache = cache_dir if cache_dir is not None else nextclade_cache(dataset.path, dataset.tag)
    tips = build_pool(dataset, cache_dir=cache, logger=logger)
    clade_of: dict[Path, str] = {}
    for tip in tips:
        header = first_header(tip)
        parts = header.split(None, 1)
        clade = parts[1].strip() if len(parts) > 1 else ""
        if clade and clade not in ("NA", "example"):
            clade_of[tip] = clade
    # one representative per clade so nearest-neighbour is not dominated by a big clade
    label_by_name = {strip_sequence_extension(p.name): c for p, c in clade_of.items()}
    reps, _ = dereplicate(list(clade_of), clade_of=label_by_name, logger=logger)
    return {p: clade_of[p] for p in reps if p in clade_of}


def assign_lineages(
    genomes: list[Path], *, user_lineage_map: Path | None = None,
    taxon: str | None = None, nextclade_dataset: str | None = None,
    ref_ani_floor: float = DEFAULT_REF_ANI_FLOOR, cluster_ani: float = DEFAULT_CLUSTER_ANI,
    email: str | None = None, cache_dir: Path | None = None, logger: logging.Logger,
) -> list[tuple[str, str, str]]:
    """Assign a lineage to every genome via the metadata -> nearest-reference ->
    de-novo ladder. Returns ``(label, lineage, source)`` rows for ``write_lineage_map``."""
    rows = build_lineage_map(
        user_tsv=user_lineage_map,
        title_by_label=titles_from_collection(genomes),
        organism=taxon,
    )
    typed = {label for label, _lin, _src in rows}
    untyped = [g for g in genomes if strip_sequence_extension(g.name) not in typed]

    if untyped:
        labeled = _reference_tips(query=untyped[0], nextclade_dataset=nextclade_dataset,
                                  email=email, cache_dir=cache_dir, logger=logger)
        if labeled:
            nn = assign_by_nearest_reference(untyped, labeled, ani_floor=ref_ani_floor,
                                             logger=logger)
            rows += [(label, clade, "nextclade-nn") for label, clade in nn.items()]
            typed |= set(nn)
            untyped = [g for g in untyped if strip_sequence_extension(g.name) not in typed]

    if len(untyped) >= 2:
        edges = _ani_edges(untyped, threshold=cluster_ani, logger=logger)
        for lineage, members in cluster_by_ani(untyped, edges).items():
            rows += [(strip_sequence_extension(g.name), lineage, "denovo") for g in members]
    elif untyped:
        rows.append((strip_sequence_extension(untyped[0].name), "denovo_1", "denovo"))

    rows.sort()
    logger.info("Typed %d/%d genome(s) (%d untyped remain).",
                len({r[0] for r in rows}), len(genomes),
                len(genomes) - len({r[0] for r in rows}))
    return rows
