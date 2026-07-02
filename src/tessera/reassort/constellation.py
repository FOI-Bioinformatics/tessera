"""The pure reassortment call: group segments by shared parent, decide the verdict.

Every judgment is pairwise and coverage-aware. For a segment pair the relation is one of
concordant (a strain is near-best for both -> same parent), discordant (a near-best parent
of one is typed in the other's tree but not near-best there -> different parents), or
uninformative (neither's near-best parents appear in the other's tree -> cannot tell). This
distinction is what keeps sparse cross-typing from masquerading as reassortment.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations

DEFAULT_MARGIN = 0.5  # ANI percentage-point window that defines a segment's "near-best" parents


@dataclass
class ParentGroup:
    segments: list[str]
    parent_strains: list[str]


@dataclass
class Constellation:
    groups: list[ParentGroup] = field(default_factory=list)
    verdict: str = "undetermined"
    pair_notes: list[str] = field(default_factory=list)


def _near_best(ranked: list[tuple[str, float]], margin: float) -> set[str]:
    """The strains within ``margin`` ANI of this segment's best match."""
    if not ranked:
        return set()
    best = ranked[0][1]
    return {strain for strain, ani in ranked if ani >= best - margin}


def call_constellation(
    candidates: dict[str, list[tuple[str, float]]],
    universes: dict[str, set[str]],
    margin: float = DEFAULT_MARGIN,
) -> Constellation:
    """Call reassortment from per-segment candidate strains and dataset strain universes."""
    segments = list(candidates)
    near = {s: _near_best(candidates[s], margin) for s in segments}

    concordant: list[tuple[str, str]] = []
    discordant = False
    notes: list[str] = []
    for i, j in combinations(segments, 2):
        shared = near[i] & near[j]
        if shared:
            concordant.append((i, j))
            notes.append(f"{i}~{j} concordant on {sorted(shared)}")
        elif (near[i] & universes[j]) or (near[j] & universes[i]):
            discordant = True
            notes.append(f"{i}!{j} discordant (cross-typed, different parents)")
        else:
            notes.append(f"{i}?{j} uninformative (no shared strain vocabulary)")

    groups = _components(segments, concordant, near)
    if discordant:
        verdict = "reassortant"
    elif concordant and len(groups) == 1:
        verdict = "clonal"
    else:
        verdict = "undetermined"
    return Constellation(groups=groups, verdict=verdict, pair_notes=notes)


def _components(
    segments: list[str], edges: list[tuple[str, str]], near: dict[str, set[str]]
) -> list[ParentGroup]:
    """Connected components over concordant edges -> one ParentGroup each."""
    parent = {s: s for s in segments}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        parent[find(a)] = find(b)

    members: dict[str, list[str]] = {}
    for s in segments:
        members.setdefault(find(s), []).append(s)

    groups = []
    for comp in members.values():
        shared: set[str] = set.intersection(*(near[s] for s in comp)) if comp else set()
        groups.append(ParentGroup(segments=sorted(comp), parent_strains=sorted(shared)))
    return sorted(groups, key=lambda g: (-len(g.segments), g.segments))
