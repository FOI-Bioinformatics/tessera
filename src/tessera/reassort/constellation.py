"""The pure reassortment call: group segments by shared parent, decide the verdict.

Every judgment is pairwise and coverage-aware. For a segment pair the relation is one of
concordant (a strain is near-best for both -> same parent), discordant (a near-best parent
of one is typed in the other's tree but not near-best there -> different parents), or
uninformative (neither's near-best parents appear in the other's tree -> cannot tell). This
distinction is what keeps sparse cross-typing from masquerading as reassortment.

A reassortant verdict needs the discordant segments in different parent groups. Groups are
the connected components of the concordant pairs, and concordance is not transitive: a
middle segment can share a strain with each of two segments that disagree with each other.
Such a discordant pair inside one group gives a verdict of ``undetermined``, with a note,
rather than a reassortant verdict beside a single group that names no parent.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations

DEFAULT_MARGIN = 0.5  # ANI percentage-point window that defines a segment's "near-best" parents
# The widest margin accepted. The window absorbs near-ties between tips and ANI estimation
# noise, a fraction of a point; once it spans the ANI spread of the candidates every pair
# shares a strain and every query reads clonal (measured: 0.5 gives reassortant on a case
# that 40 calls clonal). A fixed sanity bound, not a calibrated optimum.
MAX_MARGIN = 5.0


@dataclass
class PairCall:
    """The relation of one segment pair, with the strains it rests on, in words."""

    segment_a: str
    segment_b: str
    relation: str  # "concordant" | "discordant" | "discordant_within_group" | "uninformative"
    note: str


@dataclass
class ParentGroup:
    segments: list[str]
    parent_strains: list[str]


@dataclass
class Constellation:
    groups: list[ParentGroup] = field(default_factory=list)
    verdict: str = "undetermined"
    pair_notes: list[str] = field(default_factory=list)
    pairs: list[PairCall] = field(default_factory=list)


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
    discordant: list[tuple[str, str]] = []
    pairs: list[PairCall] = []
    for i, j in combinations(segments, 2):
        shared = near[i] & near[j]
        if not near[i] or not near[j]:
            # A segment with no candidate strain has no parent to compare. Unreachable
            # through the CLI, which passes assigned segments only, but the function is
            # public: an empty list must not read as a different parent.
            pairs.append(PairCall(i, j, "uninformative",
                                  f"{i}-{j}: uninformative (no candidate strain for one segment)"))
        elif shared:
            concordant.append((i, j))
            pairs.append(PairCall(i, j, "concordant",
                                  f"{i}-{j}: concordant on {', '.join(sorted(shared))}"))
        elif (near[i] & universes[j]) or (near[j] & universes[i]):
            discordant.append((i, j))
            pairs.append(PairCall(
                i, j, "discordant",
                f"{i}-{j}: discordant (cross-typed, different parents: "
                f"{', '.join(sorted(near[i]))} vs {', '.join(sorted(near[j]))})",
            ))
        else:
            pairs.append(PairCall(i, j, "uninformative",
                                  f"{i}-{j}: uninformative (no shared strain vocabulary)"))

    groups = _components(segments, concordant, near)
    group_of = {s: k for k, g in enumerate(groups) for s in g.segments}
    within = {(i, j) for i, j in discordant if group_of[i] == group_of[j]}
    for pair in pairs:
        if (pair.segment_a, pair.segment_b) in within:
            pair.relation = "discordant_within_group"
            pair.note += (
                "; but both are linked into one group through segments concordant with "
                "each, so the parents cannot be separated at this margin (a smaller "
                "--margin may resolve it)"
            )
    if within:
        verdict = "undetermined"
    elif discordant:
        verdict = "reassortant"
    elif concordant and len(groups) == 1:
        verdict = "clonal"
    else:
        verdict = "undetermined"
    return Constellation(
        groups=groups, verdict=verdict, pair_notes=[p.note for p in pairs], pairs=pairs,
    )


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
        shared: set[str] = set.intersection(*(near[s] for s in comp))
        groups.append(ParentGroup(segments=sorted(comp), parent_strains=sorted(shared)))
    return sorted(groups, key=lambda g: (-len(g.segments), g.segments))
