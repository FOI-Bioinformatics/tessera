"""Parent-free recombination diagnostics: the PHI test and Hudson-Kaplan Rmin.

The HMM and 3SEQ callers both need candidate *parents* in the panel and answer
"which parents, where is the breakpoint". This module answers the prior,
parent-agnostic question -- *is there recombination in this alignment at all, and
where is it concentrated* -- which matters most when the true donor is absent from
the panel and at low divergence.

Both statistics are built from the alignment's **biallelic informative columns** (a
0/1 matrix over the sequences) and the **four-gamete test**: two columns are
*incompatible* when all four gametes (00, 01, 10, 11) occur, which under the
infinite-sites model means at least one recombination happened between them.

- **PHI** (Pairwise Homoplasy Index; Bruen, Philippe & Bryant 2006): the mean
  incompatibility among columns that are *near each other* on the genome. Under
  clonal evolution incompatibility (homoplasy) is independent of genomic distance,
  so nearby columns are no more compatible than a random reordering; under
  recombination nearby columns share a local genealogy and are more compatible, so
  the windowed mean drops below the permutation null. The p-value is one-sided
  (lower = recombination), by permutation of the column order. A sliding window
  gives a per-position *profile* of where the signal concentrates.
- **Rmin** (Hudson & Kaplan 1985): the maximum number of disjoint
  four-gamete-incompatible intervals -- a lower bound on the number of recombination
  events, with the intervals themselves as breakpoint candidates.

Everything is numpy; no external tool. Panels here are small (a query plus a handful
of references) so the z*z incompatibility matrix over the z informative columns is
cheap even when z runs to a few thousand columns.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .similarity import CANONICAL_BASES

# Below this many informative columns the statistics are not meaningful.
_MIN_INFORMATIVE = 10
_PERMUTATIONS = 1000


@dataclass
class RecombinationSignal:
    """Parent-free recombination evidence for an alignment."""

    n_informative: int  # biallelic informative columns used
    phi_p: float  # PHI permutation p-value (one-sided; small = recombination)
    phi_observed: float  # the windowed mean incompatibility
    phi_window: int  # window width, in informative-column ranks
    rmin: int  # Hudson-Kaplan minimum number of recombination events
    rmin_intervals: list[tuple[int, int]] = field(default_factory=list)  # query-coord
    profile: list[tuple[int, int, float]] = field(default_factory=list)  # (msa, query, phi)


def biallelic_columns(rows: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(allele1, allele0, positions)`` for the parsimony-informative columns.

    A column is kept when exactly two canonical bases segregate, each carried by at
    least two sequences (so no third allele and no singletons). ``allele1``/``allele0``
    are ``(z, N)`` 0/1 indicator matrices (one row per kept column) marking which
    sequences carry the minor / major allele; a non-canonical sequence there is 0 in
    both (missing, excluded from the four-gamete counts). ``positions`` are the
    columns' MSA indices. The candidate columns are found vectorised, so only the few
    kept columns are looped over.
    """
    labels = list(rows)
    stack = np.array([rows[label] for label in labels], dtype=np.uint8)  # (N, width)
    bases = np.array(sorted(CANONICAL_BASES), dtype=np.uint8)
    counts = np.stack([(stack == b).sum(axis=0) for b in bases])  # (4, width)
    keep_cols = np.where(((counts >= 1).sum(axis=0) == 2) & ((counts >= 2).sum(axis=0) == 2))[0]
    if keep_cols.size == 0:
        empty = np.zeros((0, len(labels)))
        return empty, empty, np.array([], dtype=int)

    a1_cols: list[np.ndarray] = []
    a0_cols: list[np.ndarray] = []
    for col in keep_cols:
        lo, hi = np.where(counts[:, col] >= 2)[0]  # two distinct base indices, lo < hi
        # Major = more frequent; an equal-count tie goes to the lower base value
        # (deterministic). The two bases are distinct so the indicators never collide;
        # which one is "major" only labels them and the four-gamete test is symmetric.
        major, minor = (lo, hi) if counts[lo, col] >= counts[hi, col] else (hi, lo)
        column = stack[:, col]
        a1_cols.append((column == bases[minor]).astype(np.float64))
        a0_cols.append((column == bases[major]).astype(np.float64))
    return np.array(a1_cols), np.array(a0_cols), keep_cols.astype(int)


def incompatibility_matrix(allele1: np.ndarray, allele0: np.ndarray) -> np.ndarray:
    """The ``z x z`` four-gamete incompatibility matrix (True = all four gametes present).

    Gamete counts come from indicator dot products, so sequences missing at either
    column never contribute a gamete. The diagonal is False.
    """
    n11 = allele1 @ allele1.T
    n10 = allele1 @ allele0.T
    n01 = allele0 @ allele1.T
    n00 = allele0 @ allele0.T
    incompatible = (n11 > 0) & (n10 > 0) & (n01 > 0) & (n00 > 0)
    np.fill_diagonal(incompatible, False)
    return incompatible


def phi(incompatible: np.ndarray, window: int) -> float:
    """Mean incompatibility over column pairs within ``window`` ranks (the PHI statistic)."""
    z = incompatible.shape[0]
    total = count = 0
    for offset in range(1, min(window, z - 1) + 1):
        diag = np.diagonal(incompatible, offset)
        total += int(diag.sum())
        count += diag.size
    return total / count if count else float("nan")


def phi_pvalue(
    incompatible: np.ndarray, window: int, *, permutations: int = _PERMUTATIONS, seed: int = 0
) -> tuple[float, float]:
    """One-sided permutation p-value for PHI; returns ``(p_value, observed)``.

    The column order is permuted ``permutations`` times and the windowed mean
    recomputed; ``p`` is the fraction of permutations whose PHI is at most the
    observed (recombination clusters compatible sites locally, lowering PHI).
    """
    observed = phi(incompatible, window)
    z = incompatible.shape[0]
    if z < 3 or np.isnan(observed):
        return 1.0, observed
    rng = np.random.default_rng(seed)
    perms = np.array([rng.permutation(z) for _ in range(permutations)])  # (K, z)
    width = min(window, z - 1)
    total = np.zeros(permutations)
    count = 0
    for offset in range(1, width + 1):
        a = perms[:, : z - offset]
        b = perms[:, offset:]
        total += incompatible[a, b].sum(axis=1)
        count += z - offset
    perm_phi = total / count
    p = (np.count_nonzero(perm_phi <= observed) + 1) / (permutations + 1)
    return float(p), observed


def phi_profile(
    incompatible: np.ndarray, positions: np.ndarray, column_to_query, window: int
) -> list[tuple[int, int, float]]:
    """Per-position windowed PHI: ``[(msa_pos, query_pos, phi), ...]`` for each
    informative column, the mean incompatibility with its ``window`` neighbours on
    each side. A dip marks a recombination hotspot."""
    z = incompatible.shape[0]
    profile: list[tuple[int, int, float]] = []
    for i in range(z):
        lo, hi = max(0, i - window), min(z, i + window + 1)
        block = incompatible[i, lo:hi]
        denom = block.size - 1  # exclude self
        value = float((block.sum()) / denom) if denom > 0 else float("nan")
        msa_pos = int(positions[i])
        profile.append((msa_pos, int(column_to_query(msa_pos)), value))
    return profile


def hudson_kaplan_rmin(
    incompatible: np.ndarray, positions: np.ndarray
) -> tuple[int, list[tuple[int, int]]]:
    """Hudson-Kaplan ``(Rmin, intervals)`` -- the maximum set of disjoint
    four-gamete-incompatible intervals (a lower bound on recombination events).

    Each incompatible column pair ``(i, j)`` requires a recombination in the rank
    interval ``[i, j]``; the greedy maximum set of non-overlapping intervals (by right
    endpoint) is Rmin. ``intervals`` are the chosen ones in MSA column coordinates.
    """
    pairs = [(int(i), int(j)) for i, j in zip(*np.where(np.triu(incompatible, 1)), strict=True)]
    if not pairs:
        return 0, []
    pairs.sort(key=lambda p: (p[1], p[0]))  # by right endpoint
    rmin = 0
    intervals: list[tuple[int, int]] = []
    last_right = -1
    for i, j in pairs:
        if i >= last_right:  # disjoint from the last selected interval
            rmin += 1
            intervals.append((int(positions[i]), int(positions[j])))
            last_right = j
    return rmin, intervals


def recombination_signal(
    rows: dict[str, np.ndarray], query_label: str, column_to_query, *,
    window: int = 100, seed: int = 0,
) -> RecombinationSignal | None:
    """Compute the parent-free recombination signal for an alignment.

    Returns ``None`` when there are fewer than ``_MIN_INFORMATIVE`` biallelic
    informative columns (too little variation to test). The PHI p-value is reported
    as-is; the significance threshold is a reporting concern, applied downstream.
    ``query_label`` is accepted for interface symmetry; the statistics use every
    sequence in ``rows``.
    """
    allele1, allele0, positions = biallelic_columns(rows)
    z = positions.size
    if z < _MIN_INFORMATIVE:
        return None
    incompatible = incompatibility_matrix(allele1, allele0)
    p, observed = phi_pvalue(incompatible, window, seed=seed)
    rmin, intervals = hudson_kaplan_rmin(incompatible, positions)
    query_intervals = [(int(column_to_query(a)), int(column_to_query(b))) for a, b in intervals]
    profile = phi_profile(incompatible, positions, column_to_query, window)
    return RecombinationSignal(
        n_informative=int(z), phi_p=p, phi_observed=observed, phi_window=window,
        rmin=rmin, rmin_intervals=query_intervals, profile=profile,
    )
