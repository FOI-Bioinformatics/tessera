"""Native sliding-window similarity engine.

Replaces the ``recan`` dependency. For each window of alignment columns the
similarity between the query and a reference is the fraction of identical
positions among columns where both sequences carry a canonical base (``1.0`` =
identical, ``0.0`` = no shared positions match). Windows where no position is
comparable score ``nan`` (see below).

The query sequence itself is excluded from the returned per-dataset similarities
(it is the thing being compared against everything else, so it would otherwise
score a constant ``1.0`` and win every window).

Only canonical bases (A/C/G/T) are comparable: a column counts toward a window
only when both sequences carry a canonical base there. Gaps, ``N`` and IUPAC
ambiguity codes are excluded from both the numerator and the denominator, so an
``N`` never inflates identity. A window with no comparable position (e.g. an
inter-contig gap in a fragmented query) is **uninformative** and scores ``nan``
rather than ``0.0``; downstream tallies, statistics and region calling ignore
``nan`` windows.

Each window also carries a **query coordinate**: the alignment is anchored to a
backbone reference, so an MSA column does not equal a query base position. The
mapping is the running count of non-gap characters in the query row, so
``column_to_query(col)`` is the number of query bases before column ``col``
(a 0-based query offset). Query gaps -- backbone insertions relative to the
query -- leave the query coordinate unchanged across those columns; for a
multi-contig query the mapping follows MSA order.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from Bio import AlignIO

from ..core.errors import UserInputError

GAP = ord("-")
# Canonical bases that count toward a window (everything else -- gaps, N, IUPAC
# ambiguity codes -- is treated as non-comparable).
CANONICAL_BASES = frozenset(b"ACGT")
SUPPORTED_METRICS = ("pdist",)


def _canonical_mask(seq: np.ndarray) -> np.ndarray:
    """Boolean mask of positions holding a canonical base (A/C/G/T).

    Shape-preserving, so it works on a 1-D row or a 2-D ``(refs, columns)`` stack.
    """
    mask = np.zeros_like(seq, dtype=bool)
    for base in CANONICAL_BASES:
        mask |= seq == base
    return mask


@dataclass
class WindowSimilarity:
    """Per-window similarity of the query to each reference dataset.

    ``positions`` are window-center MSA columns; ``query_positions`` are the
    corresponding query coordinates (parallel list). ``similarities`` maps each
    reference label to its per-window similarity list (aligned with
    ``positions``). ``query`` is the analysed query's label.
    """

    positions: list[int]
    query_positions: list[int]
    similarities: dict[str, list[float]]
    query: str
    width: int
    query_cumulative: np.ndarray  # length width+1: non-gap counts of the query row
    # Per-window "best any reference can do": the closest reference and its
    # similarity. A low best_sim is the fingerprint of an absent reference --
    # the query's true source for that window is not in the collection.
    best_sim: list[float] = field(default_factory=list)
    best_label: list[str | None] = field(default_factory=list)
    # Comparable query bases per window (canonical A/C/G/T): how much evidence the
    # window carries, used to tell a genuinely divergent gap from a sparse one.
    informative: list[int] = field(default_factory=list)
    # The raw binomial counts behind every ratio: per reference, per window, the
    # matches (numerator) and comparable positions (denominator). These drive the
    # confidence intervals, the significance test and the HMM emissions.
    numerators: dict[str, list[int]] = field(default_factory=dict)
    denominators: dict[str, list[int]] = field(default_factory=dict)
    # The aligned rows (label -> upper-cased uint8 bytes), kept so a region can be
    # tested on its discordant sites -- columns where the query matches one
    # candidate parent but not the other (the only sites that distinguish them).
    rows: dict[str, np.ndarray] = field(default_factory=dict)
    # MSA column span ``[start, end)`` of each window. Empty in the default
    # base-pair windowing; populated by :func:`compute_similarity_informative`,
    # whose windows span a fixed number of informative columns rather than a fixed
    # base-pair width. Only :func:`recomfi.recomb.clusters._window_bounds` reads it.
    window_spans: list[tuple[int, int]] = field(default_factory=list)

    def to_dataframe(self) -> pd.DataFrame:
        """Return a DataFrame indexed by dataset label, columns = MSA positions."""
        return pd.DataFrame.from_dict(
            self.similarities, orient="index", columns=self.positions
        )

    def column_to_query(self, col: int) -> int:
        """Query base offset (0-based) at alignment column ``col``."""
        col = max(0, min(int(col), self.width))
        return int(self.query_cumulative[col])


def _read_alignment(msa_path: str) -> dict[str, np.ndarray]:
    """Read an aligned FASTA into label -> uint8 array of (upper-cased) bytes."""
    alignment = AlignIO.read(msa_path, "fasta")
    rows: dict[str, np.ndarray] = {}
    width: int | None = None
    for record in alignment:
        seq = np.frombuffer(str(record.seq).upper().encode("ascii"), dtype=np.uint8)
        if width is None:
            width = seq.size
        elif seq.size != width:
            raise UserInputError(
                f"Sequence '{record.id}' length {seq.size} != alignment width {width}; "
                "the input must be an aligned FASTA (all rows equal length)."
            )
        rows[record.id] = seq
    if not rows:
        raise UserInputError(f"No sequences found in MSA: {msa_path}")
    return rows


def _window_starts(width: int, window: int, step: int) -> range:
    return range(0, width - window + 1, step)


def _windowed_similarity(
    query_canonical: np.ndarray, query: np.ndarray, reference: np.ndarray,
    window: int, step: int,
) -> tuple[list[float], list[int], list[int]]:
    """Per-window identity fraction plus the raw binomial counts behind it.

    Returns ``(ratios, numerators, denominators)``: ``ratio = numer/denom``
    (``nan`` when ``denom == 0``), ``numer`` = matching canonical positions,
    ``denom`` = comparable canonical positions. The counts feed the confidence
    intervals, the significance test and the HMM emissions.
    """
    comparable = query_canonical & _canonical_mask(reference)
    matches = comparable & (query == reference)

    # Prefix sums let each window be evaluated in O(1).
    comp_cumsum = np.concatenate(([0], np.cumsum(comparable)))
    match_cumsum = np.concatenate(([0], np.cumsum(matches)))

    ratios: list[float] = []
    numers: list[int] = []
    denoms: list[int] = []
    for start in _window_starts(query.size, window, step):
        end = start + window
        denom = int(comp_cumsum[end] - comp_cumsum[start])
        numer = int(match_cumsum[end] - match_cumsum[start])
        numers.append(numer)
        denoms.append(denom)
        ratios.append(float("nan") if denom == 0 else numer / denom)
    return ratios, numers, denoms


def compute_similarity(
    msa_path: str,
    query_label: str,
    window_size: int = 1000,
    window_step: int = 100,
    metric: str = "pdist",
) -> WindowSimilarity:
    """Compute per-window similarity of ``query_label`` to every other sequence.

    Raises :class:`UserInputError` for an unknown metric, a missing query, or a
    window larger than the alignment.
    """
    if metric not in SUPPORTED_METRICS:
        raise UserInputError(
            f"Unknown metric '{metric}'. Supported: {', '.join(SUPPORTED_METRICS)}."
        )
    if window_size <= 0 or window_step <= 0:
        raise UserInputError("window-size and window-step must be positive.")

    rows = _read_alignment(msa_path)
    if query_label not in rows:
        available = ", ".join(sorted(rows))
        raise UserInputError(
            f"Query '{query_label}' not found in the MSA. Supply the query label "
            f"exactly as it appears in the alignment. Available: {available}"
        )

    query = rows[query_label]
    width = query.size
    if window_size > width:
        raise UserInputError(
            f"window-size ({window_size}) is larger than the alignment width ({width})."
        )

    # Cumulative count of non-gap query characters: query_cumulative[k] is the
    # number of query bases in columns [0, k).
    query_cumulative = np.concatenate(([0], np.cumsum(query != GAP))).astype(int)

    starts = list(_window_starts(width, window_size, window_step))
    positions = [start + window_size // 2 for start in starts]
    query_positions = [int(query_cumulative[c]) for c in positions]

    query_canonical = _canonical_mask(query)
    similarities: dict[str, list[float]] = {}
    numerators: dict[str, list[int]] = {}
    denominators: dict[str, list[int]] = {}
    for label, seq in rows.items():
        if label == query_label:
            continue
        ratios, numers, denoms = _windowed_similarity(
            query_canonical, query, seq, window_size, window_step
        )
        similarities[label] = ratios
        numerators[label] = numers
        denominators[label] = denoms

    best_sim, best_label = _best_per_window(similarities, len(positions))

    # Comparable query bases per window: a window may be a "gap" because the query
    # is genuinely divergent from every reference (high informative) or merely
    # because it has few comparable bases (low informative); the two differ.
    qcanon_cumsum = np.concatenate(([0], np.cumsum(query_canonical)))
    informative = [int(qcanon_cumsum[s + window_size] - qcanon_cumsum[s]) for s in starts]

    return WindowSimilarity(
        positions=positions,
        query_positions=query_positions,
        similarities=similarities,
        query=query_label,
        width=width,
        query_cumulative=query_cumulative,
        best_sim=best_sim,
        best_label=best_label,
        informative=informative,
        numerators=numerators,
        denominators=denominators,
        rows=rows,
    )


def _informative_column_mask(
    rows: dict[str, np.ndarray], ref_labels: list[str]
) -> np.ndarray:
    """Boolean mask of MSA columns polymorphic among the references.

    A column is informative when at least two references carry *different*
    canonical bases there -- the only columns that distinguish one reference from
    another, and so the only ones that carry recombination signal. A reference's
    non-canonical positions (gap/N/IUPAC) are excluded from the comparison (a gap
    is not evidence of a different base). The query row plays no part in the test.
    """
    width = next(iter(rows.values())).size
    if len(ref_labels) < 2:
        return np.zeros(width, dtype=bool)
    ref_array = np.array([rows[label] for label in ref_labels], dtype=np.uint8)  # (R, W)
    canon = _canonical_mask(ref_array)  # broadcasts over the 2-D array
    n_canon = canon.sum(axis=0)
    # Sentinels so a non-canonical cell cannot create a false agreement/disagreement:
    # canonical ASCII codes are 65-84, well inside [1, 254].
    col_min = np.where(canon, ref_array, np.uint8(255)).min(axis=0)
    col_max = np.where(canon, ref_array, np.uint8(0)).max(axis=0)
    return (n_canon >= 2) & (col_min != col_max)


def _informative_windowed_similarity(
    query_canonical: np.ndarray, query: np.ndarray, reference: np.ndarray,
    info_cols: np.ndarray, info_window: int, info_step: int,
) -> tuple[list[float], list[int], list[int]]:
    """Per-window identity restricted to ``info_cols``, windowed in informative-site
    index space (``info_window`` sites, step ``info_step``).

    The same prefix-sum trick as :func:`_windowed_similarity`, evaluated only on the
    informative-column subset, so each window is O(1). The comparable filter still
    applies: a gap/N in the query or reference at an informative column drops it from
    that window's denominator (matching the base-pair path's semantics).
    """
    q_info = query[info_cols]
    r_info = reference[info_cols]
    comparable = query_canonical[info_cols] & _canonical_mask(r_info)
    matches = comparable & (q_info == r_info)
    comp_cs = np.concatenate(([0], np.cumsum(comparable)))
    match_cs = np.concatenate(([0], np.cumsum(matches)))

    n_info = int(info_cols.size)
    ratios: list[float] = []
    numers: list[int] = []
    denoms: list[int] = []
    for start in range(0, n_info - info_window + 1, info_step):
        end = start + info_window
        denom = int(comp_cs[end] - comp_cs[start])
        numer = int(match_cs[end] - match_cs[start])
        numers.append(numer)
        denoms.append(denom)
        ratios.append(float("nan") if denom == 0 else numer / denom)
    return ratios, numers, denoms


def compute_similarity_informative(
    msa_path: str,
    query_label: str,
    info_window: int = 40,
    info_step: int = 5,
    metric: str = "pdist",
) -> WindowSimilarity:
    """Like :func:`compute_similarity`, but windows span a fixed number of
    *informative* (polymorphic) columns rather than a fixed base-pair width.

    On near-identical panels (intra-species sets, DNA viruses) a base-pair window
    holds only a handful of discriminating sites diluted by hundreds of identical
    ones, so the HMM's per-reference emission contrast collapses. Restricting each
    window to informative columns concentrates the signal: at an informative column
    the query matches its true local parent and mismatches the others, so the
    emission becomes near-binary and segmentation regains power. Positions, spans
    and coordinates remain in MSA columns, so every downstream consumer is unchanged.

    Raises :class:`UserInputError` when the alignment has fewer than ``info_window``
    informative sites (the caller should fall back to base-pair windowing).
    """
    if metric not in SUPPORTED_METRICS:
        raise UserInputError(
            f"Unknown metric '{metric}'. Supported: {', '.join(SUPPORTED_METRICS)}."
        )
    if info_window <= 0 or info_step <= 0:
        raise UserInputError("informative-window and informative-step must be positive.")

    rows = _read_alignment(msa_path)
    if query_label not in rows:
        available = ", ".join(sorted(rows))
        raise UserInputError(
            f"Query '{query_label}' not found in the MSA. Available: {available}"
        )
    query = rows[query_label]
    width = query.size
    ref_labels = [label for label in rows if label != query_label]

    info_cols = np.where(_informative_column_mask(rows, ref_labels))[0]
    n_info = int(info_cols.size)
    if n_info < info_window:
        raise UserInputError(
            f"Only {n_info} informative site(s) in the alignment (need >= {info_window} "
            "for informative-site windowing); the panel has too little variation."
        )

    query_cumulative = np.concatenate(([0], np.cumsum(query != GAP))).astype(int)
    starts = list(range(0, n_info - info_window + 1, info_step))
    positions = [int((info_cols[s] + info_cols[s + info_window - 1]) // 2) for s in starts]
    query_positions = [int(query_cumulative[c]) for c in positions]
    window_spans = [
        (int(info_cols[s]), int(info_cols[s + info_window - 1]) + 1) for s in starts
    ]

    query_canonical = _canonical_mask(query)
    similarities: dict[str, list[float]] = {}
    numerators: dict[str, list[int]] = {}
    denominators: dict[str, list[int]] = {}
    for label in ref_labels:
        ratios, numers, denoms = _informative_windowed_similarity(
            query_canonical, query, rows[label], info_cols, info_window, info_step
        )
        similarities[label] = ratios
        numerators[label] = numers
        denominators[label] = denoms

    best_sim, best_label = _best_per_window(similarities, len(positions))
    qc_cs = np.concatenate(([0], np.cumsum(query_canonical[info_cols])))
    informative = [int(qc_cs[s + info_window] - qc_cs[s]) for s in starts]

    return WindowSimilarity(
        positions=positions,
        query_positions=query_positions,
        similarities=similarities,
        query=query_label,
        width=width,
        query_cumulative=query_cumulative,
        best_sim=best_sim,
        best_label=best_label,
        informative=informative,
        numerators=numerators,
        denominators=denominators,
        rows=rows,
        window_spans=window_spans,
    )


def informative_site_count(rows: dict[str, np.ndarray], query_label: str) -> int:
    """Number of MSA columns polymorphic among the references (query excluded).

    Lets a caller decide -- from an already-computed result's ``rows`` -- whether a
    base-pair window would carry enough discriminating sites, without re-reading the
    alignment.
    """
    ref_labels = [label for label in rows if label != query_label]
    return int(_informative_column_mask(rows, ref_labels).sum())


def discordant_counts(
    rows: dict[str, np.ndarray], query: str, major: str, minor: str,
    start: int, end: int,
) -> tuple[int, int]:
    """Over columns ``[start, end)``, count sites that distinguish the two parents.

    Returns ``(favor_minor, favor_major)``: canonical sites where the query matches
    the minor but not the major, and vice versa. Sites matching both or neither are
    uninformative for telling the parents apart. Operates on alignment columns, so
    it is free of the window-overlap autocorrelation.
    """
    q = rows[query][start:end]
    big = rows[major][start:end]
    small = rows[minor][start:end]
    canon = _canonical_mask(q) & _canonical_mask(big) & _canonical_mask(small)
    q_minor = q == small
    q_major = q == big
    favor_minor = int(np.count_nonzero(canon & q_minor & ~q_major))
    favor_major = int(np.count_nonzero(canon & q_major & ~q_minor))
    return favor_minor, favor_major


def _best_per_window(
    similarities: dict[str, list[float]], n_windows: int
) -> tuple[list[float], list[str | None]]:
    """Per window, the closest reference and its similarity (nan / None if all nan)."""
    labels = list(similarities)
    if not labels:
        return [float("nan")] * n_windows, [None] * n_windows
    matrix = np.array([similarities[label] for label in labels], dtype=float)
    best_sim: list[float] = []
    best_label: list[str | None] = []
    for i in range(n_windows):
        col = matrix[:, i]
        if np.isfinite(col).any():
            j = int(np.nanargmax(col))
            best_sim.append(float(col[j]))
            best_label.append(labels[j])
        else:
            best_sim.append(float("nan"))
            best_label.append(None)
    return best_sim, best_label
