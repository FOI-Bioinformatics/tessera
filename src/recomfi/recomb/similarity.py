"""Native sliding-window similarity engine.

Replaces the ``recan`` dependency. For each window of alignment columns the
similarity between the query and a reference is the fraction of identical
positions among columns where neither sequence has a gap (``1.0`` = identical,
``0.0`` = no shared positions match). Windows where no position is comparable
score ``0.0``.

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

from dataclasses import dataclass

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
    """Boolean mask of positions holding a canonical base (A/C/G/T)."""
    mask = np.zeros(seq.size, dtype=bool)
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
) -> list[float]:
    """Per-window identity fraction over columns where both bases are canonical.

    Returns ``nan`` for a window with no comparable position (uninformative).
    """
    comparable = query_canonical & _canonical_mask(reference)
    matches = comparable & (query == reference)

    # Prefix sums let each window be evaluated in O(1).
    comp_cumsum = np.concatenate(([0], np.cumsum(comparable)))
    match_cumsum = np.concatenate(([0], np.cumsum(matches)))

    out: list[float] = []
    for start in _window_starts(query.size, window, step):
        end = start + window
        denom = comp_cumsum[end] - comp_cumsum[start]
        if denom == 0:
            out.append(float("nan"))
        else:
            numer = match_cumsum[end] - match_cumsum[start]
            out.append(float(numer) / float(denom))
    return out


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
    for label, seq in rows.items():
        if label == query_label:
            continue
        similarities[label] = _windowed_similarity(
            query_canonical, query, seq, window_size, window_step
        )

    return WindowSimilarity(
        positions=positions,
        query_positions=query_positions,
        similarities=similarities,
        query=query_label,
        width=width,
        query_cumulative=query_cumulative,
    )
