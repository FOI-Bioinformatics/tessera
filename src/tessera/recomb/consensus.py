"""Per-lineage consensus sequences -- a denoised reference for one clade.

jpHMM (Schultz et al. 2006), the method Tessera's HMM is styled on, represents each
subtype not by an arbitrary member genome but by a *profile* of the whole subtype. A
single member genome is a noisy stand-in: at low divergence it can lose the backbone /
donor competition to a near-identical genome of an adjacent clade by chance. A consensus
-- the majority base at each column across the clade's members -- is the clade's centre
and a far more stable reference. This module builds that consensus; the recruitment path
(:mod:`tessera.discover.nextclade`) can emit one per clade instead of many noisy tips.
"""

from __future__ import annotations

import numpy as np

from .similarity import CANONICAL_BASES


def consensus_sequence(sequences: list[str]) -> str:
    """The majority-base consensus of equally-long sequences (reference coordinates).

    At each column the most frequent canonical base (A/C/G/T) wins; a column with no
    canonical base in any sequence keeps the first sequence's character (a gap or N), so
    the consensus stays the same length as its inputs. A single input is returned as-is.
    """
    if not sequences:
        raise ValueError("consensus_sequence needs at least one sequence")
    if len(sequences) == 1:
        return sequences[0].upper()
    width = len(sequences[0])
    if any(len(s) != width for s in sequences):
        raise ValueError("consensus_sequence needs equal-length sequences")
    stack = np.array([np.frombuffer(s.upper().encode(), dtype=np.uint8) for s in sequences])
    bases = np.array(sorted(CANONICAL_BASES), dtype=np.uint8)  # A, C, G, T
    counts = np.stack([(stack == b).sum(axis=0) for b in bases])  # (4, width)
    majority = bases[np.argmax(counts, axis=0)]  # most frequent canonical base per column
    has_canonical = counts.sum(axis=0) > 0
    out = np.where(has_canonical, majority, stack[0])  # else keep the first row's char
    return out.astype(np.uint8).tobytes().decode()
