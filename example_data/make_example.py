"""Generate the synthetic example alignments shipped in this directory.

Two small, pre-aligned MSAs (substitutions only, equal length) that run directly with
``tessera recomb`` -- no aligner needed -- and contrast the detection methods:

- ``divergent.msa.fasta``       parents ~11 %% apart, a large insert. The default HMM
                                caller localizes the A -> B -> A mosaic confidently
                                (small q-value, sharp breakpoint); 3SEQ agrees. The
                                similarity plot shows an obvious dip. The easy case.
- ``cryptic_insert.msa.fasta``  parents ~1 %% apart, a short (800 bp) insert with only
                                ~10 discriminating sites. A fixed base-pair window
                                averages the tract away and the HMM segmentation finds
                                nothing, but ``--method 3seq`` -- which pools the
                                discriminating sites into an exact triplet test --
                                recovers it. This is why the 3SEQ caller exists.

Both queries are an A backbone carrying a B-derived insert. Deterministic (fixed seed),
standard library only. Re-run with: ``python make_example.py``.
"""

from __future__ import annotations

import random
from pathlib import Path

BASES = "ACGT"


def mutate(seq: list[str], rate: float, rng: random.Random) -> list[str]:
    """Return a copy of ``seq`` with each site substituted to a different base at ``rate``."""
    return [rng.choice([b for b in BASES if b != c]) if rng.random() < rate else c for c in seq]


def build(length: int, divergence: float, outgroup_divergence: float,
          insert: tuple[int, int], seed: int) -> dict[str, list[str]]:
    """One recombination scenario: parents A/B each ``divergence`` from a shared ancestor,
    a more distant outgroup, and a query = A backbone with B's bases over ``insert``."""
    rng = random.Random(seed)
    ancestor = [rng.choice(BASES) for _ in range(length)]
    parent_a = mutate(ancestor, divergence, rng)
    parent_b = mutate(ancestor, divergence, rng)
    outgroup = mutate(ancestor, outgroup_divergence, rng)
    lo, hi = insert
    query = parent_a[:lo] + parent_b[lo:hi] + parent_a[hi:]
    return {"query": query, "parent_A": parent_a, "parent_B": parent_b, "outgroup": outgroup}


def write_msa(path: Path, records: dict[str, list[str]]) -> None:
    with open(path, "w") as fo:
        for name, seq in records.items():
            fo.write(f">{name}\n{''.join(seq)}\n")


def main() -> None:
    here = Path(__file__).parent
    # Divergent, large insert: the easy case the HMM handles cleanly.
    write_msa(here / "divergent.msa.fasta",
              build(length=3000, divergence=0.06, outgroup_divergence=0.15,
                    insert=(1000, 2000), seed=1))
    # Low divergence, short insert: the HMM misses it; --method 3seq recovers it.
    write_msa(here / "cryptic_insert.msa.fasta",
              build(length=10000, divergence=0.006, outgroup_divergence=0.0084,
                    insert=(4600, 5400), seed=5))


if __name__ == "__main__":
    main()
