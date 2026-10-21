"""pick_pool_query.py <pool-dir>: pick an unrecombined pool tip to use as a gate query.

Prints ``<path>\t<clade>\t<clade size>``.

Plan C Task 6's end-to-end gate originally took the 200th file in the pool directory. For
the HIV-1 pool that chose AY900575, typed CRF24_BG, a designated circulating recombinant
form, so the gate's stated premise -- that each query is an unrecombined tree tip -- did
not hold, and a region reported on that query was not evidence of a false positive.

A tip qualifies here only when its label is a clade of the pool in its own right: not one
of the sentinels that mean "no clade", not a composite of two other clades of the same
pool, and not a recombinant-named lineage. Clades below ``--min-clade-size`` tips are
skipped, so the query comes from a clade the panel can also represent.
"""

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

from tessera.recomb.typing import is_recombinant_lineage, pool_labels


def choose_query(
    pool: Path, min_clade_size: int = 5
) -> tuple[Path, str, int] | None:
    """``(path, clade, clade size)`` of the first qualifying tip, or ``None``.

    Deterministic: tips are considered in sorted order, so the same pool always
    yields the same query. Shared with ``run_pool_specificity.py``, which needs the
    identical choice to keep its runs comparable between invocations.
    """
    files = sorted(pool.glob("*.fasta"))
    if not files:
        return None

    labels = pool_labels(files)  # sentinels and composites are already left out
    sizes = collections.Counter(labels.values())
    path_of = {p.stem: p for p in files}

    for stem in sorted(labels):
        clade = labels[stem]
        if is_recombinant_lineage(clade) or sizes[clade] < min_clade_size:
            continue
        return path_of[stem], clade, sizes[clade]
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("pool", type=Path, help="a cached Nextclade or NCBI Virus pool directory")
    ap.add_argument("--min-clade-size", type=int, default=5,
                    help="skip clades with fewer tips than this (default 5)")
    args = ap.parse_args()

    if not sorted(args.pool.glob("*.fasta")):
        print(f"no genomes in {args.pool}", file=sys.stderr)
        return 2

    picked = choose_query(args.pool, args.min_clade_size)
    if picked is None:
        print(f"no unrecombined tip in {args.pool}", file=sys.stderr)
        return 1
    path, clade, size = picked
    print(f"{path}\t{clade}\t{size}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
