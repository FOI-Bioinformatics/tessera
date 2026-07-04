#!/usr/bin/env python
"""Benchmark the `tessera reassort` verdict with precision / recall / F1.

Scores the shipped reassortment caller (``tessera.reassort.assign_segments``) the way the
influenza-reassortment literature does (TreeSort / TreeKnit / CoalRe use precision, recall and
F1 against known reassortment events). Flu H3N2 HA and NA Nextclade tips are labelled by strain,
and the same strain is typed in both segment datasets, so a labelled query set can be built:

- **clonal**: the HA and NA of one strain that is present in BOTH datasets (truly one parent);
- **reassortant**: the HA of strain A with the NA of strain B (two parents).

Each query is run through ``assign_segments``; the ``reassortant`` verdict is the positive class.
Needs skani and network access to the Nextclade datasets; opt-in, not part of CI. The pure
scoring (confusion matrix + F1) is unit-tested.

    export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
    python validation/run_reassort_benchmark.py            # default 15 clonal + 15 reassortant
    python validation/run_reassort_benchmark.py --n 25
"""

from __future__ import annotations

import logging
import random
import sys
import tempfile
from pathlib import Path

from tessera.core.cache import nextclade_cache
from tessera.core.io import read_fasta, strip_sequence_extension, write_fasta_record
from tessera.discover.nextclade import build_pool, resolve_dataset
from tessera.reassort import assign_segments

HA_DATASET = "nextstrain/flu/h3n2/ha/EPI1857216"
NA_DATASET = "nextstrain/flu/h3n2/na/EPI1857215"


def prf1(tp: int, fp: int, fn: int) -> tuple[float, float, float]:
    """Precision, recall, F1 from confusion counts (0.0 when undefined)."""
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return precision, recall, f1


def score(labeled: list[tuple[bool, str]]) -> dict:
    """Score ``[(is_reassortant_truth, verdict), ...]`` with reassortant as the positive class.

    ``undetermined`` on a reassortant counts as a miss (FN); on a clonal it is not a false
    call (reported separately as ``undetermined_clonal``)."""
    tp = fp = fn = tn = undet_clonal = 0
    for truth, verdict in labeled:
        called = verdict == "reassortant"
        if truth and called:
            tp += 1
        elif truth and not called:
            fn += 1
        elif not truth and called:
            fp += 1
        elif not truth and verdict == "undetermined":
            undet_clonal += 1
        else:
            tn += 1
    precision, recall, f1 = prf1(tp, fp, fn)
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "undetermined_clonal": undet_clonal,
            "precision": precision, "recall": recall, "f1": f1}


def _tip_seqs(dataset_path: str, logger: logging.Logger) -> dict[str, str]:
    """``{strain -> segment sequence}`` for every reconstructed tip of a Nextclade dataset."""
    ds = resolve_dataset(Path("/dev/null"), dataset_path, email=None, logger=logger)
    tips = build_pool(ds, cache_dir=nextclade_cache(ds.path, ds.tag), logger=logger)
    out: dict[str, str] = {}
    for tip in tips:
        records = read_fasta(tip)
        if records:
            out[strip_sequence_extension(tip.name)] = records[0][1]
    return out


def _verdict(ha_seq: str, na_seq: str, tmp: Path, logger: logging.Logger) -> str:
    q = tmp / "q.fasta"
    with open(q, "w") as fo:
        write_fasta_record(fo, "HA", ha_seq)
        write_fasta_record(fo, "NA", na_seq)
    result = assign_segments(q, dataset_overrides={"HA": HA_DATASET, "NA": NA_DATASET},
                             logger=logger)
    return result.verdict


def main(argv: list[str]) -> int:
    logging.basicConfig(level=logging.ERROR, format="%(message)s")
    logger = logging.getLogger("reassort-benchmark")
    n = int(_opt(argv, "--n", "15"))

    from tessera.discover.panel import skani_available
    if not skani_available():
        print("[SKIP] reassort benchmark needs skani (and network to the Nextclade datasets).")
        return 0

    ha = _tip_seqs(HA_DATASET, logger)
    na = _tip_seqs(NA_DATASET, logger)
    shared = sorted(set(ha) & set(na))
    if len(shared) < 2 * n:
        print(f"[SKIP] only {len(shared)} strains cross-typed in both segments; need >= {2 * n}.")
        return 0

    rng = random.Random(0)
    rng.shuffle(shared)
    clonal_strains = shared[:n]
    react_pairs = [(shared[n + i], shared[n + (i + 1) % n]) for i in range(n)]  # HA(a) + NA(b)

    labeled: list[tuple[bool, str]] = []
    with tempfile.TemporaryDirectory() as tmp:
        for s in clonal_strains:
            labeled.append((False, _verdict(ha[s], na[s], Path(tmp), logger)))
        for a, b in react_pairs:
            labeled.append((True, _verdict(ha[a], na[b], Path(tmp), logger)))

    res = score(labeled)
    print(f"\nTessera reassort verdict benchmark -- {n} clonal + {n} reassortant queries")
    print(f"(flu H3N2 HA+NA, {len(shared)} cross-typed strains; reassortant = positive class)\n")
    print(f"  TP {res['tp']}  FP {res['fp']}  FN {res['fn']}  TN {res['tn']}  "
          f"(undetermined-on-clonal {res['undetermined_clonal']})")
    print(f"  precision {res['precision']:.2f}  recall {res['recall']:.2f}  F1 {res['f1']:.2f}")
    return 0


def _opt(argv: list[str], name: str, default: str) -> str:
    return argv[argv.index(name) + 1] if name in argv else default


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
