#!/usr/bin/env python
"""Opt-in harness: what does Tessera report on clonal data that has been aligned?

``run_specificity.py`` simulates sequences that are already aligned, so it has no alignment
step and no insertions or deletions. This harness adds both. Each replicate simulates a
star-like panel (four references, two of them related to the query's lineage), gives every
lineage single-base insertions and deletions at a stated rate, **introduces no
recombination**, aligns the genomes with MAFFT through ``tessera msa`` and scans the result
with ``tessera recomb``. Every region reported is a false positive by construction.

The grid crosses pairwise divergence with the indel rate per site per lineage. The rates are
per-site probabilities on each branch, so 0.004 means roughly 24 indel events along each
6 kb lineage; that is more than most viral genomes carry, and is chosen to find where the
problem starts, not to describe a particular virus.

    python validation/run_indel_specificity.py                # 6 replicates per cell
    python validation/run_indel_specificity.py --reps 2       # quick look
    python validation/run_indel_specificity.py --keep out/    # keep MSAs and region tables

Needs ``tessera`` and ``mafft`` on ``PATH`` (append the aligner environment; prepending it
swaps the Python), and no network. Opt-in only: a full run is several minutes.

Caveat: one generator (random root, uniform substitutions, single-base indels), small
replicate counts, and a single window setting. The counts show whether and where a failure
mode exists, not its magnitude on real data.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import random
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

LENGTH = 6000
WINDOW, STEP = 600, 120
DIVERGENCES = (0.01, 0.03, 0.06, 0.12)
INDEL_RATES = (0.0, 0.001, 0.004, 0.01)
SHORT_BP = 100  # a region at most this long is counted as "short"


def mutate(seq: str, rate: float, rng: random.Random, indel: float = 0.0) -> str:
    """Substitute at ``rate`` per site; delete or insert one base at ``indel`` per site."""
    out: list[str] = []
    for base in seq:
        draw = rng.random()
        if draw < rate:
            out.append(rng.choice([b for b in "ACGT" if b != base]))
        elif draw < rate + indel:
            if rng.random() < 0.5:
                continue  # deletion
            out.append(base)
            out.append(rng.choice("ACGT"))  # insertion after this base
        else:
            out.append(base)
    return "".join(out)


def simulate(divergence: float, indel: float, seed: int) -> tuple[str, dict[str, str]]:
    """A clonal query and its four references. The query descends from ``refA``'s lineage
    and carries no tract from any other reference."""
    rng = random.Random(seed)
    root = "".join(rng.choice("ACGT") for _ in range(LENGTH))
    lineage_a = mutate(root, divergence / 2, rng, indel)
    lineage_b = mutate(root, divergence / 2, rng, indel)
    sub = divergence / 4
    refs = {
        "refA": mutate(lineage_a, sub, rng, indel),
        "refB": mutate(lineage_b, sub, rng, indel),
        "out1": mutate(root, divergence, rng, indel),
        "out2": mutate(root, divergence, rng, indel),
    }
    query = mutate(refs["refA"], divergence / 8, random.Random(seed + 99))
    return query, refs


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


def read_regions(tsv: Path) -> list[dict[str, str]]:
    """Called regions only: donor-absent coverage rows are not a caller's region."""
    with open(tsv, encoding="utf-8", newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t")
                if r["donor_absent"] == "no" and r["methods"] != "NA"]


def callers(row: dict[str, str]) -> list[str]:
    return row["methods"].split(",")


def run_replicate(work: Path, divergence: float, indel: float, rep: int) -> list[dict[str, str]]:
    seed = 7000 + round(divergence * 1000) * 100 + round(indel * 1000) * 10 + rep
    query, refs = simulate(divergence, indel, seed)
    shutil.rmtree(work, ignore_errors=True)
    (work / "coll").mkdir(parents=True)
    write_fasta(work / "q.fa", {"q": query})
    for name, seq in refs.items():
        write_fasta(work / "coll" / f"{name}.fasta", {name: seq})
    msa = subprocess.run(
        ["tessera", "msa", "-q", str(work / "q.fa"), "-c", str(work / "coll"),
         "-o", str(work / "m.fa"), "--aligner", "mafft"],
        capture_output=True, text=True, check=False)
    if msa.returncode:
        raise RuntimeError(f"msa failed: {(msa.stdout + msa.stderr)[-300:]}")
    scan = subprocess.run(
        ["tessera", "recomb", "-i", str(work / "m.fa"), "-o", str(work / "out"), "-q", "q",
         "--window-size", str(WINDOW), "--window-step", str(STEP)],
        capture_output=True, text=True, check=False)
    if scan.returncode:
        raise RuntimeError(f"recomb failed: {(scan.stdout + scan.stderr)[-300:]}")
    return read_regions(work / "out" / "recombination_regions.tsv")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--reps", type=int, default=6, help="replicates per cell (default 6)")
    parser.add_argument("--jobs", type=int, default=2, help="cells run in parallel (default 2)")
    parser.add_argument("--keep", type=Path, default=None,
                        help="keep each replicate's MSA and outputs under this directory")
    args = parser.parse_args(argv)
    for tool in ("tessera", "mafft"):
        if shutil.which(tool) is None:
            print(f"{tool} is not on PATH", file=sys.stderr)
            return 2

    base = args.keep or Path(tempfile.mkdtemp(prefix="indel_spec_"))
    base.mkdir(parents=True, exist_ok=True)
    jobs = list(itertools.product(DIVERGENCES, INDEL_RATES, range(args.reps)))

    def one(job: tuple[float, float, int]):
        div, ind, rep = job
        return job, run_replicate(base / f"d{div}_i{ind}_r{rep}", div, ind, rep)

    results: dict[tuple[float, float], list[list[dict[str, str]]]] = {}
    with ThreadPoolExecutor(args.jobs) as pool:
        for (div, ind, _), regions in pool.map(one, jobs):
            results.setdefault((div, ind), []).append(regions)

    print(f"\nClonal queries through MAFFT, {args.reps} replicates per cell, "
          f"window {WINDOW}/{STEP}, CLI defaults (--min-methods 1)")
    print(f"\n{'divergence':>10} {'indel':>7} {'runs':>5} {'any region':>11} "
          f"{'gate 2':>7} {'short':>6}")
    by_caller: Counter[str] = Counter()
    lengths: list[int] = []
    for (div, ind), reps in sorted(results.items()):
        any1 = sum(bool(r) for r in reps)
        any2 = sum(any(len(callers(x)) >= 2 for x in r) for r in reps)
        short = sum(any(int(x["length_bp"]) <= SHORT_BP for x in r) for r in reps)
        print(f"{div:>10} {ind:>7} {len(reps):>5} {any1:>11} {any2:>7} {short:>6}")
        for regions in reps:
            for row in regions:
                lengths.append(int(row["length_bp"]))
                by_caller.update(callers(row))
    print("\n'any region' and 'gate 2' count runs with at least one region, and with one "
          "called by two or more callers; 'short' counts runs with a region of "
          f"{SHORT_BP} bp or less.")
    print(f"regions by caller: {dict(sorted(by_caller.items()))}")
    print(f"region lengths (bp): {sorted(lengths)}")
    if args.keep is None:
        shutil.rmtree(base, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
