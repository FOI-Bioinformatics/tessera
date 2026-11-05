#!/usr/bin/env python
"""Opt-in harness: which short true tracts does a minimum-length floor cost? (issue #121)

``run_indel_specificity.py`` measures the regions reported on clonal data. A floor on region
length removes some of those, and the question it leaves is what it does to true tracts.
This harness splices a donor tract of known length into a backbone, gives every lineage
single-base indels (0.4 % per site), aligns with MAFFT through ``tessera msa`` and scans with
``tessera recomb`` at several ``--min-region`` values. A tract counts as **localised** when the
regions naming the donor cover at least half of it. A region that does not touch the tract
(100 bp margin) is **stray**.

    python validation/run_indel_sensitivity.py              # 6 replicates per cell
    python validation/run_indel_sensitivity.py --reps 2     # quick look

``--min-region 1`` stands for "no floor": it is the smallest value the option accepts. It is
not identical to a run that predates the default 50 bp floor, because an explicit
``--min-region`` also sets the run length that lineage clustering and sibling exclusion use.

Needs ``tessera`` and ``mafft`` on ``PATH`` (append the aligner environment) and no network.
One generator, small cells, one window setting: it shows whether a floor costs detections in
this setting, not what it costs on real genomes.
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import random
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

_PATH = Path(__file__).resolve().parent / "run_indel_specificity.py"
_SPEC = importlib.util.spec_from_file_location("run_indel_specificity", _PATH)
base = importlib.util.module_from_spec(_SPEC)
sys.modules["run_indel_specificity"] = base
_SPEC.loader.exec_module(base)

LENGTHS = (20, 40, 80, 150, 300, 600)
DIVERGENCES = (0.06, 0.12)
INDEL = 0.004
START = 2500
FLOORS = (1, 20, 50, 100)
MARGIN = 100


def simulate(divergence: float, length: int, seed: int) -> tuple[str, dict[str, str]]:
    """A backbone-lineage query carrying a ``length`` bp tract from the other lineage."""
    rng = random.Random(seed)
    root = "".join(rng.choice("ACGT") for _ in range(base.LENGTH))
    lineage_a = base.mutate(root, divergence / 2, rng, INDEL)
    lineage_b = base.mutate(root, divergence / 2, rng, INDEL)
    sub = divergence / 4
    refs = {
        "refA": base.mutate(lineage_a, sub, rng, INDEL),
        "refB": base.mutate(lineage_b, sub, rng, INDEL),
        "out1": base.mutate(root, divergence, rng, INDEL),
        "out2": base.mutate(root, divergence, rng, INDEL),
    }
    spliced = (lineage_a[:START] + lineage_b[START:START + length]
               + lineage_a[START + length:])
    return base.mutate(spliced, sub / 2, rng), refs


def score(rows: list[dict[str, str]], length: int) -> tuple[bool, int]:
    covered: set[int] = set()
    stray = 0
    for row in rows:
        lo, hi = int(row["query_start"]), int(row["query_end"])
        if row["minor_parent"] == "refB":
            covered |= set(range(max(lo, START), min(hi, START + length)))
        if hi <= START - MARGIN or lo >= START + length + MARGIN:
            stray += 1
    return len(covered) / length >= 0.5, stray


def run_replicate(
    work: Path, divergence: float, length: int, rep: int
) -> dict[int, tuple[bool, int]]:
    seed = 9000 + round(divergence * 1000) * 100 + length + rep * 7
    query, refs = simulate(divergence, length, seed)
    shutil.rmtree(work, ignore_errors=True)
    (work / "coll").mkdir(parents=True)
    base.write_fasta(work / "q.fa", {"q": query})
    for name, seq in refs.items():
        base.write_fasta(work / "coll" / f"{name}.fasta", {name: seq})
    msa = subprocess.run(
        ["tessera", "msa", "-q", str(work / "q.fa"), "-c", str(work / "coll"),
         "-o", str(work / "m.fa"), "--aligner", "mafft"],
        capture_output=True, text=True, check=False)
    if msa.returncode:
        raise RuntimeError(f"msa failed: {(msa.stdout + msa.stderr)[-300:]}")
    out: dict[int, tuple[bool, int]] = {}
    for floor in FLOORS:
        scan = subprocess.run(
            ["tessera", "recomb", "-i", str(work / "m.fa"), "-o", str(work / f"out_{floor}"),
             "-q", "q", "--window-size", str(base.WINDOW), "--window-step", str(base.STEP),
             "--min-region", str(floor)],
            capture_output=True, text=True, check=False)
        if scan.returncode:
            raise RuntimeError(f"recomb failed: {(scan.stdout + scan.stderr)[-300:]}")
        out[floor] = score(base.read_regions(work / f"out_{floor}" / "recombination_regions.tsv"),
                           length)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--reps", type=int, default=6, help="replicates per cell (default 6)")
    parser.add_argument("--jobs", type=int, default=3, help="cells run in parallel (default 3)")
    args = parser.parse_args(argv)
    for tool in ("tessera", "mafft"):
        if shutil.which(tool) is None:
            print(f"{tool} is not on PATH", file=sys.stderr)
            return 2
    root = Path(tempfile.mkdtemp(prefix="indel_sens_"))
    jobs = list(itertools.product(DIVERGENCES, LENGTHS, range(args.reps)))

    def one(job):
        div, length, rep = job
        return job, run_replicate(root / f"d{div}_L{length}_r{rep}", div, length, rep)

    cells: dict[tuple[float, int], list[dict[int, tuple[bool, int]]]] = {}
    with ThreadPoolExecutor(args.jobs) as pool:
        for (div, length, _), res in pool.map(one, jobs):
            cells.setdefault((div, length), []).append(res)

    print(f"\nShort true tracts, indels at {INDEL:.1%}, {args.reps} replicates per cell, "
          f"window {base.WINDOW}/{base.STEP}")
    print("localised tracts / stray regions, by --min-region\n")
    head = f"{'divergence':>10} {'tract bp':>9} " + " ".join(f"{f:>12}" for f in FLOORS)
    print(head)
    totals = {f: [0, 0] for f in FLOORS}
    for (div, length), reps in sorted(cells.items()):
        cols = []
        for f in FLOORS:
            hit = sum(r[f][0] for r in reps)
            stray = sum(r[f][1] for r in reps)
            totals[f][0] += hit
            totals[f][1] += stray
            cols.append(f"{hit}/{len(reps)} ({stray})".rjust(12))
        print(f"{div:>10} {length:>9} " + " ".join(cols))
    n = sum(len(r) for r in cells.values())
    print("\ntotal localised of", n, "and stray regions:",
          {f: tuple(v) for f, v in totals.items()})
    shutil.rmtree(root, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
