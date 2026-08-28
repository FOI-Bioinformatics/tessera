#!/usr/bin/env python
"""Opt-in benchmark: Tessera's *region callers* on the Jaya 2023 simulated set.

`run_benchmark.py` scores the parent-free PHI test on this data. PHI answers "does this
alignment carry recombination at all", which is the PhiPack column of that study. The
region callers -- the HMM, 3SEQ, MaxChi and Bootscan ensemble -- are the part of Tessera
comparable to the study's *sequence-level* methods (3SEQ, GENECONV, RDP, MaxChi,
Chimaera), and they were never scored against it. This harness closes that.

Jaya, Brito & Darling (2023), "Evaluation of recombination detection methods for viral
sequencing", *Virus Evolution* 9(2):vead066; data on Dryad (doi:10.5061/dryad.d7wm37q6f),
alignments named ``msa_m<mut>_rc<rec>_n<seqs>_dual<d>_rep<k>.fasta`` with the population
recombination rate in the filename.

**Read the per-query column, not the alignment-level one.** The filename gives a
population recombination *rate*, not which sequences are recombinant, and the study notes
only a small proportion of sequences per population actually are. Tessera scans one query
against a panel, so this harness samples ``--queries`` sequences per alignment.

Two rates are reported. *Alignment-level* ("any sampled query fired") compounds with
``--queries``: at a 3 % per-query false-positive rate, ten queries give
``1-(1-0.03)^10`` = 26 % at the alignment level. Measured, it behaves exactly that way --
the ``mut=0.1`` stratum reports 60/60 alignments at ten queries and 38/60 at three. So the
alignment-level number is **not** comparable to the study's PhiPack column, which runs one
test per alignment, and is only meaningful against another run at the same ``--queries``.

The *per-query* rate does not compound and is the honest figure: it is the per-scan
detection rate, closest in spirit to the study's per-sequence scoring of 3SEQ and GENECONV,
though still not identical because those score every sequence and label recombinants
individually.

**What this measured (3 queries, full grid).** Per-query detection on clonal versus
recombining alignments, by mutation rate -- the gap between the two columns is the caller's
actual discriminating power:

======== ============ ============= ==============
mut rate clonal       recombining   discrimination
======== ============ ============= ==============
0                  0%            0%          0 pt
1e-05              0%            0%          0 pt
0.0001             3%            3%         -1 pt
0.001              0%            1%          1 pt
0.01               7%           17%        +11 pt
0.1               33%           37%         +4 pt
======== ============ ============= ==============

Detection needs diversity -- below ``mut=0.001`` there is nothing to find -- but it also
**degrades at saturation**: at ``mut=0.1`` the callers fire on a third of clonal scans and
barely separate them from recombining ones. The single usable operating point on this grid
is ``mut=0.01``. A high raw detection rate at high divergence is noise, not power, which is
why the aggregate "power/specificity" line below averages three regimes that behave
differently and should not be quoted on its own.

Needs the Dryad data; SKIPs without it. No aligner and no network required -- the
alignments are already aligned. The pure sampling and scoring logic is unit-tested in CI.

    python validation/run_caller_benchmark.py                    # all alignments
    python validation/run_caller_benchmark.py --max 60 --queries 5
    python validation/run_caller_benchmark.py --method hmm,3seq,maxchi
    python validation/run_caller_benchmark.py --rows grid.tsv       # keep the per-alignment grid

A full run is tens of minutes, so progress is streamed to stderr and ``--rows`` writes each
alignment's outcome as it is scored -- an interrupted run still leaves usable results, and
the grid can be re-analysed without re-scanning.
"""

from __future__ import annotations

import logging
import random
import sys
from collections import Counter
from pathlib import Path

from tessera.recomb.analyze import analyze
from tessera.recomb.ensemble import consensus_regions, reconcile_major
from tessera.recomb.regions import DEFAULT_METHODS, RegionParams, call_regions, parse_methods
from tessera.recomb.similarity import _read_alignment, compute_similarity

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_benchmark import (  # noqa: E402 - sibling validation module
    DEFAULT_CACHE,
    parse_msa_name,
    power_specificity,
)

# The alignments are ~1.7 kb, so the shipped 1 kb default window would leave almost no
# windows to scan. These are sized for the data and reported with the results.
WINDOW, STEP = 200, 20
DEFAULT_QUERIES = 10


def pick_queries(labels: list[str], n: int, *, seed: int) -> list[str]:
    """A deterministic sample of up to ``n`` sequence labels to use as queries."""
    if not labels:
        return []
    if len(labels) <= n:
        return list(labels)
    return random.Random(seed).sample(sorted(labels), n)


def alignment_detected(region_counts: list[int]) -> bool:
    """An alignment counts as detected when any sampled query yields a region."""
    return any(count > 0 for count in region_counts)


def scan_query(path: str, query: str, methods: tuple[str, ...]) -> tuple[int, set[str]]:
    """``(number of consensus regions, callers that fired)`` for one query."""
    result = compute_similarity(path, query, window_size=WINDOW, window_step=STEP)
    analysis = analyze(result)
    per_method, majors = {}, {}
    for method in methods:
        regions, major, _ = call_regions(
            result, analysis, WINDOW, RegionParams.with_defaults(WINDOW, method=method)
        )
        per_method[method] = regions
        majors[method] = major
    major, _ = reconcile_major(majors, window_wins=analysis.winners_with_ties)
    consensus, _ = consensus_regions(per_method, major=major)
    fired = {m for region in consensus for m in region.methods}
    return len(consensus), fired


def _flag(argv: list[str], name: str, default):
    return type(default)(argv[argv.index(name) + 1]) if name in argv else default


def main(argv: list[str]) -> int:
    cap = _flag(argv, "--max", 0)
    n_queries = _flag(argv, "--queries", DEFAULT_QUERIES)
    methods = parse_methods(_flag(argv, "--method", ",".join(DEFAULT_METHODS)))

    cache = DEFAULT_CACHE
    alignments = sorted(p for p in cache.rglob("*.fasta") if parse_msa_name(p.name))
    if not alignments:
        print(f"[SKIP] no benchmark alignments in {cache}. Download performance.tar.gz from "
              f"Dryad (doi:10.5061/dryad.d7wm37q6f) and extract the msa_*.fasta files there "
              f"(or set $TESSERA_BENCHMARK_DIR).")
        return 0
    if cap:
        alignments = alignments[:cap]

    logging.getLogger("tessera").addHandler(logging.NullHandler())
    print(f"Tessera region-caller benchmark -- {len(alignments)} alignment(s), "
          f"{n_queries} query/queries each, --method {','.join(methods)}")
    print(f"(window {WINDOW}/{STEP}; detected = any sampled query yields a consensus region)\n")

    rows_path = _flag(argv, "--rows", "")
    rows_fh = open(rows_path, "w") if rows_path else None
    if rows_fh:
        rows_fh.write("alignment\tmut\trec\tdual\trep\tdetected\tqueries\tqueries_hit\tregions\tcallers\n")

    results: list[tuple[float, float]] = []
    per_rec: dict[float, list[bool]] = {}
    per_rec_q: dict[float, list[int]] = {}  # (queries_hit, queries) per alignment
    callers: Counter[str] = Counter()
    for i, path in enumerate(alignments, 1):
        params = parse_msa_name(path.name)
        rows = _read_alignment(str(path))
        counts, fired_here = [], set()
        for query in pick_queries(list(rows), n_queries, seed=0):
            n_regions, fired = scan_query(str(path), query, methods)
            counts.append(n_regions)
            fired_here |= fired
            callers.update(fired)
        hit = alignment_detected(counts)
        # Encode as a pseudo p-value so the PHI benchmark's scorer can aggregate it
        # unchanged: 0.0 = "called recombinant", 1.0 = "did not call".
        results.append((params["rec"], 0.0 if hit else 1.0))
        per_rec.setdefault(params["rec"], []).append(hit)
        per_rec_q.setdefault(params["rec"], []).append(
            (sum(1 for c in counts if c > 0), len(counts)))
        if rows_fh:
            rows_fh.write(
                f"{path.name}\t{params['mut']}\t{params['rec']}\t{params['dual']}\t"
                f"{params['rep']}\t{int(hit)}\t{len(counts)}\t"
                f"{sum(1 for c in counts if c > 0)}\t{sum(counts)}\t"
                f"{','.join(sorted(fired_here))}\n"
            )
            rows_fh.flush()
        print(f"  [{i}/{len(alignments)}] {path.name} rec={params['rec']} "
              f"-> {'DETECTED' if hit else 'none'}", file=sys.stderr, flush=True)
    if rows_fh:
        rows_fh.close()

    print(f"  {'rec-rate':<10} {'n':>4}  {'alignment-level':>16}  {'per-query':>14}")
    for rec in sorted(per_rec):
        hits = per_rec[rec]
        qh = sum(a for a, _ in per_rec_q[rec])
        qn = sum(b for _, b in per_rec_q[rec])
        label = "spec" if rec == 0 else "power"
        print(f"  {rec:<10} {len(hits):>4}  {label} {sum(hits):>3}/{len(hits):<3}"
              f"      {qh:>5}/{qn:<6} ({qh / qn:.0%})")
    print("\n  Alignment-level is 'any of the sampled queries fired', so it compounds with")
    print("  the number of queries; the per-query rate does not. Compare per-query against")
    print("  the study's per-sequence methods, and alignment-level only against itself.")

    agg = power_specificity(results, 0.5)
    power = f"{agg['power']:.2f}" if agg["power"] is not None else "n/a"
    spec = f"{agg['specificity']:.2f}" if agg["specificity"] is not None else "n/a"
    print("\n" + "=" * 60)
    print(f"  power {power} ({agg['n_recombining']} recombining) | "
          f"specificity {spec} ({agg['n_clonal']} clonal)")
    print("\n  callers firing on the regions found:")
    for method in methods:
        print(f"    {method:10s} {callers[method]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
