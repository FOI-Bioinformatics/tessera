#!/usr/bin/env python
"""Opt-in harness: what does Tessera report on data containing NO recombination?

The nine ``neg_pure`` controls in ``run_hybrids.py`` ask a pass/fail question of nine
real, curated Nextclade panels. That measures a verdict on favourable cases; it cannot
estimate a *rate*, and it does not sample the redundant-panel regime -- several
near-equidistant relatives of the query -- where the callers are most likely to mistake
a stochastic window-vote flip for an event.

This harness closes that gap. Sequences are simulated down a fixed clade tree under
JC69, so the ground truth is known by construction: **no recombination is introduced
anywhere**, and every region reported is a false positive. Four scenarios escalate the
confounders this literature blames for spurious recombination signal:

    clean         uniform site rates, uniform lineage rates (the model's own null)
    asrv          gamma(0.2) among-site rate variation
    lineage_rate  the query's own branch evolves 4x faster than the rest
    rate_shift    genome halves evolve at 0.3x and 3x -- differing constraint between
                  genes, the classic false-positive driver (Posada & Crandall 2001)

A matching positive control splices a known clade-B tract into a clade-A backbone, so
**specificity and sensitivity are reported together**: a detection claim without a
paired false-positive rate is half a result.

Unlike the other harnesses this one needs no aligner, no network and no downloaded
data -- the simulated sequences are already aligned, so it is entirely self-contained.
It is opt-in only because a full run is minutes of wall clock; the pure simulation and
scoring logic is unit-tested in CI (``tests/unit/test_specificity_scoring.py``).

    python validation/run_specificity.py                 # 10 replicates per scenario
    python validation/run_specificity.py --reps 3        # quick look
    python validation/run_specificity.py --min-methods 1 # without the agreement gate

Caveat: JC69 on a fixed topology is simpler than real viral evolution, and modest
replicate counts carry real sampling error -- hence the confidence intervals. These
numbers establish whether a failure mode exists and roughly how large it is, not its
exact magnitude on real panels.
"""

from __future__ import annotations

import logging
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np

from tessera.recomb.run import RecombParams, run_recomb
from tessera.recomb.stats import wilson_ci

QUERY = "QUERY"
SCENARIOS = ("clean", "asrv", "lineage_rate", "rate_shift")
CLADES = ("A", "B", "C", "D")
TIPS_PER_CLADE = 4
LENGTH = 12000
TRACT = (4000, 7000)  # donor tract of the positive control, in query coordinates

_BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


# --- simulation ------------------------------------------------------------------

def _evolve(seq: np.ndarray, distance: float, site_rates: np.ndarray, rng) -> np.ndarray:
    """JC69: mutate each site with probability 3/4(1 - exp(-4/3 * d * r))."""
    p = 0.75 * (1.0 - np.exp(-4.0 / 3.0 * distance * site_rates))
    hit = rng.random(seq.size) < p
    out = seq.copy()
    n = int(hit.sum())
    if n:  # a mutation always changes the base
        out[hit] = (seq[hit] + rng.integers(1, 4, size=n)) % 4
    return out


def _site_rates(scenario: str, rng) -> np.ndarray:
    if scenario == "asrv":
        return rng.gamma(0.2, 1 / 0.2, size=LENGTH)
    if scenario == "rate_shift":
        rates = np.empty(LENGTH)
        rates[: LENGTH // 2] = 0.3
        rates[LENGTH // 2 :] = 3.0
        return rates
    return np.ones(LENGTH)


def _as_text(seqs: dict[str, np.ndarray]) -> dict[str, str]:
    return {k: _BASES[v].tobytes().decode("ascii") for k, v in seqs.items()}


def simulate_clonal(scenario: str, seed: int) -> dict[str, str]:
    """A clonal panel plus a non-recombinant query. No recombination, by construction."""
    rng = np.random.default_rng(seed)
    rates = _site_rates(scenario, rng)
    root = rng.integers(0, 4, size=LENGTH)

    seqs: dict[str, np.ndarray] = {}
    ancestors: dict[str, np.ndarray] = {}
    for clade in CLADES:
        ancestor = _evolve(root, 0.08, rates, rng)
        ancestors[clade] = ancestor
        for i in range(TIPS_PER_CLADE):
            seqs[f"{clade}{i}"] = _evolve(ancestor, 0.012, rates, rng)

    # An extra tip inside clade A -- a plain descendant, nothing recombined.
    query_branch = 0.048 if scenario == "lineage_rate" else 0.012
    seqs[QUERY] = _evolve(ancestors["A"], query_branch, rates, rng)
    return _as_text(seqs)


def simulate_recombinant(seed: int) -> tuple[dict[str, str], tuple[int, int]]:
    """The positive control: a clade-A backbone carrying one clade-B tract."""
    rng = np.random.default_rng(seed)
    rates = _site_rates("clean", rng)
    root = rng.integers(0, 4, size=LENGTH)

    seqs: dict[str, np.ndarray] = {}
    ancestors: dict[str, np.ndarray] = {}
    for clade in CLADES:
        ancestor = _evolve(root, 0.08, rates, rng)
        ancestors[clade] = ancestor
        for i in range(TIPS_PER_CLADE):
            seqs[f"{clade}{i}"] = _evolve(ancestor, 0.012, rates, rng)

    backbone = _evolve(ancestors["A"], 0.012, rates, rng)
    donor = _evolve(ancestors["B"], 0.012, rates, rng)
    lo, hi = TRACT
    backbone[lo:hi] = donor[lo:hi]
    seqs[QUERY] = backbone
    return _as_text(seqs), TRACT


# --- scoring ---------------------------------------------------------------------

def score_negative(rows: list[dict]) -> tuple[int, dict[str, int]]:
    """``(false regions, per-caller counts)`` for a run that must report nothing.

    ``donor_absent`` rows are excluded: a coverage gap says "no close reference covers
    this stretch", which is a true statement about the panel rather than a claim of
    recombination. ``run_hybrids.py::_score_neg_pure`` scores a negative control the
    same way.
    """
    present = [r for r in rows if r.get("donor_absent") != "yes"]
    per_caller: Counter[str] = Counter()
    for row in present:
        for method in (row.get("methods") or "").split(","):
            if method.strip():
                per_caller[method.strip()] += 1
    return len(present), dict(per_caller)


def score_positive(
    rows: list[dict], tract: tuple[int, int], *, donor_prefix: str
) -> dict:
    """Detection, donor attribution and breakpoint error against a known tract."""
    lo, hi = tract
    present = [r for r in rows if r.get("donor_absent") != "yes"]
    overlapping = [
        r for r in present
        if int(r["query_start"]) < hi and int(r["query_end"]) > lo
    ]
    if not overlapping:
        return {"detected": False, "donor_ok": False, "breakpoint_error": None}
    # the region covering most of the true tract is the one being judged
    best = max(
        overlapping,
        key=lambda r: min(int(r["query_end"]), hi) - max(int(r["query_start"]), lo),
    )
    return {
        "detected": True,
        "donor_ok": best["minor_parent"].startswith(donor_prefix),
        "breakpoint_error": (
            abs(int(best["query_start"]) - lo) + abs(int(best["query_end"]) - hi)
        ) // 2,
    }


def false_positive_rate(n_with_false: int, n_runs: int) -> tuple[float, float, float]:
    """``(rate, lo, hi)`` -- a rate from a few replicates needs its interval shown."""
    if n_runs <= 0:
        return 0.0, 0.0, 1.0
    lo, hi = wilson_ci(n_with_false, n_runs)
    return n_with_false / n_runs, lo, hi


# --- running ---------------------------------------------------------------------

def _write_fasta(path: Path, seqs: dict[str, str]) -> None:
    with path.open("w") as fh:
        for label, seq in seqs.items():
            fh.write(f">{label}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i : i + 70] + "\n")


def _read_regions(path: Path) -> list[dict]:
    if not path.exists():
        return []
    lines = [x for x in path.read_text().splitlines() if x.strip()]
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, x.split("\t"), strict=False)) for x in lines[1:]]


def _scan(seqs: dict[str, str], logger: logging.Logger, min_methods: int) -> list[dict]:
    """Run the real shipped pipeline on one simulated alignment."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        msa = td / "aln.fasta"
        _write_fasta(msa, seqs)
        out = td / "out"
        run_recomb(
            RecombParams(msa=msa, output=out, query=QUERY, plot_format="png",
                         min_methods=min_methods),
            logger,
        )
        return _read_regions(out / "recombination_regions.tsv")


def main(argv: list[str]) -> int:
    reps = 10
    min_methods = 2
    for flag, target in (("--reps", "reps"), ("--min-methods", "min_methods")):
        if flag in argv:
            value = int(argv[argv.index(flag) + 1])
            if target == "reps":
                reps = value
            else:
                min_methods = value

    logger = logging.getLogger("tessera.specificity")
    logger.addHandler(logging.NullHandler())
    logger.propagate = False  # the scan is chatty; the table below is the output

    print(f"Tessera specificity harness -- {reps} replicate(s) x {len(SCENARIOS)} "
          f"scenario(s), --min-methods {min_methods}")
    print("No recombination is simulated: every region reported is a false positive.\n")

    print(f"{'scenario':14s} {'runs w/ a false region':>23s} {'false regions':>14s}  "
          f"{'per caller'}")
    print("-" * 88)
    total_runs = total_bad = total_regions = 0
    all_callers: Counter[str] = Counter()
    for scenario in SCENARIOS:
        bad = regions = 0
        callers: Counter[str] = Counter()
        for rep in range(reps):
            rows = _scan(simulate_clonal(scenario, seed=1000 + rep), logger, min_methods)
            n_false, per_caller = score_negative(rows)
            if n_false:
                bad += 1
            regions += n_false
            callers.update(per_caller)
        rate, lo, hi = false_positive_rate(bad, reps)
        print(f"{scenario:14s} {f'{bad}/{reps}':>10s} "
              f"{f'({rate:.0%}, CI {lo:.0%}-{hi:.0%})':>13s} {regions:>14d}  "
              + (", ".join(f"{k}={v}" for k, v in sorted(callers.items())) or "--"))
        total_runs += reps
        total_bad += bad
        total_regions += regions
        all_callers.update(callers)

    rate, lo, hi = false_positive_rate(total_bad, total_runs)
    print("-" * 88)
    print(f"{'TOTAL':14s} {f'{total_bad}/{total_runs}':>10s} "
          f"{f'({rate:.0%}, CI {lo:.0%}-{hi:.0%})':>13s} {total_regions:>14d}  "
          + (", ".join(f"{k}={v}" for k, v in sorted(all_callers.items())) or "--"))

    # Paired positive control: specificity means nothing without the matching power.
    print(f"\nPositive control -- clade-B tract spliced at query {TRACT[0]}-{TRACT[1]}:")
    detected = donor_ok = 0
    errors: list[int] = []
    for rep in range(reps):
        seqs, tract = simulate_recombinant(seed=2000 + rep)
        got = score_positive(_scan(seqs, logger, min_methods), tract, donor_prefix="B")
        detected += got["detected"]
        donor_ok += got["donor_ok"]
        if got["breakpoint_error"] is not None:
            errors.append(got["breakpoint_error"])
    median_error = f"{int(np.median(errors))} bp" if errors else "n/a"
    print(f"  detected {detected}/{reps} | correct donor {donor_ok}/{reps} | "
          f"median breakpoint error {median_error}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
