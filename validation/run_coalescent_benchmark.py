#!/usr/bin/env python
"""Posada & Crandall (2001) style coalescent-with-recombination benchmark.

Reproduces the evaluation design of Posada & Crandall (2001, "Evaluation of methods for
detecting recombination from DNA sequences: computer simulations", *PNAS* 98:13757):
simulate neutral sequences under the **coalescent with recombination** across a grid of
population mutation and recombination rates, and measure a detection method's false-positive
rate (**specificity**, at ``recombination_rate=0``) and **power** (``recombination_rate>0``).
Here the method under test is Tessera's parent-free PHI test (``recomb/diagnostics.py``).

The 2001 study's own simulated alignments were never archived, so the design is regenerated
with a modern coalescent-with-recombination simulator, **msprime** (Hudson's algorithm). The
harness is self-contained -- sequences are generated, nothing is downloaded -- and opt-in:
it needs msprime (``pip install msprime``) and SKIPs without it. The pure power/specificity
scorer is shared with ``run_benchmark.py`` and unit-tested in CI.

    python validation/run_coalescent_benchmark.py
    python validation/run_coalescent_benchmark.py --reps 20 --alpha 0.05
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_benchmark import power_specificity  # noqa: E402 - sibling validation module

# Grid, in the coalescent-with-recombination parameterisation. rec-rate 0 is the clonal
# control (specificity); the mutation rates set the diversity (informative-site count).
MUT_RATES = (1e-5, 3e-5)
REC_RATES = (0.0, 1e-6, 5e-6, 2e-5)
SAMPLES = 40
LENGTH = 2000
POP_SIZE = 1e4


def _simulate_phi(mut: float, rec: float, *, seed: int, window: int) -> float | None:
    """Coalescent-with-recombination alignment -> Tessera PHI p-value (or None if the
    alignment has too few informative columns to test)."""
    import msprime  # lazy: opt-in dependency
    import numpy as np

    from tessera.recomb.diagnostics import recombination_signal

    ts = msprime.sim_ancestry(samples=SAMPLES, sequence_length=LENGTH, recombination_rate=rec,
                              population_size=POP_SIZE, ploidy=1, random_seed=seed)
    ts = msprime.sim_mutations(ts, rate=mut, random_seed=seed)
    rows = {f"s{i}": np.frombuffer(a.encode(), dtype=np.uint8)
            for i, a in enumerate(ts.alignments(reference_sequence="A" * LENGTH))}
    sig = recombination_signal(rows, next(iter(rows)), lambda c: c, window=window)
    return None if sig is None else sig.phi_p


def main(argv: list[str]) -> int:
    reps = int(_flag(argv, "--reps", 10))
    alpha = _flag(argv, "--alpha", 0.05)
    window = int(_flag(argv, "--window", 100))

    try:
        import msprime  # noqa: F401
    except ImportError:
        print("[SKIP] coalescent benchmark needs msprime (opt-in): pip install msprime")
        return 0

    print("Posada & Crandall coalescent-with-recombination benchmark (PHI test)")
    print(f"({SAMPLES} samples x {LENGTH} bp, {reps} reps/cell; power = P[PHI p < {alpha} | "
          f"recombining], specificity = P[PHI p >= {alpha} | clonal])\n")
    print(f"  {'mut-rate':>9} {'rec-rate':>9}   power/spec")
    all_results: list[tuple[float, float | None]] = []
    for mut in MUT_RATES:
        for rec in REC_RATES:
            cell = [(rec, _simulate_phi(mut, rec, seed=1000 * i + 1, window=window))
                    for i in range(reps)]
            all_results += cell
            testable = [p for _r, p in cell if p is not None]
            if rec > 0:
                hits = sum(p < alpha for p in testable)
                lab = f"power {hits}/{len(testable)}" if testable else "power n/a"
            else:
                hits = sum(p >= alpha for p in testable)
                lab = f"spec {hits}/{len(testable)}" if testable else "spec n/a"
            print(f"  {mut:>9.0e} {rec:>9.0e}   {lab}")

    agg = power_specificity(all_results, alpha)
    pw = f"{agg['power']:.2f}" if agg["power"] is not None else "n/a"
    sp = f"{agg['specificity']:.2f}" if agg["specificity"] is not None else "n/a"
    print("\n" + "=" * 60)
    print(f"  power {pw} ({agg['n_recombining']} recombining) | "
          f"specificity {sp} ({agg['n_clonal']} clonal) | {agg['n_untestable']} untestable")
    return 0


def _flag(argv: list[str], name: str, default: float) -> float:
    return float(argv[argv.index(name) + 1]) if name in argv else default


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
