#!/usr/bin/env python
"""Opt-in benchmark: Tessera's parent-free PHI/Rmin detection on a published simulated set.

Jaya, Brito & Darling (2023), "Evaluation of recombination detection methods for viral
sequencing", *Virus Evolution* 9(2):vead066, released SANTA-SIM alignments on Dryad
(doi:10.5061/dryad.d7wm37q6f), named ``msa_m<mut>_rc<rec>_n<seqs>_dual<d>_rep<k>.fasta`` with
the recombination rate encoded in the filename (there are no explicit breakpoint coordinates
in the archive). Tessera's PHI test (``recomb/diagnostics.py``) is an *alignment-level*
recombination test, so this harness runs it on each alignment and reports detection **power**
(``rc>0``) and **specificity** (``rc=0``) across the mutation x recombination grid -- directly
comparable to the PhiPack/Profile column of that study.

Dryad gates its downloads behind a login, so this harness does not fetch the data. Download
``performance.tar.gz`` (and optionally ``scale.tar.gz``) from the Dryad page in a browser and
extract the ``.fasta`` files into the benchmark cache; the harness SKIPs when none are present.

    # default cache: validation/data/benchmark/  (override with $TESSERA_BENCHMARK_DIR)
    python validation/run_benchmark.py                    # power / specificity table
    python validation/run_benchmark.py --alpha 0.05 --window 100 --max 200
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

from tessera.recomb.diagnostics import recombination_signal
from tessera.recomb.similarity import _read_alignment

HERE = Path(__file__).resolve().parent
DEFAULT_CACHE = HERE / "data" / "benchmark"

_NAME = re.compile(
    r"^msa_m(?P<mut>[0-9.]+)_rc(?P<rec>[0-9.]+)_n(?P<n>\d+)_dual(?P<dual>[0-9.]+)_rep(?P<rep>\d+)"
)


def parse_msa_name(name: str) -> dict | None:
    """Parse a Dryad ``msa_m..._rc..._n..._dual..._rep....fasta`` name -> parameter dict.

    Returns ``None`` for a name that does not match (so unrelated files are skipped)."""
    m = _NAME.match(name)
    if not m:
        return None
    g = m.groupdict()
    return {"mut": float(g["mut"]), "rec": float(g["rec"]),
            "n": int(g["n"]), "dual": float(g["dual"]), "rep": int(g["rep"])}


def phi_pvalue_of(path: Path, *, window: int) -> float | None:
    """PHI permutation p-value for one alignment, or ``None`` when it has too few
    informative columns to test (the parent-free signal is undefined)."""
    rows = _read_alignment(str(path))
    if not rows:
        return None
    label = next(iter(rows))
    sig = recombination_signal(rows, label, lambda c: c, window=window)
    return None if sig is None else sig.phi_p


def power_specificity(results: list[tuple[float, float | None]], alpha: float) -> dict:
    """Aggregate ``[(rec_rate, phi_p_or_None), ...]`` into power/specificity.

    Power = fraction of *testable* recombining alignments (``rec>0``) with ``phi_p < alpha``;
    specificity = fraction of *testable* non-recombining alignments (``rec=0``) with
    ``phi_p >= alpha``. Untestable alignments (too few informative columns) are counted
    separately and excluded from the rates."""
    pos = [p for rec, p in results if rec > 0]
    neg = [p for rec, p in results if rec == 0]
    pos_t = [p for p in pos if p is not None]
    neg_t = [p for p in neg if p is not None]
    power = sum(p < alpha for p in pos_t) / len(pos_t) if pos_t else None
    spec = sum(p >= alpha for p in neg_t) / len(neg_t) if neg_t else None
    return {
        "power": power, "specificity": spec,
        "n_recombining": len(pos), "n_clonal": len(neg),
        "n_untestable": sum(p is None for _rec, p in results),
    }


def _cache_dir() -> Path:
    override = os.environ.get("TESSERA_BENCHMARK_DIR")
    return Path(override) if override else DEFAULT_CACHE


def main(argv: list[str]) -> int:
    alpha = _flag_float(argv, "--alpha", 0.05)
    window = int(_flag_float(argv, "--window", 100))
    cap = int(_flag_float(argv, "--max", 0))  # 0 = no cap

    cache = _cache_dir()
    alignments = sorted(p for p in cache.rglob("*.fasta") if parse_msa_name(p.name))
    if not alignments:
        print(f"[SKIP] no benchmark alignments in {cache}. Download the Jaya 2023 archives from "
              f"Dryad (doi:10.5061/dryad.d7wm37q6f) and extract the msa_*.fasta files there "
              f"(or set $TESSERA_BENCHMARK_DIR).")
        return 0
    if cap:
        alignments = alignments[:cap]

    print(f"Tessera PHI benchmark -- {len(alignments)} alignment(s) from {cache}")
    print(f"(power = P[PHI p < {alpha} | recombining]; specificity = P[PHI p >= {alpha} | clonal])")
    results: list[tuple[float, float | None]] = []
    per_rec: dict[float, list[float | None]] = {}
    for path in alignments:
        params = parse_msa_name(path.name)
        p = phi_pvalue_of(path, window=window)
        results.append((params["rec"], p))
        per_rec.setdefault(params["rec"], []).append(p)

    print("\n  rec-rate   n   power/spec   detected/testable")
    for rec in sorted(per_rec):
        ps = per_rec[rec]
        testable = [x for x in ps if x is not None]
        if rec > 0:
            hits = sum(x < alpha for x in testable)
            label = f"power {hits}/{len(testable)}" if testable else "power n/a"
        else:
            hits = sum(x >= alpha for x in testable)
            label = f"spec {hits}/{len(testable)}" if testable else "spec n/a"
        print(f"  {rec:<8} {len(ps):>4}   {label:<12} ({len(testable)} testable)")

    agg = power_specificity(results, alpha)
    print("\n" + "=" * 60)
    pw = f"{agg['power']:.2f}" if agg["power"] is not None else "n/a"
    sp = f"{agg['specificity']:.2f}" if agg["specificity"] is not None else "n/a"
    print(f"  power {pw} ({agg['n_recombining']} recombining) | "
          f"specificity {sp} ({agg['n_clonal']} clonal) | "
          f"{agg['n_untestable']} untestable")
    return 0


def _flag_float(argv: list[str], flag: str, default: float) -> float:
    if flag in argv:
        return float(argv[argv.index(flag) + 1])
    return default


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
