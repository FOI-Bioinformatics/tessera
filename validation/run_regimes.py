#!/usr/bin/env python
"""Opt-in harness: two regimes the other harnesses do not sample.

``run_specificity.py`` simulates clades ~16 % apart, so it always runs under base-pair
windowing, and its positive control has parents that differ everywhere. Two situations
found by the post-1.2.0 audit fall outside it:

    near_identical   a panel whose references are ~0.2 % apart (the mpox / VZV /
                     within-lineage SARS-CoV-2 regime, analysed under informative-site
                     windowing). Run twice: with a donor tract spliced into the query
                     (is it found, and does the HMM contribute?) and clonal (what is
                     reported when there is nothing to find?).
    sibling_tract    the panel holds a sibling recombinant whose donor tract runs 600 bp
                     further than the query's. Against that sibling the query differs
                     only over those 600 bp; downstream of them the donor and the
                     backbone are identical, so nothing marks where the region ends.

Like ``run_specificity.py`` this needs no aligner, no network and no downloaded data:
the simulated sequences are already aligned. It is opt-in because a full run is minutes
of wall clock; the simulation and scoring logic is unit-tested in CI
(``tests/unit/test_regime_scoring.py``).

    python validation/run_regimes.py                  # 10 replicates, --min-methods 1
    python validation/run_regimes.py --reps 3         # quick look
    python validation/run_regimes.py --min-methods 2  # with the agreement gate
    python validation/run_regimes.py --no-cluster-lineages

``--min-methods`` defaults to 1 here because that is the CLI default; the harness
reports what a user gets.

Caveat: a star tree under JC69 is simpler than real viral evolution, and a handful of
replicates carries real sampling error. These numbers show whether a failure mode
exists and roughly how large it is, not its magnitude on real panels.
"""

from __future__ import annotations

import logging
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np

from tessera.recomb.run import RecombParams, run_recomb

QUERY = "QUERY"
_BASES = np.frombuffer(b"ACGT", dtype=np.uint8)

# near_identical: five references on a star tree, each this far from the root.
NEAR_LENGTH = 120_000
NEAR_DISTANCE = 0.001
NEAR_REFS = ("R0", "R1", "R2", "R3", "R4")
NEAR_TRACT = (40_000, 48_000)  # donor tract (from R1) in the R0-derived query

# sibling_tract: parents 10 % apart; the sibling's tract runs 600 bp past the query's.
SIB_LENGTH = 30_000
SIB_DISTANCE = 0.05
SIB_QUERY_TRACT = (12_000, 18_000)
SIB_EXTRA = 600


# --- simulation ------------------------------------------------------------------

def _evolve(seq: np.ndarray, distance: float, rng) -> np.ndarray:
    """JC69: mutate each site with probability 3/4(1 - exp(-4/3 * d))."""
    p = 0.75 * (1.0 - np.exp(-4.0 / 3.0 * distance))
    hit = rng.random(seq.size) < p
    out = seq.copy()
    n = int(hit.sum())
    if n:  # a mutation always changes the base
        out[hit] = (seq[hit] + rng.integers(1, 4, size=n)) % 4
    return out


def _as_text(seqs: dict[str, np.ndarray]) -> dict[str, str]:
    return {k: _BASES[v].tobytes().decode("ascii") for k, v in seqs.items()}


def simulate_near_identical(
    seed: int, *, recombinant: bool, length: int = NEAR_LENGTH,
    distance: float = NEAR_DISTANCE, tract: tuple[int, int] = NEAR_TRACT,
) -> tuple[dict[str, str], tuple[int, int] | None]:
    """A near-identical panel and a query descended from ``R0``.

    With ``recombinant`` the query carries ``R1`` over ``tract``; otherwise it is a
    plain descendant of ``R0`` and nothing is recombined. Returns ``(sequences, tract)``
    with ``tract`` ``None`` for the clonal case.
    """
    rng = np.random.default_rng(seed)
    root = rng.integers(0, 4, size=length)
    seqs = {label: _evolve(root, distance, rng) for label in NEAR_REFS}
    query = _evolve(seqs["R0"], distance / 10.0, rng)
    if recombinant:
        lo, hi = tract
        query[lo:hi] = seqs["R1"][lo:hi]
    seqs[QUERY] = query
    return _as_text(seqs), (tract if recombinant else None)


def simulate_sibling_tract(
    seed: int, *, length: int = SIB_LENGTH, distance: float = SIB_DISTANCE,
    tract: tuple[int, int] = SIB_QUERY_TRACT, extra: int = SIB_EXTRA,
) -> tuple[dict[str, str], tuple[int, int]]:
    """A query and a sibling recombinant whose donor tract runs ``extra`` bp further.

    Both are a ``parent_A`` backbone carrying ``parent_B``; the query over ``tract``, the
    sibling over ``tract`` extended by ``extra``. Returns ``(sequences, span)`` where
    ``span`` is the only stretch on which the query and the sibling differ -- there the
    query matches ``parent_A``. Everywhere downstream ``parent_A`` and the sibling are
    identical.
    """
    rng = np.random.default_rng(seed)
    root = rng.integers(0, 4, size=length)
    parent_a, parent_b, other = (_evolve(root, distance, rng) for _ in range(3))
    lo, hi = tract
    query = parent_a.copy()
    query[lo:hi] = parent_b[lo:hi]
    sibling = parent_a.copy()
    sibling[lo:hi + extra] = parent_b[lo:hi + extra]
    seqs = {QUERY: query, "sibling": sibling, "parent_A": parent_a,
            "parent_B": parent_b, "other": other}
    return _as_text(seqs), (hi, hi + extra)


# --- scoring ---------------------------------------------------------------------

def _present(rows: list[dict]) -> list[dict]:
    """Rows that claim a donor; a ``donor_absent`` row is a coverage statement."""
    return [r for r in rows if r.get("donor_absent") != "yes"]


def _methods(row: dict) -> list[str]:
    return [m.strip() for m in (row.get("methods") or "").split(",") if m.strip()]


def score_tract(rows: list[dict], tract: tuple[int, int], *, donor: str) -> dict:
    """Detection of a known tract, whether the HMM was among the callers, donor
    attribution and breakpoint error."""
    lo, hi = tract
    overlapping = [
        r for r in _present(rows)
        if int(r["query_start"]) < hi and int(r["query_end"]) > lo
    ]
    if not overlapping:
        return {"detected": False, "hmm": False, "donor_ok": False,
                "breakpoint_error": None}
    best = max(
        overlapping,
        key=lambda r: min(int(r["query_end"]), hi) - max(int(r["query_start"]), lo),
    )
    return {
        "detected": True,
        "hmm": any("hmm" in _methods(r) for r in overlapping),
        "donor_ok": best["minor_parent"] == donor,
        "breakpoint_error": (
            abs(int(best["query_start"]) - lo) + abs(int(best["query_end"]) - hi)
        ) // 2,
    }


def score_clonal(rows: list[dict]) -> tuple[int, dict[str, int]]:
    """``(false regions, per-caller counts)`` for a run that must report nothing."""
    present = _present(rows)
    per_caller: Counter[str] = Counter()
    for row in present:
        per_caller.update(_methods(row))
    return len(present), dict(per_caller)


def score_span(rows: list[dict], span: tuple[int, int]) -> dict:
    """How much longer than the true difference the reported region is.

    ``reported_bp`` is the length of the longest reported region overlapping ``span``
    (``None`` when none does); ``excess_bp`` is that length minus the span's.
    """
    lo, hi = span
    overlapping = [
        r for r in _present(rows)
        if int(r["query_start"]) < hi and int(r["query_end"]) > lo
    ]
    if not overlapping:
        return {"detected": False, "reported_bp": None, "excess_bp": None}
    reported = max(int(r["query_end"]) - int(r["query_start"]) for r in overlapping)
    return {"detected": True, "reported_bp": reported,
            "excess_bp": max(0, reported - (hi - lo))}


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


def scan(
    seqs: dict[str, str], logger: logging.Logger, *, min_methods: int = 1,
    cluster_lineages: bool = True,
) -> list[dict]:
    """Run the shipped pipeline on one simulated alignment; return its region rows."""
    with tempfile.TemporaryDirectory() as td:
        msa = Path(td) / "aln.fasta"
        _write_fasta(msa, seqs)
        out = Path(td) / "out"
        run_recomb(
            RecombParams(msa=msa, output=out, query=QUERY, plot_format="png",
                         min_methods=min_methods, cluster_lineages=cluster_lineages),
            logger,
        )
        return _read_regions(out / "recombination_regions.tsv")


def _median(values: list[int]) -> str:
    return f"{int(np.median(values))} bp" if values else "n/a"


def main(argv: list[str]) -> int:
    reps, min_methods = 10, 1
    if "--reps" in argv:
        reps = int(argv[argv.index("--reps") + 1])
    if "--min-methods" in argv:
        min_methods = int(argv[argv.index("--min-methods") + 1])
    cluster = "--no-cluster-lineages" not in argv

    logger = logging.getLogger("tessera.regimes")
    logger.addHandler(logging.NullHandler())
    logger.propagate = False  # the scan is chatty; the table below is the output

    def run(seqs: dict[str, str]) -> list[dict]:
        return scan(seqs, logger, min_methods=min_methods, cluster_lineages=cluster)

    print(f"Tessera regime harness -- {reps} replicate(s), --min-methods {min_methods}, "
          f"lineage clustering {'on' if cluster else 'off'}\n")

    detected = hmm = donor_ok = 0
    errors: list[int] = []
    for rep in range(reps):
        seqs, tract = simulate_near_identical(3000 + rep, recombinant=True)
        assert tract is not None
        got = score_tract(run(seqs), tract, donor="R1")
        detected += got["detected"]
        hmm += got["hmm"]
        donor_ok += got["donor_ok"]
        if got["breakpoint_error"] is not None:
            errors.append(got["breakpoint_error"])
    print(f"near_identical, tract {NEAR_TRACT[0]}-{NEAR_TRACT[1]} from R1:")
    print(f"  detected {detected}/{reps} | called by the HMM {hmm}/{reps} | "
          f"correct donor {donor_ok}/{reps} | median breakpoint error {_median(errors)}")

    bad = regions = 0
    callers: Counter[str] = Counter()
    for rep in range(reps):
        seqs, _ = simulate_near_identical(4000 + rep, recombinant=False)
        n_false, per_caller = score_clonal(run(seqs))
        bad += bool(n_false)
        regions += n_false
        callers.update(per_caller)
    print("near_identical, clonal (every region is a false positive):")
    print(f"  runs with a false region {bad}/{reps} | false regions {regions} | "
          + (", ".join(f"{k}={v}" for k, v in sorted(callers.items())) or "--"))

    found = 0
    reported: list[int] = []
    excess: list[int] = []
    for rep in range(reps):
        seqs, span = simulate_sibling_tract(5000 + rep)
        got = score_span(run(seqs), span)
        found += got["detected"]
        if got["reported_bp"] is not None:
            reported.append(got["reported_bp"])
            excess.append(got["excess_bp"])
    print(f"sibling_tract, true difference {SIB_EXTRA} bp:")
    print(f"  detected {found}/{reps} | median reported length {_median(reported)} | "
          f"median excess {_median(excess)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
