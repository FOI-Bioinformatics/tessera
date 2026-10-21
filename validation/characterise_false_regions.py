#!/usr/bin/env python
"""Describe each false region the callers report on the specificity harness's clonal data.

``run_specificity.py`` counts false regions; this lists them, one row per region and
per caller, with the quantities needed to see *why* each was called: its length, how
many sites distinguish the two parents and which way they lean, its p- and q-value,
whether donor and major belong to the same simulated clade, and whether lineage
clustering pooled anything. It exists to inform a decision about the ``--min-methods``
default and the HMM's significance gate, not to gate a change.

Needs no aligner, network or downloaded data.

    python validation/characterise_false_regions.py             # 10 replicates
    python validation/characterise_false_regions.py --reps 3
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
import tempfile
from pathlib import Path

from tessera.recomb.analyze import analyze
from tessera.recomb.clusters import all_singletons, cluster_references
from tessera.recomb.regions import DEFAULT_METHODS, RegionParams, call_regions
from tessera.recomb.similarity import compute_similarity, discordant_counts

WINDOW, STEP = 1000, 100  # the CLI defaults, as run_specificity.py uses
COLUMNS = (
    "scenario", "seed", "caller", "minor", "major", "msa_start", "msa_end", "length",
    "favor_minor", "favor_major", "pvalue", "qvalue", "sim_minor", "sim_major",
    "same_clade", "panel_clustered",
)


def _load_specificity():
    path = Path(__file__).resolve().parent / "run_specificity.py"
    spec = importlib.util.spec_from_file_location("run_specificity", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["run_specificity"] = module
    spec.loader.exec_module(module)
    return module


def false_region_rows(scenario: str, seed: int, seqs: dict[str, str], query: str) -> list[dict]:
    """One row per region any default caller reports on a clonal alignment."""
    rows: list[dict] = []
    with tempfile.TemporaryDirectory() as td:
        msa = Path(td) / "aln.fasta"
        msa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
        result = compute_similarity(str(msa), query, WINDOW, STEP)
    analysis = analyze(result)
    clustered = not all_singletons(
        cluster_references(result, WINDOW, RegionParams.with_defaults(WINDOW))
    )
    for method in DEFAULT_METHODS:
        params = RegionParams.with_defaults(WINDOW, method=method)
        regions, _major, _siblings = call_regions(result, analysis, WINDOW, params)
        for r in regions:
            favor_minor, favor_major = discordant_counts(
                result.rows, query, r.major_parent, r.minor_parent, r.msa_start, r.msa_end
            )
            rows.append({
                "scenario": scenario, "seed": seed, "caller": method,
                "minor": r.minor_parent, "major": r.major_parent,
                "msa_start": r.msa_start, "msa_end": r.msa_end,
                "length": r.msa_end - r.msa_start,
                "favor_minor": favor_minor, "favor_major": favor_major,
                "pvalue": "" if r.pvalue is None else f"{r.pvalue:.3g}",
                "qvalue": "" if r.qvalue is None else f"{r.qvalue:.3g}",
                "sim_minor": r.mean_sim_minor, "sim_major": r.mean_sim_major,
                # tips are named <clade letter><index>, e.g. A2
                "same_clade": r.minor_parent[0] == r.major_parent[0],
                "panel_clustered": clustered,
            })
    return rows


def parse_args(argv: list[str]) -> argparse.Namespace:
    """The harness options; ``--help`` exits before any simulation."""
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        prog="characterise_false_regions.py",
        description="List each false region the callers report on the specificity "
        "harness's clonal data. Needs no aligner, network or data.",
    )
    parser.add_argument("--reps", type=int, default=10,
                        help="replicates per scenario (default 10)")
    args = parser.parse_args(argv)
    if args.reps < 1:
        parser.error("--reps must be at least 1")
    return args


def main(argv: list[str]) -> int:
    reps = parse_args(argv).reps
    rs = _load_specificity()
    print("\t".join(COLUMNS))
    total = 0
    for scenario in rs.SCENARIOS:
        for rep in range(reps):
            seed = 1000 + rep
            for row in false_region_rows(scenario, seed,
                                         rs.simulate_clonal(scenario, seed=seed), rs.QUERY):
                print("\t".join(str(row[c]) for c in COLUMNS))
                total += 1
    print(f"# {total} false region(s) over {reps * len(rs.SCENARIOS)} clonal alignment(s); "
          "single-caller rows, before the ensemble merge and any --min-methods gate",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
