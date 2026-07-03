#!/usr/bin/env python
"""Reproduce the RecombinHunt (Bianchi et al., 2024) noise-robustness benchmark design.

RecombinHunt ("Data-driven recombination detection in viral genomes", *Nat. Commun.* 2024,
s41467-024-47464-5) measured detection sensitivity as increasing numbers of *non-characteristic*
mutations were injected into simulated recombinants (SARS-CoV-2 BA.2 x AY.45, one/two breakpoints,
0-30 injected mutations). Their archived sequences are GISAID-gated, so this reproduces the
**design** on Nextclade lineage genomes: splice a one-breakpoint recombinant from the two
most-divergent clade consensuses of a dataset, inject a grid of random point mutations, run
Tessera's recombination scan against the two parents, and report the fraction of replicates in
which a called region overlaps the true breakpoint span (breakpoint-recovery **sensitivity**) at
each noise level.

Needs skani and an aligner; opt-in, not part of CI. Self-contained -- sequences are generated.

    export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
    python validation/run_recombinhunt_benchmark.py                    # default dataset
    python validation/run_recombinhunt_benchmark.py --dataset nextstrain/flu/h3n2/ha/EPI1857216
"""

from __future__ import annotations

import logging
import random
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

NOISE_GRID = (0, 3, 5, 10, 15, 20, 30)  # RecombinHunt's injected-mutation levels
# RecombinHunt used SARS-CoV-2 BA.2 x AY.45 (GISAID-gated); the design is virus-agnostic, so the
# default is a divergent, cleanly-localising Nextclade pair (override with --dataset). A
# low-divergence dataset (e.g. the XBB Omicron sublineages) degrades faster and is less informative.
DEFAULT_DATASET = "nextstrain/flu/h3n2/ha/EPI1857216"
_BASES = "ACGT"


def inject_noise(seq: str, k: int, rng: random.Random) -> str:
    """Replace ``k`` random positions of ``seq`` with a different base (point-mutation noise)."""
    if k <= 0:
        return seq
    chars = list(seq)
    for pos in rng.sample(range(len(chars)), min(k, len(chars))):
        chars[pos] = rng.choice([b for b in _BASES if b != chars[pos].upper()])
    return "".join(chars)


def sensitivity_by_noise(results: list[tuple[int, bool]]) -> dict[int, tuple[int, int]]:
    """``[(noise, recovered), ...]`` -> ``{noise: (n_recovered, n_total)}``."""
    out: dict[int, tuple[int, int]] = {}
    for noise, ok in results:
        rec, tot = out.get(noise, (0, 0))
        out[noise] = (rec + int(bool(ok)), tot + 1)
    return out


def _run(dataset: str, reps: int, logger: logging.Logger) -> int:
    import run_hybrids as rh
    from run_reassort_scan import region_overlaps_span

    from tessera.core.io import write_fasta_record
    from tessera.msa.build import MsaParams, build_msa
    from tessera.recomb.regions import DEFAULT_METHODS, parse_methods
    from tessera.recomb.run import RecombParams, run_recomb

    _g, ref, tips = rh._load_species(dataset, None, logger)
    try:
        clade_a, clade_b, src_a, src_b = rh.pick_parents(tips, ref, [], logger)
    except rh.CaseSkipped as exc:
        print(f"[SKIP] {dataset}: {exc}")
        return 0
    a_muts, b_muts = tips[src_a][1], tips[src_b][1]
    parent_a = rh.reconstruct_gapped(ref, a_muts).replace("-", "").upper()
    parent_b = rh.reconstruct_gapped(ref, b_muts).replace("-", "").upper()
    logger.info("parents: %s x %s", clade_a, clade_b)

    rng = random.Random(1)
    results: list[tuple[int, bool]] = []
    with tempfile.TemporaryDirectory() as tmp:
        collection = Path(tmp) / "collection"
        collection.mkdir()
        for label, seq in (("A", parent_a), ("B", parent_b)):
            with open(collection / f"{label}.fasta", "w") as fo:
                write_fasta_record(fo, label, seq)
        for rep in range(reps):
            clean, q0, q1 = rh.make_hybrid(ref, a_muts, b_muts)
            for k in NOISE_GRID:
                query_seq = inject_noise(clean, k, rng)
                out = Path(tmp) / f"r{rep}_k{k}"
                out.mkdir()
                qf = out / "query.fasta"
                with open(qf, "w") as fo:
                    write_fasta_record(fo, "query", query_seq)
                window, step, _ = rh.window_params(len(query_seq))
                msa = out / "msa.fasta"
                build_msa(MsaParams(query=qf, collection=collection, output=msa, aligner="mafft"),
                          logger)
                run_recomb(RecombParams(msa=msa, output=out, query="query",
                                        window_size=window, window_step=step,
                                        methods=parse_methods(",".join(DEFAULT_METHODS))), logger)
                regions = rh.parse_regions(out / "recombination_regions.tsv")
                results.append((k, region_overlaps_span(regions, q0, q1)))

    print(f"\nRecombinHunt-style noise-robustness -- {clade_a} x {clade_b} ({dataset})")
    print(f"(breakpoint-recovery sensitivity vs injected mutations; {reps} reps/level)\n")
    print("  injected-noise   sensitivity")
    for k, (rec, tot) in sorted(sensitivity_by_noise(results).items()):
        bar = "#" * rec + "." * (tot - rec)
        print(f"  {k:>13}   {rec}/{tot}  {bar}")
    return 0


def main(argv: list[str]) -> int:
    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    logger = logging.getLogger("recombinhunt-benchmark")
    dataset = _opt(argv, "--dataset", DEFAULT_DATASET)
    reps = int(_opt(argv, "--reps", "5"))
    if "--skani-missing-ok" not in argv:
        from tessera.discover.panel import skani_available
        if not skani_available():
            print("[SKIP] needs skani (and an aligner) on PATH.")
            return 0
    return _run(dataset, reps, logger)


def _opt(argv: list[str], name: str, default: str) -> str:
    return argv[argv.index(name) + 1] if name in argv else default


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
