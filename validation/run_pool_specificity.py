#!/usr/bin/env python
"""Opt-in harness: does Tessera call recombination on a clonal genome from a real panel?

``run_specificity.py`` asks the same question of simulated data, where the ground truth
is exact but the sequences are simpler than real viral evolution. This harness asks it
of real genomes instead. The query is a tip taken from a cached Nextclade reference
tree whose label is a clade in its own right -- not one of the sentinels that mean "no
clade", not a composite of two other clades, and not a recombinant-named lineage -- so
it is a single sampled genome rather than a mosaic. **Any region reported on it is a
false positive**, give or take the real recombination that exists in natural
populations, which is why a call here is reported for inspection rather than treated as
proof of a defect.

Unlike every other panel harness this one runs with **no network at all**, which is the
point. Both proxies are pointed at a closed local port in the child environment, so an
attempt to reach NCBI is refused in milliseconds instead of stalling for minutes. That
matters because seeding from a cached pool does not by itself keep a run offline:

- ``--nextclade-dataset`` fetches the dataset index over the network before it consults
  the cache, so the pool is given as ``--candidate-pool <directory>`` instead;
- ``--max-rounds 1`` does **not** prevent the fill, because the search runs inside round
  one whenever a coverage gap remains. Only a panel that leaves no gap stays local.

A pool whose panel does not cover its own query therefore cannot be scored here, and is
reported as **not testable** -- an untested absence, neither a pass nor a failure. That
is a real outcome rather than a nuisance: it says the selection step could not build a
covering panel from that pool, which is worth knowing on its own.

The run is driven through the command line rather than the Python API, unlike the other
harnesses, so that the CLI's own seed-option validation is exercised along with the
recruitment path.

Needs ``mafft`` and ``skder`` on PATH (the ``recomfi-aln`` env) and a warm pool cache;
it downloads nothing. Append the env, never prepend, or the Python is swapped:

    export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
    python validation/run_pool_specificity.py
    python validation/run_pool_specificity.py --pools rubella,ebov   # a subset
    python validation/run_pool_specificity.py --keep /tmp/poolspec   # keep the outputs

Exit status is 1 when a testable pool reports a region, 0 otherwise; pools that are not
testable do not change it.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from pick_pool_query import choose_query

# Routing both proxies at a closed port makes "offline" an enforced property of the run
# rather than an assumption about what the code will do.
BLOCKED_PROXY = "http://127.0.0.1:9"
CONVERGED = "Converged: no coverage gaps remain"


def cached_pools(cache: Path) -> list[Path]:
    """Cached pool directories, newest-looking name last, sorted for determinism."""
    root = cache / "nextclade"
    if not root.is_dir():
        return []
    return sorted(d for d in root.iterdir() if d.is_dir() and any(d.glob("*.fasta")))


def pool_label(pool: Path) -> str:
    """The dataset slug without the cache-key suffix, for display."""
    return pool.name.rsplit("_", 1)[0]


def run_detect(query: Path, pool: Path, outdir: Path, timeout: int) -> tuple[int, str]:
    """Run ``tessera detect`` offline against ``pool``; return ``(rc, log text)``."""
    env = dict(os.environ)
    env.update(
        http_proxy=BLOCKED_PROXY, https_proxy=BLOCKED_PROXY,
        HTTP_PROXY=BLOCKED_PROXY, HTTPS_PROXY=BLOCKED_PROXY,
    )
    cmd = [
        sys.executable, "-m", "tessera", "detect",
        "-q", str(query), "-o", str(outdir),
        "--candidate-pool", str(pool), "--max-rounds", "1",
    ]
    try:
        proc = subprocess.run(
            cmd, env=env, capture_output=True, text=True, timeout=timeout, check=False,
        )
    except subprocess.TimeoutExpired:
        return 124, f"timed out after {timeout}s"
    text = proc.stdout + proc.stderr
    log = outdir / "tessera.log"
    if log.exists():
        text += log.read_text(errors="replace")
    return proc.returncode, text


def read_regions(outdir: Path) -> list[str] | None:
    """Region rows (without the header), or ``None`` when no table was written."""
    path = outdir / "recombination_regions.tsv"
    if not path.exists():
        return None
    rows = [ln for ln in path.read_text(errors="replace").splitlines() if ln.strip()]
    return rows[1:] if rows else []


def score(pool: Path, workdir: Path, min_clade_size: int, timeout: int) -> dict:
    """Run and classify one pool."""
    name = pool_label(pool)
    picked = choose_query(pool, min_clade_size)
    if picked is None:
        return {"pool": name, "verdict": "SKIP",
                "note": "no unrecombined tip whose clade is large enough"}

    query, clade, clade_size = picked
    # Keyed on the full cache directory name, not the display label: a dataset can be
    # cached more than once under different cache keys, and those share a label, so
    # keying on it would make the second run overwrite the first's outputs.
    outdir = workdir / pool.name
    shutil.rmtree(outdir, ignore_errors=True)
    outdir.mkdir(parents=True, exist_ok=True)
    rc, text = run_detect(query, pool, outdir, timeout)

    row = {
        "pool": name, "query": query.stem, "clade": clade, "clade_size": clade_size,
        "rc": rc, "refs": len(list((outdir / "collection").glob("*")))
        if (outdir / "collection").is_dir() else 0,
    }
    if CONVERGED not in text:
        row.update(verdict="NOT TESTABLE",
                   note="panel left a coverage gap, so the run needed the network")
        return row

    regions = read_regions(outdir)
    if regions is None:
        row.update(verdict="NOT TESTABLE", note=f"no region table written (rc={rc})")
    elif not regions:
        row.update(verdict="PASS", note="no region called", n_regions=0)
    else:
        row.update(verdict="CALL", n_regions=len(regions),
                   note=f"{len(regions)} region(s) on a clonal query",
                   regions=regions)
    return row


def parse_args(argv: list[str]) -> argparse.Namespace:
    default_cache = Path(os.environ.get("TESSERA_CACHE", Path.home() / ".cache/tessera"))
    parser = argparse.ArgumentParser(
        description="False-positive control on clonal tips from cached Nextclade pools.",
    )
    parser.add_argument("--cache", type=Path, default=default_cache,
                        help="Tessera cache directory (default $TESSERA_CACHE)")
    parser.add_argument("--pools", default="",
                        help="comma-separated substrings; only matching pools are run")
    parser.add_argument("--min-clade-size", type=int, default=5,
                        help="skip clades with fewer tips than this (default 5)")
    parser.add_argument("--timeout", type=int, default=1800,
                        help="per-run timeout in seconds (default 1800)")
    parser.add_argument("--keep", type=Path, default=None,
                        help="keep run outputs in this directory instead of a temp one")
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    pools = cached_pools(args.cache)
    if args.pools:
        wanted = [w.strip() for w in args.pools.split(",") if w.strip()]
        pools = [p for p in pools if any(w in p.name for w in wanted)]
    if not pools:
        print(f"no cached pools under {args.cache / 'nextclade'}", file=sys.stderr)
        return 2

    workdir = args.keep or Path(tempfile.mkdtemp(prefix="pool-specificity-"))
    workdir.mkdir(parents=True, exist_ok=True)

    print(f"Tessera pool specificity -- {len(pools)} cached pool(s), no network used.")
    print("The query is a clonal tree tip: a region reported on it is a false positive.\n")
    print(f"{'pool':42s} {'query':14s} {'clade':12s} {'refs':>4s}  verdict")
    print("-" * 92)

    rows = []
    for pool in pools:
        row = score(pool, workdir, args.min_clade_size, args.timeout)
        rows.append(row)
        print(f"{row['pool'][:42]:42s} {row.get('query', '-')[:14]:14s} "
              f"{str(row.get('clade', '-'))[:12]:12s} {row.get('refs', 0):>4d}  "
              f"{row['verdict']}: {row['note']}")
        for region in row.get("regions", []):
            print(f"{'':>62s}   {region[:90]}")

    testable = [r for r in rows if r["verdict"] in ("PASS", "CALL")]
    called = [r for r in rows if r["verdict"] == "CALL"]
    untestable = [r for r in rows if r["verdict"] == "NOT TESTABLE"]
    skipped = [r for r in rows if r["verdict"] == "SKIP"]

    print("-" * 92)
    print(f"{len(testable)} pool(s) testable: {len(testable) - len(called)} clean, "
          f"{len(called)} reporting a region.")
    print(f"{len(untestable)} not testable (no covering panel without the network), "
          f"{len(skipped)} skipped (no usable tip).")
    if called:
        total = sum(r["n_regions"] for r in called)
        print(f"False-positive regions: {total} across {len(called)} pool(s) -- "
              "inspect before treating as a defect (natural recombination exists).")
    if args.keep:
        print(f"Outputs kept in {workdir}")
    else:
        shutil.rmtree(workdir, ignore_errors=True)
    return 1 if called else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
