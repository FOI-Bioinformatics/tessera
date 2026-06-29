#!/usr/bin/env python
"""Run Tessera end-to-end on published recombination datasets and check the result.

This is an opt-in validation harness, separate from the pytest suite. It reads
``datasets.json``, builds an MSA and runs the recombination scan for each enabled
dataset, then asserts the outcome against a documented expectation (the known,
published recombination event). A dataset is SKIPPED (not failed) when the
aligner binary is missing or its sequences have not been fetched.

Usage:
    # ensure the aligner is on PATH, e.g.
    #   export PATH="$HOME/miniforge3/envs/recomfi-aln/bin:$PATH"
    python validation/fetch.py            # download the fetch-based datasets
    python validation/run_validation.py   # build + scan + check, print a table
"""

from __future__ import annotations

import json
import logging
import shutil
import sys
from pathlib import Path

from tessera.core.io import strip_sequence_extension
from tessera.msa.build import MsaParams, build_msa
from tessera.recomb.run import RecombParams, run_recomb

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DATA = HERE / "data"


def _load_datasets() -> list[dict]:
    return json.loads((HERE / "datasets.json").read_text())["datasets"]


def _resolve_inputs(ds: dict) -> tuple[Path | None, Path | None, str]:
    """Return ``(query_path, collection_dir, query_label)`` or ``(None, None, "")``.

    A ``fetch`` dataset reads from ``data/<name>/``; a bundled dataset reads the
    paths given in the manifest (relative to the repo root).
    """
    if "fetch" in ds:
        base = DATA / ds["name"]
        q_label = ds["fetch"]["query"]["label"]
        query = next((base / "query").glob(f"{q_label}.fasta*"), None)
        collection = base / "collection"
        if query is None or not collection.is_dir() or not any(collection.iterdir()):
            return None, None, q_label
        return query, collection, q_label

    query = REPO / ds["query"]
    collection = REPO / ds["collection"]
    if not query.exists() or not collection.is_dir():
        return None, None, ds.get("query_label", "")
    label = ds.get("query_label") or strip_sequence_extension(query.name)
    return query, collection, label


def _parse_regions(path: Path) -> list[dict]:
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"), strict=False)) for ln in lines[1:]]


def _parse_winners(path: Path) -> list[tuple[str, int]]:
    lines = path.read_text().splitlines()[1:]
    out = []
    for ln in lines:
        name, num = ln.split("\t")
        out.append((name, int(num)))
    return out


def _check(ds: dict, out_dir: Path) -> tuple[bool, list[str]]:
    """Assert the dataset's ``expected`` block against the recomb output."""
    exp = ds["expected"]
    msgs: list[str] = []
    ok = True

    winners = _parse_winners(out_dir / "window_winners.tsv")
    regions = _parse_regions(out_dir / "recombination_regions.tsv")
    top_winner = winners[0][0] if winners else None
    region_minors = {r["minor_parent"] for r in regions}
    region_majors = {r["major_parent"] for r in regions}
    present = {w for w, _ in winners} | region_minors | region_majors

    if "major_parent" in exp:
        good = top_winner == exp["major_parent"]
        ok &= good
        msgs.append(f"major={top_winner} (want {exp['major_parent']}) {'OK' if good else 'FAIL'}")

    if "region_minor_any" in exp:
        good = bool(region_minors & set(exp["region_minor_any"]))
        ok &= good
        msgs.append(
            f"region minors={sorted(region_minors)} (want any of "
            f"{exp['region_minor_any']}) {'OK' if good else 'FAIL'}"
        )

    if "parents_present" in exp:
        missing = set(exp["parents_present"]) - present
        good = not missing
        ok &= good
        msgs.append(f"parents present: missing={sorted(missing)} {'OK' if good else 'FAIL'}")

    if "min_regions" in exp:
        good = len(regions) >= exp["min_regions"]
        ok &= good
        msgs.append(f"#regions={len(regions)} (>= {exp['min_regions']}) {'OK' if good else 'FAIL'}")

    if "breakpoint" in exp:
        bp = exp["breakpoint"]
        lo_key = "msa_start" if bp["coord"] == "msa" else "query_start"
        hi_key = "msa_end" if bp["coord"] == "msa" else "query_end"
        overlap = any(
            int(r[lo_key]) <= bp["hi"] and int(r[hi_key]) >= bp["lo"] for r in regions
        )
        ok &= overlap
        msgs.append(
            f"breakpoint {bp['coord']} [{bp['lo']},{bp['hi']}] overlapped: "
            f"{'OK' if overlap else 'FAIL'}"
        )

    return ok, msgs


def _run_one(ds: dict, logger: logging.Logger) -> tuple[str, list[str]]:
    query, collection, label = _resolve_inputs(ds)
    if query is None:
        return "SKIP", ["sequences not present (run validation/fetch.py)"]

    work = DATA / ds["name"] / "_run"
    work.mkdir(parents=True, exist_ok=True)
    msa = work / "msa.fasta"
    recomb_out = work / "recomb_out"

    rec = ds.get("recomb", {})
    build_msa(
        MsaParams(
            query=query, collection=collection, output=msa,
            aligner=ds.get("aligner", "sibeliaz"),
            reference=ds.get("reference"),
            threads=ds.get("threads", 4),
        ),
        logger,
    )
    run_recomb(
        RecombParams(
            msa=msa, output=recomb_out, query=label,
            window_size=rec.get("window_size", 1000),
            window_step=rec.get("window_step", 100),
        ),
        logger,
    )
    ok, msgs = _check(ds, recomb_out)
    return ("PASS" if ok else "FAIL"), msgs


def main() -> int:
    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    logger = logging.getLogger("tessera.validation")

    datasets = [d for d in _load_datasets() if d.get("enabled", True)]
    have_aligner = shutil.which("sibeliaz") is not None

    print(f"\nTessera validation -- {len(datasets)} enabled dataset(s)\n" + "=" * 60)
    results: list[tuple[str, str]] = []
    for ds in datasets:
        if not have_aligner:
            status, msgs = "SKIP", ["sibeliaz not on PATH"]
        else:
            try:
                status, msgs = _run_one(ds, logger)
            except Exception as exc:  # noqa: BLE001 - report, don't abort the batch
                status, msgs = "ERROR", [f"{type(exc).__name__}: {exc}"]
        results.append((ds["name"], status))
        print(f"\n[{status}] {ds['name']} -- {ds.get('description', '')}")
        for m in msgs:
            print(f"    {m}")

    print("\n" + "=" * 60)
    for name, status in results:
        print(f"  {status:6} {name}")
    failed = [n for n, s in results if s in ("FAIL", "ERROR")]
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
