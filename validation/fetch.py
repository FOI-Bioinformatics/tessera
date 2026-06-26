#!/usr/bin/env python
"""Download the fetch-based validation datasets from NCBI.

Reads ``datasets.json`` and, for every dataset that carries a ``fetch`` block,
downloads each accession to ``data/<name>/query/<label>.fasta`` or
``data/<name>/collection/<label>.fasta`` using NCBI Entrez ``efetch`` (part of
``entrez-direct``; ``conda install -c bioconda entrez-direct``).

Accessions whose value starts with ``TODO`` are skipped with a warning -- fill
them in ``datasets.json`` (see the per-dataset ``provenance`` note) and re-run.

Usage:
    python validation/fetch.py                 # all datasets with a fetch block
    python validation/fetch.py sarscov2_xbb    # only the named dataset(s)
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"


def _load_datasets() -> list[dict]:
    return json.loads((HERE / "datasets.json").read_text())["datasets"]


def _efetch(accession: str, dest: Path) -> bool:
    """Fetch one nucleotide accession to ``dest`` as FASTA. Returns success."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        out = subprocess.run(
            ["efetch", "-db", "nuccore", "-id", accession, "-format", "fasta"],
            capture_output=True, text=True, timeout=300, check=True,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        print(f"    ! efetch failed for {accession}: {exc}")
        return False
    if not out.stdout.startswith(">"):
        print(f"    ! {accession}: no FASTA returned")
        return False
    dest.write_text(out.stdout)
    return True


def _fetch_dataset(ds: dict) -> None:
    base = DATA / ds["name"]
    print(f"\n{ds['name']}: {ds.get('description', '')}")

    items: list[tuple[str, dict]] = [("query", ds["fetch"]["query"])]
    items += [("collection", c) for c in ds["fetch"]["collection"]]

    for kind, rec in items:
        acc, label = rec["accession"], rec["label"]
        if acc.startswith("TODO"):
            print(f"    - skip {label}: accession not set ({acc})")
            continue
        dest = base / kind / f"{label}.fasta"
        if dest.exists():
            print(f"    = {label}: already present")
            continue
        if _efetch(acc, dest):
            print(f"    + {label} <- {acc}")


def main(argv: list[str]) -> int:
    if shutil.which("efetch") is None:
        print("efetch not found. Install with: conda install -c bioconda entrez-direct")
        return 1

    wanted = set(argv)
    datasets = [d for d in _load_datasets() if "fetch" in d]
    if wanted:
        datasets = [d for d in datasets if d["name"] in wanted]
    if not datasets:
        print("No matching fetch-based datasets.")
        return 0

    for ds in datasets:
        _fetch_dataset(ds)

    print("\nDone. Fill any TODO accessions in datasets.json, set enabled=true, "
          "then run validation/run_validation.py.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
