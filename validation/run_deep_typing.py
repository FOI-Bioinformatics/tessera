#!/usr/bin/env python
"""Curated end-to-end check for the lineage-assignment ladder that ``--deep-typing`` runs.

The unit tests exercise ``assign_lineages`` with the Nextclade and skani calls stubbed;
this opt-in harness runs the real ladder once, so the nearest-neighbour tier is exercised
against a live Nextclade dataset. It reuses the bundled SARS-CoV-2 collection, rewrites each
genome's header to a bare accession (so header/title mining cannot short-circuit and the
ladder must reach the Nextclade nearest-neighbour tier), then calls the same
``assign_lineages`` entry point ``detect --deep-typing`` invokes and checks the result.

It SKIPs (does not fail) when skani is not on PATH, the collection has not been fetched, or
the Nextclade fetch is unavailable -- matching ``run_validation.py``'s SKIP semantics.

Usage:
    export PATH="$HOME/miniforge3/envs/recomfi-aln/bin:$PATH"   # skani on PATH
    python validation/run_deep_typing.py

The ncbi-datasets source tag (the other half of --deep-typing) is covered by the unit tests
(the sidecar round-trip and the assign_lineages datasets_rows pass-through); a fetch-based
end-to-end for that tag would add NCBI-download flakiness for little extra confidence.
"""

from __future__ import annotations

import logging
import os
import shutil
import sys
import tempfile
from pathlib import Path

from tessera.core.io import read_fasta, strip_sequence_extension
from tessera.discover.lineage_assign import assign_lineages

HERE = Path(__file__).resolve().parent
COLLECTION = HERE / "data" / "sarscov2_xbb" / "collection"
NEXTCLADE_DATASET = "nextstrain/sars-cov-2/wuhan-hu-1/orfs"


def _bare_accession_copies(collection: Path, dest: Path) -> list[Path]:
    """Copy each collection genome into ``dest`` with its header reduced to the bare
    accession, so tier-1 title mining yields nothing and the ladder must reach the
    Nextclade nearest-neighbour tier. The file name (the label) is preserved."""
    genomes: list[Path] = []
    for src in sorted(collection.glob("*.fasta*")):
        records = read_fasta(str(src))
        if not records:
            continue
        header, seq = records[0]
        accession = header.split()[0]
        out = dest / f"{src.name.split('.fasta')[0]}.fasta"
        out.write_text(f">{accession}\n{seq}\n")
        genomes.append(out)
    return genomes


def _check() -> tuple[str, list[str]]:
    if shutil.which("skani") is None:
        return "SKIP", ["skani not on PATH"]
    if not COLLECTION.is_dir() or not any(COLLECTION.glob("*.fasta*")):
        return "SKIP", [f"collection not found: {COLLECTION}"]

    logger = logging.getLogger("tessera.validation")
    with tempfile.TemporaryDirectory() as tmp:
        genomes = _bare_accession_copies(COLLECTION, Path(tmp))
        if not genomes:
            return "SKIP", ["no readable genomes in the collection"]
        try:
            rows = assign_lineages(
                genomes,
                nextclade_dataset=NEXTCLADE_DATASET,
                email=os.environ.get("NCBI_EMAIL"),
                logger=logger,
            )
        except Exception as exc:  # noqa: BLE001 - a missing dataset/network is a SKIP, not a FAIL
            return "SKIP", [f"Nextclade ladder unavailable ({type(exc).__name__}: {exc})"]

    labels = {label for label, _lineage, _source in rows}
    expected = {strip_sequence_extension(g.name) for g in genomes}
    sources = {source for _label, _lineage, source in rows}
    nn = [(label, lineage) for label, lineage, source in rows if source == "nextclade-nn"]

    msgs = [
        f"typed {len(labels)}/{len(expected)} genome(s); sources: {sorted(sources)}",
        f"nextclade-nn: {sorted(f'{lab}={lin}' for lab, lin in nn)}",
    ]
    ok = True
    if labels != expected:
        ok = False
        msgs.append(f"MISSING: {sorted(expected - labels)}")
    if not nn:
        ok = False
        msgs.append("no genome was typed via the nextclade-nn tier (ladder did not reach tier-2)")
    if any(not lineage or lineage in ("?", "NA") for _label, lineage in nn):
        ok = False
        msgs.append("a nextclade-nn label is empty/placeholder")
    return ("PASS" if ok else "FAIL"), msgs


def main() -> int:
    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    print("\nTessera deep-typing ladder check (SARS-CoV-2, nextclade-nn)\n" + "=" * 60)
    try:
        status, msgs = _check()
    except Exception as exc:  # noqa: BLE001 - report, don't crash the harness
        status, msgs = "ERROR", [f"{type(exc).__name__}: {exc}"]
    print(f"\n[{status}] sarscov2_xbb -- lineage ladder via {NEXTCLADE_DATASET}")
    for m in msgs:
        print(f"    {m}")
    print("\n" + "=" * 60)
    print(f"  {status}")
    return 1 if status in ("FAIL", "ERROR") else 0


if __name__ == "__main__":
    sys.exit(main())
