#!/usr/bin/env python
"""Measure Option B (issue #90): when is a donor sub-clade not identifiable from the panel?

The scan names one donor. When the donor's true lineage is missing, the named reference
can be a neighbour the query does not actually prefer. ``flag_subclade_identifiability``
(src/tessera/recomb/coverage.py) annotates a region as sub-clade-unresolved when the query
does not significantly prefer the named donor's lineage over some *other* present lineage.
This measures that flag's reach on the hybrid harness, by running the shipped function over
the per-case artifacts an earlier ``run_hybrids.py`` left on disk -- so the number in the
changelog is reproducible and measures the real code, not a reimplementation.

It reports, per case, whether the primary region (most callers agreeing, longest breaking
the tie) is flagged, split by whether the donor's lineage is present (the ordinary hybrid
cases) or absent (``donorabsent_rsv`` and the ``absent_*`` probe runs from
``probe_donor_absent.py``). A useful flag catches donor-absent cases without caveating
donor-present ones.

Needs no aligner, network or downloaded data: it reads ``panel.msa.fasta``,
``lineages.tsv`` and ``recombination_regions.tsv`` under ``validation/data/hybrids`` and
``validation/data/donor_absent_probe``. Run ``run_hybrids.py`` (and, for the absent side,
``probe_donor_absent.py``) first if those directories are empty.

    python validation/measure_subclade_identifiability.py
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from tessera.recomb.coverage import flag_subclade_identifiability
from tessera.recomb.regions import Region
from tessera.recomb.similarity import _read_alignment

HERE = Path(__file__).resolve().parent
ROOTS = (HERE / "data" / "hybrids", HERE / "data" / "donor_absent_probe")
ALPHA = 0.05
# probe_donor_absent.py rebuilds the harness's own donor-absent case; counting both would
# double it, and it is kept only for the side-by-side check in measure_pendant.py.
SKIP = {"absent_rsv_a"}


def _primary(regions: list[dict[str, str]]) -> dict[str, str]:
    return max(regions, key=lambda r: (len(r["methods"].split(",")), int(r["length_bp"])))


def _region_from_row(row: dict[str, str]) -> Region:
    # Only the fields flag_subclade_identifiability reads need to be real.
    return Region(
        minor_parent=row["minor_parent"], major_parent=row["major_parent"],
        msa_start=int(row["msa_start"]), msa_end=int(row["msa_end"]),
        query_start=int(row["query_start"]), query_end=int(row["query_end"]),
        n_windows=int(row["n_windows"]), mean_sim_minor=0.0, mean_sim_major=0.0, margin=0.0,
    )


def _case_dirs(roots: tuple[Path, ...]):
    for root in roots:
        if root.is_dir():
            yield from sorted(p for p in root.iterdir() if p.is_dir())


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, nargs="+", default=list(ROOTS))
    args = ap.parse_args(argv)

    present: list[tuple[str, bool | None, tuple[str, ...]]] = []
    absent: list[tuple[str, bool | None, tuple[str, ...]]] = []
    all_regions: list[tuple[str, int, bool | None]] = []  # (class, n callers, identifiable)
    for case in _case_dirs(tuple(args.data)):
        if case.name in SKIP:
            continue
        reg_tsv = case / "recombination_regions.tsv"
        lin_tsv = case / "lineages.tsv"
        msa = case / "panel.msa.fasta"
        if not (reg_tsv.exists() and lin_tsv.exists() and msa.exists()):
            continue
        # Negatives have no donor to attribute; skip them and the un-scanned oropouche.
        if case.name.startswith(("neg_", "negwithin_")) or case.name == "oropouche":
            continue
        lineage: dict[str, str] = {}
        for row in csv.DictReader(open(lin_tsv, encoding="utf-8"), delimiter="\t"):
            if row["source"] != "query":
                lineage[row["label"]] = row["genotype"]
        regions = [r for r in csv.DictReader(open(reg_tsv, encoding="utf-8"), delimiter="\t")
                   if r.get("donor_absent") != "yes"]
        if not regions:
            continue
        region = _region_from_row(_primary(regions))
        rows = _read_alignment(str(msa))
        flag_subclade_identifiability(
            [region], rows, "hybrid", lineage.get, alpha=ALPHA
        )
        is_absent = case.name.startswith("absent_") or case.name == "donorabsent_rsv"
        record = (case.name, region.subclade_identifiable, region.tied_lineages)
        (absent if is_absent else present).append(record)
        # Every called region, as a real run reports them: the primary region above is one
        # per case, but a run usually reports more, and the flag behaves differently on them.
        every = [_region_from_row(r) for r in regions]
        flag_subclade_identifiability(every, rows, "hybrid", lineage.get, alpha=ALPHA)
        for row, reg in zip(regions, every, strict=True):
            all_regions.append((
                "absent" if is_absent else "present",
                len([m for m in row["methods"].split(",") if m]),
                reg.subclade_identifiable,
            ))

    for label, group in (("donor-present", present), ("donor-absent", absent)):
        judged = [r for r in group if r[1] is not None]
        flagged = [r for r in judged if r[1] is False]
        print(f"\n=== {label}: {len(group)} cases, {len(judged)} judgeable, "
              f"{len(flagged)} flagged 'sub-clade unresolved' ===")
        for name, identifiable, tied in group:
            if identifiable is None:
                verdict = "n/a (no competing lineage present)"
            elif identifiable:
                verdict = "resolved"
            else:
                verdict = "FLAG tied=" + ",".join(tied)
            print(f"  {name:<24} {verdict}")

    if not present and not absent:
        print("No case artifacts found. Run validation/run_hybrids.py first.", file=sys.stderr)
        return 2
    _report_all_regions(all_regions)
    return 0


def _report_all_regions(all_regions: list[tuple[str, int, bool | None]]) -> None:
    """The flag over every called region, split by how many callers found the region.

    The per-case figures above use each case's primary region. A run reports more than that,
    and a region only one caller found is the usual case where the sub-clade is not resolved.
    """
    print("\n=== every called region (not only each case's primary), by number of callers ===")
    for label in ("present", "absent"):
        evaluated = [r for r in all_regions if r[0] == label and r[2] is not None]
        if not evaluated:
            continue
        no = [r for r in evaluated if r[2] is False]
        one = [r for r in evaluated if r[1] == 1]
        many = [r for r in evaluated if r[1] >= 2]
        print(f"  donor-{label}: {len(no)} of {len(evaluated)} evaluated regions marked no | "
              f"found by 1 caller: {sum(r[2] is False for r in one)} of {len(one)} | "
              f"found by 2+ callers: {sum(r[2] is False for r in many)} of {len(many)}")
    evaluated = [r for r in all_regions if r[2] is not None]
    no = [r for r in evaluated if r[2] is False]
    one = [r for r in evaluated if r[1] == 1]
    many = [r for r in evaluated if r[1] >= 2]
    print(f"  both: {len(no)} of {len(evaluated)} evaluated regions marked no | "
          f"found by 1 caller: {sum(r[2] is False for r in one)} of {len(one)} | "
          f"found by 2+ callers: {sum(r[2] is False for r in many)} of {len(many)}")


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
