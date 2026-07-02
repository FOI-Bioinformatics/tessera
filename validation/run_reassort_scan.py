"""Opt-in probe for ``reassort --scan-segments``: does the intragenic scan localize a real
recombinant?

Builds an HA query from two divergent H3N2 HA clade tips (a clade-B insert in a clade-A
backbone, real genomes with an exact known breakpoint) plus a clonal NA control, runs the
shipped ``assign_segments(scan_segments=True)``, and reports whether a called region overlaps
the true insert span (localization, the XPASS gate), whether the donor clade is named
(attribution, reported), and whether the clonal control is flagged (specificity, reported).

Needs skani, an aligner, and network access to the Nextclade datasets. Not part of CI. Run:

    export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
    python validation/run_reassort_scan.py
"""

from __future__ import annotations

import logging
import sys
import tempfile
from pathlib import Path

HA_DATASET = "nextstrain/flu/h3n2/ha/EPI1857216"
NA_DATASET = "nextstrain/flu/h3n2/na/EPI1857215"


def region_overlaps_span(regions: list[dict], q_start: int, q_end: int) -> bool:
    """True iff any present (non-``donor_absent``) region overlaps ``[q_start, q_end]``."""
    for r in regions:
        if r.get("donor_absent") == "yes":
            continue
        try:
            r_start, r_end = int(r["query_start"]), int(r["query_end"])
        except (KeyError, ValueError):
            continue
        if r_start <= q_end and r_end >= q_start:
            return True
    return False


def _clonal_na(tips_na, reference_na, reconstruct_gapped, largest_clade):
    """A real single-clade NA genome (the largest clade's central tip) -> (clade, sequence)."""
    clade, src = largest_clade(tips_na)
    seq = reconstruct_gapped(reference_na, tips_na[src][1]).replace("-", "").upper()
    return clade, seq


def _run_probe(out_dir: Path, logger: logging.Logger) -> int:
    """Build the query, run the shipped scan, print the report. Always returns 0 (an opt-in
    probe never fails the caller); a SKIP is printed when no divergent HA pair is available."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import run_hybrids as rh

    from tessera.core.io import write_fasta_record
    from tessera.reassort import assign_segments

    _g, ref_ha, tips_ha = rh._load_species(HA_DATASET, None, logger)
    try:
        clade_a, clade_b, src_a, src_b = rh.pick_parents(tips_ha, ref_ha, [], logger)
    except rh.CaseSkipped as exc:
        print(f"[SKIP  ] HA localization  {exc}")
        return 0
    ha_query, q_start, q_end = rh.make_hybrid(ref_ha, tips_ha[src_a][1], tips_ha[src_b][1])

    _g2, ref_na, tips_na = rh._load_species(NA_DATASET, None, logger)
    na_clade, na_query = _clonal_na(tips_na, ref_na, rh.reconstruct_gapped, rh._largest_clade)

    query = out_dir / "query.fasta"
    with open(query, "w") as fo:
        write_fasta_record(fo, "HA", ha_query)
        write_fasta_record(fo, "NA", na_query)
    logger.info("HA insert: clade %s in clade %s backbone (true span %d-%d); NA clonal %s",
                clade_b, clade_a, q_start, q_end, na_clade)

    assign_segments(query, dataset_overrides={"HA": HA_DATASET, "NA": NA_DATASET},
                    scan_segments=True, output=out_dir, logger=logger)

    ha_regions = rh.parse_regions(out_dir / "HA" / "recombination_regions.tsv")
    na_regions = rh.parse_regions(out_dir / "NA" / "recombination_regions.tsv")
    localized = region_overlaps_span(ha_regions, q_start, q_end)
    attributed = any(clade_b in r.get("minor_parent", "") or clade_b in r.get("major_parent", "")
                     for r in ha_regions)
    na_present = [r for r in na_regions if r.get("donor_absent") != "yes"]

    print("\nTessera reassort --scan-segments probe (opt-in; needs the aligner env)")
    print("=" * 70)
    verdict = "XPASS " if localized else "XFAIL "
    print(f"[{verdict}] HA localization  insert clade {clade_b} (true span {q_start}-{q_end}); "
          f"{'a region overlaps it' if localized else 'no region overlaps it'}")
    print(f"[report] HA attribution   donor clade {clade_b} named in a region: "
          f"{'yes' if attributed else 'no'}")
    print(f"[report] NA specificity   clonal control ({na_clade}): "
          f"{len(na_present)} region(s) called")
    print("\n(opt-in probe; localization is the gate, attribution/specificity are reported)")
    return 0


def main(argv: list[str]) -> int:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    logger = logging.getLogger("reassort-scan-probe")
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(argv[0]) if argv else Path(tmp)
        out.mkdir(parents=True, exist_ok=True)
        return _run_probe(out, logger)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
