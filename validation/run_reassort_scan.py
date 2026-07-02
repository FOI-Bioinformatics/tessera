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
