# A cleaner `reassort --scan-segments` validation probe -- design

## Context

The `reassort --scan-segments` intragenic scan (v2, merged) was first exercised with an ad-hoc
scratch script that spliced two per-clade **consensus** sequences into a synthetic HA query. That
probe was inconclusive: a consensus-vs-consensus query against a broad panel produced whole-genome,
mis-attributed regions, and it could not cleanly show localization or specificity
(`validation/attribution-results.md`). It also lived only in scratch, so it was not reproducible.

This spec replaces it with a durable, opt-in probe that builds the HA query from **real tip
genomes** (not consensuses) with an exact known breakpoint, runs the shipped scan, and reports
whether the scan **localizes** the intragenic recombination -- plus a clonal control for specificity.

## Decisions (from brainstorming)

- **Mosaic shape:** a central **insert** -- a clade-B segment (fraction 0.35-0.65) spliced into a
  clade-A backbone -- the localizable recombinant shape the caller is tuned for (not a half/half
  breakpoint).
- **Gating:** XPASS is gated on **localization only** (a called region overlaps the true insert
  span). Donor-clade attribution and the clonal-control specificity are **reported, not gating**.
- **Home:** a dedicated `validation/run_reassort_scan.py`, opt-in (needs the aligner env, not CI),
  importing helpers from `run_hybrids` rather than duplicating them.

## Architecture

New opt-in harness `validation/run_reassort_scan.py`. It reuses, by import from `run_hybrids`:
`_load_species` (dataset -> genomes, reference, tips), `pick_parents` (most-divergent clade pair +
their central-genome accessions), `make_hybrid` (splice an insert-shaped query, returning
`(query, q_start, q_end)` in query coordinates), `reconstruct_gapped`, `parse_regions`, and the
module constants (`INSERT`, `MIN_MEMBERS`). It drives the **shipped** code path via
`tessera.reassort.assign_segments(..., scan_segments=True, output=...)` -- not a reimplementation.

Two flu H3N2 datasets are used (the same ones the reassort probes already use):
`nextstrain/flu/h3n2/ha/EPI1857216` and `nextstrain/flu/h3n2/na/EPI1857215`.

### Data flow

```
HA dataset -> _load_species -> (genomes, ref_ha, tips_ha)
  pick_parents(tips_ha, ref_ha) -> (clade_a, clade_b, src_a, src_b)   # most-divergent pair
  make_hybrid(ref_ha, tips_ha[src_a][1], tips_ha[src_b][1])           # real-tip insert
    -> (ha_query, q_start, q_end)                                     # true donor span
NA dataset -> _load_species -> (_, ref_na, tips_na)
  a real single-clade NA tip -> reconstruct_gapped(ref_na, tips_na[src][1]).replace("-","")
write multi-FASTA: >HA <ha_query>, >NA <na_query>
assign_segments(query, dataset_overrides={HA:..., NA:...}, scan_segments=True, output=out)
  -> reads out/HA/recombination_regions.tsv and out/NA/recombination_regions.tsv
score:
  HA localization: XPASS iff a called region overlaps [q_start, q_end]      (gating)
  HA attribution : does any region name clade_b as a parent?                (reported)
  NA specificity : any region called on the clonal control?                 (reported)
```

### The scorer (pure, unit-tested)

```
def region_overlaps_span(regions: list[dict], q_start: int, q_end: int) -> bool
```

`regions` are `parse_regions` rows (each has `query_start`, `query_end`, `donor_absent`). Returns
True iff any present (non-`donor_absent`) region `[query_start, query_end]` overlaps
`[q_start, q_end]` (standard interval overlap: `r_start <= q_end and r_end >= q_start`). This is the
only nontrivial pure logic and is CI-unit-tested with a synthetic TSV.

## Output

Prints a small report (mirroring the frontier table style), for example:

```
Tessera reassort --scan-segments probe (opt-in; needs the aligner env)
======================================================================
[XPASS ] HA localization  insert clade B (0.35-0.65) recovered; region 612-1190 overlaps 601-1122
         HA attribution   donor clade B named in a region: yes
[report] NA specificity   clonal control: 1 region called (a false positive)
```

The end-to-end result is appended to `validation/attribution-results.md` and the harness is listed
in `validation/README.md`. A clonal-control false positive or an HA miss is stated plainly, never
retuned to force a pass.

## Error handling

- Needs skani + an aligner + network (Nextclade). Missing skani/aligner -> the shipped
  `assign_segments`/`require_aligner` raises a clear error up front.
- Too few eligible clades for a divergent pair -> `pick_parents` raises `CaseSkipped`; the probe
  reports SKIP, not failure.
- A per-segment scan failure is already non-fatal in `assign_segments` (recorded as not-scanned);
  the probe reports it faithfully.

## Testing / verification

- **CI unit test (no binaries):** `region_overlaps_span` on synthetic `parse_regions` rows --
  overlap, no-overlap, touching-edge, and a `donor_absent` row excluded.
- **Opt-in end-to-end (aligner env, not CI):** run `python validation/run_reassort_scan.py`; expect
  HA localization XPASS (a region overlapping the insert) and the NA specificity line, recorded as
  measured.
- `ruff check src tests validation` clean; CI coverage floor (80%) holds.

## Scope / non-goals

- One HA-insert + one NA-clonal probe, on the two flu H3N2 datasets. Not a broad multi-pathogen
  sweep and not part of CI.
- Localization is the only gate. Improving low-divergence specificity or attribution of the scan is
  the single-backbone caller's concern, out of scope here.
- No change to `src/` (the probe drives the shipped `assign_segments`); the only non-validation file
  touched is a new CI unit test for the pure scorer.
