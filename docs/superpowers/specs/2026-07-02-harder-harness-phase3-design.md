# Harder harness -- Phase 3 (frontier: inter-species + reassortment) -- design

## Context

Phases 1-2 hardened the hybrid harness on cases the tool is designed to handle. Phase 3 probes the
**edges of the design envelope**: cross-species recombination (which stresses the HMM's identity
floor and the aligner's divergence ceiling) and segment reassortment (which the pipeline has no
concept of). These are not "should-pass" cases; they are **measured known-limitations**. The
load-bearing requirement is faithful reporting: frontier cases run only under `--frontier`, are
scored XFAIL/XPASS/KNOWN-LIMIT in a separate block, and can **never** count as a regression against
the must-pass headline. This is Phase 3 of the approved harder-harness plan; it is harness-only.

## Architecture

- **`tier` gate.** Cases gain `tier` (default `"must_pass"`; frontier cases set `"frontier"`).
  `main` runs `must_pass` only unless `--frontier` is passed; the frontier set prints in its own
  table and its accounting lives in a separate dict that never touches `passed`/`errored` or the
  return code.
- **Two frontier case types**, each needing a **second dataset** (`case["second_dataset"]`):
  `inter_species` (Category 6) and `reassortant` (Category 5). A `_prepare_frontier_case` resolves
  both datasets and builds the two-source query + panel; a `_score_frontier` emits the verdict.
- Reuse: `resolve_dataset`, `build_pool` (`discover/nextclade.py`); `reconstruct_gapped`,
  `consensus_sequence`, `pct_identity`, `make_hybrid`, `clade_members`, `collect_tips`,
  `_download_text`, `read_reference` (harness); `select_regional`/`representative_panel` for the
  panel; `run_recomb` + `parse_regions` for the scan.

## Category 6 -- inter-species (RSV-A x RSV-B)

**Build.** Resolve dataset A (`nextstrain/rsv/a/...`) and B (`nextstrain/rsv/b/...`); build both
tip pools. Because the two clade systems are disjoint, prefix each pool's clade labels with the
species (`A:<clade>`, `B:<clade>`) to avoid collision. Pick a backbone clade from A and a donor
clade from B (their central genomes). Splice a cross-species hybrid with `make_hybrid` (A backbone,
B insert); the combined panel is the union of both species' source-removed pools.

**Envelope + scoring (`_score_frontier`, "inter_species").** Compute the A-B pair identity via
`pct_identity` at setup.
- **>= 80% identity (in-envelope):** run and score like a positive (detection + backbone clade +
  donor-region recovery). Record **XPASS** if it detects and attributes, else **XFAIL** -- either
  way in the frontier block, never a must-pass FAIL.
- **< 80% identity (out-of-envelope):** the HMM clamps identity to `[0.80, 0.999]` (`hmm.py:54-57`)
  and the default aligners top out ~91%, so attribution is not meaningful. Mark **KNOWN-LIMIT**
  with the measured identity -- a documented boundary, not a bug.

RSV-A x RSV-B is ~81% genome-wide -- right at the floor -- so this case measures exactly where the
envelope ends.

## Category 5 -- reassortment (flu H3N2 HA + NA), feasibility probe

**Build (`_prepare_frontier_case`, "reassortant").** Resolve the HA dataset and the NA dataset
(`case["second_dataset"]`). For each segment, group tips by clade and take a per-clade majority
consensus (`consensus_sequence`, as `consensus_panel` does). Choose a backbone clade X and a donor
clade Y (per segment; the two segments' clade systems differ, so X/Y are chosen independently as
the two most-divergent clades within each segment's own consensus set). The panel references are the
**concatenated** per-clade consensuses `HA_consensus[c] ++ NA_consensus[c]` (compound label
`HA:<c>|NA:<c>`); the reassortant **query** is `HA_consensus[X] ++ NA_consensus[Y]` -- a mosaic whose
single breakpoint is the HA|NA junction.

**Scoring (`_score_frontier`, "reassortant") -- detection-gated only.** PASS/XPASS = a recombinant
region is called overlapping the junction (query coordinate = `len(HA_consensus[X])`). The
whole-segment-swap **semantic labeling** (which "parent" the query is) is a documented
**KNOWN-LIMIT**: the pipeline has no reassortment concept and no genome-wide clade exists for a
reassortant, so we assert only that a breakpoint is *found*, not that it is labeled a reassortment.

This is a probe. If the two segment datasets cannot be assembled into a consistent panel (e.g. too
few clades per segment, or an empty consensus), the case is `CaseSkipped` -- reported as skipped in
the frontier block, not failed.

## Faithful-reporting invariant

- Without `--frontier`: only `must_pass` runs; the headline (`sensitivity`, `specificity`) is
  unchanged from Phase 2 (frontier cases are absent).
- With `--frontier`: a separate table lists each frontier case as `XPASS` / `XFAIL` / `KNOWN-LIMIT`
  / `SKIP` with the measured identity (inter-species) or junction outcome (reassortant). The
  frontier tally is printed apart from the must-pass tally and does not affect the process return
  code. An out-of-envelope case is *structurally incapable* of appearing as a must-pass FAIL.
- Record the measured envelope in `validation/attribution-results.md`; document `--frontier` in
  `validation/README.md`.

## Edge cases
- A frontier dataset that fails to resolve / build -> `CaseSkipped` (reported skipped).
- Inter-species pair >= 91% (unusually close species) -> still scored in-envelope (XPASS/XFAIL).
- Reassortant where a segment has < 2 clades -> `CaseSkipped`.
- Running without `--frontier` -> frontier cases are filtered out before scoring; zero effect.

## Testing / verification
- **CI unit tests (no binaries):** the `tier` filter (must_pass-only vs `--frontier`); the
  identity-envelope branch of `_score_frontier` (>=80% -> XPASS/XFAIL; <80% -> KNOWN-LIMIT) with a
  stubbed identity; the reassortant junction-overlap check with a synthetic regions TSV; the
  frontier tally is separate from the must-pass tally and never changes the return code.
- **Opt-in harness gate** (`recomfi-aln`): `python validation/run_hybrids.py --frontier` runs the
  two frontier cases; record the measured RSV-A x RSV-B identity + verdict and the flu reassortant
  junction outcome in `attribution-results.md`. `python validation/run_hybrids.py` (no flag)
  confirms the must-pass headline is unchanged.
- `ruff check src tests validation` clean; CI coverage floor (80%) holds.

## Scope / non-goals
- Only the harness: the `tier`/`--frontier` gate, the two frontier case types + their prepare/score
  functions, and the frontier reporting block. No product/pipeline or data-layer (src/) change.
- No combinatorial HA x NA panels; no attempt to make the pipeline *understand* reassortment. The
  reassortant case is a detection probe with the labeling limitation documented.
- No new runtime dependency.
