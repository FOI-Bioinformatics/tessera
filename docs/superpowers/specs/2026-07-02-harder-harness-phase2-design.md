# Harder harness -- Phase 2 (hard topologies + masking sibling) -- design

## Context

Phase 1 turned the hybrid harness into a precision/recall instrument (specificity, low-divergence,
panel-adversarial) and closed the `donorabsent` finding. Every positive case is still a **single**
centered insert. Phase 2 adds harder *topologies* -- multi-breakpoint mosaics, short tracts,
asymmetric splits, terminal breakpoints -- and a **masking-sibling** attribution case, so the
harness exercises segmentation with many switches and near-neighbour donor discrimination. This is
Phase 2 of the approved harder-harness plan (`~/.claude/plans/jazzy-squishing-quasar.md`); it is
harness-only (no product/data-layer change). Phase 3 (frontier inter-species / reassortment, XFAIL)
remains a separate cycle.

## Architecture

Extend `validation/run_hybrids.py`. The centerpiece is generalising the single splice to a mosaic,
done **additively** to keep the Phase-1 scorers untouched:

- **`make_mosaic(reference, segments) -> (query, true_spans)`** where `segments` is
  `list[(donor_muts, frac_start, frac_end)]` (each segment reconstructed from a clade source's
  mutations over `[frac_start, frac_end)` of the reference). `true_spans` is
  `list[(q_start, q_end, donor_clade)]` for the **non-backbone** segments only. `make_hybrid`
  becomes a thin wrapper: `make_mosaic(ref, [(A,0,.35),(B,.35,.65),(A,.65,1)])` -> one B span, so
  existing single-insert cases are byte-identical.
- **`CaseSetup` gains `true_spans: list[tuple[int, int, str]]`** (default `[]`). The existing
  scalar `q_start`/`q_end` are **kept** as the overall envelope (min start / max end of
  `true_spans`, or `0` when empty), so the Phase-1 scorers (`_score_single_insert`, `_score_low_div`,
  `_score_panel_*`) read them unchanged. Only the new mosaic scorer loops `true_spans`.
- **`pick_parents_n(tips, reference, n, *, floor) -> list[(clade, src)]`** -- greedily pick the `n`
  most mutually-divergent clades (pick the most-divergent pair, then repeatedly add the clade
  maximising the minimum distance to those chosen). Reuses `pct_identity` / `clade_representative`.

`_prepare_case` branches query construction on `case_type`; a new `mosaic` type reads a `pattern`
(and `n_parents`) and builds the `segments` + `true_spans`. `_score_regions` gains a `_score_mosaic`
scorer (registered in the dispatch dict).

## Components

### Category 3 -- hard topologies (`case_type="mosaic"`)
`_prepare_case` builds segments per `pattern`:
- **`ABAC`** (`n_parents=3`, `pick_parents_n(...,3)`): segments `[(A,0,.3),(B,.3,.45),(A,.45,.7),
  (C,.7,1)]` -> two non-backbone spans (B, C). Backbone A = the largest/most-central clade.
- **`AB_9010`**: `[(A,0,.45),(B,.45,.55),(A,.55,1)]` -- a 10% donor tract.
- **`AB_short`**: donor width `< one recomb window` (`window_params` gives `window`); set the B
  fraction to `~0.5 * window / genome_len`.
- **`AB_terminal`**: `[(B,0,.15),(A,.15,1)]` -- donor at the genome start, no left flank.

**`_score_mosaic`** (new): `major_clade = clade_of(regions[0].major_parent)`. For each
`(q_start, q_end, donor_clade)` in `true_spans`, a span is *recovered* when some span-overlapping
present region has `donor_match(clade_of(minor), donor_clade, clade_a)`. Verdict by pattern:
- `ABAC`, `AB_9010`: PASS = detection + backbone `clade_match(major, clade_a)` + **all** spans
  recovered (report `spans_hit/spans_total`).
- `AB_terminal`: PASS as above **and** the recovered region's `query_start` is within one window of
  0 (guards a whole-genome mis-call spuriously "overlapping").
- `AB_short`: **detection-gated only** -- PASS = a region is detected; span recovery is *reported*,
  not gating (a miss is visible, never a regression, because a sub-window tract may be windowed
  away). This is the one intentionally-lenient pattern.

Candidate taxa (need >=3 well-separated clades under the one-dataset-per-genome constraint):
`dengue` (4 serotypes), `hiv1` (subtypes A/B/C/D) for `ABAC`; `measles`/`mumps` for short/terminal.

### Category 4b -- masking sibling (`case_type="mask_sibling"`)
Pin a case where the donor clade B has a **sibling sub-clade** present in the panel (e.g. RSV
`A.D.1.8` donor with `A.D.5.2` present -- the diagnosed rsv_a contender). Construction is a standard
single splice from B (reuse `make_hybrid`), asserting (via `CaseSkipped` if not) that a sibling of B
remains in `members_by_clade`. **Scorer** tightens the donor requirement to **exact**:
`_score_mask_sibling` = detection + backbone `clade_match` + a span-overlapping region with
`attribution_tier(clade_of(minor), clade_b) == "exact"` (not merely `sibling`). This is a direct
regression guard for the plurality-major fix (a sibling winning the donor is a FAIL). Report the
sibling-confusion rate.

## Data flow
```
pattern -> segments -> make_mosaic -> (query, true_spans)   [+ derived q_start/q_end envelope]
run scan -> _score_mosaic loops true_spans (per-span donor recovery) / _score_mask_sibling (exact)
```
Phase-1 case types and their scorers are unchanged; `single_insert` now routes through
`make_mosaic` via the `make_hybrid` wrapper (byte-identical output).

## Edge cases
- A dataset with `< n_parents` eligible clades -> `pick_parents_n` raises `CaseSkipped`.
- `mask_sibling` where B has no sibling sub-clade present -> `CaseSkipped`.
- `AB_short` where the tract is below `MIN_GENOME`-scale resolution -> detection-gated, so a miss is
  reported, not failed.
- Single-insert / neg / low_div / panel cases: `true_spans` has 0 or 1 entries; the derived
  `q_start`/`q_end` reproduce today's values exactly.

## Testing / verification
- **CI unit tests (no binaries)**: `make_mosaic` span math (segment boundaries, `true_spans` for a
  3-parent pattern, the `make_hybrid` wrapper equals the old single-insert output); `pick_parents_n`
  (n most-divergent, `CaseSkipped` when too few); `_score_mosaic` (all-spans-recovered PASS, a
  missed span FAIL for ABAC, detection-gated PASS for AB_short) and `_score_mask_sibling`
  (exact-donor PASS, sibling-donor FAIL) with synthetic regions TSVs; the derived-envelope
  `q_start`/`q_end` from `true_spans`.
- **Opt-in harness gate** (`recomfi-aln`): `python validation/run_hybrids.py` -- the Phase-1 set
  holds (single-insert byte-identical), and the new `mosaic`/`mask_sibling` cases are scored;
  record the outcome in `validation/attribution-results.md`. A genuine FAIL is a finding (record it,
  as with `donorabsent`), not something to tune away.
- `ruff check src tests validation` clean; CI coverage floor (80%) holds.

## Scope / non-goals
- Only the harness: `make_mosaic`/`true_spans`, `pick_parents_n`, the `mosaic` + `mask_sibling` case
  types and scorers, and their case entries. No product/pipeline or data-layer change.
- No frontier inter-species / reassortment (Phase 3). No new runtime dependency.
