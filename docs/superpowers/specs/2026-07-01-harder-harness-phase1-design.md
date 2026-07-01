# Harder harness -- Phase 1 (specificity + low-divergence + panel-adversarial) -- design

## Context

The synthetic-hybrid harness (`validation/run_hybrids.py`) reached a clean **18/18** and is now
too easy to discriminate: every case is a single insert spliced from the **most-divergent** viable
clade pair, it measures **sensitivity only** (no false-positive concept anywhere), and below 4%
divergence it drops the backbone-attribution requirement. This is Phase 1 of a phased hardening
(full design: `~/.claude/plans/jazzy-squishing-quasar.md`, user-approved). Phase 1 is **harness-only,
no data-layer changes**, and delivers the biggest honesty gain: it turns "18/18 detection" into a
**precision + recall** harness and adds genuinely hard low-divergence and adversarial-panel cases.

Phases 2 (topology refactor: `make_mosaic`/`true_spans`, multi-breakpoint, masking sibling) and 3
(frontier: inter-species + reassortment, XFAIL known-limitations) are separate later cycles.

## Architecture

Extend `run_hybrids.py` in place (no fork). Add two discriminator keys to each `HYBRIDS` dict
(`:74-105`), defaulted so all 24 existing cases are byte-unchanged:

- `case_type: str = "single_insert"` -- selects query construction + scoring branch.
- `tier: str = "must_pass"` -- Phase 1 introduces the field; all Phase-1 cases are `must_pass`
  (the `frontier` value is used in Phase 3, but the key is added now so the summary split is stable).

Branch on `case.get("case_type")` in `_prepare_case` (`:472`) and `_score_regions` (`:643`). Reuse
without reimplementation: `pick_parents` (`:337-380`), `make_hybrid` (`:383-392`),
`reconstruct_gapped` (`:127-141`), `clade_match`/`donor_match`/`attribution_tier`/
`shared_clade_depth` (`:203-262`), the source-removed pool build (`:507-528`), and the pipeline's
existing donor-absent regions (`coverage.gaps_as_regions`, `coverage.py:151-181`; TSV field
`donor_absent`, already filtered at `_score_regions:651`).

## Components

### C2 -- negative controls / specificity (highest value)
- **`case_type="neg_pure"`**: `_prepare_case` writes the clade-A central genome *unspliced* as the
  query (`reconstruct_gapped` + strip gaps); its source is removed from the pool as today;
  `q_start=q_end=0` (no true span). Scoring branch: `passed = (len(present) == 0)`; record
  `n_false_regions = len(present)`.
- **`case_type="neg_within"`**: a new helper (beside `pick_parents`) picks the two most-divergent
  *members of the same clade* as src_a/src_b; splice via `make_hybrid`. Scoring:
  `passed = no called region whose minor-parent top-level clade differs from the major's`. A
  cross-clade call is the false positive (guards `plurality_major` manufacturing a donor from noise).
- **Summary (`_run_default` `:778-782`)**: split counters by case_type; print two lines --
  `sensitivity P/N` (positive cases) and `specificity Q/M (K false calls)`. Never merge.
- Candidate taxa: `measles`, `mumps`, `dengue`, `hiv1` (HIV subtype B internal diversity for
  `neg_within`).

### C1 -- low-divergence attribution ("less difference")
- **`pick_parents` (`:372-377`)**: add `case.get("pair_objective")` (`"max"` default keeps current
  behaviour; `"min"` selects the closest pair, guarded `div >= min_divergence` so <0.5% pairs that
  collapse under dereplication `CaseSkipped` rather than score).
- Add per-case `divergence_band: (lo, hi)`; checked in `_prepare_case` after `divergence` is known
  (`:488-491`) -> `CaseSkipped` if outside, so the case asserts it tests the hard regime.
- **Honest scoring floor** (`case_type="low_div"`, at `_score_regions:665`): replace today's
  "drop backbone < 4%" with **PASS = detection + donor-region overlaps span + donor top-level clade
  correct + backbone same top-level clade** (`attribution_tier != "mismatch"` on both). The
  exact-vs-sibling rate is reported via existing `backbone_tier`/`donor_tier` (non-gating), so
  tightening later is a data change. Forbids the "any backbone" free pass without claiming
  sub-clade resolution the informative-site path (auto-triggers <8% polymorphic cols) cannot give.
- Candidate taxa: `dengue` within-serotype (~3-4%), `rsv_a` `A.D.1.8` vs `A.D.5.2` (~1-3%), `wnv`
  1a/1b.

### C4a + C4c -- panel-adversarial (pool manipulation only)
- **`case_type="panel_donor_absent"`**: after `members_by_clade` (`:518-528`) drop *all* clade-B
  members from pool + members. Scoring: `passed = the true span overlaps a donor_absent region AND
  no confident cross-clade present-region mis-attributes the span` (read `donor_absent=="yes"`,
  the inverse of the `:651` filter). Exercises `coverage.gaps_as_regions` (the honest missing-donor
  path).
- **`case_type="panel_equidistant"`**: pin clade-B plus a clade-C with
  `divergence(A,B) ~= divergence(A,C)` (within ~0.5%), both present; splice donor from B only.
  Scoring: `passed = donor attributed to B, not C` -- guards the `plurality_major` tie-break
  (`hmm.py:78-83`). Report which won.
- Candidate taxa: `rsv_a` (donor-absent), `dengue`/`hiv1` (equidistant).

## Data flow

Unchanged pipeline. Each case: `_prepare_case` builds the query (spliced / unspliced / same-clade)
and the pool (standard / donor-clade-removed), runs the scan, and `_score_regions` selects the
verdict branch by `case_type`. The default and `--compare` runners partition the summary into
sensitivity (positive cases) and specificity (negative cases).

## Edge cases
- All 24 existing cases default to `single_insert` / `must_pass` -> identical behaviour and the
  existing 18/18 must hold.
- `neg_within` where the clade has <2 divergent members -> `CaseSkipped`.
- `low_div` where no pair clears `min_divergence` or the band -> `CaseSkipped` (not a FAIL).
- `panel_donor_absent` where removing clade-B empties the pool -> `CaseSkipped`.

## Testing / verification
- **CI unit tests** (no binaries, `tests/unit/`): the `neg_pure` / `neg_within` / `low_div` /
  `panel_donor_absent` / `panel_equidistant` scoring branches driven by a synthetic regions TSV;
  the `pick_parents` `"min"` objective and the same-clade member picker; the summary counter split
  (sensitivity vs specificity). These cover the harness *logic* independent of any aligner.
- **Opt-in harness run** (aligner env `recomfi-aln`, not CI):
  `python validation/run_hybrids.py --compare` -- the existing 18/18 positives hold, the new
  must-pass hard cases are scored, and a specificity line is reported. Record in
  `validation/attribution-results.md` faithfully (a new hard case that fails is a real finding).
- `ruff check src tests validation` clean; CI coverage floor (80%) holds.

## Scope / non-goals
- Phase 1 only: `case_type`/`tier` plumbing, the five new case types above, and the summary split.
- No `make_mosaic`/`true_spans` refactor (Phase 2), no multi-breakpoint/short/terminal topologies
  (Phase 2), no frontier inter-species/reassortment (Phase 3), no data-layer or pipeline changes.
- No new runtime dependency.
