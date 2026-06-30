# Donor re-attribution -- design

## Context

The gap analysis (`docs/gap-analysis.md`, G1) found that Tessera's two synthetic-harness FAILs
are *attribution* misses, not detection misses: the recombination is found but a parent is
named as the wrong adjacent sub-clade. The G1 measurement
(`validation/attribution-results.md`) then ruled out the obvious panel-composition fixes:

- the **barcode** caller is inert (no verdict or attribution change on 18 cases);
- a **per-clade consensus** panel is a net regression (16 -> 14): it sharpens the *donor*
  (`rsv_a` sibling -> exact, `marburg` sibling -> exact) but degrades the *backbone* on
  `measles` and `rubella`;
- a **mixed** panel (tips + consensus) is inconsistent: it removed the `measles` backbone
  regression but not `rubella`'s, and did not close `rsv_a` (the consensus genome wins the
  donor region *and* steals the backbone slot).

The root cause is structural: **one reference panel competes for both parent roles**, so any
reference good enough to fix the donor is also good enough to displace the backbone. The fix
has to separate the two representations -- score the **major parent (backbone)** against
individual genomes (as today) and the **minor parent (donor)** against denoised per-clade
consensuses. The codebase already learned (and recorded) that attribution fixes must be
**post-hoc**, never by changing what competes for detection (an earlier lineage-aware
competition change regressed the panel and was reverted). This design honours that: detection
and the backbone are untouched; only the donor *label* is refined after the regions are called.

The intended outcome: close `rsv_a` (donor recovered exactly as `A.D.1.8`, backbone still
`A.1`) without regressing the 16 currently-passing cases -- proven on the `--compare` harness
before the behaviour is ever a default.

## Goal

An **opt-in** post-merge step that re-labels each called region's donor to the clade whose
consensus best matches the query over that region, guarded by a similarity **margin** so it
only acts on a clear improvement. The backbone is never changed; detection, support, q-values,
breakpoints and method agreement are never changed. Off by default; validated, then promoted
to default-when-typed in a separate follow-up only if the harness shows it is clean.

## Components

### `src/tessera/recomb/reattribute.py` (new -- the substance)

A single pure entry point, dependency-free numpy:

```
def reattribute_donors(
    regions: list[Region],
    result: WindowSimilarity,
    lineage_map: LineageMap,
    major_clade: str | None,
    *,
    margin: float = 0.03,
    min_sites: int = 20,
) -> list[Region]
```

1. **Per-clade consensus, once per run.** Group the panel's reference rows
   (`result.rows`, minus the query) by clade via `lineage_of(label, lineage_map)`. For each
   clade build a column-wise majority canonical base array (the `barcode.clade_markers` /
   `consensus.consensus_sequence` primitive: count A/C/G/T per column, take argmax). Record a
   representative member label per clade (`by_clade[clade][0]`, as the barcode caller does), so
   a re-attributed region can name an actual genome the report types by lineage.
2. **Per region, score each clade over the span.** Restrict to the region's MSA columns
   (`region.msa_start:region.msa_end`). For each clade, similarity = fraction of columns where
   the query and that clade's consensus carry the *same* canonical base, over columns where
   both are canonical; require at least `min_sites` comparable columns (else skip that clade).
3. **Margin-guarded full re-attribution.** Let `current = lineage_of(region.minor_parent)`.
   Among all clades **except `major_clade`**, find the best-scoring `best`. If
   `score[best] - score[current] >= margin` (treating a missing `current` score as 0), re-label:
   `minor_parent := rep[best]`, and recompute `mean_sim_minor` over the span for transparency.
   Otherwise return the region unchanged. Never re-attribute to `major_clade`; never act when a
   clade lacks the markers/sites to be scored.
4. Return a new region list (regions are replaced, not mutated in place); log each change at
   INFO (`re-attributed donor <old> -> <new> over query <span>`).

The function is a no-op (returns `regions` unchanged) when `lineage_map` is falsy or fewer than
two clades can be scored.

### `src/tessera/recomb/run.py` (integration -- one call site)

After `consensus_regions(...)` produces `regions` (`run.py:243-247`), and only when
`params.reattribute_donors` and a `lineage_map` are present:

```
regions = reattribute_donors(regions, result, lineage_map, major_parent,
                             margin=params.reattribute_margin)
```

Add `reattribute_donors: bool = False` and `reattribute_margin: float = 0.03` to
`RecombParams`. The step runs before the report is written, so `method_breakdown` (built from
the pre-re-attribution merge) keeps recording which callers found each region; the re-attributed
label is what the report and TSV show. (The breakdown's `minor_parent` is cosmetic; if it
diverges confusingly we re-point it to the new label in the same step -- a one-line follow-up,
noted, not required for correctness.)

### CLI surface

`--reattribute-donors/--no-reattribute-donors` (default off) and `--reattribute-margin`
(default 0.03) on `cli/cmd_recomb.py`, threaded through `RecombParams`. `cli/cmd_detect.py` and
`cli/cmd_fill_references.py` (+ `FillParams`) gain the same flag so the recruit-then-detect
paths can opt in. Help text notes it needs a typed panel (Nextclade or `--lineage-map`) and is
a donor-label refinement only.

## Data flow

```
callers -> consensus_regions -> regions (donor from tips)
                                   |
              params.reattribute_donors & lineage_map?
                                   | yes
   per-clade consensus over result.rows, score each region's span,
   re-label minor_parent when best non-backbone clade beats current by >= margin
                                   |
                         regions (donor refined) -> report / TSV
```

Backbone (`major_parent` / `reconcile_major`) and detection are entirely upstream and unchanged.

## Validation

- **Harness gate.** Add a `reattribute` config to the `--compare` grid in
  `validation/run_hybrids.py` (tip panel, default ensemble, `reattribute_donors=True`); this
  needs `_build_and_score` to pass the flag into `RecombParams`. Acceptance: the `reattribute`
  column flips `rsv_a` to PASS (donor `exact`, backbone still `exact`) with **zero regressions**
  on the 16 baseline PASSes and no new donor `mismatch`. Run with the `recomfi-aln` env on PATH.
- **Unit tests (CI, no binaries)** in `tests/unit/test_reattribute.py` on a small synthetic
  typed `WindowSimilarity`: (a) a region whose span matches clade C's consensus better than its
  current clade B by > margin is re-labelled to C's representative; (b) below margin it is
  unchanged; (c) the backbone/major clade is never chosen even if it scores best; (d) no-op when
  `lineage_map` is empty or < 2 scorable clades; (e) a clade with < `min_sites` comparable
  columns in the span is skipped.
- `ruff check src tests validation` clean; the CI coverage floor (80%) holds.

## Edge cases / error handling

- **Untyped panel** -> no-op (the flag has no effect without a lineage map; documented).
- **Region span with too few informative columns** -> no clade clears `min_sites`, region
  unchanged (avoids re-labelling on noise).
- **Tie within margin** -> unchanged (the margin is the guard against churn).
- **Re-attributing to the backbone** -> forbidden by construction (`major_clade` excluded).
- **Donor-absent regions** (`donor_absent=True`) -> still scored; if every clade is a poor
  match the margin will not be met and the region stays flagged absent.

## Scope / non-goals

- This does **not** change detection, the backbone, support, q-values, breakpoints, or
  agreement -- only the donor label and its `mean_sim_minor`.
- It does **not** become a default in this cycle; default-when-typed is a separate follow-up
  gated on a clean harness, mirroring how MaxChi/Bootscan were promoted.
- It does not address `flu_h3n2_ha` (a short-segment *backbone* mismatch -- a different gap).

## Verification summary

`pytest -m "not requires_binary"` green incl. the new `test_reattribute.py`; `ruff` clean;
coverage floor holds; and `python validation/run_hybrids.py --compare` shows the `reattribute`
column closing `rsv_a` with zero regressions. Only then is promotion to default considered.
