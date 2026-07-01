# Plurality best-match major parent for the HMM caller -- design

## Context

On the attribution harness, `rsv_a` (A.1 backbone x A.D.1.8 donor, ~6.6% divergence) FAILs, and a
diagnostic run showed the failure is **not** in panel selection: the true backbone clade A.1
(genome `PP_001CC3R`) is in the panel and is the strict best match in the most windows (203),
yet the HMM caller labels it a *minor/donor* and picks a uniformly-mediocre sibling
(`PP_001NMAF`, clade A.D.5.2, ~0.92 similarity *everywhere*) as the **major parent** across the
whole genome. The model inverted the attribution: it anchored on the equidistant genome and
recast both true parents (A.1 backbone + A.D.1.8 donor) as "donations."

Root cause: `_major_state` (hmm.py:205) derives the major parent as the **most-occupied Viterbi
state** -- a majority vote over the decoded path. A genome that is *best nowhere but decent
everywhere* can occupy the most path positions (its emission never collapses), so it wins the
"major" label over a genome that is *excellent over a large span but poor elsewhere* (a real
dominant parent). This contradicts the module's own definition (`regions.py`: "the query is most
similar to one dataset across the genome -- the major parent") and the heuristic caller, which
already sets `major = most windows won` (`rank_datasets(analysis, 1)[0]`). Only the HMM diverges.

## Goal

Make the HMM caller's major parent the **plurality best-match** genome -- the one that is the
strict best match in the most informative windows -- so a uniformly-mediocre genome can no longer
become the backbone. This corrects both the label and the region topology (the A.1 backbone stops
being called a donation; the A.D.1.8 span becomes the single donor region). Target: `rsv_a` PASS
with no regression on the other cases (a clean 17/18).

## Approach

Plurality best-match (the raw window-wins plurality), not most-occupied HMM state. Chosen over a
guard/repair heuristic because it is the principled definition already stated by the module and
used by the heuristic caller; the trade-off is that the major changes on every HMM call, so a full
harness re-run is a mandatory gate.

## Components

### `plurality_major(result, labels) -> str | None` (new, `recomb/hmm.py`)
Over the informative windows, count strict-argmax wins per label restricted to `labels`; return
the label with the most wins, tie-broken by the highest mean similarity over informative windows;
`None` when there is no informative window. Reuses `_best_per_window` (`recomb/similarity.py:421`)
for the per-window winner. The tiebreak is deterministic (mean similarity, then first label), so
the same panel always yields the same major.

### `segment_query` (`recomb/hmm.py:86`)
Replace `major = _major_state(path, labels)` with `major = plurality_major(result, labels)`. It is
computed over the same `result`/`labels` the HMM competes (the clustered / sibling-excluded `work`
view passed by `_call_regions_hmm`), so the major stays consistent with the segments and the
downstream donor test. `_major_state` becomes unused and is removed (with any test that pins it).

Everything downstream is unchanged: `_call_regions_hmm` defines a donor region as a non-major
segment whose state beats the major on discordant sites, and `reconcile_major` still takes the
HMM's major when the HMM ran.

## Data flow

```
result -> segment_query -> major = plurality_major(work, labels)
       -> _call_regions_hmm (donors tested vs this major) -> reconcile_major -> verdict
```

For `rsv_a`: A.1 (203 window wins) becomes the major; the A.D.1.8 span is the single donor region;
the spurious "A.1-donates-the-backbone" regions over a A.D.5.2 major disappear.

## Testing / verification

### CI unit tests (no binaries)
- `plurality_major`: genome X is the strict best match in most windows while a uniform genome Y
  would occupy the most HMM path -> returns X. Tie -> the higher-mean-similarity label. No
  informative window -> `None`. Single label -> that label.
- `segment_query` on the rsv_a-shaped inversion: X excellent over the first ~70% and poor over the
  last ~30%; Y uniformly mediocre throughout -> major is X (today it would be Y). This test pins
  the fix.
- The existing HMM / region unit tests must still pass unchanged (in the common case, occupancy
  and plurality agree -- a genome that dominates the path also wins the most windows).

### Harness gate (opt-in, not CI) -- mandatory before merge
- Full `python validation/run_hybrids.py --compare`. Accept only if `rsv_a` flips to PASS and the
  other 15 runnable cases hold -> a clean **17/18**. Because the major changes on every HMM call,
  this is the real safety check.
- Confirm `flu_h3n2_ha` is unaffected (its FAIL is a separate short-segment backbone issue), and
  `hiv1` and the low-divergence `mpox`/`vzv` cases do not shift.
- Record the before/after in `validation/attribution-results.md` faithfully -- if any case
  regresses, that is a reported finding, not something to hide, and the merge decision is revisited.

### Static
- `ruff check src tests validation` clean; CI coverage floor (80%) holds.

## Edge cases
- Non-recombinant query closest to one genome everywhere -> plurality = that genome (unchanged).
- All-`nan` window similarities / a single label -> `None` / that label, as today.
- Ties -> mean-similarity tiebreak, then deterministic first-label.

## Scope / non-goals
- Only the HMM caller's major derivation changes. No change to `select_regional` / panel selection,
  the HMM segmentation model, sibling exclusion, the 3SEQ / MaxChi / Bootscan / barcode callers, or
  `reconcile_major`.
- No new runtime dependency.
