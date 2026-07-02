# Caveat a confident donor when its region overlaps a coverage gap -- design

## Context

The harder-harness Phase 1 `donorabsent_rsv` case surfaced a genuine attribution limitation.
When the true donor clade is removed from the panel, at moderate divergence (~6.6%) the caller
attributes the donated region to a present, close sibling with high confidence -- and the
"donor may be missing" signal is lost. Diagnosis (`run.py:290-300`): the coverage scan *does*
detect the divergent stretch (a `coverage_gaps` entry, best-similarity 0.938 below the 0.950
adaptive floor), but the bridge that turns coverage gaps into donor-absent regions **drops** a
gap whenever a confident donor-present region overlaps it:

```python
absent = [ r for r in gaps_as_regions(coverage_gaps, result, major_parent)
           if not any(r.msa_start < p.msa_end and p.msa_start < r.msa_end for p in regions) ]
```

The overlapping donor-present region is also *not* flagged `donor_undercovered`, because
`flag_undercovered_regions` compares the region's **mean** similarity (0.979 over 4698 bp) to
the floor -- the 650 bp gap is diluted by the well-matched flanks. So the tool has the signal
(it detected the gap at the window level) and discards it in favour of a confident-but-wrong
sibling attribution.

## Goal

When a coverage gap overlaps a confident donor-present region, **caveat that region**
(`donor_undercovered = True`) instead of silently dropping the gap -- keeping the confident call
but honestly annotating "the true donor for part of this region may be absent; this is the
closest available reference." Non-overlapping gaps still become donor-absent regions, unchanged.

## Approach (chosen)

Caveat-the-region, keep-attribution. Preferred over emitting a second overlapping donor-absent
row (noisier; every legitimate close-donor case would also show a gap) and over redefining
`flag_undercovered_regions` to a window-level test (a broader change to flagging semantics
affecting all cases). This reuses the existing `donor_undercovered` field, which is already
serialized to the regions TSV (`report_text.py:201,217`) and surfaced in the HTML report
(`report_html.py:82,225`).

## Components

### `src/tessera/recomb/run.py` (the reconciliation, ~lines 293-297)
Replace the drop-on-overlap filter with a loop: for each coverage-gap region, find the called
regions it overlaps; if any, set `p.donor_undercovered = True` on each and do not emit the gap;
if none overlap, keep the gap as a donor-absent region (as today).

```python
absent = []
for gap in gaps_as_regions(coverage_gaps, result, major_parent):
    overlapping = [p for p in regions
                   if gap.msa_start < p.msa_end and p.msa_start < gap.msa_end]
    if overlapping:
        for p in overlapping:
            p.donor_undercovered = True  # the true donor for part of this span may be absent
    else:
        absent.append(gap)
regions = sorted(regions + absent, key=lambda r: r.msa_start)
```

`flag_undercovered_regions` (line 282) still runs first, so a region already below the floor by
its mean is flagged as before; this only *adds* the caveat where a sub-region gap overlaps a
region whose mean stayed above the floor.

### `validation/run_hybrids.py` -- `_score_panel_donor_absent`
Accept the caveat as the honest signal: PASS when the true span overlaps a region flagged
`donor_absent=="yes"` **or** `donor_undercovered=="yes"`; FAIL only on a *clean* confident
cross-clade region over the span (present, multi-method, not undercovered, wrong top-level).
This keeps the case a real test (a silent confident mis-attribution still fails) while crediting
the fixed, honest behaviour.

## Data flow

```
coverage_gaps --gaps_as_regions--> per gap:
    overlaps a called region?  yes -> mark that region donor_undercovered=True (caveat)
                               no  -> emit as a donor-absent region
regions TSV donor_undercovered column -> report + harness scorer
```

## Edge cases
- A gap overlapping several regions caveats all of them.
- A region already `donor_undercovered` by its mean (flag_undercovered_regions) stays flagged.
- The 18 positive harness cases have the true donor present -> best-sim high -> no coverage gap
  over the donor region -> no new caveats; their verdicts are unchanged (checked on the gate).
- No coverage gaps (`coverage_gaps` empty) -> loop is a no-op; behaviour identical to today.

## Testing / verification
- **CI unit test (recomb, no binaries)**: a synthetic `regions` + `coverage_gaps` where a gap
  overlaps a called region -> after reconciliation the region has `donor_undercovered=True` and
  no duplicate donor-absent region is emitted; a non-overlapping gap still yields a donor-absent
  region. (Test the reconciliation logic directly, factoring it into a small pure helper if that
  keeps the test clean.)
- **CI unit test (harness)**: `_score_panel_donor_absent` PASSes when the span overlaps a
  `donor_undercovered=="yes"` region; still FAILs on a clean confident cross-clade region.
- **Opt-in harness gate** (`recomfi-aln`): `python validation/run_hybrids.py donorabsent_rsv`
  now PASSes; a full run confirms the 18 positives + `lowdiv_rsv` + `neg_measles` hold (no
  regression). Record in `validation/attribution-results.md`.
- `ruff check src tests validation` clean; CI coverage floor (80%) holds.

## Scope / non-goals
- Only the gap/region reconciliation and the harness scorer change. No change to the coverage
  threshold, the callers, `flag_undercovered_regions`'s mean-based rule, or the report layout
  (the caveat column already exists).
- No new runtime dependency.
