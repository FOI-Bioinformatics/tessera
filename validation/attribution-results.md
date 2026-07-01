# G1 attribution measurement -- results

The gap analysis (`docs/gap-analysis.md`) prioritised parental **attribution** and recommended
measuring, before any product change, whether the **barcode** caller and a **per-clade
consensus** panel close the two attribution FAILs (`rsv_a`, `flu_h3n2_ha`). This records the
result of running `python validation/run_hybrids.py --compare` over all 24 datasets (18 run,
6 SKIP), each scored under the 2x2 grid {tip, consensus} panel x {default, +barcode}.

The numbers below are a one-off measurement on the current Nextclade trees, not a contract.
The raw per-case table is `validation/data/attribution_compare.tsv` (git-ignored, regenerate
with `--compare`).

## Result

| config | PASS | donor attribution improved vs baseline |
|--------|------|----------------------------------------|
| baseline (tip, default ensemble) | 16/18 | -- |
| +barcode (tip) | 16/18 | 0 |
| consensus (default ensemble) | **14/18** | 2 |
| consensus+barcode | **14/18** | 2 |

`baseline` reproduces the headline 16 PASS / 2 FAIL exactly, so the comparison is anchored.

## What the data shows

**Barcode is inert in this harness.** Adding `barcode` to the ensemble (tip panel) changed no
verdict and improved no case's donor-attribution depth (0 of 18). Where the other four callers
already recover the donor region, the barcode caller's region merges without changing the
reported donor clade; where they miss the sub-clade (`rsv_a`), barcode misses it too. On these
datasets it adds neither a fix nor a regression.

**The consensus panel trades backbone accuracy for donor accuracy, and nets negative.** It
improved donor attribution on two cases -- `rsv_a` (donor `sibling` depth 3 -> `exact` depth 4,
i.e. it recovers `A.D.1.8` exactly) and `marburg` (donor `sibling` -> `exact`, and that case
stays PASS) -- and fixed the `vzv` backbone (`mismatch` -> `exact`). But it **broke the
backbone on two cases that previously passed**: `measles` (7.5%) and `rubella` (9.0%), both
`exact` -> `mismatch` on the major parent, dropping the PASS count 16 -> 14. Collapsing a clade
to a single consensus genome makes the backbone (major-parent) competition lose to an adjacent
clade's consensus at moderate divergence.

**Neither feature closes the two FAILs.**
- `rsv_a` (A.1 x A.D.1.8, 6.6%): the consensus panel makes the **donor** exact, but the
  **backbone** degrades to a sibling A-sub-clade, so it still FAILs -- the failure moves from
  donor to backbone rather than resolving.
- `flu_h3n2_ha` (C.1 x K, 6.8%): the donor (`K`) is already `exact` in every config; the
  backbone (`C.1`) is `mismatch` in every config. Neither feature moves it.

## Decision (data-gated)

The Phase 2 gates from the implementation plan are **not met**:

- **Do not promote barcode to a default** (Gate A): it shows no measurable benefit here.
- **Do not default `--pool-consensus`** (Gate B): it is a net regression (-2 verdicts) on this
  set, because it degrades backbone attribution at moderate divergence.

This is the value of measuring first: both promotions would have been net-neutral or net-harmful.

## Refined hypothesis (for a separate future cycle, not built here)

The consensus panel's effect is asymmetric: it helps the **donor** (minor parent) and hurts the
**backbone** (major parent). A targeted variant -- use denoised per-clade consensus only for
*donor* attribution while keeping individual tip genomes for the major parent -- might capture
the `rsv_a` / `marburg` donor gains without the `measles` / `rubella` backbone regression. That
is a new design and would get its own brainstorm -> spec -> plan, with the same `--compare`
measurement as its gate. `flu_h3n2_ha`'s backbone mismatch on a short HA segment is a separate
problem and likely needs segment-aware handling, not a panel change.

## Caller-level donor re-attribution (`--reattribute-donors`)

The asymmetry above motivated a caller-level fix (design + plan under
`docs/superpowers/`): a post-merge step that re-labels each region's **donor** to the clade
whose denoised consensus best matches the query over that region, margin-guarded, leaving
detection and the **backbone** untouched. It is opt-in (`--reattribute-donors`, default off)
and a no-op on untyped panels. Measured as a fifth `--compare` column:

| config | PASS | donor attribution improved vs baseline |
|--------|------|----------------------------------------|
| baseline | 16/18 | -- |
| reattribute | 16/18 | 0 |

**Safe but inert on this harness.** Unlike the consensus panel, re-attribution causes **no
regression** -- `measles` and `rubella` stay PASS, because it competes only the donor and
excludes the backbone clade (a backbone-exclusion bug, caught in final review, was fixed
before this run: the call site now maps the backbone *genome* to its *clade*). But it does
**not** close `rsv_a`: the donor stays `A.D.1.6` (sibling), unchanged from baseline.

**Why -- a representation artifact, not an attribution-mechanism limit.** The true donor
clade `A.D.1.8` is **absent from the panel re-attribution sees**. The harness's
`select_regional` dereplication reduces the tip panel to `A.D.1.4 / A.D.1.6 / A.D.1.9`
(among others) and drops `A.D.1.8`. Re-attribution scores against the clades present in the
aligned panel (`result.rows`), so it cannot recover a clade that was dereplicated away. The
`consensus`-panel column reached `A.D.1.8` only because it builds one consensus per clade
from the *full* source-removed pool (all ~36 `A.D` sub-clades), bypassing the dereplication
-- which is also why it pays the backbone regression elsewhere.

**Conclusion.** The feature is correct, safe, and opt-in-merge-worthy (it would sharpen a
donor when the true clade *is* represented but a noisy adjacent genome wins a few windows),
but on this harness `rsv_a` is fundamentally a **G2 representation** problem: the true donor
sub-clade is not in the panel. Per the plan's gate, do **not** promote to default-when-typed
(it did not flip `rsv_a`). The real levers for `rsv_a` are G2: keep/recruit the true donor
sub-clade in the panel (don't dereplicate it away), or score the donor against per-clade
consensuses built from a fuller clade set than the dereplicated panel exposes. Both are
separate cycles. `flu_h3n2_ha` remains a separate short-segment backbone problem.

## Caveat on the product `--pool-consensus`

The harness builds the consensus with the two source genomes removed (to keep the test honest),
so it is a proxy for, not identical to, the shipped `--pool-consensus` (which includes all
tips). The backbone-degradation-at-moderate-divergence effect is a property of one-genome-per-
clade competition and may apply to the product feature too; `--pool-consensus` is best treated
as a tool for low-divergence, near-identical panels rather than an always-on default.

## Lineage-aware panel selection (default-when-typed)

Follow-up measurement of the G2 lever the section above called for: reduce the panel **by
lineage** instead of by clade-blind ANI. When the references carry lineage labels,
`select_regional` now keeps one query-closest representative per non-recombinant lineage
(recombinant CRF/URF/X lineages excluded unless `--keep-recombinant-lineages`); untyped
genomes still fall back to ANI dereplication. The harness types its tip pool by tree clade
and passes it in. Re-ran `--compare` over all 24 datasets (18 run, 6 SKIP).

| config | PASS | vs prior baseline |
|--------|------|-------------------|
| baseline (tip, lineage-selected) | 16/18 | unchanged (16/18) |

**No change to the headline count, but the `rsv_a` failure moved.** The prior baseline failed
`rsv_a` on the **donor** (dereplication dropped the true donor sub-clade `A.D.1.8`, so it was
scored `A.D.1.6`, a sibling). Lineage selection keeps `A.D.1.8` as its own lineage, so the
**donor is now exact** -- the G2 representation gap this document diagnosed is closed on the
donor side. But the same one-representative-per-lineage rule, applied across RSV-A's ~45 fine
sub-lineages, preserves a near-sibling of the query backbone that ANI dereplication used to
collapse, so the **backbone** attribution regresses `exact -> sibling`. Net: `rsv_a` still
FAILs, with the failure relocated from donor to backbone (`bb sibling  don exact`).

The other cases are unchanged: `flu_h3n2_ha` still FAILs on its short-segment backbone
mismatch; `hiv1` was already PASS on the main baseline (the "masking CRF" regression belonged
to the abandoned clade-aware-dereplication branch, not to main), and stays PASS; the remaining
15 hold. Spot-checked `rsv_a` / `hiv1` / `dengue` again after the `select_regional`
simplification (single up-front query-ANI, one reduction switch) -- verdicts identical.

**Conclusion.** Lineage selection is the correct primitive for donor **representation**: it
recovers a true parent clade that ANI dereplication discards, exactly the fix the section above
asked for. It is merged **default-when-typed** because it is net-neutral on the pass count,
closes the donor-representation gap, and is the consumer of the `type-lineages` output plus the
requested `--keep-recombinant-lineages` toggle. It does **not** by itself close `rsv_a`: on a
fine-grained sub-lineage system, keeping one representative per lineage must be paired with
suppression of the query's whole-genome backbone siblings, or the recovered donor is offset by
a backbone-sibling regression. That pairing is a separate cycle.
