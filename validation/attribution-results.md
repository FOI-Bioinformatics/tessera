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

## Follow-up: the mixed panel (tips + per-clade consensus)

The consensus result above is asymmetric -- it helps the donor and hurts the backbone -- so a
natural refinement is a **mixed** panel: keep the individual tips (so the backbone keeps its
genomes) *and* add one consensus genome per clade (so the donor region can be won by a clade
consensus). The harness supports it as a fifth `--compare` config (`--method` default ensemble,
`panel_mode="mixed"`). Measured on the decisive cases:

| case | baseline | consensus | mixed |
|------|----------|-----------|-------|
| `measles` (7.5%) | PASS | FAIL (backbone) | **PASS** (backbone recovered) |
| `rubella` (9.0%) | PASS | FAIL (backbone) | **FAIL** (backbone still mismatched) |
| `rsv_a` (6.6%) | FAIL (donor) | FAIL (backbone) | **FAIL** (backbone sibling, donor exact) |
| `marburg` (21.5%) | PASS | PASS | not completed (slow full genome; not decisive) |

**The mixed panel is not a clean win either.** It removed the `measles` backbone regression but
not `rubella`'s, and it did not close `rsv_a`: the consensus genome wins the donor region
(donor becomes `exact`) but *also* outcompetes the backbone's own tips for the major-parent
slot, so the backbone still degrades to a sibling sub-clade. Adding a denoised consensus to the
panel makes it a stronger competitor for **both** parent roles, not just the donor.

### What this closes, and the open direction

The panel-composition lever is exhausted for G1: neither barcode, a consensus panel, nor a
mixed panel cleanly improves attribution. The root cause is structural -- a single reference
panel competes for both the major and the minor parent, so any reference good enough to fix the
donor is also good enough to steal the backbone. A real fix has to *separate the
representations*: score the **minor parent (donor)** against denoised per-clade consensuses
while scoring the **major parent (backbone)** against individual genomes. That is a
**caller-level** change (inside the attribution scoring), not a panel-composition change, and it
is the open G1 direction -- deferred to its own brainstorm -> spec -> plan cycle. `flu_h3n2_ha`
remains a separate short-segment backbone problem.

## Caveat on the product `--pool-consensus`

The harness builds the consensus with the two source genomes removed (to keep the test honest),
so it is a proxy for, not identical to, the shipped `--pool-consensus` (which includes all
tips). The backbone-degradation-at-moderate-divergence effect is a property of one-genome-per-
clade competition and may apply to the product feature too; `--pool-consensus` is best treated
as a tool for low-divergence, near-identical panels rather than an always-on default.
