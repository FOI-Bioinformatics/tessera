# Tessera gap analysis

> **Status (2026-07-01, updated after this analysis).** This is the original analysis that set
> the attribution/representation agenda; its headline below ("16 PASS, 2 FAIL") is the state at
> the time of writing. The two attribution FAILs it identifies (`rsv_a`, `flu_h3n2_ha`) have
> since been **closed** -- the harness baseline is now **18/18** -- by lineage-aware panel
> selection and, decisively, deriving the HMM major parent by plurality best-match. See
> `validation/attribution-results.md` for the current, measured state. The gaps and reasoning
> below are retained as the record of how that agenda was framed.

A review of where Tessera underperforms, grounded in its testing results and the published
methods it draws on. The aim is to separate what the tool does well from what it does not,
and to order the remaining work by impact. The language is deliberately cautious: each gap
is tied to specific evidence, and recommendations are starting points, not commitments.

Inputs:

- The synthetic-hybrid harness (`validation/run_hybrids.py`), 24 pathogens, one
  representative dataset each.
- The four curated published-recombination datasets (`validation/datasets.json`): orthopox,
  SARS-CoV-2 XBB, HIV-1 CRF01_AE, norovirus GII.
- The 218-test unit and integration suite.
- The published-method landscape Tessera cites (`docs/detection-methods.md`).

## Headline

On the synthetic harness with the default ensemble (`hmm,3seq,maxchi,bootscan`):

```
16 PASS, 2 FAIL, 6 SKIP (0 errors); donor agreement 17 of 18 running cases
```

The result that organises this whole analysis: **there are no detection failures.** Every
recombinant that runs is detected. The Hudson-Kaplan Rmin lower bound is non-zero for all of
them (23-785 events) across the full divergence range tested, from 33% (dengue x DENV4) down
to 0.5% (mpox) and 0.2% (vzv). Seventeen of the eighteen running cases reach a recombinant
region whose donor is corroborated by more than one independent caller.

Every non-PASS is one of two things, neither of which is a missed event:

- **Attribution** -- the recombination is found but a parent is named as the wrong adjacent
  sub-clade (2 FAILs).
- **Representation** -- the true parent lineage is not in the panel, so attribution cannot be
  tested (6 SKIPs).

So the useful question is not "does Tessera miss recombination?" (it does not, on these
tests) but "does it name the parents correctly, is the panel adequate, and do the numbers
cover enough?"

## What the testing shows: detection, attribution, representation

**Detection is the strength to protect.** Across a 150-fold spread of divergence the parent-
free PHI/Rmin track and the four-caller ensemble both fire. Below roughly 8% column
divergence, base-pair windowing loses contrast and Tessera switches automatically to
informative-site windowing (`docs/detection-methods.md`), which is why mpox (0.5%) and vzv
(0.2%) still detect. Ensemble agreement -- a region independently called by two or more of
HMM, 3SEQ, MaxChi, Bootscan -- is treated as a first-class confidence signal, and restricting
the ensemble to `hmm,3seq` alone drops corroboration from 17/18 to 12/18 without changing the
verdict set. That is direct evidence the added callers buy confidence, not just verdicts.

**Attribution is the visible weakness** -- see G1.

**Representation is the SKIP cause** -- see G2.

## Gap families, prioritized

The families below are ordered by priority. Attribution accuracy is first because it accounts
for both genuine failures and is the area where the relevant machinery already exists but is
unmeasured.

### G1 -- Attribution accuracy at low divergence (top priority)

**Evidence.** The only two genuine FAILs are wrong-sibling-clade attribution at low
divergence:

- `rsv_a` (A.1 backbone x A.D.1.8 donor, 6.6% divergence): the donor region is detected but
  recovered as a different A-sublineage than A.D.1.8.
- `flu_h3n2_ha` (C.1 x K on the HA segment, 6.8%): near-tied subclades; the backbone is
  recovered as G.1.3.

Both *detect* the recombinant region. Both *misname* a parent.

**Root cause.** The four default callers compete individual genomes on pairwise discriminating
sites. When adjacent clades are near-identical, a single genome from the wrong clade can win
windows by chance. This is the exact failure mode that two existing features were built to
resist:

- the **barcode caller** (`src/tessera/recomb/barcode.py`), which attributes by denoised
  per-clade characteristic markers rather than single-genome competition (the RecombinHunt /
  rebar principle); and
- the **per-clade consensus pool** (`--pool-consensus`), which replaces every tree tip with
  one denoised consensus genome per clade, harder to lose to an adjacent-clade sibling.

Neither is in the default ensemble, and -- importantly -- the harness cannot currently show
their benefit. Its design removes the two source genomes and relies on other genomes
representing their clades; where the donor clade has no members left, the barcode caller is
correctly silent and contributes nothing to score. So barcode and pool-consensus are
**unmeasured on the cases they were built for** (typed panels where the donor lineage is
represented).

**Sub-gaps.**

1. Barcode + pool-consensus are not evaluated in the harness on typed, represented-donor
   panels, so their real-world benefit is unquantified and absent from the headline numbers.
2. There is no attribution scoring that credits a correct *sibling-clade* call when divergence
   is below the resolution at which the exact sub-clade is statistically identifiable. (The
   harness already concedes this below 4% with its attribution floor; the 6-7% FAILs sit just
   above it.)
3. Reference-coordinate barcodes from Nextclade internal nodes -- the most tree-faithful form
   -- remain deferred; they compose from existing parts via `--pool-consensus` + `--method
   barcode` without coordinate reconciliation.

**Recommendation (evaluate before building).** Add a harness configuration that keeps the
donor lineage represented and runs `--method hmm,3seq,maxchi,bootscan,barcode` with
`--pool-consensus`, and report those numbers next to the default ensemble. If barcode closes
`rsv_a` / `flu_h3n2_ha` on represented-donor panels, consider promoting it to the default
*when the panel is typed*. Consider a lineage-aware attribution tolerance that scores a
sibling-clade donor as correct below the resolution floor. RIPPLES (parsimony placement)
stays out: it needs the UShER/matUtils toolchain, against the dependency-free design.

### G2 -- Panel composition and missing-donor recovery

**Evidence.** All six SKIPs are representation failures, not detection failures:

- The backbone clade has no panel stand-in once the source is removed (`iav_h5_ha`,
  Am-nonGsGD).
- The donor lineage is absent from the panel (`prrsv2`, L1C on ORF5).
- Fewer than two clades have at least three genomes (`hantavirus`, `oropouche`, `cchfv`).

Beyond the harness, the HIV-1 CRF01_AE A/E mosaic cannot be recovered from the Nextclade
community pool: that pool is saturated with AE-derived second-generation CRFs and contains no
pure-E (env) endpoint, so the true env donor is reported DONOR-ABSENT. The same query
resolves cleanly against a LANL pure-subtype pool.

**Root cause.** Nextclade pools are clade-typing reference trees, not collections of clean
recombinant-parent endpoints, and recruitment optimises query coverage rather than
parental-endpoint representation. Tessera already handles a missing donor honestly --
`donor_absent` in `recombination_regions.tsv`, a coverage caveat banner, `coverage_gaps.tsv`,
and `find-references` / `fill-references` -- but this is reactive: it reports the gap after the
fact rather than recruiting to avoid it.

**Recommendation.** Document composition-over-size explicitly (a small pure-lineage panel
beats a large clade-typing tree for resolving a specific recombinant's parents). Bias
recruitment toward pure-lineage / endpoint references when the goal is attribution. Surface
DONOR-ABSENT and coverage gaps more prominently in the verdict so a representation gap is not
mistaken for a clean call.

### G3 -- Evaluation rigor (trust the numbers)

The harness is a good sensitivity probe but has blind spots that bound how far its numbers can
be trusted:

- **One representative dataset per pathogen.** A single dataset per organism cannot
  characterise within-organism variability.
- **Single-breakpoint synthetic hybrids only.** No multi-breakpoint mosaics, and no real
  labelled benchmark beyond the four curated events.
- **No specificity measurement.** There are no pure non-recombinant negative controls, so the
  harness measures sensitivity (find the planted recombinant) but not the false-positive rate
  (leave a non-recombinant alone). Specificity is the other half of a detection claim.
- **SKIPs hide organism classes.** The six SKIPs silently remove segmented / reassortant
  taxa from the reported denominator.
- **Barcode and pool-consensus are absent from the headline counts** (see G1).

**Recommendation.** Add non-recombinant negative controls and report specificity alongside
sensitivity. Add multi-breakpoint hybrids. Cover the segmented / reassortant organism classes
explicitly rather than letting them fall to SKIP. Publish barcode / pool-consensus columns
next to the default ensemble so their contribution is visible.

### G4 -- Method breadth versus published knowledge

- **No in-tool tree-based confirmation (GARD).** This is by design -- GARD infers its own
  trees and would breach the dependency-free constraint. It is positioned as external
  confirmation (with RDP), not a gap to close inside Tessera.
- **No genome-wide multiple-testing correction across regions.** Benjamini-Hochberg is applied
  within a caller's candidate segments, not across the whole genome's regions; a cross-region
  correction would tighten the false-discovery claim.
- **Reassortment is not modelled distinctly from recombination.** The harness already
  special-cases it (the ebola `/r2021b` re-emergence suffix), which is a sign the underlying
  distinction -- segment exchange in segmented genomes versus within-segment crossover -- is
  conceptual, not implemented. For segmented viruses this matters.
- **PHI is conservative on clean parental panels.** When the panel is tidy parental clades
  around a single hybrid, the genome-wide average dilutes the few incompatibilities, so PHI
  can sit near its permutation floor even when Rmin is clearly non-zero.

**Recommendation.** Document the reassortment-versus-recombination distinction and when each
applies. Consider a cross-region FDR. Keep GARD and RDP as named external confirmation rather
than reimplementing them.

## Strengths to protect

Any change in service of the gaps above should preserve what the testing shows already works:

- Detection robustness across two orders of magnitude of divergence.
- The dependency-free numpy reimplementations (PHI, Rmin, 3SEQ, MaxChi, Bootscan).
- Ensemble agreement as an explicit confidence signal.
- Faithful PASS / FAIL / SKIP reporting, including the scoring rules that make the harness
  well-posed.
- The parent-free PHI/Rmin track, which works precisely where attribution is hardest.

## Prioritized roadmap

Ordered by priority, with each item tagged by the kind of work it is.

1. **G1, evaluate-first.** Add a represented-donor harness configuration and measure
   barcode + `--pool-consensus`. Report against the default ensemble. *(no new dependency)*
2. **G1, build (conditional).** If the evaluation supports it: a lineage-aware attribution
   tolerance, and promotion of barcode to the default ensemble when the panel is typed.
   *(no new dependency)*
3. **G3, build.** Add non-recombinant negative controls and a specificity number; add
   multi-breakpoint hybrids. *(no new dependency)*
4. **G2, document then build.** Composition-over-size guidance; recruitment biased toward
   parental endpoints; more prominent DONOR-ABSENT surfacing. *(no new dependency)*
5. **G4, document.** The reassortment-versus-recombination distinction and GARD/RDP
   positioning; a cross-region FDR as a later build. *(no new dependency for the docs;
   FDR is numpy-only)*

G1 leads because it owns both failures and the machinery already exists. G2 and G3 can proceed
in parallel after the G1 evaluation. G4 is mostly documentation plus one optional numeric
addition. None of the recommendations requires a new runtime dependency; each promoted item
would get its own design and plan.

## Deferred / out of scope

- **RIPPLES** (parsimony placement) -- needs the UShER / matUtils toolchain; against the
  dependency-free design.
- **Reference-coordinate barcodes from Nextclade internal nodes** -- the most tree-faithful
  attribution, but it composes from existing parts (`--pool-consensus` + `--method barcode`)
  without reconciling Nextclade reference coordinates to Tessera's backbone-anchored
  alignment.
- **GARD / RDP** -- kept as external confirmation, not reimplemented.
