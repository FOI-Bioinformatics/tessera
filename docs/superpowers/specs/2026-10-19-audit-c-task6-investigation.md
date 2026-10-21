# Plan C Task 6: investigation of the failed end-to-end gate

Date: 2026-10-19. Investigated the two pools that failed the Task 6 acceptance
criterion ("this branch must not report a present-donor region that main does not
report"). Evidence is the preserved gate output in `e2e2-before/`, `e2e2-after/` and
`e2e-before/`, plus the cached pools. No code was changed and nothing was re-run.

## Summary

The gate failure is largely an artifact of two broken premises, but it also exposed one
genuine defect in the patch. Recommendation: fix that defect, rebase, rewrite the gate,
then re-run. Do not ship as-is and do not revert.

## Finding 1: the mpox baseline never scanned

`e2e-before/nextstrain_mpox_all_clades.log` ends with:

```
[WARNING] Curation left 1 reference(s) in the panel; detection needs at least 2.
```

Selection gave the baseline one genome, the 82-minute NCBI fill added nothing, and the
scan refused. "0 regions before" is therefore an untested absence, not a negative, and
the acceptance criterion cannot be evaluated for mpox as written.

## Finding 2: the mpox locus is a pre-existing alignment artifact

All three mpox regions are nested at one locus, MSA 180512-183402, query
177938-180659, in the right terminal region of a ~197 kb genome. Gap density there
against each sequence's genome-wide rate:

| sequence | role | gap% genome | gap% locus |
|---|---|---|---|
| PP_000T3GC | major parent | 1.64 | 8.65 |
| PP_000TMJP | query | 1.59 | 5.85 |
| PP_004DYJ3 | minor parent (HMM) | 3.00 | 3.36 |

The major parent is five times gappier at the locus than genome-wide while the minor
parent is unchanged. Where the backbone is locally gap-degraded, any gap-stable donor
matches the query better and the HMM switches state. The three method families report
three different minor parents at the same locus, which is the signature of a locus
problem rather than one bad donor. The patch did not create this; it supplied enough
references for the scan to run and so made it visible.

The mpox query PP_000TMJP is labelled `IIb`, a single clade, so it is an unrecombined
tip and these regions are false positives.

## Finding 3: the coverage guard is blind to a degraded major parent

`coverage_gaps.tsv` holds no rows. With a 0.95 best-similarity threshold and roughly
1-2% inter-clade mpox divergence, a locally gap-degraded major parent never trips it.
The `donor_undercovered` flag covers the minor parent only, so all three mpox rows read
`donor_undercovered = no`. Independent of plan C.

## Finding 4: the patch admits two non-clade labels as panel lineages (REAL DEFECT)

Label vocabulary of the cached mpox pool (1728 genomes, 7 distinct labels):

| label | count |
|---|---|
| IIb | 1054 |
| Ia | 315 |
| Ib | 313 |
| IIa | 33 |
| example | 11 |
| Ib/IIb | 1 |
| unassigned | 1 |

`_NO_POOL_LABEL` covers only `""`, `na` and `example`. So the two singletons each become
their own lineage and each contributes a reference:

- `Ib/IIb` is a composite label spanning two clades. Its genome PP_004DYJ3 is the minor
  parent of the 2721 bp HMM region. `is_recombinant_lineage` matches only `^(CRF|URF)\d`,
  `^X[A-Z]` and the literal `recombinant`, so a slash-composite label slips past the
  recombinant-lineage exclusion.
- `unassigned` is the dataset's `reconstructed_ancestral_mpox`, a synthetic
  maximum-likelihood reconstruction, admitted as a candidate donor.

This is the patch's own contribution and must be fixed before it ships. Reading labels
verbatim is correct; promoting every distinct label to a selection bucket is not.

## Finding 5: the HIV-1 gate failure rests on a false premise

The criterion states that each query is an unrecombined tree tip. The HIV-1 query
AY900575 is typed `CRF24_BG` in both runs, a designated circulating recombinant form
(B/G). Regions reported on a recombinant query are not automatically false positives, so
the criterion is inapplicable to this case.

## Finding 6: the HIV-1 baseline is not clean either

| run | regions | detail |
|---|---|---|
| before | 2 | 4 bp (q=1.4e-2), 964 bp (q=3.4e-12); minor parents typed only `HIV-1` |
| after | 3 | 28 bp, 6 bp, 50 bp; the 50 bp one names subtype G but is flagged `donor_undercovered` |

The baseline reports 968 bp of present-donor region on this query with no subtype
resolution; the patched run reports 84 bp, and the one region naming a subtype carries
its own undercovered-donor caveat. The region count rises by one while the extent falls
about 90%. Neither output is a clean false-positive count, because the query is a
recombinant.

## Finding 7: the patch no longer applies cleanly

`git apply --check --3way` applies five of six files cleanly and conflicts on
`src/tessera/discover/iterate.py`, which plan E rewrote. The three target functions
(`_select_from`, `_seed_from_pool`, `_type_panel`) still exist with the same signatures;
the conflict is the variable rename `genomes` -> `kept` at the `_copy_into(_select_from(...))`
call. The rebase is mechanical.

## Recommended next steps

1. Fix finding 4: extend the sentinel set past `""`/`na`/`example` to cover `unassigned`
   (and likely `unclassified`, `unknown`), and stop a composite label containing `/`
   from becoming a distinct selection bucket.
2. Rebase onto current HEAD; only `iterate.py` needs hand work.
3. Rewrite the gate before re-running. Choose queries that are genuinely unrecombined
   (the HIV-1 one was not) and require the baseline to scan (mpox did not), or the
   criterion is untestable.
4. File findings 2 and 3 separately. They are pre-existing defects in low-divergence
   coverage handling, not plan C's to fix.

## Outcome

Acted on the same day, in two commits on `fix-audit-c-task6`.

Task 6 was rebased onto the plan D-G work; the one conflict was in `discover/iterate.py`,
where plan E had renamed the post-exclusion pool variable, and plan E's exclusion logic is
unchanged. The label-vocabulary defect of finding 4 was then fixed: `unassigned` and
`unclassified` join the sentinels that mean "no clade", and a label whose slash-separated
parts are all clades of the same pool is not treated as a lineage. Over the 33 cached
pools the two rules drop 221 of 34027 typed genomes and no legitimate clade; enterovirus
D68's `A2/D` and Ebola's `Ebov-2013/r2021` both survive, because neither names its parts
as standalone clades.

Measured on the mpox dataset with the same query the failing gate used, PP_000TMJP of
clade IIb:

| | panel | present-donor regions |
|---|---|---|
| Task 6 as written | 6 labels, including `Ib/IIb` and `unassigned` | 3 |
| Task 6 plus the fix | 4 labels, the four mpox clades | 2 |

The region the fix removes is exactly the one the investigation attributed to the label
defect: the 2721 bp HMM call whose minor parent was the `Ib/IIb` genome, which is no
longer selected. The two that remain are the locus artifact of findings 2 and 3, with
different minor parents, as predicted.

Findings 2 and 3 are filed as issues #81 and #82 and are not addressed here. The gate's
acceptance criterion in the plan has been corrected on both premises, and the query
chooser it now uses is `validation/pick_pool_query.py`.

## The HIV-1 gate was abandoned

The paired HIV-1 run, using a pure subtype A1 tip in place of the recombinant the original
gate picked, was killed after about 90 minutes while still on the baseline side. Round 1
found a 1947 bp coverage gap at query 4990-6937, the fill split it into nine remote
`blastn/nt` submissions, and NCBI returned its ten-minute warning on each of the three
issued. Neither side reached the scan, so there is no region table to compare.

The fix is unaffected: all 1052 genomes of that pool type identically with and without the
two new rules, because the pool carries no sentinel and no slash label, so the mpox
comparison above is the fix's whole gate. What remains untested is Task 6's panel change
on an unrecombined HIV-1 query. That is an untested absence rather than a negative, the
same distinction the corrected criterion now requires, and it must not be read as a pass.

One observation survives from the baseline side: `Lineage selection: 0 representative(s)
from 0 lineage(s); 412 untyped genome(s) -> 166 after ANI dereplication`. Zero lineages
typed from 1052 labelled tips is Task 6's premise, confirmed by the CLI rather than by
reading the pool offline.
