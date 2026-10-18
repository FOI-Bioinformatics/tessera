# Decision record: false regions at the default agreement gate

**Status:** open -- for the maintainer. **Date opened:** 2026-10-01. Numbers below were
measured on 2026-10-02 on the audit plan C branch unless a row says otherwise.
**Question:** at the CLI default (`--min-methods 1`) the specificity harness reports a
false region in roughly four of ten clonal runs. Should the default change, should the
HMM's significance gate change, or should the default stay and be documented?

## What was measured

All on simulated clonal data (`validation/run_specificity.py`: four clades about 16 %
apart, four tips each about 2.4 % apart, a clonal query in clade A). Every region is a
false positive by construction. Replicate counts are small; the intervals are wide.

| setting (10 replicates, 40 clonal runs) | runs with a false region | false regions by caller | positive control (3 kb tract, 10 replicates) |
|---|---|---|---|
| `--min-methods 1` (CLI default) | 17/40 (42 %, CI 29-58 %) | hmm = 16, maxchi = 4, 3seq = 2 | 10/10 detected, 10/10 donor, 9 bp |
| `--min-methods 2` (harness default) | 2/40 (5 %, CI 1-17 %) | 3seq = 2, maxchi = 2 (the two agree) | 9/10 detected, 9/10 donor, 10 bp |

Source files: `validation/data/audit-c/task9/spec_mm1.txt`, `spec_mm2.txt`.

**These figures predate the fourth audit's caller-statistics branch**, which lowered them:
putting the shallow donors into the triplet callers' Benjamini-Hochberg families moved
`--min-methods 1` to 16/40 (19 regions; hmm = 16, maxchi = 3, 3seq = 1) and
`--min-methods 2` to 1/40, measured the same way at 10 replicates
(`validation/data/audit-g/final/spec_mm1.txt`, `spec_mm2.txt`). The decision below was
taken on the rows above; the direction of every argument is unchanged, and the gate is
still not free of false regions.

The agreement gate is therefore not free of false regions at this replicate count (2/40),
and it lost one of ten positive-control detections.

On near-identical panels (`validation/run_regimes.py`, references about 0.2 % apart),
after lineage clustering stopped pooling the whole panel (audit item C1):

| setting | clonal runs with a false region (20 replicates) | tract detected (10 replicates) | HMM among the callers |
|---|---|---|---|
| `--min-methods 1` | 1/20 (hmm; CI 1-24 %) | 10/10 | 10/10 |
| `--min-methods 2` | 0/20 (CI 0-16 %) | 10/10 | 10/10 |

Source files: `validation/data/audit-c/task9/near_clonal.txt`, `regimes_mm1.txt`,
`regimes_mm2.txt`.

On the hybrid harness (`validation/run_hybrids.py`, `--min-methods 1`, same branch):
39/41 cases pass, sensitivity 29/30, specificity 10/10, with no negative control
reporting a region. One case changed with item C1: in `sars_cov_2` (0.4 % divergent) the
HMM no longer joins the 3SEQ call of the tract, which is now a single-method region.
Source file: `validation/data/audit-c/task9/hybrids.txt`.

## What the false regions are

From `validation/characterise_false_regions.py` (10 replicates, 40 clonal alignments,
22 single-caller rows in 17 alignments, listed before the ensemble merge;
`validation/data/audit-c/task10/false_regions.tsv`):

- 16 of 22 are HMM calls, 4 are MaxChi and 2 are 3SEQ. At 3 replicates all eight were
  HMM calls; the site callers' false regions appear at 10.
- 19 of 22 are switches between two tips of the query's own clade, about 2.4 % apart:
  all 16 HMM rows, 2 of 4 MaxChi rows and 1 of 2 3SEQ rows. The three cross-clade rows
  are short site-caller calls (45 and 47 bp, on 6 distinguishing sites).
- none of the panels was pooled by lineage clustering (the tips sit below its 98.5 %
  floor), so the tips competed as individual genomes;
- HMM rows: lengths 0.9-4.3 kb, resting on 17-101 distinguishing sites; over all rows
  the median length is 1543 bp and the median number of distinguishing sites is 42;
- 14 of 22 have a q-value above 0.01 (all between 0.01 and 0.05): 11 of the 16 HMM rows,
  2 of 4 MaxChi rows and 1 of 2 3SEQ rows.

The likely mechanism for the HMM rows: the HMM picks a segment because the donor wins
there, and the sign test is then run on the same sites. It is not aware that a segment
was searched for. 3SEQ and MaxChi test against nulls that account for the scan, and
produce fewer false regions, though not none.

## Options

### A. Change the default to `--min-methods 2`

- For: 2/40 and 0/20 runs with a false region here, against 17/40 and 1/20.
- Against: this was tried and reverted (PR #47). On the hybrid harness it lost three
  true detections (`rsv_a`, `mpox`, `masksib_rsv`) and gained nothing, because every
  negative control already passed. At very low divergence only the HMM has power, so a
  gate of two is not reachable there. On this branch the `sars_cov_2` hybrid tract is
  also a single-method region, and the positive control fell from 10/10 to 9/10.
- Evidence still needed: `run_hybrids.py` at `--min-methods 2` on the current branch.

### B. Keep `--min-methods 1`; make the HMM's gate aware of the scan

One concrete form was measured, not implemented: an HMM segment must pass both the sign
test and the 3SEQ max-descent test for its (major, donor) pair, taking the larger
p-value. These numbers were measured on 2026-10-01 on a scratch copy that is not part of
this branch, and were not repeated.

| measurement (3 or 10 replicates as stated) | current gate | scan-aware gate |
|---|---|---|
| `run_specificity.py --reps 3 --min-methods 1`: runs with a false region | 7/12 | 2/12 |
| positive control | 3/3, 3/3, 13 bp | 3/3, 3/3, 13 bp |
| near-identical, 10 reps: HMM among the callers | 10/10 | 10/10 |
| near-identical clonal, 10 reps: runs with a false region | 1/10 | 1/10 |
| `run_validation.py` verdicts | 6 PASS, 1 FAIL, 1 SKIP | the same |
| `hiv1_crf`, AE_env region 5569-8177 | `hmm,bootscan` | `bootscan` (the HMM vote is lost) |
| p-value of the four strongest regions | 1e-32 to 1e-67 | 5e-05 (permutation floor) |

- For: removes about two thirds of the false regions without an agreement gate.
- Against: it removed a real HMM call on HIV-1 CRF01_AE; it replaces exact small
  p-values by the permutation floor wherever the triplet is too large for the exact
  test; it was not run on the hybrid harness, where the HMM-only detections live.
- Evidence still needed: `run_hybrids.py`; a design that keeps the exact p-value.

### C. Keep everything; document the rate and recommend `--min-methods 2` for redundant panels

- For: no behaviour change; the README and `--min-methods` help already describe the
  trade.
- Against: a default that reports a false region in about four of ten clonal runs on a
  redundant panel is a weak default for a first-time user, and `hcv_clonal_1b` (a real
  clonal control) currently fails on a single-caller region.
- Needed: state the measured rate in `README.md` and `validation/README.md`, and make
  the harness and the CLI use the same default or say plainly that they do not.

### D. Lower `--alpha` for the HMM only

14 of the 22 false regions have q between 0.01 and 0.05 (11 of the 16 HMM rows).

Measured with `--alpha 0.01` for the whole scan at `--min-methods 1`, 10 replicates
(`validation/data/audit-c/task10/alpha_001.txt`): 10/40 clonal runs with a false region
(25 %, CI 14-40 %), against 17/40 at the default alpha; positive control 10/10 detected.
It was not measured on real data or on the hybrid harness, where weak true regions
(several real HMM calls have q between 0.01 and 0.05) would be lost, and an alpha for
one caller that differs from the others needs a reason beyond "it removes these rows".

## Recommendation from the audit

Not A on its own -- the hybrid harness already answered that, and at 10 replicates the
agreement gate is neither free of false regions nor free of cost on the positive
control. D removes fewer than half of the affected runs. B is the option that addresses
the mechanism, but the form measured here is not ready: it cost a real HMM call and
coarsened p-values. The audit's suggestion is to take C now (it is documentation only)
and treat B as a design task of its own, gated on the hybrid harness.

## Addendum (2026-10-12, audit 4): what a vote is worth

Two findings of the fourth audit bear on this decision. Neither changes a default.

**3SEQ and MaxChi are one signal, not two (audit item A12).** MaxChi computes its
chi-square over the tract 3SEQ's maximum descent finds. On 300 simulated clonal walks
the two p-values have a Spearman correlation of 0.95 (0.88 over the 143 walks MaxChi
permutes, excluding donors with a descent too shallow to permute, which receive p = 1);
at alpha 0.05 both reject in 0.030 of walks against 0.0012 expected under independence,
and MaxChi rejects in 0.90 of the walks where 3SEQ does. Both are individually
calibrated (3SEQ 0.033; MaxChi 0.037 without its depth floor and 0.033 with it; nominal
0.05). The two false regions that remained at
`--min-methods 2` in the table above were exactly such a pair. For option A this means
the agreement gate of two is weaker than it reads whenever the two agreeing callers are
3SEQ and MaxChi; a gate that counted them as one vote is a further option (A'), not
measured here.

**The barcode caller votes without a significance test (audit item A14).** A barcode
region has a support and no p-value, q-value or test, yet counts toward `--min-methods`.
The options, for the maintainer:

- keep it counting (status quo): a typed panel's barcode call corroborates like any
  other caller;
- record barcode in `methods` but do not count it toward the gate until it has a null
  (for example a permutation of marker positions), stated in the `--min-methods` help;
- give it a null first, then count it.

The barcode caller is exercised in the opt-in attribution comparison of
`validation/run_hybrids.py` (the `+barcode` and `consensus+barcode` entries of
`COMPARE_CONFIGS`; see `validation/README.md`), where it is recorded as inert on those
datasets. What no harness measures is its effect on the `--min-methods` agreement gate,
which is what this decision needs.

## Decision

_To be filled in by the maintainer: option chosen, date, and the harness numbers it was
taken on._
