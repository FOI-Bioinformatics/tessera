# Issue #90, Option A: the pendant-branch statistic -- results

Date: 2026-10-29. Measured on `main` at `fa6ed28` (Tessera 1.3.0).

The design note
[`docs/superpowers/specs/2026-10-29-issue-90-donor-absent-design.md`](../docs/superpowers/specs/2026-10-29-issue-90-donor-absent-design.md)
recommends exactly one thing and schedules nothing: **measure Option A before building any
of it.** This is that measurement. No package code changed, and nothing here is proposed for
merge.

Reproduce with:

```
export PATH="$PATH:/path/to/aligner-env/bin"    # append; prepending swaps the Python
python validation/probe_donor_absent.py                      # the real donor-absent runs
python validation/measure_pendant.py --out validation/data/pendant.tsv
```

Every number below is printed by `measure_pendant.py`. The per-region table is
`validation/data/pendant.tsv` (git-ignored, regenerate as above).

## Verdict

**Option A is not falsified, and it is not a fix.** It is the first signal measured on this
issue that puts donor-absent cases materially above donor-present ones, and it clearly beats
the identity discriminator the issue already eliminated, on the same cases. But no threshold
on it both caveats every donor-absent case and leaves every correct attribution alone, and
the design note's central claim for it -- that calibrating against the panel's own pendant
lengths makes the statistic "scale-free across pathogens and divergences" -- **is not
supported by the measurement.**

So it is recorded here as a **surviving candidate with a measured ceiling**, not as a fifth
eliminated discriminator and not as work to schedule. The design note's stated prior ("the
prior is not favourable", because RSV's sub-clades are nearly equidistant) was wrong: the
statistic does carry information that identity does not.

## What was measured, and against what

| | |
|---|---|
| donor-present | 20 cases, full panels, **every one attributed correctly** |
| donor-absent | **17 real runs** -- panel re-recruited from the surviving lineages and re-aligned |
| control | the issue's discriminator 1, identity over the tract, re-run on the same 37 cases |
| statistic | the query tract's NJ pendant branch over the region's columns, against the panel's own pendant lengths |

The harness has one donor-absent case, and one case calibrates nothing. The other 16 were
built for real by `validation/probe_donor_absent.py`, which reuses each positive case's
dataset and parent pair with `case_type: panel_donor_absent`, removing the donor's parent
sub-clade where clade names are hierarchical (`A.D.1.8` -> remove `A.D.1`, as the existing
case does) and the donor's own clade where naming is flat (`DENV4`, `B3`, `denovo_3`). It is a
probe: `run_hybrids.HYBRIDS` is untouched and the must-pass headline is unchanged.

**The probe is validated by reproducing the case that already existed**, and that replicate is
then excluded rather than counted twice. `absent_rsv_a` lands on the same panel (9 members),
the same region (4504 bp), the same named donor (`A.D.5`) and the same statistic (3.2231) as
`donorabsent_rsv`. Only `donorabsent_rsv` is counted in the 17.

**A cheap simulation was rejected, with the reason measured.** Dropping the donor's rows from
an existing alignment needs no realignment, but it is biased rather than merely approximate:
a real run re-recruits a panel from the lineages that remain, while dropping rows just leaves
a *sparser* panel, whose own pendant lengths are longer -- and those are the statistic's
denominator. On `rsv_a` the simulation reads **1.02** where the real run reads **3.22**. Both
simulated variants are still reported by the script, separately and labelled, so the proxy is
never read as the measurement.

The neighbour joining is plain Saitou-Nei in numpy -- no tree library, no binary, per the
dependency-free convention. It recovers the generating pendant lengths on 400 random additive
trees to 2e-16 and reproduces the published worked example exactly.

## Result

Separation of the 17 real donor-absent runs from the 20 donor-present cases. AUC is the
chance a random donor-absent value exceeds a random donor-present one.

| scaling | donor-present median | donor-absent median | AUC | best cut | errors |
|---|---|---|---|---|---|
| identity over the tract **(control)** | -0.9917 | -0.9615 | 0.791 | -0.9697 | 7/37 |
| `quantile` | 0.2071 | 0.8750 | 0.882 | 0.6250 | 6/37 |
| `ratio_max` | 0.1752 | 0.7009 | 0.926 | 0.4770 | 5/37 |
| **`ratio_median`** | 0.4684 | 1.5853 | **0.932** | 0.8777 | **4/37** |
| `robust_z` | -0.9271 | 1.3688 | 0.912 | -0.1030 | 5/37 |

The two operating points a decision actually turns on, for `ratio_median`:

* to caveat **every** donor-absent case: **9 of 20** donor-present cases falsely caveated;
* to caveat **none** wrongly: **47%** of donor-absent cases caught.

9 of the 17 donor-absent values fall inside the donor-present range. The highest
donor-present `ratio_median` is 1.6436, and the nine runs at or below it cannot be flagged
without caveating a correct attribution:

| flagged (8) | missed (9) |
|---|---|
| `absent_measles` 10.86, `absent_flu_h3n2_ha` 4.03, `absent_sars_cov_2` 4.00, `donorabsent_rsv` 3.22, `absent_enterovirus_d68` 2.93, `absent_mpox` 2.52, `absent_ebola` 2.51, `absent_iav_h5_ha` 2.40 | `absent_prrsv2` 1.59, `absent_yellow_fever` 1.55, `absent_mumps` 1.16, `absent_rubella` 1.15, `absent_hantavirus` 0.98, `absent_vzv` 0.97, `absent_hepatitis_a` 0.89, `absent_hiv1` 0.88, `absent_hmpv` 0.59 |

**The control is the point.** Identity sits at AUC 0.791 with **12 of 17** donor-absent values
inside the donor-present range; catching them all would falsely caveat 15 of 20 correct
attributions, and catching none wrongly reaches 29%. `absent_vzv` (identity 0.9993) and
`absent_sars_cov_2` (identity 0.9994) carry the same identity and opposite pendant verdicts.
So the pendant branch is not discriminator 2 in another form -- it is evidence of a different
kind, as the issue required. It is simply not a sufficient one.

## Why it is not scale-free, which is the substantive finding

The design note's mechanism for the statistic being comparable across pathogens is that the
panel calibrates itself. That is testable: a scaling that self-calibrates carries no residual
correlation with the panel's geometry *within* either group. Spearman rho against the panel's
spread share (its own median pendant over its median pairwise distance):

| scaling | donor-present | donor-absent | \|sum\| |
|---|---|---|---|
| `ratio_max` | +0.041 | -0.309 | **0.349** |
| `ratio_median` | -0.105 | **-0.654** | 0.760 |
| `robust_z` | -0.364 | -0.483 | 0.847 |
| `quantile` | -0.400 | -0.458 | 0.858 |
| identity (control) | +0.526 | +0.583 | 1.110 |

Panel geometry leaks straight into the best-separating scaling, and no single global
threshold can absorb that. Splitting the cases at the median spread share -- chosen from panel
geometry alone, before any verdict -- confirms it, and AUC needs no threshold, so neither
column can be overfitted by picking one:

| scaling | pooled | dense panels | sparse panels |
|---|---|---|---|
| `quantile` | 0.882 | 0.869 | 0.925 |
| `ratio_max` | 0.926 | 0.917 | 0.938 |
| `ratio_median` | 0.932 | **0.940** | **0.975** |
| `robust_z` | 0.912 | 0.917 | 0.950 |
| identity (control) | 0.791 | 0.750 | 0.825 |

`ratio_median` separates better within **either** density stratum than it does pooled -- one
error in 19 on dense panels at a cut of 2.51, one in 18 on sparse ones at a cut of 0.88.
Pooling is what costs the separation, because the two strata want cuts a factor of three
apart. Only `ratio_median` and `robust_z` improve in both strata; `quantile`, `ratio_max` and
identity all get *worse* on dense panels than pooled, so the sparse stratum being easier for
everything is part of the story and should not be read as a calibration win.

**The real tension, which the de-duplicated numbers sharpen rather than resolve:**
`ratio_max` separates almost as well as `ratio_median` (0.926 against 0.932) and is far more
geometry-stable (0.349 against 0.760) -- but it **misses the one case the issue is actually
about.** `donorabsent_rsv` sits at `ratio_max` 0.2967 inside a donor-present range of
[0.0000, 0.6806]. The scaling that generalises best across these 37 cases is the one that
fails on the motivating case, which is a warning that the motivating case may not be
representative of donor-absence in general.

## Where the statistic does not exist at all

Option A calibrates against the panel's own pendant lengths, so a panel of two or three
members has no distribution to calibrate against and the statistic is undefined. That is a
limit on where it can be used, not a number that came out badly.

* **15 cases** carry no statistic on at least one called region, all for this reason: panels
  of 2--3 members (the dengue, WNV, zika, chikv, CCHFV, hantavirus and marburg families).
* `absent_zika` could not be **scanned** at all: removing the donor clade left one reference
  besides the query, and detection needs two donors to switch between. The same limitation,
  one step earlier.
* `lowdiv_rsv` admits no donor-absent variant: removing the donor's parent level (`A.D`)
  would take its backbone (`A.D.1`) with it.
* So of 25 eligible positives, 17 yielded a measurable donor-absent statistic.

On roughly a third of this harness, removing a donor lineage leaves either no panel to detect
against or none to calibrate against. Any Option A implementation has to say "not applicable"
there, and that is a large fraction of the small-panel, high-divergence cases where an
unsampled parent is most likely in the first place.

## What this changes for the design note's three options

1. **Option A: do not build it on this evidence, and do not bury it either.** It is a real
   signal -- AUC 0.932 against identity's 0.791 on the same cases -- with a measured ceiling
   of 47% sensitivity at zero false caveats, undefined on a third of the harness, and a
   geometry dependence the note assumed away. If anyone takes it further, the measurement
   names the next question precisely: a normalisation that removes the panel-spread-share
   dependence, measured the same way, on the same 37 cases. Not a threshold tuned on the
   RSV case -- which, as the `ratio_max` result shows, is exactly how one would pick the
   scaling that generalises worst.
2. **Option B is strengthened, and is now the thing the measurement supports.** A signal good
   for 47% sensitivity at zero false alarms cannot honestly set a binary `donor_absent` flag,
   but it is about the right strength for the weaker, true statement Option B proposes --
   that the sub-clade is not identifiable and which lineages are tied. The note already
   recommended Option B regardless of Option A's outcome; this result does not change that
   and gives it a second, independent reason.
3. **Option C remains the maintainer's call**, unchanged. The measurement removes one
   argument for it: the information is not simply absent from the alignment, because a
   statistic that reads it separates these cases at AUC 0.932. What is missing is a
   calibration, not the signal.

## Honest notes on strength of evidence

* `n = 17` donor-absent and `n = 20` donor-present. The strata in the stratified table are
  `n = 19` and `n = 18`. These are small; the AUCs carry wide intervals, and the per-stratum
  cuts are fitted and must not be read as thresholds.
* An early partial run, with only 4 donor-absent cases finished, showed *perfect* separation
  (AUC 1.000, zero errors). That did not survive the full set. It is recorded here because it
  is exactly the small-sample result that a threshold tuned on one or two cases would have
  been built on.
* The first pass of this measurement counted `absent_rsv_a` and `donorabsent_rsv` as two
  cases. They are one. The figures above exclude the replicate; the earlier, inflated pass
  read AUC 0.936 on `ratio_median` and 0.908 on `ratio_max`, so de-duplicating moved
  `ratio_max` by nearly two points and narrowed the gap between the two scalings
  considerably.
* The donor-absent runs re-recruit their panels, so the two groups' panels are not matched
  pairs. Panel sizes span 4--14 in both groups, and the geometry analysis above addresses the
  confound directly rather than assuming it away.
* The probe removes a donor lineage by clade label. Where a true parent's lineage is
  unsampled in reality, it is unsampled by accident of surveillance, not by label -- a
  cleaner removal than reality provides.
