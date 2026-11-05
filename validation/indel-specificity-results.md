# Issue #121: stray regions on aligned, indel-bearing clonal data -- results

Date: 2026-11-05. Measured on `main` at `45d771c` (Tessera 1.3.0), MAFFT 7.526 through
`tessera msa --aligner mafft`, `tessera recomb` with the CLI defaults except
`--window-size 600 --window-step 120`. No package code changed.

Reproduce with `python validation/run_indel_specificity.py` (6 replicates per cell, about
four and a half minutes). The queries are clonal, so every region is a false positive.

## What was run

A 6 kb genome and four references; divergence 1, 3, 6 and 12 %; single-base indels at 0,
0.1, 0.4 and 1 % per site per lineage; 6 replicates per cell, 96 runs. The indel rates are
high: 0.4 % is about 24 events along each lineage. They were chosen to find where the
problem starts.

## Result

Runs with at least one region (of 6), at the CLI default (`--min-methods 1`). "Gate 2"
counts runs that also have a region called by two or more callers.

| indel rate | runs | any region | gate 2 |
|---|---|---|---|
| 0 | 24 | 1 | 0 |
| 0.1 % | 24 | 0 | 0 |
| 0.4 % | 24 | 8 | 5 |
| 1 % | 24 | 14 | 8 |

By caller: 3SEQ 18 regions, HMM 15, MaxChi 14. Lengths of the 30 regions: eleven of 2-22 bp
and nineteen of 391-2548 bp; none in between.

An earlier, separate set of 96 runs from the same generator (a different draw of the
sequences) gave the same pattern with lower counts: no region at 0 or 0.1 %, 2 of 24 runs at
0.4 % and 13 of 24 at 1 %; 3 of 96 runs had a multi-caller region, and all 15 short regions
(3-17 bp, plus one of 86 bp) were called by a single caller.

## What this supports

- **The rate depends on the indel rate.** No region in 48 runs at 0 or 0.1 % in the earlier
  set; one in the 24 indel-free runs of the set above, which this harness did not retain, so
  its cause is not known. At 0.4 % and above, a fifth to a half of the runs have one.
- **There are two kinds of region, not one.** Short ones (2-22 bp) come from 3SEQ and
  MaxChi, at 0.4 % and above; in the earlier set all were single-caller. Long ones
  (hundreds of bp to 2.5 kb) come from the HMM, alone or with the other callers.
- **A minimum length does not address the long ones.** Re-scanning the earlier 96 alignments
  with `--min-region 20`, `50` and `100`: runs with a region fell from 15 to 9 for all three
  floors, the 9 being the long regions. A floor removes the short regions and nothing else.
- **`--min-methods 2` removes short regions** (they were single-caller) but not the
  long multi-caller ones. Its cost on true detections is documented elsewhere and is not
  re-measured here.
- Local gap density is higher around the short regions than the genome average in just over
  half of them (4-7 gap columns in a 65-75 column neighbourhood, against 1-2
  expected), which is consistent with clustered indels contributing. It is not the whole
  explanation: the rest sit in ordinary gap density.

## What it does not support

- No rate for real data. One generator (random root, uniform substitutions, single-base
  indels, star-like panel), small cells, one window setting. Real indels are longer and
  rarer in coding sequence.
- No change to a default. Whether a floor is worth its cost on true short tracts was
  checked only against the lengths of the cases in `run_hybrids.py`, where the shortest
  gated tract is one window; the adversarial tier (0.15 to 1 window) is not part of the
  gate and was not run.
- No mechanism. The gap-density observation is a pattern, not a test.

## For a reader

If an aligned set is indel-rich and a region of a few bp to a few tens of bp is reported by
one caller only, it is likely to be of this kind. A region of hundreds of bp or more from
the HMM on such data deserves a look at the alignment around it before it is read as an
event. Neither observation justifies a change to a default on this evidence.

## Addendum: a 50 bp default floor (issue #121)

Measured with `validation/run_indel_sensitivity.py` (6
replicates per cell, 72 recombinants with a donor tract of 20-600 bp, indels at 0.4 %, same
generator and window settings). A tract counts as localised when the regions naming the donor
cover at least half of it; a region that does not come within 100 bp of the tract is stray.
`--min-region 1` stands for "no floor".

| `--min-region` | tracts localised (of 72) | stray regions |
|---|---|---|
| 1 | 38 | 73 |
| 20 | 38 | 38 |
| 50 | 37 | 12 |
| 100 | 33 | 7 |

At 50 bp the stray regions fall by 84 % and one tract is lost: a 40 bp tract at 12 %
divergence. Tracts of 20 bp are not localised at any setting, and of 40 bp in one run of six.
At 100 bp four more tracts are lost (the 80 bp ones), so 50 bp is the larger of the values
tried that cost almost nothing here; 20 bp is nearly free but leaves half the stray regions.
A first run of the same design with a throwaway script (a different draw of the sequences)
gave 33 localised at every floor and 77, 32 and 5 stray regions at none, 20 and 50 bp.

**Side effect.** Over the 72 runs, donor-absent rows went from 4 to 7 as the floor rose to
50 bp (scratch run). Each new row was a real coverage gap of 595-1070 bp in which a stray
short region had sat: with the region removed, the gap is reported as what it is, a stretch
where the query is far from every reference.

**Not shown.** What the floor costs on real genomes. The curated cases of `run_hybrids.py`
(shortest gated tract one window) are unchanged: 39 of 39 gating, sensitivity 29/29,
specificity 10/10, and all 42 case rows equal to `main`.
