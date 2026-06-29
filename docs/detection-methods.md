# Detection methods and output

Tessera computes, in sliding windows across the MSA, the similarity of the query to
each reference (1 = identical, 0 = no similarity). The reference winning the most
windows is the **major parent** (the backbone donor). On top of that windowed
similarity it runs one or more parent-attributed callers, plus a parent-free
diagnostic that runs for every method.

## Ensemble (default)

The callers are complementary: 3SEQ recovers short low-divergence tracts the HMM
segmentation dilutes (the `cryptic_insert` example), while the HMM localizes breakpoints
better on divergent parents (the `divergent` example); MaxChi and Bootscan add two more
independent votes. By default Tessera runs **four** of them
(`--method hmm,3seq,maxchi,bootscan`) and merges their regions into one consensus,
sharing the single similarity scan so the extra callers are cheap. On the synthetic
harness this leaves the verdict set unchanged but corroborates the donor in 17 of 18
cases (vs 12 with hmm,3seq alone), raising confidence without over-calling.

The merge is transparent -- no combined score is invented. Two regions are the same
event when they overlap in query coordinates and name the same donor -- the same minor
parent genome, or, when the references are typed (a lineage map), the same **lineage**,
so two callers that pick different representative genomes of one lineage still count as
the same event and as agreement. The consensus region records exactly **which** callers
found it (the report's *Method(s)* column and a *Method comparison* table) and whether
the parent-free Hudson-Kaplan Rmin signal corroborates it. A region called by more than one method is treated as **higher
confidence** -- agreement is the point of running an ensemble: the union of the callers
raises recall, their agreement raises precision.

The default ensemble's two further callers are **MaxChi** (a chi-square triplet test,
complementary to 3SEQ) and **Bootscan** (a distance + bootstrap method that yields a
support for the closest parent). Select callers with a comma list or `all`: `--method
hmm,3seq,maxchi,bootscan` (default), `--method hmm` (one caller, reproducing the
single-method report), `--method hmm,3seq` (the two fastest), or `--method all` (every
caller, including the legacy heuristic). The same `--method` is available on `tessera
detect` and `tessera
fill-references`. The individual callers are described below.

## HMM caller (`--method hmm`)

The default caller segments the query against the reference panel with a hidden
Markov model (jpHMM-style): each window emits a binomial copying likelihood per
reference, and a single jump rate (`--jump-rate`) penalises switching reference, so
near-identical references do not flip and thin windows cannot drive a call. A segment
is reported as recombinant only when its donor beats the major parent on the
**discordant sites** -- positions where the query matches one candidate parent but
not the other -- by a sign test at level `--alpha`. This is far more discriminating
than an all-sites margin (it recovers subtle breakpoints between near-identical
parents) and does not invent regions from noise.

Each region carries a **support** (the share of distinguishing sites favouring the
donor), a sign-test **p-value** with a Benjamini-Hochberg **q-value**
(false-discovery-rate across the segments tested), and a **breakpoint uncertainty
interval**. A high support with a small q-value is a confident call; a region with
strong directional support but a large q-value (few distinguishing sites, e.g. a
recombination between near-identical lineages) is flagged as marginal rather than
dropped. The legacy `--method heuristic` (margin / merge-gap / min-region) is kept
for comparison.

```
Recombination regions (major parent: cowpox_KC813504):
  Minor parent  Major parent     Query start  Query end  Length(bp)  Sim minor  Sim major  Support  q-value  Breakpoint
  --------------------------------------------------------------------------------------------------------------------
  variola       cowpox_KC813504  66268        147150     80882       0.999      0.977      0.97     2e-300   66768
```

### Low-divergence panels (intra-species sets, DNA viruses)

When the references are nearly identical -- e.g. mpox clades (~0.5 %), VZV (~0.2 %),
within-species ebola -- a fixed base-pair window holds only a handful of
discriminating columns diluted by hundreds of identical ones, so the per-reference
emission contrast collapses and the segmentation loses power. The signal is still
there (low *percentage* divergence over a large genome is still hundreds of variable
sites), just buried. Tessera therefore switches to **informative-site windowing**:
windows span a fixed number of polymorphic columns rather than a fixed base-pair
width, concentrating the signal so the HMM regains contrast. This is automatic when
the references differ at less than ~8 % of columns, and controllable with
`--informative-sites` / `--no-informative-sites` (and `--informative-window` /
`--informative-step`). Breakpoint intervals are necessarily coarser at low divergence
-- you cannot localise a switch more finely than the spacing of the discriminating
sites.

## 3SEQ caller (`--method 3seq`, scan-aware triplet test)

A second caller, complementary to the HMM, after Boni, Posada & Feldman (2007). For
the query against the major parent and each candidate donor it looks only at the
**discriminating sites** (where the two parents differ and the query matches one of
them) and measures the **maximum drawdown** of the resulting +1/-1 walk -- a
sustained run of donor matches inside the backbone, i.e. a mosaic. Its p-value is the
*exact* probability that a random arrangement of the same matches reaches that
drawdown (a dynamic program over the walk depth, no dependency; a vectorised
permutation falls back for very large inputs), Benjamini-Hochberg-corrected across
the donors tested. Because it is purely informative-site based it keeps full power on
near-identical panels -- it detects the mpox clade-I/II recombination where base-pair
windowing finds nothing -- and because the null accounts for scanning every
breakpoint it does not over-call.

## MaxChi caller (`--method maxchi`, chi-square triplet test)

A third triplet caller after Maynard Smith (1992) (the RDP family), on the same
discriminating-site +1/-1 walk as 3SEQ but with an independent statistic: rather than
the maximum *drawdown* it tests whether the **proportion** of donor matches differs
sharply inside a tract versus outside it -- a 2x2 chi-square -- with a scan-aware
permutation null. Drawdown and proportion are sensitive to different mosaics (a
shallow-but-wide tract is weak under one and strong under the other), so MaxChi and 3SEQ
are genuinely independent votes; a region significant under both is corroborated by two
tests, exactly as RDP treats them as separate methods.

## Bootscan caller (`--method bootscan`, distance + bootstrap)

The published method of choice for identifying *which parent* a region came from, after
Salminen (1995) (SimPlot++ 2022). Per window it measures the query's identity to each
candidate parent, then resamples the window's alignment columns with replacement; the
**bootstrap support** of a parent is the fraction of resamples in which it is the query's
closest match. A run of windows where a non-major parent's support clears 70 % and beats
the backbone is a region, carrying that support as a confidence the other callers express
only as a p-value.

## Parent-free recombination signal (PHI + Rmin)

Alongside the parent-attributed regions, every run reports a parent-agnostic
diagnostic that asks only whether the alignment carries recombination at all, with no
candidate parents required -- the regime where the HMM and 3SEQ callers have least to
work with (low divergence, or the true donor absent from the panel). It is built from
the alignment's biallelic informative columns and the four-gamete test (two columns
are incompatible when all four gametes are present, which implies a recombination
between them).

The **PHI test** (Pairwise Homoplasy Index; Bruen, Philippe & Bryant 2006) is
significant when columns near each other on the genome are more compatible than a
random reordering -- the signature of shared local genealogy under recombination --
with a one-sided permutation p-value. **Rmin** (Hudson & Kaplan 1985) is the minimum
number of recombination events the incompatibilities force, with the intervals as
breakpoint candidates. Both are dependency-free.

On the synthetic-hybrid harness, Rmin is non-zero for every recombinant across the
full divergence range (dengue serotypes at 33 % down to the mpox clade-I/II
recombination at 0.5 %), and the PHI test reaches the permutation floor in most
cases. PHI's genome-wide p-value is conservative in Tessera's typical setup, though:
when the panel is clean parental clades around a single hybrid query, the many
clade-defining sites are mutually compatible and dilute the few incompatibilities the
one query introduces, so PHI can stay non-significant even where detection succeeds
(e.g. yellow fever). In that regime Rmin and the per-site PHI **profile** -- which
localizes the signal rather than averaging it away -- are the more informative
parent-free outputs. The diagnostic runs for every `--method`; disable with
`--no-phi`, or widen its window with `--phi-window`.

This all remains an **indicative screen**: the built-in HMM and 3SEQ tests are fast
triplet/segmentation screens, not a full tree-based analysis (such as GARD). Treat
regions as candidates to confirm.

## Output files

| File | Contents |
|---|---|
| `report.html` | Self-contained report: run provenance, the region table, the per-dataset stats, and an embedded interactive plot |
| `recombination_regions.tsv` | Called regions: minor/major parent, start/end in **both MSA columns and query bases**, length, support, mean similarities, the calling `methods`, and `parent_free_support` |
| `recombination_methods.tsv` | Ensemble breakdown (only when several methods run): one row per region with a Y/n per method and the parent-free flag |
| `recombination_profile.tsv` | Parent-free signal: header with the PHI p-value and Rmin, then per-informative-site local incompatibility (the PHI profile) |
| `similarity_windows.tsv` | Full per-window matrix: `msa_position`, `query_position`, `winner`, and one similarity column per dataset |
| `similarity_stats.tsv` | Per-dataset similarity statistics (median, windows above identity thresholds) |
| `window_winners.tsv` | Per-dataset count of windows won (ties included) |
| `coverage_gaps.tsv` | Stretches where even the closest reference is a poor match -- possible missing references |
| `similarity_top{N}.{fmt}` | Static plot of the nearest `--top-n` datasets, called regions shaded |
| `similarity_pair.{fmt}` | Static plot of the major vs leading minor parent, region shaded |

Similarity is computed only over columns where both sequences carry a canonical base
(A/C/G/T); gaps, `N` and IUPAC ambiguity codes are ignored, so an `N` never counts as
a match. A window with no comparable position -- for example an inter-contig gap in a
fragmented query -- is uninformative and reported as `NA` in `similarity_windows.tsv`
(and excluded from the winners, statistics and region calling). Query coordinates are
reported alongside MSA coordinates, so regions need not be mapped back to the query by
hand.

![divergent recombinant: a clean similarity crossover](figures/divergent.png)

*`example_data/divergent.msa.fasta` (parents ~11 % apart). The query tracks `parent_A`
except over the shaded called region, where it switches to `parent_B`. With divergent
parents the HMM caller localizes the event confidently (q ~1e-29); 3SEQ agrees.*

![cryptic insert: the 3SEQ caller recovers a short low-divergence tract](figures/cryptic_3seq.png)

*`example_data/cryptic_insert.msa.fasta` (parents ~1 % apart, an 800 bp insert with
~10 discriminating sites). The dip is shallow and narrow, so the HMM segmentation
finds nothing across a base-pair window; `--method 3seq` pools the discriminating
sites and recovers the shaded region (q ~1e-12).*

## References

The callers and statistics reimplemented here (dependency-free numpy) follow these
sources:

- **HMM segmentation + per-clade consensus (jpHMM-style).** Schultz A-K, Zhang M,
  Leitner T, et al. (2006). A jumping profile hidden Markov model and applications to
  recombination sites in HIV and HCV genomes. *BMC Bioinformatics* 7:265.
- **MaxChi.** Maynard Smith J (1992). Analyzing the mosaic structure of genes. *Journal
  of Molecular Evolution* 34(2):126-129.
- **Bootscan / SimPlot.** Salminen MO, Carr JK, Burke DS, McCutchan FE (1995).
  Identification of breakpoints in intergenotypic recombinants of HIV type 1 by
  bootscanning. *AIDS Research and Human Retroviruses* 11(11):1423-1425. SimPlot++:
  Samson S, Lord E, Makarenkov V (2022). *Bioinformatics* 38(11):3118-3120.
- **3SEQ triplet test.** Boni MF, Posada D, Feldman MW (2007). An exact nonparametric
  method for inferring mosaic structure in sequence triplets. *Genetics* 176(2):1035-1047.
- **PHI test (Pairwise Homoplasy Index).** Bruen TC, Philippe H, Bryant D (2006). A
  simple and robust statistical test for detecting the presence of recombination.
  *Genetics* 172(4):2665-2681.
- **Four-gamete test and Rmin.** Hudson RR, Kaplan NL (1985). Statistical properties of
  the number of recombination events in the history of a sample of DNA sequences.
  *Genetics* 111(1):147-164.
- **Benjamini-Hochberg false-discovery-rate control.** Benjamini Y, Hochberg Y (1995).
  Controlling the false discovery rate: a practical and powerful approach to multiple
  testing. *Journal of the Royal Statistical Society B* 57(1):289-300.
