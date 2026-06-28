# RecomFi validation on published data

An opt-in harness that runs the full RecomFi pipeline (`build_msa` then
`run_recomb`) on real recombination datasets and checks the result against a
documented expectation. It is intentionally **separate from the pytest suite**:
CI stays fast and binary-free, while this exercises the real aligner path on
demand.

## Layout

```
validation/
  datasets.json        one entry per dataset (inputs, aligner, expected outcome)
  fetch.py             download the fetch-based datasets from NCBI (efetch)
  run_validation.py    build MSA + run recomb per dataset, check, print a table
  data/                downloaded sequences + run artifacts (gitignored)
```

## Prerequisites

An aligner binary on `PATH`. None ship with conda for every platform, so the
practical setup is a dedicated env:

```
conda create -n recomfi-aln -c bioconda -c conda-forge minimap2 sibeliaz mafft
export PATH="$HOME/miniforge3/envs/recomfi-aln/bin:$PATH"
```

`fetch.py` additionally needs Entrez Direct: `conda install -c bioconda
entrez-direct`.

## Run

```
python validation/fetch.py            # populate data/ for fetch-based datasets
python validation/run_validation.py   # build + scan + check, PASS/FAIL/SKIP table
```

A dataset is `SKIP`ped (not failed) when the aligner is absent or its sequences
have not been fetched, so a partial setup still reports cleanly.

## Datasets

| name | genome | event | aligner | observed result |
|------|--------|-------|---------|-----------------|
| `orthopox_example` | orthopoxvirus ~200 kb | cowpox backbone + variola insert | minimap2 | major cowpox_KC813504; variola region ~MSA 68-150 kb |
| `sarscov2_xbb` | SARS-CoV-2 ~30 kb | XBB = BA.2.10.1 x BA.2.75 | minimap2 | major BA.2.10.1; BA.2.75 region over spike (query ~21.8-26.6 kb) |
| `hiv1_crf` | HIV-1 ~10 kb | CRF01_AE (CM240): A backbone + E env | mafft | major A1; donor-present **AE_env** region over gp120 (~query 5.7-8.2 kb, breakpoint ~5.8 kb). Pure subtype E exists only as the CRF01_AE env, so an env-only E reference (JN388230) is in the panel; without it the env is a donor-absent region instead. |
| `norovirus_gii` | norovirus ~7.5 kb | GII.P16-GII.1, ORF1/ORF2 junction | mafft | major GII.P16-GII.4 (polymerase); GII.1 capsid region from ~nt 4.9 kb |

All four reproduce their published recombination event end-to-end (verified
with `minimap2`/`mafft` installed). Accessions are listed per dataset in
`datasets.json` (`provenance` field) and were confirmed against NCBI nuccore.

### Aligner note (orthopoxvirus)

`sibeliaz` (the default backend) OOMs on the large collinear orthopox blocks on
a typical laptop -- its spoa step is killed and it writes an empty MAF, with a
clear error directing you to more memory or another backend. `orthopox_example`
therefore uses `minimap2`, the assembly-to-reference backend, which suits the
fragmented 8-contig query and runs in seconds.

### One caveat to confirm

Every accession is confirmed on NCBI nuccore except the SARS-CoV-2 **BA.2.10.1**
parent (`OR324670`), which was matched by an NCBI text search rather than a
verified Pango lineage assignment. Confirm its lineage (e.g. with Nextclade or
the NCBI Virus lineage filter) before relying on that dataset for anything
beyond a smoke test; the recombination call (major BA.2.10.1, BA.2.75 over
spike) is consistent with the expected XBB topology regardless.

## Synthetic hybrid harness (`run_hybrids.py`)

A second, generative harness that stress-tests detection on **synthetic
recombinants** built from Nextclade datasets, so the inputs are not limited to
the handful of curated published events above. It carries one representative
dataset per distinct Nextclade pathogen -- **14 core (Nextstrain-maintained) plus
10 community** -- covering SARS-CoV-2, dengue, measles, mumps, rubella, RSV, hMPV,
WNV, yellow fever, ebola, hantavirus, influenza (flu HA), mpox, VZV, HIV-1,
Marburg, Oropouche, Zika, hepatitis A, H5 influenza, CCHFV, chikungunya,
enterovirus D68 and PRRSV. For each dataset it:

1. builds the Nextclade reference pool (the shipped `build_pool`, cached per
   dataset version);
2. picks the most-divergent pair of well-represented, **non-recombinant** clades
   (A and B), each represented by its central genome;
3. splices an A-backbone genome with a B insert over the middle 35-65 % of the
   genome, recording the true donor span in query coordinates;
4. runs RecomFi pool-only with the two exact source genomes removed (their clades
   stay represented), so the query is not a trivial self-match;
5. checks the call: recombination detected, backbone (major parent) clade == A
   (hierarchical labels match, e.g. `A` == `A.1`), a donor region recovered for
   clade B overlapping the true span. Wall-clock runtime is recorded per case.

```
export PATH="$HOME/miniforge3/envs/recomfi-aln/bin:$PATH"
python validation/run_hybrids.py            # all cases
python validation/run_hybrids.py hiv1 dengue   # only named cases
```

Needs MAFFT/skani/skDER on PATH and contacts the Nextclade dataset server on the
first run (pools are cached afterwards under `~/.cache/recomfi/nextclade`).

A dataset is **SKIP**ped (not failed) when it cannot supply a valid test: a
gene-fragment / short-segment dataset that skani rejects (flu HA, yellow-fever
prM-E, H5 HA, PRRSV ORF5), a dataset with fewer than two clades carrying at least
three genomes (hantavirus, Oropouche, CCHFV), or one that is entirely sub-lineages
of a single recombinant (SARS-CoV-2 XBB).

### Observed performance (24 pathogens; one representative dataset each)

| case | backbone x donor | divergence | result |
|------|------------------|-----------:|--------|
| `dengue` | DENV1 x DENV4 | 33.1 % | PASS |
| `marburg` | MARV.B.2 x RAVV.2 | 21.5 % | FAIL -- donor (RAVV) region not recovered |
| `wnv` | 2 x 1B | 20.2 % | FAIL -- backbone clade unlabelled in the tree |
| `hmpv` | B1 x A2.2.1 | 19.0 % | PASS |
| `hepatitis_a` | IIIA x IIA | 16.6 % | FAIL -- backbone/donor inverted |
| `chikv` | III-Asian x I-WestAfrica | 15.5 % | PASS |
| `hiv1` | A1 x B | 15.1 % | PASS |
| `enterovirus_d68` | B3 x A2/D | 11.3 % | PASS |
| `zika` | Asian x African | 10.9 % | FAIL -- backbone clade unlabelled |
| `rubella` | 2B x 1G | 9.0 % | PASS |
| `measles` | H1 x B3 | 7.5 % | FAIL -- a neighbouring genotype wins the backbone |
| `mumps` | A x K | 6.9 % | FAIL -- neighbouring genotype |
| `rsv_a` | A.1 x A.D.1.8 | 6.6 % | FAIL -- backbone recovered, donor sub-clade not |
| `ebola` | Ebov x Ebov | 3.7 % | FAIL -- intra-species, too similar |
| `mpox` | Ib x IIa | 0.5 % | FAIL -- too similar (minimap2, ~200 kb) |
| `vzv` | clade 2 x clade ... | 0.2 % | FAIL -- highly conserved genome |
| `sars_cov_2` | -- | -- | SKIP -- XBB dataset is all one recombinant's sub-lineages |
| `yellow_fever` | -- | -- | SKIP -- prM-E gene region too short for skani |
| `hantavirus` | -- | -- | SKIP -- < 2 clades with >= 3 genomes |
| `flu_h3n2_ha` | -- | -- | SKIP -- HA segment too short for skani |
| `iav_h5_ha` | -- | -- | SKIP -- HA segment too short for skani |
| `cchfv` | -- | -- | SKIP -- < 2 clades with >= 3 genomes |
| `oropouche` | -- | -- | SKIP -- < 2 clades with >= 3 genomes |
| `prrsv2` | -- | -- | SKIP -- ORF5 gene region too short for skani |

**6 PASS, 10 FAIL, 8 SKIP.** RecomFi recovers the recombinant cleanly when the
parental clades are clearly divergent and well represented -- dengue serotypes,
HIV subtypes, chikungunya/hMPV/EV-D68/rubella genotypes (all >= 9 % divergent).
It degrades in two regimes: (1) **low divergence** -- closely-related genotypes or
fine-grained sub-clades (~5-7 %: measles, mumps, RSV-A) and highly conserved
genomes (mpox 0.5 %, VZV 0.2 %, intra-species ebola 3.7 %), where the backbone or
the short donor tract is not separable from a neighbouring clade; and (2) a few
high-divergence **labelling artifacts** (WNV, Zika), where the winning backbone
genome carries no clade label in the Nextclade tree so the exact-match check scores
a miss. The pass set is a performance characterisation, not a fixed contract -- the
clades chosen follow from each dataset's current Nextclade tree.

## Expectation schema (`expected` block)

- `major_parent` -- the top window-winner must equal this label.
- `region_minor_any` -- at least one called region's minor parent is in this set.
- `parents_present` -- every label must appear as a winner or a region parent.
- `min_regions` -- at least this many recombinant regions are called.
- `breakpoint` -- `{coord: "msa"|"query", lo, hi}`; at least one region overlaps
  the band (query coordinates map to the query's own bases).
