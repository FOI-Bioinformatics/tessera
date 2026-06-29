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
   dataset version), keeping only the clade-labelled reference-tree tips (the
   dataset's example sequences carry no clade and are dropped);
2. picks the most-divergent pair of well-represented, **non-recombinant** clades
   (A and B), each represented by its central genome; a coarser clade attribute
   can be pinned per case (`clade_key`, e.g. SARS-CoV-2 `clade_nextstrain`);
3. splices an A-backbone genome with a B insert over the middle 35-65 % of the
   genome, recording the true donor span in query coordinates;
4. runs RecomFi pool-only with the two exact source genomes removed (their clades
   stay represented), so the query is not a trivial self-match; window sizes adapt
   to the genome length and the aligner is per-case (minimap2 for ~200 kb mpox/VZV).
   Sibling-dropping is **off**: the synthetic pool has no recombinant twin of the
   query, and for a close-parent hybrid the backbone parent is >95 % genome-wide ANI
   to the query (its 70 % backbone dominates the average) and would otherwise be
   discarded as a masking twin -- the documented `--seed-keep-siblings` case. On a
   near-identical panel (mpox/VZV ~0.5 %) dereplication collapses the parent clades
   into one representative, so the panel is rebuilt from one central genome per
   clade; RecomFi then auto-switches to informative-site windowing (the `mode`
   column);
5. checks the call: recombination detected, backbone (major parent) clade == A,
   and a donor region recovered for clade B overlapping the true span. Clade labels
   match hierarchically (`A` == `A.1`); a donor region is also credited to a sibling
   sub-clade of the donor (e.g. Marburg `RAVV.1` for a `RAVV.2` donor) **only** when
   the donor lineage is distinct from the backbone's, so a shared-top-level pair
   (RSV `A.1` / `A.D.1.8`) still requires resolving the exact donor. Runtime is
   recorded per case.

The summary table also carries the **parent-free** recombination signal per case
(`PHI p` and `Rmin`), read from `recombination_profile.tsv`. This is independent of the
backbone/donor verdict: it confirms the dependency-free PHI test and Hudson-Kaplan Rmin
fire on the recombinant hybrids across the whole divergence range -- including the
low-divergence mpox (0.5%) and ebola (3.7%) cases where the parent-attributed callers
have least to work with (PHI p at the permutation floor, Rmin well above zero).

```
export PATH="$HOME/miniforge3/envs/recomfi-aln/bin:$PATH"
python validation/run_hybrids.py            # all cases
python validation/run_hybrids.py hiv1 dengue   # only named cases
```

Needs MAFFT/minimap2/skani/skDER on PATH and contacts the Nextclade dataset server
on the first run (pools are cached afterwards under `~/.cache/recomfi/nextclade`).
For a short gene/segment dataset that skani rejects the panel falls back to one
central genome per clade.

A dataset is **SKIP**ped (not failed) when it cannot supply a valid test: the
most-divergent clade pair is below ~4 % divergence (too few discriminating sites:
mpox, VZV, ebola, SARS-CoV-2 within Omicron), or it has fewer than two clades with
at least three genomes -- including datasets with no clade attribute at all
(hantavirus, Oropouche, CCHFV).

### Observed performance (24 pathogens; one representative dataset each)

| case | backbone x donor | divergence | result |
|------|------------------|-----------:|--------|
| `dengue` | DENV1 x DENV4 | 33.1 % | PASS |
| `marburg` | MARV.B.2 x RAVV.2 | 21.5 % | PASS |
| `yellow_fever` | Clade VII x Clade III | 21.0 % | PASS |
| `iav_h5_ha` | Am-nonGsGD x 2.2.1.1a | 20.3 % | FAIL -- HA segment; no recombination called |
| `wnv` | 2 x 1B | 20.2 % | PASS |
| `hmpv` | B1 x A2.2.1 | 19.0 % | PASS |
| `prrsv2` | L8D x L1C.2 | 18.6 % | FAIL -- ORF5 (~600 bp); donor tract too short |
| `hepatitis_a` | IIIA x IIA | 16.6 % | PASS |
| `chikv` | III-Asian x I-WestAfrica | 15.5 % | PASS |
| `hiv1` | A1 x B | 15.1 % | PASS |
| `enterovirus_d68` | B3 x A2/D | 11.3 % | PASS |
| `zika` | Asian x African | 10.9 % | PASS |
| `rubella` | 2B x 1G | 9.0 % | PASS |
| `flu_h3n2_ha` | K x unassigned | 7.7 % | FAIL -- HA segment; fine subclades |
| `measles` | H1 x B3 | 7.5 % | PASS |
| `mumps` | A x K | 6.9 % | PASS |
| `rsv_a` | A.1 x A.D.1.8 | 6.6 % | FAIL -- backbone recovered, deep donor sub-clade not |
| `ebola` | Ebov x Ebov | 3.7 % | FAIL -- info-site; backbone recovered, deep donor sub-lineage not |
| `mpox` | Ib x IIa | 0.5 % | FAIL -- info-site; recombination + donor recovered, backbone labelled `Ib/IIb` |
| `sars_cov_2` | 22B x ... | 0.4 % | SKIP -- Omicron clades too similar |
| `vzv` | clade 2 x ... | 0.2 % | FAIL -- info-site; genome too conserved to call |
| `hantavirus` | -- | -- | SKIP -- no clade attribute in the tree |
| `oropouche` | -- | -- | SKIP -- no clade attribute in the tree |
| `cchfv` | -- | -- | SKIP -- < 2 clades with >= 3 genomes |

**13 PASS, 7 FAIL, 4 SKIP** (0 errors). RecomFi recovers the recombinant cleanly
across the full divergence range that has both parents represented -- from dengue
serotypes (33 %) down to measles and mumps genotypes (~7 %).

Low-divergence DNA-virus / intra-species panels (mpox 0.5 %, VZV 0.2 %, ebola
3.7 %) auto-trigger **informative-site windowing** (`mode` column). This is the real
win: on mpox, base-pair windowing detects *nothing* (every window is ~997/1000
identical to both parents), while informative-site windowing detects the
recombination, recovers the IIa donor across the true insert, and places the
breakpoint. Their strict PASS still fails on edges that are not detection misses:
mpox's backbone genome carries the coarser tree label `Ib/IIb` (the recombination
itself is called correctly), ebola's donor is a deep Ebov sub-lineage 3.7 % from the
backbone, and VZV at 0.2 % is simply too conserved for any caller. Only within-Omicron
SARS-CoV-2 stays SKIP (< 4 % and the dataset offers no more-divergent pair); three
segmented viruses carry no usable clade attribute.

The moderate-divergence cases (measles, mumps, and the earlier hepatitis-A
inversion) were not a caller limitation but a **panel-recruitment artifact**: for a
hybrid of close parents the backbone parent is >95 % genome-wide ANI to the query,
so the default sibling-drop discarded it and the caller crowned a neighbouring or
the donor clade. Keeping siblings (the documented `--seed-keep-siblings` setting for
close parents with no masking twin) restores the backbone parent and the calls come
out correct. The remaining single-gene / segment failures (flu HA, H5 HA, PRRSV
ORF5) are short and finely subdivided. The pass set is a performance characterisation,
not a fixed contract -- the clades chosen follow from each dataset's current
Nextclade tree.

## Expectation schema (`expected` block)

- `major_parent` -- the top window-winner must equal this label.
- `region_minor_any` -- at least one called region's minor parent is in this set.
- `parents_present` -- every label must appear as a winner or a region parent.
- `min_regions` -- at least this many recombinant regions are called.
- `breakpoint` -- `{coord: "msa"|"query", lo, hi}`; at least one region overlaps
  the band (query coordinates map to the query's own bases).
