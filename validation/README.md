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
| `hiv1_crf` | HIV-1 ~10 kb | CRF01_AE (CM240) vs pure subtypes | mafft | major A1; recombinant signal in the env region (~query 5.7-9.0 kb) |
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

## Expectation schema (`expected` block)

- `major_parent` -- the top window-winner must equal this label.
- `region_minor_any` -- at least one called region's minor parent is in this set.
- `parents_present` -- every label must appear as a winner or a region parent.
- `min_regions` -- at least this many recombinant regions are called.
- `breakpoint` -- `{coord: "msa"|"query", lo, hi}`; at least one region overlaps
  the band (query coordinates map to the query's own bases).
