# Tessera validation on published data

An opt-in harness that runs the full Tessera pipeline (`build_msa` then
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
  run_benchmark.py     PHI/Rmin power+specificity on a published simulated set
  run_coalescent_benchmark.py   Posada & Crandall coalescent design (msprime)
  run_recombinhunt_benchmark.py RecombinHunt noise-robustness of breakpoint detection
  run_reassort_benchmark.py     reassort verdict precision/recall/F1 (flu HA+NA)
  run_method_comparison.py      per-caller detection on the real recombinants
  run_method_comparison_hybrids.py  GENECONV gap probe on the hard synthetic cases
  run_deep_typing.py   run the real --deep-typing lineage ladder, check nextclade-nn
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
| `enterovirus_e11` | enterovirus ~7.3 kb | Echovirus-11 x Coxsackievirus-B1, breakpoint in P2 | mafft | recombination detected; both parents named (checks parents-present, not the ambiguous backbone direction) |
| `hiv_crf02ag` | HIV-1 ~9.2 kb | CRF02_AG (IbNG): A backbone + subtype-G segments | mafft | major A; G donor region(s) over the pol/vif and vpu-env inserts |
| `hcv_2k1b` | HCV ~9.4 kb | RF1_2k/1b: genotype-2k 5' + 1b 3', breakpoint in NS2/NS3 | mafft | major 1b; genotype-2 donor over the 5'; **precise breakpoint ~nt 3187 recovered** |
| `hcv_clonal_1b` | HCV ~9.4 kb | pure genotype-1b (non-recombinant control) | mafft | resolves to 1b throughout; **0 regions** (real-data specificity) |

Each reproduces its published event (or, for the clonal control, its *absence* of
recombination) end-to-end; the current run is **7 PASS, 0 FAIL**
(`orthopox_example` SKIPs until its 7-genome collection is built). `hcv_2k1b` is the
first real **precise-breakpoint** check and `hcv_clonal_1b` the first real
**false-positive** (specificity) check. Accessions are
listed per dataset in `datasets.json` (`provenance` field) and were confirmed
against NCBI nuccore. A `sarscov2_xe` (BA.1 x BA.2) stub is present but
`enabled:false` until a public GenBank XE genome is confirmed (most are
GISAID-only). Each dataset is gated on the aligner it actually uses, so the
mafft/minimap2 cases run even when the whole-genome aligners are absent.

### Scope boundary

Bacterial **homologous** recombination (as detected by Gubbins, ClonalFrameML, and similar on
core-genome alignments) is a different problem -- phylogenetic incongruence within a clonal frame,
not inter-lineage mosaic recombination -- and is out of scope for Tessera and these benchmarks.

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
4. runs Tessera pool-only with the two exact source genomes removed (their clades
   stay represented), so the query is not a trivial self-match; window sizes adapt
   to the genome length and the aligner is per-case (minimap2 for ~200 kb mpox/VZV).
   Sibling-dropping is **off**: the synthetic pool has no recombinant twin of the
   query, and for a close-parent hybrid the backbone parent is >95 % genome-wide ANI
   to the query (its 70 % backbone dominates the average) and would otherwise be
   discarded as a masking twin -- the documented `--seed-keep-siblings` case. On a
   near-identical panel (mpox/VZV ~0.5 %) dereplication collapses the parent clades
   into one representative, so the panel is rebuilt from one central genome per
   clade; Tessera then auto-switches to informative-site windowing (the `mode`
   column);
5. checks the call: recombination detected, backbone (major parent) clade == A,
   and a donor region recovered for clade B overlapping the true span. Clade labels
   match hierarchically (`A` == `A.1`); a donor region is also credited to a sibling
   sub-clade of the donor (e.g. Marburg `RAVV.1` for a `RAVV.2` donor) **only** when
   the donor lineage is distinct from the backbone's, so a shared-top-level pair
   (RSV `A.1` / `A.D.1.8`) still requires resolving the exact donor. Runtime is
   recorded per case.

The summary table also carries the **parent-free** recombination signal per case
(`PHI p` and `Rmin`), read from `recombination_profile.tsv`, independent of the
backbone/donor verdict. Hudson-Kaplan **Rmin** is non-zero for every recombinant that
runs (23-785 events), across the whole divergence range -- including the low-divergence
mpox (0.5%) and ebola (3.7%) cases where the parent-attributed callers have least to
work with. The **PHI** p-value reaches the permutation floor in most cases but is
conservative when the panel is clean parental clades around a single hybrid (the
clade-defining sites are mutually compatible and dilute the genome-wide average), so it
can stay non-significant even where detection succeeds (e.g. yellow fever); there the
Rmin count and the per-site PHI profile carry the parent-free signal.

Detection runs the default **ensemble** (`hmm,3seq`), so the table's `agr` column marks
cases whose recovered donor region was called by more than one method (agreement, hence
higher confidence). The ensemble's union does not lower the documented PASS/FAIL set --
it can only add 3SEQ-recovered regions the HMM alone would miss.

```
export PATH="$HOME/miniforge3/envs/recomfi-aln/bin:$PATH"
python validation/run_hybrids.py            # all cases
python validation/run_hybrids.py hiv1 dengue   # only named cases
```

Needs MAFFT/minimap2/skani/skDER on PATH and contacts the Nextclade dataset server
on the first run (pools are cached afterwards under `~/.cache/tessera/nextclade`).
For a short gene/segment dataset that skani rejects the panel falls back to one
central genome per clade.

## Deep-typing ladder check (`run_deep_typing.py`)

The unit tests exercise the lineage-assignment ladder with the Nextclade and skani
calls stubbed. This opt-in check runs the real ladder once -- the path `detect
--deep-typing` uses -- so the Nextclade nearest-neighbour tier is exercised against a
live dataset. It reuses the bundled SARS-CoV-2 collection, rewrites each genome header
to a bare accession (so header mining cannot short-circuit and the ladder must reach
the nearest-neighbour tier), calls `assign_lineages`, and checks that every genome is
typed via source `nextclade-nn` with a non-empty label.

```
export PATH="$HOME/miniforge3/envs/recomfi-aln/bin:$PATH"   # skani on PATH
python validation/run_deep_typing.py
```

It SKIPs (does not fail) when skani is absent, the collection has not been fetched, or
the Nextclade fetch is unavailable. The `ncbi-datasets` source tag -- the other half of
`--deep-typing` -- is covered by the unit tests rather than a fetch-based end-to-end.

## Intragenic scan probe (`run_reassort_scan.py`)

An opt-in probe for `reassort --scan-segments` (the per-segment intragenic scan). It builds a
real-tip insert-shaped HA recombinant -- a clade-B segment (fraction 0.35-0.65) spliced into a
clade-A backbone from two divergent H3N2 HA tree-tip genomes, with an exact known breakpoint --
plus a clonal NA control, and drives the shipped `assign_segments(scan_segments=True)`. It reports
`XPASS`/`XFAIL` gated on **localization** (a called HA region overlaps the true insert span); donor
attribution and the clonal-control specificity are reported, not gating. Needs skani and an aligner
and is not part of CI.

```
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"   # skani + mafft on PATH
python validation/run_reassort_scan.py            # ephemeral (artifacts in a temp dir)
python validation/run_reassort_scan.py scan-out/  # keep the query + per-segment TSVs
```

With no argument the artifacts are written to a temporary directory removed on exit (only the
printed report survives); pass an output path to keep the query and per-segment tables. The pure
localization scorer (`region_overlaps_span`) is unit-tested in CI; the end-to-end run is opt-in. The
measured result is recorded in `attribution-results.md`.

## Simulated detection benchmark (`run_benchmark.py`)

A published, quantitative benchmark using the SANTA-SIM alignments from Jaya, Brito & Darling
(2023, *Virus Evolution* 9(2):vead066; Dryad doi:10.5061/dryad.d7wm37q6f). The archive files are
named `msa_m<mut>_rc<rec>_n<seqs>_dual<d>_rep<k>.fasta` with the recombination rate in the filename
but **no explicit breakpoint coordinates**, so this harness scores *alignment-level* detection --
which maps directly onto Tessera's parent-free **PHI + Hudson-Kaplan Rmin** test (no query, panel,
or aligner needed, since the files are already aligned). It reports **power** (fraction of
recombining alignments, `rc>0`, with PHI p < alpha) and **specificity** (fraction of non-recombining
alignments, `rc=0`, with PHI p >= alpha) across the mutation x recombination grid -- directly
comparable to the PhiPack/Profile column of that study.

Dryad gates its downloads behind a login, so fetch `performance.tar.gz` (and optionally
`scale.tar.gz`) from the Dryad page in a browser and extract the `.fasta` files into the cache
(`validation/data/benchmark/`, or set `$TESSERA_BENCHMARK_DIR`); the harness SKIPs when none are
present.

```
python validation/run_benchmark.py                    # power / specificity table
python validation/run_benchmark.py --alpha 0.05 --window 100 --max 200
```

The pure filename parser and power/specificity aggregation are unit-tested in CI; the end-to-end
run is opt-in and needs no aligner. Breakpoint-position accuracy is out of scope here (the archive
has no breakpoint annotations); it would need our own SANTA-SIM runs with breakpoint logging.

## Coalescent benchmark (`run_coalescent_benchmark.py`)

Reproduces the evaluation design of **Posada & Crandall (2001, *PNAS* 98:13757)**: neutral
sequences simulated under the **coalescent with recombination** across a mutation x recombination
grid, scoring a method's false-positive rate (specificity, at `rec=0`) and power (`rec>0`). The
2001 alignments were never archived, so the design is regenerated with `msprime` (Hudson's
algorithm); the method under test is Tessera's parent-free PHI test. Self-contained (nothing is
downloaded) and opt-in -- it needs msprime (`pip install msprime`) and SKIPs without it.

```
python validation/run_coalescent_benchmark.py            # power / specificity grid
python validation/run_coalescent_benchmark.py --reps 20
```

Measured (40 samples x 2 kb, PHI window 100): **power 1.0, specificity 1.0** across the grid --
PHI fires on every recombining coalescent sample and never on the clonal (`rec=0`) controls. The
pure power/specificity scorer is shared with `run_benchmark.py` and unit-tested in CI.

## Noise-robustness benchmark (`run_recombinhunt_benchmark.py`)

Reproduces the noise-robustness design of **RecombinHunt (Bianchi et al., *Nat. Commun.* 2024,
s41467-024-47464-5)**, which measured detection sensitivity as increasing numbers of injected
(non-characteristic) mutations were added to simulated recombinants. Their SARS-CoV-2 BA.2 x AY.45
sequences are GISAID-gated, so the (virus-agnostic) design is reproduced on Nextclade lineage
genomes: splice a one-breakpoint recombinant from the two most-divergent clade consensuses, inject
the mutation grid (0, 3, 5, 10, 15, 20, 30), run the recomb scan, and report the fraction of
replicates where a called region overlaps the true breakpoint. Needs skani and an aligner; opt-in.

```
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
python validation/run_recombinhunt_benchmark.py
python validation/run_recombinhunt_benchmark.py --dataset nextstrain/rsv/a/EPI_ISL_412866
```

Measured on the default divergent pair (flu H3N2 HA, C.1 x K): breakpoint recovery **5/5 at every
noise level 0-30** -- robust to the full injected-mutation range on a divergent pair. A
low-divergence dataset degrades faster; the pure noise-injection and per-level aggregation are
unit-tested in CI.

## Reassortment benchmark (`run_reassort_benchmark.py`)

Scores the shipped `tessera reassort` verdict the way the influenza-reassortment literature does
(TreeSort / TreeKnit / CoalRe report **precision / recall / F1** against known reassortment events).
Flu H3N2 HA and NA tips are strain-labelled and many strains are typed in both segment datasets, so a
labelled query set is built: **clonal** = the HA and NA of one cross-typed strain (one parent);
**reassortant** = the HA of strain A with the NA of strain B. Each query runs through the shipped
`assign_segments`; `reassortant` is the positive class. Needs skani + network; opt-in.

```
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
python validation/run_reassort_benchmark.py            # 15 clonal + 15 reassortant
python validation/run_reassort_benchmark.py --n 25
```

Measured (479 cross-typed strains; 10 clonal + 10 reassortant): **precision 0.83, recall 1.00,
F1 0.91** -- every reassortant is caught, and the two false positives quantify the documented
cross-typing-coverage limit (a clonal isolate whose exact strain is not in both segments' nearest
top-k can read `reassortant`; see `docs/reference-panels.md`). The pure confusion/F1 scorer is
unit-tested in CI. Recorded as measured, not tuned.

## Method comparison (`run_method_comparison.py`)

Two comparisons on the same published-recombinant alignments (the `datasets.json` positives): the
first internal, the second against a standard external tool.

1. **Tessera's own callers** -- HMM, 3SEQ, MaxChi, Bootscan -- each run alone, plus the default
   four-caller ensemble, tabulating per-caller detection of the known recombination.
2. **OpenRDP** (the maintained RDP5 reimplementation, [PoonLab/OpenRDP](https://github.com/PoonLab/OpenRDP))
   -- its seven methods (GENECONV, Bootscan, MaxChi, SiScan, Chimaera, 3SEQ, RDP) run on the *same*
   alignment, so Tessera sits beside the standard methods on identical data.

Needs an aligner + the fetched data; opt-in. OpenRDP is pure Python but pins `numpy<2`/`h5py<3.11`,
so it is easiest to keep in its own environment; the harness resolves it via `$OPENRDP_CMD`, then an
`openrdp` on PATH, then `conda run -n $OPENRDP_CONDA_ENV openrdp` (default env `openrdp`), and skips
that half with a note if none resolves. (An earlier attempt used the `raspberryhusky` fork, which
pins `numpy==1.17.3` and needs `cblas.h`/`pybind11` and does not build; the PoonLab implementation
builds cleanly.)

```
# once, in a dedicated env:
conda create -n openrdp python=3.10 -y
conda run -n openrdp pip install "numpy<2.0.0" "scipy>=1.5.0" "h5py<3.11.0" \
    git+https://github.com/PoonLab/OpenRDP

export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
python validation/fetch.py && python validation/run_method_comparison.py
python validation/run_method_comparison.py --no-openrdp   # Tessera callers only
```

The OpenRDP detection metric: a method "detects" when the known recombinant sequence appears in any
recombination triplet it reports (Recombinant / Parent1 / Parent2). Using the whole triplet is robust
to the parent/child direction ambiguity that triplet methods -- and Tessera's single backbone --
share. The parsing and detection helpers are unit-tested in CI.

Measured (six positives; orthopox `SKIP`s until its query is fetched, ~200 kb). One table, Tessera's
callers on the left of the `|` divider and OpenRDP's methods on the right, both on the same alignment
(OpenRDP's default set; SiScan is opt-in via `--with-siscan` -- it is permutation-heavy, timed out at
>15 min on the full-genome HIV alignment, and detected nothing where it completed):

| dataset | hmm | 3seq | maxchi | bootscan | ensemble | `\|` | geneconv | bootscan | maxchi | chimaera | threeseq | rdp |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sarscov2_xbb    | yes | yes | yes | yes | yes | `\|` | yes | yes | yes | no  | no  | yes |
| hiv1_crf        | yes | yes | yes | yes | yes | `\|` | no  | yes | yes | no  | yes | yes |
| norovirus_gii   | yes | yes | yes | yes | yes | `\|` | yes | yes | no  | no  | yes | yes |
| enterovirus_e11 | yes | yes | yes | yes | yes | `\|` | yes | no  | yes | no  | yes | yes |
| hiv_crf02ag     | yes | yes | yes | yes | yes | `\|` | yes | yes | yes | no  | yes | yes |
| hcv_2k1b        | yes | yes | yes | yes | yes | `\|` | yes | yes | yes | no  | yes | yes |

Columns left of the `|` are Tessera's callers; right of it are OpenRDP's methods (`maxchi`/`bootscan`
appear on both sides -- the same algorithm family in each tool).

On these well-characterised real recombinants Tessera's four callers each fire independently, so no
single method carries the ensemble here; the ensemble earns its keep on the harder cases -- short
tracts, low divergence, multi-breakpoint mosaics -- where the callers disagree (see the synthetic
hybrid harness `--compare` table below). OpenRDP concurs on the recombination in every case -- at
least four of its six default methods fire per dataset -- though Chimaera (and SiScan, when run) add
nothing on these small taxon sets. On `hcv_2k1b` its 3SEQ/GENECONV/Bootscan place the breakpoint at
nt ~3186-3187 -- an independent, external confirmation of the breakpoint Tessera recovers (~3187) and
of the published NS2/NS3 junction. Recorded as measured.

### Would an OpenRDP method improve Tessera? (`run_method_comparison_hybrids.py`)

Of OpenRDP's seven methods, MaxChi and Bootscan are already Tessera callers, RDP is the
similarity-switch signal Tessera's HMM already carries, and Chimaera/SiScan add nothing on these
data -- leaving **GENECONV** (runs of matching polymorphic sites) as the one method with a genuinely
different signal. On the easy real recombinants Tessera already detects everything, so a benefit
could only appear on the *hard* synthetic cases (short tracts, low divergence, mosaics) or in
specificity. This probe reuses the synthetic-hybrid builder to make those cases with known ground
truth, runs Tessera's ensemble and OpenRDP on the same `panel.msa.fasta`, and counts, per OpenRDP
method, **rescues** (a positive Tessera's ensemble missed) and **false-calls** on negatives.

```
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
python validation/run_method_comparison_hybrids.py     # needs the openrdp conda env + network
```

Measured over the eight hard cases (7 positives: `low_div`, `panel_donor_absent`, the four mosaics
ABAC/AB_9010/AB_short/AB_terminal, `mask_sibling`; 1 negative: `neg_pure`):

| | positives detected | negatives false-called | GENECONV rescues | GENECONV false-calls |
|---|---|---|---|---|
| Tessera ensemble | 7/7 | 0/1 | -- | -- |
| OpenRDP (incl. GENECONV) | 7/7 | 0/1 | **0** | **0** |

Tessera's four callers each fired on all seven hard positives and the ensemble rejected the negative;
OpenRDP's GENECONV also detected all seven and false-called none, so it rescued nothing Tessera
missed. **Verdict: on this envelope, adding a GENECONV caller would not improve detection or
specificity, so none was written.** On this envelope GENECONV is detecting the same events, just
cases Tessera already covers. But the must-pass cases are ones Tessera is built to pass, so to find a
real gap the probe carries an **adversarial tier** of cases Tessera fails.

#### Adversarial sub-window tracts (`--adversarial`)

The `AB_micro` mosaic pattern splices a donor tract of `tract_windows` x the sliding window (default
0.3, i.e. below one window) and scores it strictly, to push past the windowed callers' floor. A
tract-length sweep (`python validation/run_method_comparison_hybrids.py --adversarial`):

| tract (x window) | dataset | Tessera | OpenRDP methods that detect |
|---|---|---|---|
| 1.00 | wnv   | PASS (hmm, 3seq, maxchi, bootscan) | geneconv, bootscan, threeseq, rdp |
| 0.50 | wnv   | PASS (3seq, maxchi, bootscan)      | geneconv, bootscan, threeseq, rdp |
| 0.30 | wnv   | **FAIL**                           | geneconv, bootscan, threeseq, rdp |
| 0.15 | wnv   | **FAIL**                           | geneconv, threeseq, rdp |
| 0.30 | zika  | **FAIL**                           | geneconv, threeseq |
| 0.30 | mumps | **FAIL**                           | geneconv, bootscan, threeseq, rdp |

There is a **real gap**: below ~0.3 x window Tessera's windowed scan misses the tract (HMM drops out
first, by 0.5x), while OpenRDP still flags the query. But it is **not** a GENECONV-specific gap:
across the four failures OpenRDP's GENECONV rescued 4/4 and its 3SEQ rescued 4/4 too -- and Tessera
already ships a 3SEQ caller.

To settle it, a query-vs-panel **GENECONV caller was built** (`src/tessera/recomb/geneconv.py`, opt-in
via `--method`) and the sweep re-run with it in the ensemble
(`--tessera-methods hmm,3seq,maxchi,bootscan,geneconv`):

| tract (x window) | dataset | Tessera + geneconv | did native geneconv fire? |
|---|---|---|---|
| 1.00 | wnv   | PASS | yes (with the ensemble) |
| 0.50 | wnv   | PASS | yes (with the ensemble) |
| 0.30 | wnv   | **FAIL** | no |
| 0.15 | wnv   | **FAIL** | no |
| 0.30 | zika  | **FAIL** | no |
| 0.30 | mumps | **FAIL** | no |

Tessera's own GENECONV fires on the tracts it can already detect (1.0x, 0.5x) but **rescues none of
the four sub-window failures** -- exactly where OpenRDP's GENECONV still detects. This confirms the
diagnosis directly: a query-vs-panel GENECONV hits the same few-discriminating-sites wall as 3SEQ;
OpenRDP catches these tracts because it scans the whole population alignment (many taxa, more power
per short region), not because of GENECONV's signal. **Verdict: the short-tract miss is a limit of
Tessera's query-vs-panel scan configuration, not a missing method -- adding GENECONV does not close
it.** The caller is kept opt-in (not in the default ensemble): it is correct and available for method
comparison, but does not extend Tessera's detection envelope. The lever for short tracts is Tessera's
windowing / scan design. Recorded as measured; probe and tier stay in the suite. (Several 2-parent
`AB_micro` cases SKIP when a donor clade loses its only panel representative after source removal -- a
harness invariant, not a detection result.)

A dataset is **SKIP**ped (not failed) when it cannot supply a valid test: the
most-divergent clade pair is below ~4 % divergence (too few discriminating sites:
mpox, VZV, ebola, SARS-CoV-2 within Omicron), or it has fewer than two clades with
at least three genomes -- including datasets with no clade attribute at all
(hantavirus, Oropouche, CCHFV).

### Harder cases: precision as well as recall

Beyond the single-insert positives above, the harness carries additional `case_type`s so it
measures **specificity** (false positives), not only sensitivity. The default run prints both
`sensitivity P/N` and `specificity Q/M` lines. Phase-1 case types:

- `neg_pure` -- a pure, non-recombinant genome that must yield **zero** recombinant regions.
- `low_div` -- pins the *closest* viable clade pair (`pair_objective: "min"`, within a
  `divergence_band`) and **requires** backbone + donor at the correct top-level clade, with no
  sub-4% relaxation, testing attribution near the informative-site floor.
- `panel_donor_absent` -- removes the true donor clade from the panel; the run must report a
  coverage gap rather than mis-attribute the region to a present sibling.
- `panel_equidistant` -- two candidate donors equally close to the query; the caller must
  attribute the true one (a tie-break guard).

An out-of-envelope case is `SKIP`ped, and a genuine attribution shortfall is reported as a FAIL,
not hidden -- see `attribution-results.md` for the measured result.

Phase 2 adds hard *topologies* via `make_mosaic` / `true_spans`: a `mosaic` case type with
`pattern` in `{ABAC, AB_9010, AB_short, AB_terminal}` (multi-breakpoint 3-parent, asymmetric,
sub-window short tract detection-gated, terminal breakpoint), and a `mask_sibling` case that
requires the donor attributed *exactly* when a sibling sub-clade is present (a regression guard
for the plurality-major major-parent fix).

Phase 3 adds **frontier** probes behind `--frontier` (`python validation/run_hybrids.py
--frontier`): `inter_species` (a cross-species recombinant, envelope measured by skani ANI vs the
0.80 HMM floor) and `reassortant` (a segment swap, detection-gated at the segment junction). They
are scored `XPASS`/`XFAIL`/`KNOWN-LIMIT` in a separate table that can never count as a regression;
without the flag the must-pass headline is unchanged.

### Observed performance (24 pathogens; one representative dataset each)

| case | backbone x donor | divergence | result |
|------|------------------|-----------:|--------|
| `dengue` | DENV1 x DENV4 | 33.1 % | PASS |
| `marburg` | MARV.B.2 x RAVV.2 | 21.5 % | PASS |
| `yellow_fever` | Clade VII x Clade III | 21.0 % | PASS |
| `iav_h5_ha` | 2.3.4.4h x 2.3.2.1f | 11.0 % | PASS -- pinned a well-represented divergent pair (the basal Am-nonGsGD auto-pick had no stand-in) |
| `wnv` | 2 x 1B | 20.2 % | PASS |
| `hmpv` | B1 x A2.2.1 | 19.0 % | PASS |
| `prrsv2` | L1H x L8D | 17.6 % | PASS -- pinned two large lineages (the L1C.2 auto-pick had no stand-in) |
| `hepatitis_a` | IIIA x IIA | 16.6 % | PASS |
| `chikv` | III-Asian x I-WestAfrica | 15.5 % | PASS |
| `hiv1` | A1 x B | 15.1 % | PASS |
| `enterovirus_d68` | B3 x A2/D | 11.3 % | PASS |
| `zika` | Asian x African | 10.9 % | PASS |
| `rubella` | 2B x 1G | 9.0 % | PASS |
| `measles` | H1 x B3 | 7.5 % | PASS |
| `mumps` | A x K | 6.9 % | PASS |
| `flu_h3n2_ha` | C.1 x K | 6.8 % | PASS |
| `rsv_a` | A.1 x A.D.1.8 | 6.6 % | PASS |
| `ebola` | Ebov-2013 x Ebov-2018b | 3.7 % | PASS |
| `mpox` | Ib x IIa | 0.5 % | PASS -- below the 4 % floor: scored on detection + donor |
| `sars_cov_2` | 24H x outgroup | 0.4 % | PASS -- below the 4 % floor: scored on detection + donor |
| `vzv` | clade 2 x clade 9 | 0.2 % | PASS -- below the 4 % floor: scored on detection + donor |
| `hantavirus` | denovo_1 x denovo_4 | 20.3 % | PASS -- de-novo ANI lineages (no clade attribute) |
| `cchfv` | denovo_1 x denovo_3 | 22.6 % | PASS -- de-novo ANI lineages (no clade attribute) |
| `oropouche` | -- | -- | SKIP -- L-segment genomes all > 98 % ANI; no lineage split possible (data limit) |

**23 PASS, 0 FAIL, 1 SKIP** (0 errors). All 23 cases that run detect the recombinant and
recover both parents; agreement is lineage-aware (two callers that pick different
representative genomes of one lineage still count). Tessera recovers the recombinant across
the full divergence range that has both parents represented, from dengue serotypes (33 %)
down to the mpox clade-I/II recombination at 0.5 %. The five previously-skipped datasets now
made runnable are fixed in ways that keep scoring honest: `sars_cov_2` joins `mpox`/`vzv` under the
sub-4 % detection+donor rule; `iav_h5_ha` and `prrsv2` pin a well-represented divergent
clade pair (the auto-selected most-divergent clade had no panel stand-in after source
removal); and `hantavirus` and `cchfv`, which carry no clade attribute, are typed de-novo
by ANI clustering their tips (`cluster_ani`, reusing the `type-lineages` clustering).
`oropouche` remains a documented data limit -- its 264 L-segment genomes are all > 98 %
ANI, so no >= 2 lineages of >= 3 members can be formed at any threshold; it is reported
SKIP rather than forced.

#### Scoring rules

A synthetic case is scored fairly only when it is well-posed. Four rules -- each
motivated by the harness's own design or the underlying biology, applied uniformly --
decide PASS / FAIL / SKIP:

1. **Reassortment suffix.** A `/rYYYY` tag (e.g. ebola `Ebov-2018b/r2021b`) marks the
   same parental lineage re-emerging in a later year, so `clade_match` compares base
   lineages (`Ebov-2018b/r2021b` == `Ebov-2018b`). Ebola's donor call is correct and now
   scores so.
2. **Attribution floor (4 %).** Below the harness's own meaningfulness floor the exact
   backbone clade is statistical noise -- near-identical genomes from adjacent clades win
   windows by chance -- so sub-4 % cases (`mpox` 0.5 %, `vzv` 0.2 %) are scored on
   **detection + donor-region recovery** (which still requires the donor clade and span),
   not the noisy backbone label. The detection itself is the real win: on mpox base-pair
   windowing finds *nothing*, while informative-site windowing recovers the IIa donor
   across the true insert.
3. **Representation invariant.** The harness removes the two source genomes; its design
   keeps "their clades represented by other genomes." When that fails -- the backbone
   clade (`iav_h5_ha` Am-nonGsGD) or the donor lineage (`prrsv2` L1C, even via a sibling)
   has no stand-in left -- attribution cannot be tested regardless of detection quality,
   so the case is **SKIP**ped rather than scored a detection FAIL (the app correctly flags
   those as coverage gaps).
4. **Junk labels.** `unassigned` / `unclassified` genomes are not a clean parental
   lineage and are excluded from parent selection.

MaxChi and Bootscan are part of the default ensemble because they raise confidence at no
cost: on the earlier 18-run set, restricting to `hmm,3seq` (`HARNESS_METHODS=hmm,3seq`) gave
the same verdicts but corroborated the donor in only 12 of 18 cases, whereas the four-caller
default reached 17 -- the two extra independent callers add agreement with no regression and
no new false positives. Adding the legacy heuristic (`HARNESS_METHODS=all`) does not change
either number. `--method hmm,3seq` is the lighter option for large genomes (Bootscan's
bootstrap is the main added cost).

There are no attribution FAILs on the current Nextclade trees. `flu_h3n2_ha` and `rsv_a`
were previously the two attribution misses (a fine-subclade near-tie on the ~1.7 kb HA
segment, and a donor attributed to a neighbouring `A`-sublineage rather than `A.D.1.8`); on
the current trees both recover the correct backbone and donor. The pass set is a performance
characterisation, not a fixed contract -- the clades chosen follow from each dataset's
current Nextclade tree, so a tree update can move a fine-subclade case either way. The
`--compare` analysis below was measured earlier, when those two were the FAILs; its numbers
are anchored to that run.

#### Attribution comparison (`--compare`)

Because both FAILs are *attribution* (wrong adjacent sub-clade), not detection, the harness
has a measurement mode that asks whether the two features built for attribution help: the
**barcode** caller (denoised per-clade markers) and a **source-removed per-clade consensus**
panel (the `--pool-consensus` idea, with the two spliced source genomes excluded so the test
stays honest). Run:

```
python validation/run_hybrids.py --compare            # all cases
python validation/run_hybrids.py --compare rsv_a flu_h3n2_ha
```

For each case it builds the setup once and scores a 2x2 grid -- {tip, consensus} panel x
{default ensemble, +barcode} -- writing `data/attribution_compare.tsv` and a per-config
table. Alongside the PASS/FAIL verdict it reports a finer **attribution tier**
(`exact` / `sibling` / `mismatch`) and a **clade depth** (shared leading clade components
between the recovered and the true donor), so an improvement shows even when it does not flip
the verdict. The `baseline` column reproduces the headline 16/2/6; the other columns are
this measurement, **not** the headline numbers. Whether any configuration is then
promoted to a product default is a separate, data-gated decision -- no PASS rule is loosened
here to manufacture a win.

A fifth column, `reattribute`, runs the tip panel with `--reattribute-donors` enabled: after
detection, the minor-parent label is re-evaluated by scoring per-clade consensus sequences
against the query over the recombinant insert, with the backbone unchanged. This is expected to
recover `rsv_a`'s donor (the insert is a closer match to `A.D.1.8` than the initially-called
neighbouring sublineage) in cases where the consensus panel alone could not, because it
operates on the alignment already produced rather than replacing the panel itself.

The first run of this measurement is written up in
[`attribution-results.md`](attribution-results.md): on these datasets barcode is inert (no
verdict or attribution change) and the consensus panel is a net regression (it improves donor
attribution on two cases but degrades backbone attribution on two others), so neither was
promoted -- the value of measuring before committing.

## Expectation schema (`expected` block)

- `major_parent` -- the top window-winner must equal this label.
- `region_minor_any` -- at least one called region's minor parent is in this set.
- `parents_present` -- every label must appear as a winner or a region parent.
- `min_regions` -- at least this many recombinant regions are called.
- `breakpoint` -- `{coord: "msa"|"query", lo, hi}`; at least one region overlaps
  the band (query coordinates map to the query's own bases).
