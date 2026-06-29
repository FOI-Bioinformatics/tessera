# RecomFi
RecomFi (Recombination Finder) identifies recombination events in a query sequence, contigs or genome, against a collection of reference sequences.

# Description
RecomFi is developed to identify recombination in relatively similar datasets, such as between (sub)species of a genus or family. It generates a reference-anchored "pseudo-MSA (multiple sequence alignment)" by using one sequence as a backbone. This makes RecomFi fast but limits the resolution. With the pseudo-MSA strategy the query may be a fragmented genome, for example a set of contigs, which RecomFi organizes relative to the backbone.

Recombination events are detected by sliding a window over the MSA and computing, in each window, the similarity between the query and each reference. A recombination event is indicated where the query is near reference A across most of its length but is near another reference B over a region.

RecomFi is organized as a small, modular Python package:
- Alignment is delegated to a pluggable **aligner backend** (`progressivemauve`, `sibeliaz`, or `cactus`), discovered through entry points so additional backends can be added without changing the core.
- The recombination scan uses a built-in sliding-window distance engine (no external alignment-analysis dependency).

# Installation
RecomFi needs Python (>= 3.11) and at least one aligner backend. The aligner binaries are most easily installed with conda.

```
# create an environment with an aligner backend, then install RecomFi
conda create -n recomfi -c conda-forge -c bioconda python">=3.11" mauve "boost-cpp=1.74.0"
conda activate recomfi
pip install .
```

Optional aligner backends can be added to the same environment:
```
conda install -c bioconda sibeliaz   # or: cactus
```

An `environment.yml` is provided that installs Python, all three backends, and RecomFi in one step:
```
conda env create -f environment.yml
conda activate recomfi
```

# Aligner backends
All backends produce a reference-anchored alignment, which is what the recombination scan assumes. Choose with `--aligner` and tune with repeatable `--aligner-arg key=value`:

| Backend | Best for | Notes |
|---|---|---|
| `sibeliaz` (default) | Moderately divergent genomes, including rearrangements | Installs cleanly via conda; `kmer`, `abundance`, `bubble`, `filtermemory` |
| `mafft` | Similar, largely collinear genomes | True base-level alignment, the canonical input for the window method; adds a fragmented query with `--addfragments`. `maxiterate`, `retree`, `op`, `ep`, `sixmerpair` |
| `minimap2` | Speed and assembly/contig queries | Fast assembly-to-reference projection; `preset` (default `asm20`, e.g. `asm10` for closer genomes) |
| `progressivemauve` | Genomes with large rearrangements/inversions | Tolerant but slow, heavy, and not available as a conda build on all platforms; `seed_weight`, `single` |
| `cactus` | Same-species pangenomes | Resource heavy (Toil/containers) |

`sibeliaz` is the default: it installs cleanly across platforms and, on the example data, reproduces `progressivemauve`'s recombination coordinates. For very similar, collinear genomes `mafft` gives the most faithful base-level signal and `minimap2` the fastest run (and the best fit for a fragmented query); `progressivemauve` remains an option for genomes with large rearrangements. Reference-anchored backends drop material inserted relative to the backbone; `mafft` keeps it as a true alignment.

Example: `recomfi msa ... --aligner minimap2 --aligner-arg preset=asm10` or `recomfi msa ... --aligner mafft --aligner-arg maxiterate=1000`.

# Example dataset
Find an example dataset of orthopoxvirus in `example_data/`. The query is a short-read assembly (8 contigs) of a synthetic cowpox sample with a variola segment. The collection are reference-labelled orthopoxvirus sequences from `BV-BRC.org`.

Example folder structure (query is `cowpox_with_variolaInsert.fasta.gz`):
```
.
├── collection
│   ├── camelpox.fasta.gz
│   ├── cowpox.fasta.gz
│   ├── cowpox_KC813504.fasta.gz
│   ├── monkeypox.fasta.gz
│   ├── taterapox.fasta.gz
│   ├── vaccinia.fasta.gz
│   └── variola.fasta.gz
└── cowpox_with_variolaInsert.fasta.gz
```

# Usage

## One-shot detection (no genomes needed)
Give RecomFi only a query and it detects the taxon, recruits a diverse reference
panel from NCBI, aligns, and calls recombination:
```
recomfi detect --query CRF01_AE.fasta --output out/ --email you@example.org
```
It recruits the parental lineages organism-agnostically -- a negative-lineage BLAST
(exclude the query's own over-represented lineage so the divergent parents surface)
plus NCBI Virus taxonomy diversity -- drops the query's siblings, and competes
*lineages* rather than individual genomes, so a recombinant whose own lineage is
common in NCBI is not masked. The report carries a plain-language verdict with a
confidence label. Fetched panels are cached per taxon (`--cache-dir`), so a repeat
run is fast. Needs an aligner and Entrez Direct (skani/skDER and the `datasets` CLI
improve recruitment). For a heavily sequenced taxon (e.g. SARS-CoV-2) supply a local
panel with `--candidate-pool`. The steps below show the underlying `msa` / `recomb`
commands when you want to drive them manually.

### Build the panel and detect separately
To recruit the donor panel without running detection -- so you can inspect the panel,
or re-run detection with different window/HMM settings without re-fetching or
re-aligning -- use `build-panel`, then `recomb`:
```
recomfi build-panel --query CRF01_AE.fasta --output panel/ --email you@example.org
recomfi recomb --msa panel/panel.msa.fasta --query CRF01_AE --output panel/
```
`build-panel` uses the same parent-recruiting, sibling-dropping recruitment as
`detect` but stops at the panel, writing `panel/collection/` (the donor genomes), a
stable `panel/panel.msa.fasta` alignment, and `panel/panel_lineages.tsv`. The exact
`recomfi recomb` follow-up command is logged on completion. The panel is recruited
for this query (donors are its regional best matches), so it is not a reusable
organism-wide reference set. Pass `--collection` to grow an existing collection
instead of recruiting from scratch. The full-featured `fill-references` command takes
`--no-report` for the same effect.

### Typed lineage names in the report
When `detect`, `build-panel`, or `fill-references` recruit a panel, each reference is
typed with a genotype mined from its genome header -- the NCBI lineage note
(`pangolinClassification` / `isolate.lineage`) or a designation token in the GenBank
title (e.g. `GII.P16-GII.4`, `CRF01_AE`, `B.1.1.7`). The mapping is written to
`<output>/lineages.tsv`, and the report then names parents by lineage, e.g.
"GII.P16-GII.4 (MK573073)" instead of the bare accession. A standalone
`recomfi recomb` picks up a `lineages.tsv` sitting beside the output or the MSA, so the
two-step workflow above carries the typed names through. To override or supply names
(for schemes that titles do not carry well, such as HIV pure subtypes or a curated CRF
reference set), pass `--lineage-map accession_to_genotype.tsv` (a `accession<TAB>genotype`
table) to any of these commands; user-supplied names take precedence over mined ones.
With no source, references fall back to bare accessions and the report is unchanged.

The query is also typed from its own header; when its lineage differs from the file
name the verdict states it ("The query is typed as GII.P16-GII.4"). For a SARS-CoV-2
query whose lineage is a Pango recombinant (an `X` lineage such as `XBB.1.5`), RecomFi
looks up the designated parents in the Pango `alias_key.json` (fetched once into the
cache) and adds a cross-check block, so the recruited parents can be compared with the
designated ones. The lookup is best-effort: a network failure simply skips it.

Generate a multiple sequence alignment:
```
recomfi msa --query cowpox_with_variolaInsert.fasta.gz --collection collection/ --output msa.fasta

# Choose an aligner backend (default: sibeliaz) and pass tuning options:
#   recomfi msa ... --aligner sibeliaz --aligner-arg kmer=15
#   recomfi msa ... --aligner progressivemauve --aligner-arg seed_weight=11
#
# If you have a single-contig query you can use it as the backbone instead of a
# reference from the collection:
#   recomfi msa ... --query-as-backbone
```

Identify recombination events (state the query label as it appears in the MSA, i.e. the query file name without extension):
```
recomfi recomb --msa msa.fasta --query cowpox_with_variolaInsert --output recomfi_out

# The window, step, metric, number of plotted datasets and plot format are
# configurable:
#   recomfi recomb ... --window-size 1000 --window-step 100 --top-n 5 --plot-format png
#
# Region calling can be tuned (defaults derive from the window size):
#   recomfi recomb ... --min-region 1000 --margin 0.0 --merge-gap 1000
```

Run `recomfi --help`, `recomfi msa --help` or `recomfi recomb --help` for the full set of options.

# Output
RecomFi computes, in sliding windows across the MSA, the similarity of the query to each reference (1 = identical, 0 = no similarity). The reference winning the most windows is the **major parent** (the backbone donor).

To call recombinant regions, the default caller (`--method hmm`) segments the query against the reference panel with a hidden Markov model (jpHMM-style): each window emits a binomial copying likelihood per reference, and a single jump rate (`--jump-rate`) penalises switching reference, so near-identical references do not flip and thin windows cannot drive a call. A segment is reported as recombinant only when its donor beats the major parent on the **discordant sites** — positions where the query matches one candidate parent but not the other — by a sign test at level `--alpha`. This is far more discriminating than an all-sites margin (it recovers subtle breakpoints between near-identical parents) and does not invent regions from noise. Each region carries a **support** (the share of distinguishing sites favouring the donor), a sign-test **p-value** with a Benjamini-Hochberg **q-value** (false-discovery-rate across the segments tested), and a **breakpoint uncertainty interval**. The legacy `--method heuristic` (margin / merge-gap / min-region) is kept for comparison.

```
Recombination regions (major parent: cowpox_KC813504):
  Minor parent  Major parent     Query start  Query end  Length(bp)  Sim minor  Sim major  Support  q-value  Breakpoint
  --------------------------------------------------------------------------------------------------------------------
  variola       cowpox_KC813504  66268        147150     80882       0.999      0.977      0.97     2e-300   66768
```
A high support with a small q-value is a confident call; a region with strong directional support but a large q-value (few distinguishing sites, e.g. a recombination between near-identical lineages) is flagged as marginal rather than dropped.

**Low-divergence panels (intra-species sets, DNA viruses).** When the references are nearly identical — e.g. mpox clades (~0.5 %), VZV (~0.2 %), within-species ebola — a fixed base-pair window holds only a handful of discriminating columns diluted by hundreds of identical ones, so the per-reference emission contrast collapses and the segmentation loses power. The signal is still there (low *percentage* divergence over a large genome is still hundreds of variable sites), just buried. RecomFi therefore switches to **informative-site windowing**: windows span a fixed number of polymorphic columns rather than a fixed base-pair width, concentrating the signal so the HMM regains contrast. This is automatic when the references differ at less than ~8 % of columns, and controllable with `--informative-sites` / `--no-informative-sites` (and `--informative-window` / `--informative-step`). Breakpoint intervals are necessarily coarser at low divergence — you cannot localise a switch more finely than the spacing of the discriminating sites.

**`--method 3seq` (scan-aware triplet test).** A second caller, complementary to the HMM, after Boni, Posada & Feldman (2007). For the query against the major parent and each candidate donor it looks only at the **discriminating sites** (where the two parents differ and the query matches one of them) and measures the **maximum drawdown** of the resulting +1/−1 walk — a sustained run of donor matches inside the backbone, i.e. a mosaic. Its p-value is the *exact* probability that a random arrangement of the same matches reaches that drawdown (a dynamic program over the walk depth, no dependency; a vectorised permutation falls back for very large inputs), Benjamini-Hochberg-corrected across the donors tested. Because it is purely informative-site based it keeps full power on near-identical panels — it detects the mpox clade-I/II recombination where base-pair windowing finds nothing — and because the null accounts for scanning every breakpoint it does not over-call. Run it with `recomfi recomb ... --method 3seq`.

**Parent-free recombination signal (PHI + Rmin).** Alongside the parent-attributed regions, every run reports a parent-agnostic diagnostic that asks only whether the alignment carries recombination at all, with no candidate parents required — the regime where the HMM and 3SEQ callers have least to work with (low divergence, or the true donor absent from the panel). It is built from the alignment's biallelic informative columns and the four-gamete test (two columns are incompatible when all four gametes are present, which implies a recombination between them). The **PHI test** (Pairwise Homoplasy Index; Bruen, Philippe & Bryant 2006) is significant when columns near each other on the genome are more compatible than a random reordering — the signature of shared local genealogy under recombination — with a one-sided permutation p-value. **Rmin** (Hudson & Kaplan 1985) is the minimum number of recombination events the incompatibilities force, with the intervals as breakpoint candidates. Both are dependency-free; the synthetic-hybrid harness fires the PHI test (p at the permutation floor) and a non-zero Rmin across the full divergence range, from dengue serotypes (33%) down to the mpox clade-I/II recombination at 0.5%. It runs for every `--method`; disable with `--no-phi`, or widen its window with `--phi-window`.

It remains an **indicative screen**: the built-in HMM and 3SEQ tests are fast triplet/segmentation screens, not a full tree-based analysis (such as GARD). Treat regions as candidates to confirm.

Output files in the chosen directory:

| File | Contents |
|---|---|
| `recombination_regions.tsv` | Called regions: minor/major parent, start/end in **both MSA columns and query bases**, length, support, mean similarities |
| `recombination_profile.tsv` | Parent-free signal: header with the PHI p-value and Rmin, then per-informative-site local incompatibility (the PHI profile) |
| `similarity_windows.tsv` | Full per-window matrix: `msa_position`, `query_position`, `winner`, and one similarity column per dataset |
| `similarity_stats.tsv` | Per-dataset similarity statistics (median, windows above identity thresholds) |
| `window_winners.tsv` | Per-dataset count of windows won (ties included) |
| `coverage_gaps.tsv` | Stretches where even the closest reference is a poor match — possible missing references |
| `similarity_top{N}.{fmt}` | Static plot of the nearest `--top-n` datasets, called regions shaded |
| `similarity_pair.{fmt}` | Static plot of the major vs leading minor parent, region shaded |
| `report.html` | Self-contained report: run provenance, the region table, the per-dataset stats, and an embedded interactive plot |

Similarity is computed only over columns where both sequences carry a canonical base (A/C/G/T); gaps, `N` and IUPAC ambiguity codes are ignored, so an `N` never counts as a match. A window with no comparable position — for example an inter-contig gap in a fragmented query — is uninformative and reported as `NA` in `similarity_windows.tsv` (and excluded from the winners, statistics and region calling).

Query coordinates are reported alongside MSA coordinates, so regions no longer need to be mapped back to the query by hand.

A plot is generated for the nearest datasets (`--top-n`, default 5), showing similarity across the MSA with called regions shaded.
![image](wiki/plot_x5.png) \
**Similarity in each window of the nearest five sequences to the query. Values towards 1 indicate high similarity. Here the query is most similar to a Cowpox sequence but has a region in the middle similar to a Variola sequence — a putative recombination event, called automatically and reported in `recombination_regions.tsv` in both MSA and query coordinates.**

The pairwise plot shows the major parent against the leading minor parent:
![image](wiki/plot_x2.png) \
**The two sequences most likely involved in the recombination, with the called region shaded.**

# Is a reference missing?
RecomFi always reports the *closest* reference, even when every reference is far
from the query — so a recombinant whose true donor is not in the collection is
still assigned to the least-bad reference. To catch this, the scan also reports
**reference coverage**: stretches where even the closest reference is below a
best-similarity threshold (adaptive by default, or set with `--coverage-floor`).
These are written to `coverage_gaps.tsv`, shaded on the report, and flagged with a
caveat banner; a region labelled `divergent` (ample comparable bases, yet a poor
match) is the signature of an absent reference, as opposed to `low_information`
(too few comparable bases to judge).

When a gap is found, `recomfi find-references` searches NCBI for the missing
genome: it BLASTs the under-covered query subsequence against `nt`, reports
candidate references (accession, identity, whether already in your collection),
and — with `--download` — fetches the best new one into a collection directory to
re-run with.
```
recomfi find-references --msa msa.fasta --query cowpox_with_variolaInsert \
    --collection collection/ --email you@example.org --download collection/
```
This contacts NCBI over the network; `--download` needs Entrez Direct
(`conda install -c bioconda entrez-direct`). BLAST often returns the query's own
GenBank record (the MSA labels the query by name, not accession), so a
near-identical, near-full-length hit is auto-skipped; use `--keep-self-hits` to
keep it, or `--exclude <accession>` to drop specific records.

A gap longer than `--subtile` (default 400 bp) is searched in overlapping
sub-intervals, not as one subsequence. A short divergent tract inside a longer gap
is otherwise diluted by the gap's high-matching flanks: one whole-gap BLAST returns
the flanks' donor and hides the tract's. (On a norovirus capsid gap this is the
difference between recruiting only the flank reference and also recruiting the
divergent core's distinct donor.) Resolution has two levers: `--window-size`
controls whether a short region is *detected* as a gap at all (a wide window
averages it away), and `--subtile` controls whether its donor is *found* once
detected. Set `--subtile 0` to disable; very short intervals (< ~150 bp) lose BLAST
specificity.

To do this repeatedly until coverage stops improving, `recomfi fill-references`
runs the whole cycle — build MSA, scan, BLAST, download — for several rounds,
growing a copy of the collection each time:
```
recomfi fill-references --query cowpox_with_variolaInsert.fasta.gz \
    --collection collection/ --output filled/ --aligner mafft --email you@example.org
```
It stops when the gaps close, when no new reference can be found, or when a round
no longer improves the worst gap (so a genuinely hypervariable region is reported,
not chased forever). Each round is recorded in `filled/fill_summary.tsv`, and the
query's own record is auto-excluded from its FASTA header. This needs an aligner
and Entrez Direct, and rebuilds the alignment every round.

Omit `--collection` to **start fresh** with no suggested references: the first
round seeds the collection from an NCBI search, then the loop fills the remaining
gaps as usual.
```
recomfi fill-references --query CRF01_AE.fasta --output filled/ --aligner mafft \
    --seed-mode parents --curate --email you@example.org
```
How the seed is chosen (`--seed-mode`) matters for a recombinant query. A whole-query
search returns the query's closest *whole-genome* relatives — which, for a
recombinant, are its own-lineage **siblings**; they cover the whole query and mask
the recombination. The parents that donated each region are *regional* best matches,
so by default seeding searches per region:

- `whole` — one whole-query search; the closest whole-genome relatives. Best for a
  non-recombinant query.
- `windowed` (default) — search each window (`--seed-window`, default 1500 bp) and
  seed the per-window best hits. Safe for any query; surfaces regional matches.
- `parents` — like windowed, but drop near-identical full-coverage hits (the
  siblings) so each region contributes its best *divergent* source. This recruits
  the parental lineages of a recombinant query; pair it with `--curate`. If a region
  has only siblings in NCBI, it falls back to the best hit so the panel is never
  empty.

`parents` mode helps when the parents are discoverable by BLAST — the common case
of a newly observed recombinant whose own form is still rare. When the query's own
lineage saturates `nt` (for example HIV-1 CRF01_AE, with thousands of genomes),
every region's top hits are siblings and the parents never surface; the run reports
this ("no distinct parental lineage found in NCBI") rather than masking it. BLAST
ranks by similarity, so it cannot escape an over-represented lineage. For that case
seed from a **finite pool** instead (`--seed-source`):

- `local` — a directory of genomes you control (`--candidate-pool dir/`: a RefSeq
  dump, a subtype reference set, or a lab collection). Because the pool is finite it
  is dereplicated and the query's siblings (near-identical genome-wide) are dropped,
  so the divergent parents are no longer out-ranked. The query is then matched region
  by region and the best-covering genomes are kept — a focused parental panel, offline
  and reproducible.
- `ncbi-virus` — a taxon-scoped set fetched with the `datasets` CLI
  (`--taxon "HIV-1"`, auto-detected if omitted). By default it pulls the RefSeq
  representative set (one genome per lineage, diverse by construction); add
  `--source-complete` for all complete genomes, which are then dereplicated.
- `nextclade` — a pool reconstructed from a Nextclade dataset
  (https://docs.nextstrain.org/projects/nextclade/en/stable/user/datasets.html).
  recomfi auto-detects the dataset from the query (`nextclade sort` when the CLI is
  installed, otherwise BLAST taxon detection mapped to a dataset), or you pass it
  explicitly with `--nextclade-dataset <path>` (e.g. `nextstrain/sars-cov-2/XBB`,
  `community/neherlab/hiv-1/hxb2`). Every reference-tree tip is reconstructed from
  the dataset reference plus its mutations and labelled by clade, so the report
  names parents by clade. Fetched pools are cached per dataset version.

```
recomfi fill-references --query CRF01_AE.fasta --output filled/ --aligner mafft \
    --seed-source local --candidate-pool subtype_refs/ --curate
```

```
recomfi detect --query CRF01_AE.fasta --output out/ --nextclade
recomfi fill-references --query q.fasta --output out/ --seed-source nextclade \
    --nextclade-dataset nextstrain/sars-cov-2/XBB
```
A Nextclade dataset is a clade-typing reference tree, so the pool spans clade
representatives rather than clean recombinant-parental endpoints; it suits a broad
screen. Single-segment datasets (most flu segments) yield segment-length pool
genomes, which is logged. The `nextclade` CLI is optional (it only sharpens dataset
auto-detection).

The sibling drop uses an absolute identity cutoff, which suits divergent-lineage
recombinants (HIV subtypes ~12 % apart). For closely related parents that have no
masking sibling (e.g. SARS-CoV-2 sublineages), add `--seed-keep-siblings`. NCBI
Virus / `datasets` needs `conda install -c conda-forge ncbi-datasets-cli`.

By default (`--auto-diversify`), BLAST seeding will **switch to the `ncbi-virus`
diversity path automatically** when it finds only siblings — i.e. when the query's
lineage saturates `nt` and no parental lineage can be recruited by similarity. A
broad fetch is capped (`--fetch-limit`, default 2000) and dereplicated; for a
heavily sequenced taxon the capped sample may miss lineages, so a curated
`--candidate-pool` is recommended (and the run says so). Whether the diversity
panel actually contains the parents depends on the taxon: it works when they are
genotype/lineage representatives, less so for fine genotype-specific recombinants.
Disable with `--no-auto-diversify`. This complements the caller-side defence:
`recomfi recomb` excludes whole-genome siblings from the competition
(`--exclude-siblings`, on by default), so a sibling that still slips into the panel
cannot win every window and mask the event.

## Curating the panel for detection
Auto-filling maximises *coverage* by recruiting the query's closest genomes — but
for a recombinant query those closest genomes are its own-lineage relatives, which
match it everywhere and **mask** the recombination. (Filling an HIV-1 CRF01_AE
query recruited CRF01_AE relatives that hid the underlying subtype-A/E mosaic.)
Detection needs *diversity*: distinct parental lineages, and none of the query's
siblings.

`recomfi curate-panel` curates a collection for detection. It uses
[skani](https://github.com/bluenote-1577/skani) to measure each reference's
genome-wide identity (ANI) and how much of the query it covers, then drops the
query's siblings and dereplicates near-duplicates with
[skDER](https://github.com/raufs/skDER):
```
recomfi curate-panel --query CRF01_AE.fasta --collection collection/ --output curated/
```
A *sibling* is a reference whose ANI to the query exceeds the backbone's by a
margin **and** that covers most of the query — a whole-genome relative. A regional
parent (for example an env-only reference) has high ANI but low coverage and is
kept. The rule is relative to the backbone (the query's closest whole-genome
match), so it needs no per-organism identity cutoff: it flags HIV subtype
relatives (~12 % apart) and SARS-CoV-2 sublineages (<1 % apart) alike. The curated
`curated/collection/` and a `panel_lineages.tsv` (each reference's role and
ANI/coverage) are written; rebuild with `recomfi msa` then `recomfi recomb`.

The same curation runs inside the fill loop with `fill-references --curate`, which
keeps the growing panel diverse and sibling-free each round and adds a "Reference
panel" section to the report. Both need skani (and skDER for dereplication):
`conda install -c bioconda skani skder`.

# Known limitations
The HMM caller's segmentation and the discordant-site sign test address the main
weaknesses of the old heuristic (window autocorrelation, the `--margin 0.0`
over-calling), but RecomFi is still an indicative screen, not a full phylogenetic
recombination test. It compares the query to a fixed reference panel rather than
inferring trees (so it cannot resolve which lineage is ancestral), uses a single
substitution model, and applies no genome-wide multiple-testing correction across
regions. Confirm strong candidates with a dedicated method (3SEQ, GARD, RDP). When
a region's donor is itself a poor match the result flags a possible missing
reference (see above) rather than a confident event.

# Development
```
pip install -e ".[dev]"
ruff check src tests
pytest                       # add -m "not requires_binary" to skip aligner-dependent tests
```
