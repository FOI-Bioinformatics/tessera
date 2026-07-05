# Changelog

All notable changes to Tessera are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **`tessera reassort` -- per-segment reassortment detection for segmented viruses**
  (influenza, bunyaviruses, reoviruses). Reassortment -- a whole-segment swap between
  parents -- cannot be represented by the single-backbone intragenic scan; `reassort`
  types each segment of a multi-FASTA query against its own Nextclade dataset and calls
  reassortment when the segments trace to different parents. Writes `reassortment.tsv`
  and `constellation.tsv` and prints a genotype mosaic plus a
  `clonal` / `reassortant` / `undetermined` verdict.
- **Coverage-aware reassortment call with a parent constellation.** The verdict is a
  pairwise, coverage-aware comparison (concordant / discordant / uninformative over an
  ANI margin) rather than a single strain intersection, so a clonal isolate whose parent
  strain is not cross-typed in every segment's tree is not miscalled a reassortant as the
  segment count grows. Segments are grouped into parent constellations. Adds a `--margin`
  option; ranking uses skani ANI and alignment fraction.
- **`reassort --scan-segments` -- intragenic per-segment scan.** After assignment, each
  assigned segment is aligned to its own per-clade-consensus panel and run through the
  ordinary recombination scan (`run_recomb`), writing the full per-segment output to
  `out/<segment>/` and a rollup to `out/segment_scan.tsv`. Needs an aligner (`--aligner`,
  default `mafft`), checked up front; a per-segment scan failure is non-fatal.
- **Opt-in validation probe** `validation/run_reassort_scan.py` that builds a real-tip
  insert-shaped recombinant plus a clonal control and checks that `--scan-segments`
  localizes the intragenic recombination (localization-gated; attribution and specificity
  reported). Not part of CI.
- **`--method geneconv` -- an opt-in GENECONV-style caller** (Sawyer 1989): over the
  triplet discriminating sites it scores the longest uninterrupted run of consecutive
  donor-matches (a clean gene-conversion fragment), with a permutation p-value and
  Benjamini-Hochberg across donors. Complementary to 3SEQ's drawdown and MaxChi's boundary
  chi-square; not in the default ensemble.
- **External method comparison and short-tract probe (opt-in validation).**
  `validation/run_method_comparison.py` runs Tessera's callers beside OpenRDP (the
  maintained RDP5 reimplementation) on the same published-recombinant alignments, in one
  combined table; `validation/run_method_comparison_hybrids.py` measures short-tract
  sensitivity against ground truth, including an adversarial sub-window-tract tier. Both
  need external tools and are not part of CI.
- **Benchmark scoring against published datasets (opt-in validation).** A set of harnesses
  that score detection on real and simulated recombination, recording measured results:
  - **Real published recombinants** (`validation/run_validation.py`, `datasets.json`),
    extended with the HCV 2k/1b inter-genotype recombinant (a new virus family,
    Flaviviridae, and the first precise real-breakpoint check, ~nt 3187) and a real
    clonal-negative control, scored by a new `max_regions` specificity key.
  - **Reassortment precision / recall / F1** (`validation/run_reassort_benchmark.py`) for
    the shipped `reassort` verdict, the way the influenza-reassortment literature scores
    it (a labelled clonal-vs-reassortant flu H3N2 HA+NA query set).
  - **Simulated benchmarks:** a SANTA-SIM power/specificity set (Jaya 2023,
    `run_benchmark.py`), a coalescent-with-recombination design (Posada & Crandall 2001,
    `run_coalescent_benchmark.py`, msprime), and a breakpoint noise-robustness design
    (RecombinHunt 2024, `run_recombinhunt_benchmark.py`). The pure scorers are unit-tested;
    the harnesses need aligner / simulator environments and are not part of CI.
- **Broadened specificity panel (opt-in validation).** The synthetic-hybrid suite now carries
  nine non-recombinant `neg_pure` controls across the divergence and panel-size axes
  (was one), so the false-positive rate is measured over many pathogens rather than a
  single case.

### Fixed

- **`build_pool(per_clade_consensus=True)` crashed on within-clade indels.** The
  per-clade consensus reconstructed tips with deletion gaps stripped, so unequal-length
  tips reached `consensus_sequence` and raised. It now reconstructs in gapped reference
  coordinates before the consensus. This also repairs `detect`'s consensus-pool mode on
  datasets with indels (e.g. influenza).
- **Example-genome clade label.** Nextclade example sequences were written with a mined
  clade prefix that could masquerade as a real clade; they now carry an `example` marker
  and are treated as untyped by the pool consumers (`reassort`, `type-lineages`).
- **`reassort` now rejects duplicate segment names** up front, rather than silently
  dropping a colliding segment (which could manufacture a false reassortant).
- **Transient typing failures surface.** A download / network error during segment typing
  now propagates instead of being reported as a biological "unassigned"; a genuine
  no-dataset resolution or a skani rejection of a short segment stays a non-fatal
  per-segment skip.
- **Short recombinant tracts (below one sliding window) were missed.** The site-based
  callers (3SEQ, MaxChi, GENECONV) drew candidate donors only from window winners, so a
  sub-window tract -- whose donor is diluted in the windowed vote and wins no window --
  was never tested, despite an often overwhelming discriminating-site signal. They now
  draw candidates from all references (`rank_datasets`), with the existing significance
  floors dropping the rest, so short tracts are recovered at no measured specificity cost.
- **Backbone mis-identified when the query is nearly identical to its parent.** For a very
  short tract the query is almost identical to its true backbone, so the HMM's whole-genome
  sibling exclusion could drop that backbone as a twin and pick a genotype that wins zero
  windows as the major -- which `reconcile_major` then adopted, mislabelling every region's
  backbone. `reconcile_major` now ignores a HMM major that wins zero windows when the
  windowed vote has a clear winner; a HMM major that wins any window (including the
  masking-twin case the exclusion exists for) is left untouched.
- **Over-calling from the widened donor search (FDR control).** After the candidate-selection
  change above, every caller still gated on the raw p-value while the Benjamini-Hochberg
  q-value was computed but unused -- so testing more candidate donors raised the false-call
  rate (the broadened negative panel measured 10 false calls across the controls). All four
  callers (3SEQ, MaxChi, GENECONV, and the HMM segment test) now gate on the q-value, so the
  false-call rate no longer scales with the number of candidates scanned. Measured: detection
  sensitivity unchanged (30/30 on the hybrid suite, including the low-divergence case), with
  fewer false positives.

## [2.0.0]

The current declared version (`pyproject.toml`); the rebuild of the tool as the `tessera`
package. Not yet tagged as a dated release, so the entries below are grouped here rather than
in `[Unreleased]` for readability. Renamed the package and CLI from **`recomfi`** to
**`tessera`** (`tessera` entry point, `src/tessera/` layout).

### Added

**Core engine.** Rebuilt as a dependency-free Python package with a Typer CLI, an entry-point
aligner registry, and a native numpy sliding-window similarity engine (no external `recan`
dependency). Builds a reference-anchored pseudo-MSA from XMFA / MAF / HAL, excludes ambiguous
bases (N / IUPAC) from comparability, and reports uninformative windows as `NA`.

**Commands.**

- `detect` -- one-shot, zero-config entry point: taxon detection, diverse panel recruitment,
  alignment, and recombination calling, with a plain-language confidence verdict and a run log.
- `build-panel`, `recomb`, `msa` -- panel building, detection, and alignment as separate steps
  (re-run detection with different settings without re-fetching or re-aligning).
- `find-references`, `fill-references` -- BLAST-based reference discovery and the iterative
  build-scan-BLAST-download gap-filling loop (`fill_summary.tsv` per round).
- `curate-panel` -- skani sibling-drop (backbone-relative, no per-organism cutoff) plus skDER
  dereplication; also available as `fill-references --curate`.
- `type-lineages` -- assign a lineage to every genome in a collection via a ladder
  (header / metadata designation, nearest Nextclade-dataset tip by skani ANI, de-novo ANI
  clustering); writes a `lineages.tsv` sidecar the other commands read.

**Detection methods.** A default four-caller **ensemble** (HMM segmentation, 3SEQ, MaxChi,
Bootscan) with a lineage-aware consensus merge; each caller is also selectable via `--method`.

- HMM segmentation (jpHMM-style) with a discordant-site sign test at `--alpha`, support scores,
  Benjamini-Hochberg q-values, and breakpoint-uncertainty intervals.
- 3SEQ maximum-drawdown triplet test, MaxChi chi-square triplet test, and Bootscan
  bootstrap-support caller.
- Parent-free **PHI** (Pairwise Homoplasy Index) and Hudson-Kaplan **Rmin** diagnostics that
  fire even when the true donor is absent from the panel (`recombination_profile.tsv`).
- Opt-in clade-barcode (lineage-marker) caller for typed panels.
- Informative-site windowing for low-divergence panels (auto-enabled below ~8% between-reference
  divergence), so windows span polymorphic columns rather than base pairs.

**Reference panels.**

- Seed sources: BLAST (`--seed-mode whole|windowed|parents`), a finite local pool (skDER
  dereplication + skani region matching), NCBI Virus (taxon-scoped `datasets` download, cached),
  and Nextclade (reference-tree tips reconstructed and clade-labelled, dataset auto-detected or
  pinned). Auto-diversify and negative-lineage BLAST recruit divergent parents when a lineage
  saturates NCBI nt.
- Sibling exclusion and lineage-cluster competition in the caller so near-identical relatives
  cannot mask an event.
- Lineage-aware panel selection; typed lineage names in the report; query self-typing; a
  SARS-CoV-2 Pango recombinant cross-check; `--deep-typing` (full typing ladder on the recruited
  panel); `--reattribute-donors` (opt-in post-hoc donor re-labelling); `--pool-consensus`
  (one denoised consensus genome per clade).

**Aligner backends.** mafft (base-level, `--addfragments` for fragmented queries), minimap2,
SibeliaZ (default), progressiveMauve, and cactus, via a pluggable entry-point registry.

**Output.** A self-contained HTML report (query-mosaic track, called-regions table, window-winner
chart, per-dataset stats, embedded similarity plot, methods glossary); `recombination_regions.tsv`
(query + MSA coordinates, support, p / q-values, breakpoint interval, calling methods),
`recombination_methods.tsv`, `recombination_profile.tsv`, `coverage_gaps.tsv`,
`panel_lineages.tsv`, and a persisted run log.

### Fixed

- MAF reverse-strand projection scrambled SibeliaZ / cactus coordinates (~43% of blocks place the
  reference on the `-` strand); blocks are now placed at their true forward-reference coordinate.
- `find-references` returned the query's own record as a self-hit; near-identical near-full-length
  hits are auto-skipped (`--keep-self-hits` to override, `--exclude` to drop accessions).
- HMM emission produced degenerate likelihoods on some near-identical panels (the CRF01_AE
  env-region case).
- The Nextclade examples fetch was fatal when a dataset had no examples or the network failed; it
  is now non-fatal, and clade labelling prefers `LANL_subtype` where present.
- A coverage gap overlapping a low-donor-similarity region is now flagged as donor-undercovered
  rather than reported as a confident call.
- NCBI Virus full-download for mega-taxa (HIV cold-start): a spurious `--limit`, an `ARG_MAX`
  overflow from inline skDER paths, and skDER representative symlinks resolved after their temp
  directory was removed.
