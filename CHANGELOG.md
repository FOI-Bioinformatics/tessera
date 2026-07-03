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
