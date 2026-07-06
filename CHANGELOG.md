# Changelog

All notable changes to Tessera are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Nothing yet -- the initial release below is not yet tagged, so all current changes are recorded
under [1.0.0].

## [1.0.0] - 2026-07-06

The initial `tessera` release: the rebuild of the `recomfi` prototype as the `tessera` package --
package and CLI renamed from **`recomfi`** to **`tessera`** (`tessera` entry point, `src/tessera/`
layout). `recomfi` was never released, so the first public version is **1.0.0** (a stable first
release), not a continuation of a prior series.

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
- `reassort` -- per-segment reassortment detection for segmented viruses (see **Reassortment**).

**Detection methods.** A default four-caller **ensemble** (HMM segmentation, 3SEQ, MaxChi,
Bootscan) with a lineage-aware consensus merge; each caller is also selectable via `--method`.

- HMM segmentation (jpHMM-style) with a discordant-site sign test at `--alpha`, support scores,
  Benjamini-Hochberg q-values, and breakpoint-uncertainty intervals.
- 3SEQ maximum-drawdown triplet test, MaxChi chi-square triplet test, and Bootscan
  bootstrap-support caller.
- Opt-in GENECONV-style clean-fragment caller (`--method geneconv`, Sawyer 1989): the longest
  uninterrupted run of donor-matches over the triplet discriminating sites, with a permutation
  p-value and Benjamini-Hochberg across donors. Not in the default ensemble.
- Parent-free **PHI** (Pairwise Homoplasy Index) and Hudson-Kaplan **Rmin** diagnostics that
  fire even when the true donor is absent from the panel (`recombination_profile.tsv`).
- Opt-in clade-barcode (lineage-marker) caller for typed panels.
- Informative-site windowing for low-divergence panels (auto-enabled below ~8% between-reference
  divergence), so windows span polymorphic columns rather than base pairs.
- Site callers draw candidate donors from all references (not only window winners) and gate on the
  Benjamini-Hochberg q-value, so short sub-window tracts are recovered while the false-call rate
  does not scale with the number of candidates scanned.

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

**Reassortment.** For segmented viruses (influenza, bunyaviruses, reoviruses); a whole-segment swap
between parents cannot be represented by the single-backbone intragenic scan.

- `reassort` types each segment of a multi-FASTA query against its own Nextclade dataset and calls
  reassortment when the segments trace to different parents; writes `reassortment.tsv` and
  `constellation.tsv` and prints a genotype mosaic plus a `clonal` / `reassortant` / `undetermined`
  verdict. The verdict is a pairwise, coverage-aware comparison (concordant / discordant /
  uninformative over an ANI margin, `--margin`), so a clonal isolate whose parent strain is not
  cross-typed in every segment's tree is not miscalled a reassortant as the segment count grows;
  segments are grouped into parent constellations (ranking uses skani ANI and alignment fraction).
- `reassort --scan-segments` aligns each assigned segment to its own per-clade-consensus panel and
  runs the ordinary recombination scan (`run_recomb`), writing the full per-segment output to
  `out/<segment>/` and a rollup to `out/segment_scan.tsv` (needs an aligner, `--aligner` default
  `mafft`; a per-segment scan failure is non-fatal).

**Validation harnesses (opt-in, not part of CI).** Harnesses that score detection on real and
simulated recombination against a documented expectation, recording measured results.

- `validation/run_reassort_scan.py` -- builds a real-tip insert-shaped recombinant plus a clonal
  control and checks that `reassort --scan-segments` localizes the intragenic recombination
  (localization-gated; attribution and specificity reported).
- `validation/run_method_comparison.py` -- runs Tessera's callers beside OpenRDP (the maintained
  RDP5 reimplementation) on the same published-recombinant alignments, in one combined table;
  `validation/run_method_comparison_hybrids.py` measures short-tract sensitivity against ground
  truth, including an adversarial sub-window-tract tier.
- Benchmark scoring against published datasets: real published recombinants
  (`validation/run_validation.py`, `datasets.json`, including the HCV 2k/1b inter-genotype
  recombinant -- a Flaviviridae case with a precise ~nt 3187 breakpoint -- and a real clonal-negative
  control scored by a `max_regions` key); reassortment precision / recall / F1
  (`validation/run_reassort_benchmark.py`); and simulated benchmarks -- SANTA-SIM power/specificity
  (Jaya 2023, `run_benchmark.py`), a coalescent design (Posada & Crandall 2001,
  `run_coalescent_benchmark.py`, msprime), and breakpoint noise-robustness (RecombinHunt 2024,
  `run_recombinhunt_benchmark.py`). The pure scorers are unit-tested.
- Synthetic-hybrid harness (`validation/run_hybrids.py`): a broadened specificity panel of nine
  non-recombinant `neg_pure` controls across the divergence and panel-size axes, plus two
  `neg_within` intra-clade-splice controls (the panel keeps both same-clade sources so the splice
  is credited intra-clade). `neg_hiv1` is reported as an investigated `KNOWN-LIMIT` (not gating): a
  pure HIV subtype genome is statistically indistinguishable from a real recombinant against a
  source-removed subtype panel, so it cannot be suppressed without breaking genuine detection.
  Datasets with no Nextclade clade attribute are typed de-novo by ANI-clustering their reconstructed
  tips, making five of six previously-skipped datasets runnable.

### Changed

- CI installs an aligner (mafft) and runs the `requires_binary` tests, so the MSA-build and
  end-to-end recombination paths are exercised in CI rather than only locally.

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
- `build_pool(per_clade_consensus=True)` crashed on within-clade indels: the per-clade consensus
  reconstructed tips with deletion gaps stripped, so unequal-length tips reached
  `consensus_sequence` and raised. It now reconstructs in gapped reference coordinates first, which
  also repairs `detect`'s consensus-pool mode on datasets with indels (e.g. influenza).
- Nextclade example sequences were written with a mined clade prefix that could masquerade as a real
  clade; they now carry an `example` marker and are treated as untyped by the pool consumers.
- `reassort` rejects duplicate segment names up front, rather than silently dropping a colliding
  segment (which could manufacture a false reassortant).
- A download / network error during segment typing now propagates instead of being reported as a
  biological "unassigned"; a genuine no-dataset resolution or a skani rejection of a short segment
  stays a non-fatal per-segment skip.
- Short recombinant tracts (below one sliding window) were missed: the site-based callers drew
  candidate donors only from window winners, so a sub-window tract's donor -- diluted in the windowed
  vote and winning no window -- was never tested despite an often overwhelming discriminating-site
  signal. They now draw candidates from all references, with the significance floors dropping the
  rest, recovering short tracts at no measured specificity cost.
- Backbone mis-identified when the query is nearly identical to its parent: for a very short tract
  the HMM's whole-genome sibling exclusion could drop the true backbone as a twin and pick a genotype
  that wins zero windows as the major, which `reconcile_major` then adopted. `reconcile_major` now
  ignores a HMM major that wins zero windows when the windowed vote has a clear winner; a HMM major
  that wins any window (including the masking-twin case the exclusion exists for) is left untouched.
- Over-calling from the widened donor search: every caller gated on the raw p-value while the
  Benjamini-Hochberg q-value was computed but unused, so testing more candidate donors raised the
  false-call rate. All four callers (3SEQ, MaxChi, GENECONV, and the HMM segment test) now gate on
  the q-value; detection sensitivity is unchanged (30/30 on the hybrid suite, including the
  low-divergence case) with fewer false positives.
