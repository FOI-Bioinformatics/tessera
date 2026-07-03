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
