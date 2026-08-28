# Changelog

All notable changes to Tessera are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **`validation/run_caller_benchmark.py`** -- scores the *region callers* on the Jaya et al. (2023)
  published simulated grid (Dryad doi:10.5061/dryad.d7wm37q6f), the part of Tessera comparable to
  that study's sequence-level methods and previously never benchmarked. Reports a per-query rate
  (which does not compound with the number of sampled queries) alongside the alignment-level one.
  It measured a previously undocumented **upper** divergence bound: at the grid's highest mutation
  rate the callers fire on 33 % of scans over clonal data and barely separate them from recombining
  ones, so a high call rate on a saturated panel is noise rather than power. Documented in
  `docs/detection-methods.md`, "The divergence window".
- **`validation/run_specificity.py`** -- a harness that measures the false-positive *rate* on
  simulated clonal data, where the ground truth is known by construction: no recombination is
  introduced anywhere, so every region reported is a false positive. Four scenarios escalate the
  confounders the literature blames for spurious signal (uniform rates, among-site rate variation,
  a fast query lineage, and a rate shift between genome halves), and a paired positive control
  reports sensitivity beside specificity. Rates carry Wilson confidence intervals. Unlike the other
  harnesses it needs no aligner, no network and no downloaded data; the simulation and scoring
  logic is unit-tested in CI. This closes gap **G3** ("no specificity measurement") from
  `docs/gap-analysis.md`: the nine `neg_pure` controls measure a verdict on nine favourable real
  panels, which can neither estimate a rate nor sample the redundant-panel regime where the
  callers are most likely to over-call.

### Fixed

- **`recombination_regions.tsv` rows now name the test behind their own `pvalue`.** The four
  callers do not share a statistic -- a one-sided sign test on discordant sites, an exact or
  permuted max-descent test, a scan-aware chi-square permutation, and a block-permutation run
  length all wrote into one `pvalue` column, with `support` variously a discordant-site
  fraction, a tract-match fraction or mean bootstrap support. Only `methods` hinted at which,
  and for a merged ensemble row that hint was wrong: the merge keeps the single most
  significant member's p-value, so a region listing four agreeing callers carried a p-value
  from one of them with no way to tell which. New `test` and `statistic` columns state it per
  row, and the ensemble merge carries the retained member's names along with its p-value. For
  3SEQ the column also distinguishes the exact test from the permutation fallback.
- **`length_bp` had two definitions in one column.** Every caller set it from the MSA span
  (`msa_end - msa_start`) and the ensemble merge from the query span (`query_end -
  query_start`), so its meaning depended on which caller produced the row. It is now always
  query bases, matching its name and the adjacent query coordinates, with the alignment width
  reported separately as `length_msa`. Both are derived from the coordinate columns rather
  than stored, so they cannot drift apart again.
- **p- and q-values are written at full precision.** They were rounded to two significant
  figures on the way out, which discarded precision a reader cannot recover and made an
  underflow indistinguishable from an exact zero. The gate that decides which regions are
  reported always used the untruncated value, so no detection behaviour changes -- only what
  the file records (`1.4e-12` becomes `1.422084672242363e-12`). The report still formats them
  for display.
- **The scope of the FDR correction is now recorded.** Benjamini-Hochberg is applied within
  one caller's own candidate segments, not across callers or across the genome, so an ensemble
  run has no single genome-wide false-discovery rate. `run_provenance.json` states this under
  `multiple testing`, and `docs/detection-methods.md` explains how to read the column. This
  documents the existing behaviour; it does not change the statistics (see gap G4).
- **An interrupted NCBI Virus fetch no longer becomes a permanent partial panel.** Genomes were
  written one file at a time straight into the shared cache directory, so a run killed part way
  through -- Ctrl-C, OOM, a dropped connection -- left a truncated set behind. The next run saw a
  non-empty directory, logged "Using the cached NCBI Virus panel" and carried on. Panel
  composition decides which donors can be found at all, so this did not produce a slow run but a
  wrong one, silently narrowing what detection could report. The fetch is now staged in a sibling
  directory and installed with `os.replace` only once complete, and each entry carries a
  `tessera_cache.json` manifest; a cache directory without one is ignored rather than trusted,
  which also rejects partial caches left by earlier versions. Two concurrent runs fetching the
  same panel no longer race: the first to finish wins and the second reuses it.
  **This invalidates existing NCBI Virus caches** -- they carry no manifest and will be refetched
  once. The Nextclade pool cache already worked this way and is unaffected.
- **The NCBI Virus cache key now covers the fetch scope.** It was `sha1(taxon)` alone, so a
  `--source-refseq` run and a default run shared one slot and whichever ran second silently
  received the other's differently-scoped panel.
- **A result can now name the aligner that produced it.** `core/binaries.py` exists to resolve
  tool versions for provenance, and `build_msa` called `preflight()` and then only logged the
  answer, so the scan recorded the alignment as a bare path: given a `report.html` you could not
  tell which aligner, at which version, with which arguments had built the alignment it
  described. `tessera msa` now writes a `<msa stem>.provenance.json` sidecar (aligner, resolved
  versions, `--aligner-arg` values, backbone, genome count) and `tessera recomb` picks it up
  automatically, showing it in the report and writing the whole run record to
  `run_provenance.json` beside the TSVs. An alignment you built yourself has no sidecar and the
  record says so rather than guessing. The run also states that the permutation callers are
  deterministically seeded, so a reader can tell a reproducible p-value from a lucky one.
- **A panel too small to test is no longer reported as "no recombination".** Recombination is a
  switch between donors, so it takes two references before there is anything to switch between.
  Below that every caller early-returned on `len(labels) < 2`, `reconcile_major` yielded no major
  parent and the coverage threshold collapsed to zero, so the scan finished successfully and
  reported nothing found -- a result indistinguishable from a genuine clean negative on an
  adequate panel. `run_recomb` now refuses an alignment carrying fewer than two references and
  says why. The similarity engine itself is unchanged: computing per-window identity against a
  single reference is still a meaningful operation, it is *calling recombination* on it that is
  not.
- **Duplicate sequence identifiers are rejected instead of silently discarding data.**
  `Bio.AlignIO` preserves repeated headers, but the alignment was stored as a label to sequence
  mapping, so only the last record with a given name survived: panel members vanished without a
  message, the reported dataset count disagreed with the input FASTA, and a duplicated *query*
  label meant the row actually analysed was whichever came last. Reading now collects duplicates
  and raises, naming them.
- **A ragged (unaligned) FASTA now names the offending record.** The check that does so was
  unreachable: `AlignIO.read` raised `Sequences must all be the same length` first, which
  identifies no record and surfaced through the CLI as an "Unexpected error". Records are read
  with `SeqIO.parse` so Tessera's own message -- which gives the record and both lengths -- is
  the one that fires. Non-ASCII input is likewise reported against the record rather than as a
  raw `UnicodeEncodeError`.
- **Bootscan is now gated on a null model.** It previously carried no significance test at all:
  `benjamini_hochberg` was applied in the HMM, 3SEQ and MaxChi callers but absent from
  `bootscan.py`, whose regions were built with neither `pvalue` nor `qvalue` and gated only on a
  fixed 0.70 bootstrap-support threshold. Bootstrap support answers *which reference is closest
  here*, not *whether the query recombined*. Each candidate run of won windows is now tested
  against a permutation null holding the win count fixed, permuted in blocks one window wide so
  the autocorrelation of overlapping windows survives into the null, then BH-corrected and gated
  on the q-value like the other callers.
- **`parent_free_support` no longer fires on recurrent mutation.** The flag was set from the
  Hudson-Kaplan Rmin intervals alone, but Rmin bounds recombination only under the
  infinite-sites model -- under finite sites, recurrent mutation produces four-gamete violations
  on strictly clonal data (measured: Rmin 13-991 across a 2-24 % divergence ladder with no
  recombination anywhere, versus 0 throughout under infinite sites). The flag consequently
  corroborated almost every false positive. It now requires the **PHI test** to be significant
  as well -- PHI establishes that the alignment carries recombination, the Rmin intervals say
  where. On a clonal panel (PHI p = 0.16, Rmin = 292) the flag is now `no`; on a true recombinant
  (PHI p = 0.001, Rmin = 312) it remains `yes`.
- **Reassortment: the candidate cap no longer truncates the near-best ANI window.** `TOP_K` cut
  each segment's candidate list by *rank* before `constellation._near_best` applied the ANI
  `margin`, so when more than `TOP_K` strains sat inside that window a genuinely shared parent
  could be dropped and two segments of one clonal strain read as having different parents -- a
  clonal isolate reported as a reassortant. `validation/README.md` recorded this as an inherent
  cross-typing-coverage limit; it was not. On the seeded flu H3N2 benchmark the false positives
  go from 3 to 0 (precision 0.82 -> 1.00, F1 0.87 -> 0.89), at the cost of recall 0.93 -> 0.80 --
  the old behaviour gained sensitivity by accident, from an arbitrarily narrowed near-best set.
  `--margin` is the principled lever for that trade.
- **Window spacing is no longer read off a single sampled gap.** `clusters.py` and `siblings.py`
  converted a base-pair threshold (`--min-region`) into a window count using
  `positions[1] - positions[0]`. That is exact under base-pair windowing, where every gap equals
  `window_step`, but wrong under **informative-site windowing** -- automatic below ~8 % divergence,
  so precisely for the low-divergence pathogens (mpox, VZV, ebola, within-Omicron) -- where
  positions are midpoints of informative-site windows and their spacing tracks how densely
  polymorphic sites fall. On a constructed alignment the true spacing ranged 9-765 columns while
  the code used 13 throughout. Both sites now share `clusters.median_window_step`. Since the
  affected code decides lineage clustering and sibling exclusion, the change was validated on the
  real 24-pathogen harness: sensitivity 30/30, specificity 10/10, 0 false calls.
- **RNA-alphabet input is no longer silently mis-read.** `U` was not a canonical base and no
  `U`->`T` normalisation existed anywhere, in a tool whose target domain is RNA viruses. A single
  U-alphabet record silently changed the result for *every* reference -- in one check it cut a
  neighbouring reference's usable windows from 25 to 9 and pushed a third's median similarity from
  0.5 to 0.0 -- with no warning, because it also perturbs the automatic informative-site switch.
  `U` is now folded to `T` on read, and a record whose non-gap positions are mostly not A/C/G/T is
  reported as a warning rather than quietly contributing nothing.
- **Donor re-attribution no longer bypasses its own margin guard.** `best_score - (cur or 0.0)`
  collapsed an *unscorable* current donor -- one absent from the lineage map, or in a span with
  too few comparable sites -- to a similarity of 0.0, making the margin test vacuous: a donor
  matching the query at 98.2 % was replaced by a clade matching at 58.7 %. Re-attribution now
  declines when the current donor cannot be scored.

### Changed

- **Ensemble agreement can now be used as a gate** (`--min-methods`), the control RDP5 exposes as
  "list events detected by more than N methods". It is **off by default** (`1`, the previous union
  behaviour). Measured both ways: on a redundant panel -- several near-equidistant relatives of the
  query -- a gate of 2 cuts false regions on simulated clonal data from 21 to 2; but on the curated
  panels the hybrid harness builds it loses three true detections (`rsv_a`, `mpox`, `masksib_rsv`)
  and gains nothing, since every negative control passes without it. The callers have different
  applicability domains -- at mpox's 0.5 % divergence only the HMM has power -- so requiring two is
  structurally impossible in the regime where detection is hardest. Raise it when your panel
  carries near-duplicates; leave it at 1 for a curated panel. Clamped to the number of callers
  actually run, and any suppressed count is logged rather than silently dropped.

Measured on simulated clonal panels (4 clades, 16 tips, JC69, no recombination anywhere), 40 runs
of the default ensemble, counting donor-present regions as `validation/run_hybrids.py` scores a
negative control:

| | before | bootscan null (default) | + `--min-methods 2` |
|---|---|---|---|
| runs reporting a false region | 39/40 | 17/40 | 2/40 |
| false regions | 139 | 21 | 2 |
| of which bootscan | 132 | 0 | 0 |

On the real 24-pathogen hybrid harness the bootscan null alone carries the specificity:
**sensitivity 30/30, specificity 10/10 with 0 false calls**. Adding `--min-methods 2` there
costs three true detections and gains nothing, which is why it is not the default.

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
