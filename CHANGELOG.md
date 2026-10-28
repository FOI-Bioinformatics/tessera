# Changelog

All notable changes to Tessera are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### A draft query could not be aligned at all by MAFFT 7.525

- **The defect.** MAFFT 7.525 segfaults inside `pairlocalalign` whenever
  `--addfragments` is given more than one record and the pairwise-local path is taken.
  Tessera adds every genome with `--addfragments`, and only reaches for `--6merpair`
  above 20 kb, so on that release **no** multi-contig query of a smaller genome could be
  aligned: the run died with the aligner's exit 1 and a segfault line.
- **It is the release, not the input.** Measured against 7.525 and 7.526 here: two
  300 bp fragments crash as readily as three 2 kb ones; it crashes with and without
  `--adjustdirection` and with and without `--keeplength`, so neither the strand
  handling nor the length rule is involved; a single-record file is fine; and 7.526
  aligns every one of those inputs cleanly.
- **The fix.** Where the reported MAFFT version is older than 7.526, a genome holding
  more than one record is added with `--6merpair`, which derives its guide distances
  from k-mer counts and never enters the crashing routine. It is the same escape the
  memory rule already uses, and on anything above 20 kb -- most real genomes -- Tessera
  asks for it anyway. The run logs which genomes were affected and that upgrading MAFFT
  restores the pairwise alignment. A version string that cannot be parsed is left on the
  pairwise path: a clear failure from the tool beats quietly coarsening every alignment.
- **Why it mattered beyond a draft query.** `ubuntu-latest` installs 7.525 from apt, so
  this had been failing CI's `requires_binary` step on the default branch since
  2026-10-25, on `test_draft_query_with_one_reversed_contig_aligns`. The same test now
  passes under 7.525 and under 7.526, verified against both.

### Text IO is UTF-8 whatever the locale says

- **The defect.** `report.html` declares UTF-8 in its own meta tag, but it and every
  TSV were written in the locale's encoding, and the alignment was handed to Biopython
  as a path, which opens it the same way. Under `LC_ALL=C` -- an ordinary setting on a
  compute node -- an alignment carrying a non-ASCII label could not be read at all, and
  the report could not be written even without one, because the plotly bundle inlined
  in it is not ASCII. The run died at its last step.
- **The fix.** Every text read and write in the package names its encoding, the
  alignment is opened by Tessera rather than by Biopython, and ruff's
  `unspecified-encoding` rule is switched on so a new call site cannot reintroduce it
  (`explicit-preview-rules` keeps the rest of the preview set out). A test runs a whole
  analysis under `LC_ALL=C` on an alignment whose donor is labelled `Ebola/Guéckédou`,
  in a child interpreter because the locale is fixed at startup; it fails on the
  previous commit and passes here.

### Four ways a broken or unexpected input ended a run

- **A Pango alias key that is not JSON was cached before it was parsed**, so a proxy
  error page or a truncated read stayed in the cache and every later call -- online or
  not -- warned that the cache was unreadable and returned nothing. One bad response
  disabled the lineage cross-check until somebody deleted the file by hand. The body is
  parsed first and only then cached.
- **A Nextclade dataset that publishes no tree raised a bare `KeyError`**, which
  reached the user as a traceback and, through `reassort`, ended a whole multi-segment
  run on one tree-less segment dataset. It is refused by name, naming the dataset and
  the file it lacks. A cut transfer (`IncompleteRead`) is an `HTTPException` rather than
  an `OSError`, so it escaped the handler that exists for it; both download paths now
  catch it.
- **A cached BLAST result whose row fields no longer match `Hit` raised `TypeError`**,
  which is not the `BlastError` the callers catch to skip a gap and carry on, so the run
  ended -- for the fortnight the entry stays valid. An unreadable entry is now a cache
  miss: the search runs again and replaces it.
- **A FASTA header or an accession named a file without being sanitised.** An empty
  header raised `IndexError` while splitting a downloaded archive, one beginning `../`
  wrote outside the destination directory, and `efetch` took the same route. Both go
  through `safe_filename_stem`, which `reassort` already used for this. The inputs come
  from NCBI rather than from a user, so this is defence in depth.
- **An unnamed record could take a named segment's place.** `reassort` named an unnamed
  record `segment_N` by position, which could equal a real header: the two then shared
  one key, so one was never scanned and the other was typed twice under the same name.
  The synthetic name now skips names the query already uses.

### Inputs that were accepted and then ignored, or accepted and then failed

- **`msa --output <directory>`** warned "Overwriting the existing alignment", ran the
  aligner to completion and only then failed moving the result into place, with
  `Is a directory` under "Unexpected error". Every other command takes a directory
  there, so this is an easy slip, and it cost the whole alignment. It is rejected
  before any work starts, as are a `--collection` that is not a directory and a
  missing `--query`.
- **`find-references --collection <typo>`** was silently ignored, so every BLAST hit
  read as a reference the collection lacked -- the opposite of what the command
  reports. A `--download` path that is an existing file is checked too.
- **`recomb --min-region 0`, `--min-region -5`, `--merge-gap -500`** all ran to a clean
  exit 0 having quietly meant something else. They are range-checked like every other
  numeric option on the command.
- **A taxon was mapped to a Nextclade dataset by bare substring**, so on the
  BLAST-detection fallback (no `--nextclade-dataset`, no `nextclade` binary) "human
  immunodeficiency virus 2" resolved to the HIV-1 dataset and "human respiratory
  syncytial virus b" to the RSV-A one: the wrong subtype's clade tree became the
  recruitment pool without a word. Where a keyword names one member of a family whose
  siblings have their own datasets, a serotype or subtype token following it must be one
  that dataset covers; an unqualified name keeps its old mapping, and a longer trailing
  word ("strain", "isolate") is not a qualifier. An unmappable taxon is refused by name
  rather than guessed at.
- **The BLAST client's `Entrez.email` assignment no longer depends on the mypy
  release.** Biopython declares the attribute as `None`, so assigning an address is a
  type error -- which mypy 2.4 reports and 2.3 does not. With `warn_unused_ignores` on,
  a `type: ignore` fails under whichever of the two does not need it, so the check's
  result depended on which patch release `mypy>=1.11,<3.0` happened to resolve to. It
  goes through `setattr` instead, which suppresses nothing and cannot go stale with a
  release. (An earlier revision of this branch removed the suppression on the strength
  of a local mypy 2.3.1 and reddened CI, which runs 2.4.0; the default branch was green
  throughout.)

### Four output columns now say what their headers say

- **`posterior`** in `recombination_regions.tsv` is an HMM measurement, but the field's
  default was 1.0, so every 3SEQ, MaxChi, Bootscan, GENECONV and donor-absent row
  reported a posterior of 1.0 it never had, indistinguishable from a measured 1.0. When
  a merged row's most significant member was a site caller the HMM member's real value
  was discarded for that default. The column is now `NA` on a row no HMM segment
  contributed to, and a merged row carries its HMM member's value, as it already did
  for the breakpoint.
- **`margin` and `minor_cluster_size` after donor re-attribution.** Re-labelling a
  region's donor to a better clade replaced `minor_parent` and `mean_sim_minor` and left
  `margin` and the cluster size from the previous donor, so a row could read
  `mean_sim_minor 1.0, mean_sim_major 0.9, margin -0.3`. The new `mean_sim_minor` was
  also the clade consensus score, on a different scale from `mean_sim_major`. The row
  now describes the genome it names: identity over the span, the margin as their
  difference, and the size of that genome's clade.
- **`length_bp` in `coverage_gaps.tsv`** was the alignment-column span, gaps included,
  under a header that says bases, while `Region.length_bp` is query bases. It is now
  `query_end - query_start`. The `min_gap` filter, which compares alignment columns,
  says so in its comment instead of claiming bp; its behaviour is unchanged.
- **`nan` similarities** on a donor-absent row built with no major parent printed as
  `nan` beside the `NA` every other optional column uses. They print as `NA`.

### Bootscan's run-length null no longer depends on the block grid's phase

- **The defect.** The null permutes windows in blocks one window wide so that the
  autocorrelation of overlapping windows survives, but the grid was fixed at window 0. A
  donor run that happened to sit inside one block was carried intact into every
  permutation and scored p = 1.0 whatever its length, while the same run shifted by three
  windows was cut in two and scored about 0.1. Measured on 111 windows with block 10, the
  same 6-window run scored 0.0755 at one offset and 1.0000 at another, and the same
  10-window run likewise. The answer came from where the tract happened to fall.
- **The fix.** The phase is a nuisance parameter, so the null averages over it instead of
  fixing one: each permutation rotates the win vector circularly by a random amount
  before the blocks are cut and shuffled. Every row still keeps exactly its own number of
  wins. The spread across offsets falls from 0.925 to under 0.01, and the all-wins and
  single-block guards are unchanged.
- **It was wrong in both directions, and the cost is real.** A favourable phase made a
  short run look significant as surely as an unfavourable one made it look impossible: a
  6-window run (shorter than one window width) scored 0.0755 at one offset and now scores
  about 0.52 wherever it sits, which is what a run shorter than the correlation length is
  worth. The same averaging makes longer runs somewhat harder to call: an 18-window run
  that scored 0.005 at its lucky phase now scores 0.03. One test case changed verdict
  because of this -- on a panel of eight references 0.1 % apart, Bootscan's 18-window run
  reached q = 0.06 and no longer clears FDR, so it stops naming a donor there. On a panel
  0.25 % apart all four callers still agree on the planted donor at q = 2e-14.
- **Measured against the harness that can see this.** `validation/run_specificity.py
  --reps 3` is unchanged by the fix: 0 of 12 runs report a false region before and after,
  and the positive control is detected 3/3 with the correct donor and a median breakpoint
  error of 13 bp in both. A sensitivity re-measurement on `run_hybrids.py` needs an
  aligner and has not been run.

### The HMM breakpoint interval always contains the switch

- **The defect.** The band was searched outward from the first window of the new state,
  so when the posterior was already confident on both sides it collapsed to that one
  window centre, which is after the breakpoint. On a clean 10 kb A/B fixture with the tract
  at 4000-7000 the HMM reported (4100, 4100) while 3SEQ and MaxChi reported (3999, 4004).
  Its own test checked only `lo <= hi`.
- **The fix.** The band runs from the last window before the switch at which the outgoing
  parent is confident to the first window at or after it at which the incoming parent
  is, so it spans the two window centres the path switched between and is never
  narrower than one window step. A test now asserts the planted breakpoint lies inside
  it. Issue #96 (an interval measured on a different basis from the span it belongs
  to) is a separate symptom and stays open.

### The HMM's major parent no longer depends on reference order

- **The defect.** `plurality_major` counted window wins with `nanargmax`, which hands a
  shared best similarity to whichever label comes first, while the window vote the other
  callers rank their major by counts a tie for every label that shares it. On a panel
  where two references tie in most windows the HMM and the site callers could name
  opposite majors, and swapping the order of the references in the alignment swapped
  the HMM's answer.
- **The fix.** A shared best counts for every label that shares it, the mean-similarity
  tiebreak is unchanged, and a test asserts the HMM major equals the window vote's
  winner under both label orders.

### The report says when the backbone is degraded (issue #82, second half)

- **The defect.** `donor_undercovered` reported on the minor parent alone. A backbone
  that was itself locally degraded was never flagged, so on the mpox run behind the
  issue all three regions read `donor_undercovered = no` while the backbone was the
  degraded sequence, and a reader had nothing in the output saying so.
- **A new column, `major_degraded`,** appended after `donor_undercovered` rather than
  widening it. Which parent is at fault changes what the reader should conclude: a poor
  donor suggests the true donor is missing from the panel, a degraded backbone means the
  comparison behind the region was not between comparably aligned sequences.
- **Measured on the gap rate, not on identity.** The obvious implementation -- flag a low
  `mean_sim_major` -- is wrong, and was tried and rejected here: over a true donor region
  the backbone is the poorer match *by definition*, so it fires on every genuine
  recombinant. It dropped a clean four-caller recombinant in the test suite from high to
  low confidence, which is how it was caught. The flag instead uses issue #81's measured
  discriminator, the backbone's gap rate over the region against its own genome-wide
  rate, so a uniformly gappy assembly is not penalised for being gappy.
- **The bar, and why there.** At least 3.5 % gaps locally and at least 1.5x the
  genome-wide rate. #81 withholds a region outright above 5 % and 2.0x, so this covers
  the band below that: degraded enough to say so, not enough to remove the region.
  Checked against all 72 region rows of the 43 hybrid validation cases, whose regions are
  real recombinants: the flag fires on **none** of them, and the highest local gap rate
  any of them shows is 2.21 %. A donor-absent row is exempt, as it is for #81.
- **Where it shows.** The new column; a `backbone degraded` badge on the region's row in
  the report; and a note under the console table naming the spans. The run-level
  confidence verdict is deliberately left alone -- it already accounts for a poor donor,
  and folding this in was what produced the regression described above.

### A region resting on a locally misaligned backbone is withheld (issue #81)

- **The defect.** Region calling asks whether some reference beats the backbone over a
  stretch of columns. Where the backbone is locally gap-degraded it loses that comparison
  to any reference that stayed intact, and the callers read the switch as a minor parent.
  On the mpox panel behind the issue the backbone carried 31.2 % gaps over the called
  locus against 1.64 % genome-wide, and three method families converged on the same locus
  naming three different donors, which points at the locus rather than at any one donor.
- **The test, and where its numbers come from.** A region is withheld when the backbone's
  gap rate over it is at least `--backbone-gap-ratio` times its own genome-wide rate (2.0
  by default) and at least 5 % in absolute terms. Measured across the five seeds of
  `make_gap_degraded_msa`, the mpox panel and every region of the eight curated datasets
  in `validation/`: a backbone degraded enough to lose the comparison sits at 2.27x to
  19.06x its genome-wide rate with a local rate of 5.02 % to 31.21 %, while the backbone
  under a real recombinant region stays at or below 0.74x and 3.09 %. Both conditions are
  required, so a backbone with no gaps elsewhere is not withheld for a single gap column
  and a draft assembly gappy throughout is not withheld for being uniformly gappy. The
  test is applied per caller, since the defect drew several callers to the same locus.
- **A donor-absent coverage row is exempt.** It states that the true donor is likely
  missing, which is not a comparison the backbone can lose; one curated dataset carries
  such a row at twice its backbone's genome-wide rate.
- **Nothing is dropped silently.** Each withheld locus is logged, the count goes into
  `run_provenance.json`, and the report carries a note naming the columns.
  `--backbone-gap-ratio 0` reports the regions anyway.
- **What this does not fix, issue #86.** On one fixture seed the flanks are reported as
  a mosaic with the query's own sister cast as the donor -- the mosaic inverted. The
  reported backbone is aligned normally over those regions, so no region is withheld
  here. The same two regions are reported on the previous release, so this is not new.
  The strict `xfail` for issue #81's own symptom comes off; a strict `xfail` for #86
  goes in. (Its cause is **not** the ranking statistic this entry first named -- see the
  #86 section below, which corrects it with the measurement.)

### Reference coverage on a panel that is nearly identical to the query (issue #82)

- **The best-similarity threshold could not resolve a low-divergence panel.** It is the
  query's typical best match less `--coverage-rel-drop`, so a panel whose members all
  sit within a percent of the query put it below every similarity the panel shows: no
  window could reach it and no stretch was ever reported. The mpox run behind issues #81
  and #82 recorded a threshold of 0.9500 against a minimum best similarity of 0.9950 and
  an empty `coverage_gaps.tsv`, although two loci there are locally misaligned. When the
  fixed drop is vacuous in this way the threshold now falls back to the panel's own
  spread of best similarities, and never finer than two mismatches over the comparable
  bases a window carries -- the smallest identity difference a window can resolve. The
  fallback is bounded by `--coverage-rel-drop`, so it only ever tightens a threshold that
  could not resolve anything; a panel the fixed standard already resolves is untouched.
  Being unreachable is not sufficient: on any well-covered panel the worst window may sit
  above the standard, which is a legitimate negative rather than a vacuous threshold, so
  the fallback also requires the closest reference to be within one percent of a perfect
  match across the genome. Measured on five cached panels, the three where the fixed drop
  carries no information sit 0.0000 to 0.0010 from a perfect match, and the two where it
  is merely never reached sit 0.0408 and 0.0670; without this condition the rubella and
  enterovirus pools each gained a stretch that main does not report.
- **A new `kind`, `local_dip`, for what the fallback finds.** Such a stretch sits above
  the fixed standard: the closest reference matches it less well than it does elsewhere
  but still closely, which is a local identity dip and usually a locally poor alignment,
  not an absent reference. At 99.7 % identity nothing is missing, so a `local_dip` is
  reported in `coverage_gaps.tsv` and given its own note in the report and the console,
  is never counted in the missing-reference headline, is never bridged into a
  donor-absent region, and is labelled "local dip" rather than "low coverage" on the
  plot.
- **Recruitment leaves a local dip alone.** `find-references` and the `fill-references`
  loop treated every gap as a genome to go and find, so the tighter threshold would have
  sent them searching NCBI for a reference that is already in the panel. Both now skip
  this kind, and the fill loop converges on a near-identical panel as before.
- **`donor_undercovered` is unchanged.** It states that a region's donor may be a
  stand-in for an absent one, which is an absolute question, so it stays judged against
  the fixed standard rather than the fallback. No region call, coordinate, p-value or
  flag changes; a donor 99.7 % identical to the query is not reported as a poor match.
- **Not addressed here:** the second half of issue #82, a coverage flag that describes
  the major parent rather than the donor, and the false region of issue #81 itself. The
  strict `xfail` for the first half comes off; the one for #81 stays red.

### The mosaic is no longer reported the wrong way round (issue #86)

- **The defect.** `reconcile_major` already guards one artefact: a caller segments on a
  view that can disagree with the genome-wide evidence, and a major parent that is the
  closest reference in not one base-pair window is an artefact of that view. The guard
  corrected the run's backbone and stopped there, leaving the regions that caller had
  already built on the rejected major in the table. The report then named one backbone in
  its headline and the opposite in its rows, casting the genome-wide closest reference as
  a donor into a reference that is closest nowhere.
- **The cause, which is not the one #86 was filed on.** The issue recorded the backbone
  ranking as the culprit, with the median as the remedy. Measurement does not support it:
  `plurality_major` ranks by window wins and uses the mean only to break ties, and on the
  base-pair view it picks the correct reference on all five seeds of
  `make_gap_degraded_msa`. Mean-alone would pick the wrong one on all five, but it never
  reaches the decision, so switching to the median would change nothing. What
  distinguishes the one failing seed is that it alone crosses the informative-site
  auto-trigger; on that view the backbone inverts, and its major won zero base-pair
  windows against the true backbone's 129.
- **The fix.** A region is withheld when the reconciled backbone is its *minor* parent
  **and** its major parent won no base-pair window -- the same condition the existing
  guard used to reject that major. Both are required, so a caller that disagrees about
  the backbone with window evidence of its own is left alone rather than overruled by the
  headline, and an empty vote withholds nothing. The breakdown rows are filtered in step
  with the regions, so no row is left describing a different region than before.
- **Nothing is dropped silently.** Each withheld row is logged with the backbone that was
  rejected, and the count and reason go into `run_provenance.json`.
- The strict `xfail` on `test_a_clonal_query_yields_no_region_at_all` comes off.

### Python 3.13 is the floor, and CI tests one interpreter

- **`requires-python` is now `>=3.13`** (it was `>=3.11`), with the older classifiers
  removed and `environment.yml`, the README install line and `CONTRIBUTING.md` raised to
  match. The container already pinned 3.13. The practical effect is narrower than the
  version jump suggests: the code uses 3.12 syntax, so 3.11 and older fail at import
  with a `SyntaxError`, while 3.12 still imports and is simply no longer supported or
  tested. Declaring the floor is what keeps that honest.
- **PEP 695 type parameters.** `ruff`'s `target-version` moves to `py313`, which keeps
  the PEP 695 rewrites: `core/ncbi.py`, `core/plugins.py` and `core/executors.py`
  declare their type parameters inline (`def f[T](...)`, `class Registry[T]`) rather
  than through module-level `TypeVar`s. Raising the target from `py312` to `py313`
  surfaced no further rewrites.
- **CI is four jobs and no matrix.** With one supported interpreter the version matrix
  has nothing left to vary, so `test` runs lint, types and both test steps in a single
  job on Linux; `test-macos`, `container` and `build` are unchanged in what they check.
  Lint and types had been split into a job of their own only to stop the matrix paying
  for the same answer in every entry, and lint is the first step of the job that
  replaces it, so a style error still surfaces about as early. `mypy` no longer needs
  its note about guarding the oldest supported version by a separate job.

### Validation: an offline false-positive control on real genomes

- **`validation/run_pool_specificity.py`.** `run_specificity.py` has exact ground truth
  but simulated sequences; this asks the same question of real ones. The query is a tip
  from a cached Nextclade reference tree whose label is a clade in its own right, so it
  is a single sampled genome and a region reported on it is a false positive. It needs no
  network, and that is enforced rather than assumed: both proxies are pointed at a closed
  local port, so an attempt to reach NCBI is refused in milliseconds instead of stalling.
  Two findings made that necessary. `--nextclade-dataset` fetches the dataset index over
  the network *before* consulting the cache, so the pool is passed as `--candidate-pool`
  instead; and `--max-rounds 1` does not prevent the fill, because the search runs inside
  round one whenever a coverage gap remains. A pool whose panel cannot cover its own query
  is therefore reported as **not testable** -- an untested absence, neither a pass nor a
  failure -- rather than being counted or quietly dropped.
- **`pick_pool_query.py` exposes `choose_query`**, so the chooser is shared with the new
  harness instead of duplicated; its command-line behaviour is unchanged.

### Test data: the gap-degraded major parent

- **A fixture for the false positive behind issues #81 and #82.**
  `tests.factories.make_gap_degraded_msa` builds a clonal query on a near-identical,
  clade-structured panel whose nearest reference is misaligned over one locus
  (frame-shifted and gap-riddled at the rates measured on the mpox locus: ~9% gaps there
  against ~2% genome-wide). Nothing is copied into the query, so every region a caller
  reports on it is a false positive by construction. The suite covers what the data is
  (clonal, locally degraded, the dip confined to the locus) and records the two open
  defects as strict `xfail`: a region is called over the locus on all five seeds tried
  (#81), and neither the coverage gaps nor `donor_undercovered` flag it (#82). Fixing
  either issue turns its test red until the marker comes off.

### Caller behaviour, plan C Task 6: pool clade labels

Completes the one task of plan C
(`docs/superpowers/plans/2026-10-01-audit-c-caller-behaviour.md`) that was held back when
its end-to-end gate was not met, rebased onto the plan D-G work. The investigation behind
the decision is `docs/superpowers/specs/2026-10-19-audit-c-task6-investigation.md`.

- **Pool clade labels are read as written instead of mined.** Tessera writes pool genomes
  with the defline `accession label`, then passed that defline to a miner built for
  free-text GenBank titles, which wants a digit-bearing token of four or more characters.
  Short and multi-word clades did not survive: on the cached pools every measles tip
  (846), all but two mpox tips (1556) and 413 of 1052 HIV-1 tips came back untyped, and
  128 recombinant-named SARS-CoV-2 tips were never excluded because they were never typed.
  Those genomes fell back to clade-blind ANI dereplication. The label is now read verbatim
  for pools Tessera built and carried to the panel's lineage table through a
  `pool_labels.tsv` sidecar; a user's own `--candidate-pool` is still mined, since its
  headers are free text. This changes panel composition: selection changed for 15 of the
  33 cached pools.
- **A label that is not a clade no longer seeds a lineage.** Selection takes one
  representative per lineage, so every distinct label won a panel slot. Two kinds of label
  do not deserve one. `unassigned` is what Nextclade writes for a tip it places in no
  clade (mumps writes `unclassified`); it appears in nine of the 33 cached pools, 126 tips
  in the H3N2 HA dataset, and on the mpox dataset it is the reconstructed ancestral
  sequence, which reading it as a clade admitted to the panel as a donor. A label such as
  `Ib/IIb` spans two clades the same pool also carries standalone; the one mpox tip
  carrying it became the minor parent of a 2721 bp region called on a clonal clade-IIb
  query. A composite is only recognised when every slash-separated part is itself a clade
  of the same pool, because a slash also occurs inside ordinary clade names: enterovirus
  D68 names a clade `A2/D` across 364 cached genomes while naming neither `A2` nor `D`,
  and Ebola names `Ebov-2013/r2021`. Over the 33 cached pools the two rules drop 221 of
  34027 typed genomes and no legitimate clade.

The gate that originally failed rested on two premises that did not hold, and its
criterion has been corrected. Its mpox baseline never scanned, because selection left it
one reference and detection needs two, so "no regions before" was an untested absence
rather than a negative. Its HIV-1 query was AY900575, typed CRF24_BG, a designated
circulating recombinant form, so a region reported on it was not evidence of a false
positive. The gate now chooses a query whose label is a clade in its own right and
reports a dataset whose baseline could not scan as not testable rather than as a pass.

Two defects the investigation uncovered are independent of this task and are filed rather
than fixed here: a false positive where the major parent is locally gap-degraded on a
low-divergence panel (#81), and the coverage guard that cannot see it, because its 0.95
best-similarity threshold has no resolution on such a panel and its undercover flag
describes the minor parent only (#82).

### Fourth audit, plan F: reporting, CLI and reassortment

From `docs/superpowers/specs/2026-10-12-audit-swot-gap.md` (B1, B3, B5, B7-B9, the smaller
reporting items, E1-E4, E5 in part, E14, F1-F5; help text and documentation for B2, B4,
B6, A4). None of these changes which regions are called or their coordinates. The
specificity harness at the default gate (`run_specificity.py --reps 3`, `--min-methods 2`:
0/12 runs with a false region) gives output byte-identical to the run before the branch,
apart from its two header lines. With `--min-methods 1` it gave 7/12 runs with a false
region, 8 regions, all HMM, and positive control 3/3 detected, 3/3 correct donor, 13 bp
median breakpoint error; no byte-baseline was captured for that setting.

- A test or caller that gave no result is reported as `not run`, with the reason, in the
  report, the TSVs and `run_provenance.json`: the PHI test under `--no-phi` or with too few
  informative sites, and an HMM that lineage clustering left with one state. A PHI test
  that is not testable at the window used keeps `not testable` in `run_provenance.json`.
  `parent_free_support` reads `not run` instead of `no` in all three PHI cases.
- A `low_information` coverage gap is no longer announced as a possible missing
  reference.
- `run_provenance.json` gains `coordinates`; under informative-site windowing
  `window / step` says it describes the base-pair scan.
- The report headline names the run's reconciled major parent; the pair plot uses the
  report's colours; `recombination_methods.tsv` has one row per row of
  `recombination_regions.tsv`, donor-absent regions included.
- Every command rejects an output path that is a file. The discover commands range-check
  the coverage fractions and search counts (`--top-gaps 0` is an error). `msa` refuses
  `--reference` with `--query-as-backbone` and warns before overwriting. `--threads` is
  range-checked on `detect`, `build-panel`, `fill-references`, `msa` and `reassort`.
  `find-references -q` is kept verbatim. `build-panel` takes `--pool-consensus`,
  `--keep-recombinant-lineages` and `--deep-typing`.
- `reassort --margin` must be above 0 and at most 5 ANI points, and an `--ani-floor`
  above 0 and at most 1 is rejected as a likely fraction. A discordant pair inside one
  parent group gives `undetermined` instead of `reassortant`, unless another discordant
  pair crosses parent groups, which still gives `reassortant`; every pair's relation is
  written to `segment_pairs.tsv`. `reassort` takes `--threads` for `--scan-segments`.
- The help of `--min-region`, `--margin` and `--merge-gap` names the callers that read
  them. `validation/run_specificity.py`, `run_regimes.py` and `characterise_false_regions.py`
  parse their arguments with `argparse`; the other harnesses, including `run_hybrids.py`,
  still parse them by hand, so `--help` there is not handled.

Fixes from the post-1.2.0 audit
(`docs/superpowers/specs/2026-10-01-post-1.2.0-audit-design.md`, items A1-A10, B1-B11 and
C1-C5, C7, C8). The data-safety and reporting fixes (A, B) leave the callers alone: for a
given alignment they find the same regions, with the same coordinates and p-values, as
before. The caller-behaviour changes (C) do not: HMM region coordinates, the HMM's part on
near-identical panels, 3SEQ/MaxChi tract starts, very small 3SEQ p-values and the PHI
p-value can differ, as listed under "Callers" and "Changed" below. Three further things
around the callers change: the `donor_undercovered` flag and
donor-absent rows where a coverage gap is now recognised as a breakpoint artefact; runs
that select the barcode caller on a panel where it cannot run (it no longer counts toward
`--min-methods`, and a barcode-only run is refused); and what reaches the scan -- an
alignment built with the MAFFT or MAF-based backends, an alignment built with
progressiveMauve (the backbone row is now the reference itself), a sibeliaz alignment on
macOS (which previously came out empty), or a panel built with `--curate`, may now be a
different (corrected) one.

### Fixed

Alignment integrity, from the fourth audit
(`docs/superpowers/specs/2026-10-12-audit-swot-gap.md`, section 1.3). These change what
alignment reaches the scan, not how the scan calls regions:

- **The progressiveMauve backbone row could lose columns.** It was taken from the first
  query's projection, which is gapped wherever that query did not align; with queries
  covering different parts of the backbone, the rest of the backbone came out as gaps,
  and which part depended on the order of the queries. The backbone row is now the
  backbone's own sequence. Checked against a fake of the tool only: progressiveMauve is
  not installed in CI or on the development machine.
- **The `sibeliaz` backend wrote no alignment on macOS.** Its wrapper uses
  `paste --delimiters`, which BSD `paste` reads as a file name, and an `xargs -I` command
  longer than BSD's 255-byte default; either left the MAF without rows, and the adapter
  reported this as spoa running out of memory. The patched wrapper now rewrites both.
  Reproduced and checked with SibeliaZ 1.2.7 on macOS; Linux is unaffected.
- **Two records with one sequence ID in a genome were silently merged.** The `sibeliaz`,
  `cactus` and `minimap2` layouts are keyed on the ID, so the second record's blocks were
  written onto the first one's columns (or the backbone row lost a record) in an
  otherwise well-formed alignment. Such a file is now refused at staging, for every
  backend, naming the file, the ID and both lines.
- **`minimap2` multiplied its thread budget.** `--threads N` started N concurrent
  `minimap2` processes with N threads each. The budget is now shared between them, by the
  same rule the `mafft` backend uses.
- **A genome the `minimap2`, `mafft` or `progressivemauve` backend could not place was
  dropped in silence.** It became an all-gap row with no warning (only the MAF-based
  backends named it), and the scan later reported it as an alphabet problem. These
  backends now name such genomes, and the scan reports an all-gap row as unaligned,
  separately from a record whose letters it could not read.
- **Two MAF inputs were absorbed instead of reported.** A MAF row cut off before its
  sequence text was dropped, so its genome became an all-gap row that read as divergence;
  and a backbone contig whose length in the MAF differed from the staged FASTA had its
  blocks written onto the next contig's columns. Both now stop the run with a message
  naming the file.
- **The Nextclade pool cache could not be installed over an existing directory.** Two runs
  building the same pool at once, or a cache directory holding only unrecognised files,
  made the install fail with "Unexpected error". The pool is now installed the way the
  NCBI Virus panel is: a pool another run finished first is kept, a directory without a
  completeness record is replaced, and the pool carries such a record. Pools cached by
  earlier releases are still used.
- **A file or symlink at `<output>/collection` gave a raw error.** When a previous run's
  marker was present, clearing the working collection failed with `NotADirectoryError`,
  an `OSError` for a symlink, or `FileExistsError` for a dangling one. Each is now a user
  error that names the path; a symlink is refused rather than followed.
- **External tools inherited standard input.** A tool that read stdin or prompted would
  have waited on the terminal, and alignment has no default timeout. Every external tool,
  and every version probe, now runs with standard input closed. No current tool is known
  to have done this. Reconfiguring logging also now closes the previous run log.
- **`cactus`: two labels differing only in `.` and `_` would share one row.** Cactus
  sample names cannot contain `.`, and the adapter maps it to `_`; such inputs are now
  refused before Cactus starts. Checked without the binary, which is not installed here.

Data safety and alignment:

- **The MAFFT backend ignored strand.** A genome, or one contig of a draft assembly, on the
  opposite strand to the backbone was aligned as given and came out at chance-level identity
  (about 0.40 against 0.97 for the same genome in forward orientation), which the scan then
  read as a divergent region. MAFFT now runs with `--adjustdirection`.
- **The MAFFT backend could exhaust memory on long genomes.** `--addfragments` ran a
  pairwise local alignment whose memory grows with the product of the two lengths, and
  `--threads N` started N such jobs with N threads each: a 200 kb mpox pair passed 100 GB.
  When the backbone or a genome to add is longer than 20 kb the alignment now uses
  `--6merpair` (0.09 GB for the same pair), `--threads` is shared between concurrent jobs,
  and below 20 kb the number of concurrent jobs is limited by physical memory. Alignments
  above 20 kb may differ from before: the one pair compared in that range (SARS-CoV-2)
  was identical, and nothing could be compared at mpox size. Below the threshold the mode
  is unchanged. `--aligner-arg sixmerpair=false` restores the earlier mode.
- **Pool seeding ignored the exclusion list.** With `--seed-source nextclade`, `ncbi-virus`
  or `local`, the query's own record and any `--exclude` accession were still offered to
  panel selection. A query that is itself in the pool could be seeded as its own
  reference, which stopped the run ("Two inputs share the label") or gave the scan a
  reference identical to the query. Pool genomes whose file name is the query's label,
  the record in its header, or an `--exclude` accession are now left out before
  selection. Matching is by file name, so a copy of the query stored under an unrelated
  name is still not recognised.
- **`reassort` never applied its alignment-fraction filter.** The threshold was written as a
  fraction (0.5) and compared with skani's percentage, so a tip aligning over a fifth of a
  segment could outrank a full-length match and become the segment's nearest strain. The
  threshold is now 50 %.
- **`fill-references` did not align the last round's downloads.** On a `--max-rounds` exit the
  final downloads were counted in `fill_summary.tsv`, the report and `lineages.tsv` but were
  absent from `panel.msa.fasta`, so detection ran without them. One more alignment
  (`final.msa.fasta`) is now built from them, without a further search.
- **`fill-references --curate` did nothing unless a round both found a gap and downloaded a
  reference.** A supplied collection that contains a whole-genome sibling of the query has
  no coverage gap, so the loop converged before curation ran. Curation now runs before the
  first alignment for a supplied collection and before each later build for that round's
  downloads. A panel seeded from scratch is unchanged: seeding has its own sibling filter.
  One limit remains: the sibling test is anchored on the query's closest whole-genome
  relative, so when a sibling is itself the closest genome in the collection it becomes
  that anchor and is kept (further siblings are removed). The log names the anchor chosen.
- **`fill-references --curate --reference X` could delete X** and fail the next round with
  "Reference 'X' not found". The given reference is now never removed by curation (it is
  listed as `sibling-kept` / `redundant-kept` where it would have been). The sibling test
  stays anchored on the query's closest relative, not on the reference.
- **A run could delete an existing `<output>/collection` directory that was not its own.**
  `detect`, `fill-references`, `build-panel` and `curate-panel` clear that directory to hold
  their working copy of the references. If it already exists, is not empty, and no earlier
  run left its files beside it, the run now stops with an error and nothing is deleted.
- **`find-references --download <collection> --curate` deleted genomes that were already in
  the collection.** Curation now removes only what the run downloaded; genomes already
  present are reported as `sibling-kept` / `redundant-kept` in `panel_lineages.tsv`. A
  downloaded sibling can no longer be picked as the curation backbone either.
- **`reassort --scan-segments` could write outside the output directory.** A segment record
  named `..` resolved the scan directory to the parent of the output and replaced its
  `collection/`. Segment names are sanitised by one shared helper, and two names that
  sanitise alike no longer share a scan directory.
- **MAF-based backends (`sibeliaz`, `cactus`) could drop rows and columns.** A genome placed in
  no alignment block had no row at all; it is now an all-gap row. It is named in a warning,
  as is a genome aligned only in blocks that do not include the backbone. With
  `sibeliaz`, a multi-contig backbone was laid out in sorted-name order and lost any contig
  that no block covered, shifting later coordinates; it is now laid out as its FASTA file is,
  with uncovered contigs kept as gap columns.
- **`progressivemauve` named rows after the resolved input path.** Symlinked collections,
  unrecognised extensions and whitespace in a path gave wrong or duplicate row names, and a
  genome whose link target shared the backbone's name was dropped. Rows are now named from
  the staged genome label. (Checked against a simulated XMFA; the tool itself was not
  available.)
- **Whitespace in sequence lines was read as sequence.** A reference saved with trailing
  spaces, tabs or Windows line endings produced a ragged alignment. Such a genome is now
  aligned from a cleaned copy (a warning names it), and the FASTA reader ignores whitespace
  in sequence lines.
- A header line `> ` (no name), a truncated MAF block, an XMFA that does not list the
  reference and an empty MAFFT result are reported as input/output errors that name the
  file, instead of "Unexpected error".

Reporting:

- **The run provenance named the wrong callers.** MaxChi, Bootscan, GENECONV and the barcode
  caller were each described as `heuristic (min ... / margin ... / merge ...)` in
  `run_provenance.json` and the report, so a default run recorded
  `hmm + 3seq + heuristic + heuristic`. Each caller is now described under its own name, and
  the record gains the settings that change what is reported: the agreement gate
  (`--min-methods`), sibling exclusion, lineage clustering and donor re-attribution.
- **A clean recombinant between divergent parents was reported as a possible missing
  reference.** A window straddling a breakpoint matches neither parent well on its own, so
  its best similarity fell below the coverage threshold, the stretch was called a `divergent`
  coverage gap, and the region was marked `donor_undercovered` -- on the shipped
  `divergent` example the headline read "low confidence" for a donor identical to the query.
  A gap within one window of a called region boundary is now labelled `breakpoint` in
  `coverage_gaps.tsv` and the report when every under-threshold window in it is matched,
  to the coverage threshold, by a single switch between the region's two parents. Such a
  gap does not caveat the region, is not turned into a donor-absent region, and is left out
  of the headline. A stretch from a source outside the panel stays `divergent` even when it
  is short and sits beside a breakpoint. **`donor_undercovered` and the confidence wording
  change for regions whose only gaps were breakpoint artefacts.** One consequence is not yet
  measured on the hybrid harness: a region whose donor is a close stand-in for an absent
  lineage (above the coverage threshold) was previously caveated only through its
  breakpoint gaps, and is now reported without a caveat.
  Reference recruitment (`fill-references`, `find-references`) is unaffected.
- The report's "covering N kb (P %) of the query" added the lengths of overlapping regions
  that name different donors; it now reports their union.
- After `--reattribute-donors`, `recombination_methods.tsv` and the report's method table
  kept the donor from before re-attribution.
- `--lineage-map` pointing at a file that does not exist was ignored (exit 0, untyped
  report) by `recomb`, `type-lineages`, `detect`, `fill-references` and `build-panel`. It
  is now an error.
- **The barcode caller on an untyped panel read as a negative.** It named the first record
  of the alignment as major parent and its column in the method table said `no`. It now
  reports no major parent; in an ensemble it is logged and shown as `not run`, the
  agreement gate counts only callers that ran (a gate lowered for that reason is logged and
  recorded as, for example, `1 (requested 2; 1 caller ran)`), and a run that selected only
  `barcode` is refused. The report and provenance give the actual reason, which on a typed
  panel is that fewer than two clades carry enough markers.
- The report judged the PHI p-value at alpha 0.05 whatever `--alpha` was.
- **A PHI test that could not reject was reported as "no signal".** With no more
  informative sites than the window can hold (`--phi-window`, default 100) the permutation
  p-value is 1 for any data. This is now reported as `not testable` (`NA` in the
  `recombination_profile.tsv` header, `phi_p = None` in the API). With the cap on the
  PHI window (see "Changed") a scan no longer reaches this state; the report wording and
  the `None` value remain for it.
- Plots labelled a donor-absent region "recombinant: <backbone>"; `similarity_pair` showed
  the two leading window winners rather than the major parent and the leading donor; with
  `--top-n 1` the donor was drawn grey.
- The report footer listed `.pdf` plots under `--plot-format png` and files that were not
  written; the methods text and references described a two-caller ensemble; `--method`
  help omitted `geneconv`.
- `sibeliaz` was probed with `-v`, which it rejects, and the error line was recorded as the
  aligner version. A failed probe is no longer recorded as a version.
- Input checks: `find-references --msa <missing>` and `recomb -o <existing file>` failed
  with "Unexpected error"; `--max-rounds 0` exited 0 having built nothing; `reassort`
  accepted `--ani-floor 500`, `--margin -3` and a `--dataset` key matching no segment;
  `type-lineages` rejected `.fasta.gz` / `.fas` collections (it now reads what the other
  commands read, and like them rejects a file in the collection that is not FASTA, where it
  used to skip it); a multi-line aligner error broke `segment_scan.tsv`.

Callers:

- A 3SEQ or MaxChi tract started at the first maximum of the discriminating-site walk
  before its trough, not the last, so it could begin a few sites early and include a
  stretch that is not donor tract. MaxChi's permutation null follows the same rule.
- Exact 3SEQ p-values below about 1e-16 were written as `0.0`. They now keep their
  magnitude, and the donor choice among such p-values no longer falls to a tie.
- The ensemble merge could label a region with the same genome as donor and major when
  the callers disagreed on the backbone, and could report one donor lineage as two
  overlapping regions on a typed panel.

Panel recruitment, from the fourth audit
(`docs/superpowers/specs/2026-10-12-audit-swot-gap.md`, section 1.4). These change which
references are recruited, not how the scan calls regions:

- `fill-references` stopped after round 2 whenever there were more coverage gaps than
  `--top-gaps`: the stop test watched the worst gap, while the search reached only the
  longest ones, so an unsearched short gap held the measure still.
- The fill loop re-downloaded, every round, the genomes curation had just removed, because
  they were never added to the exclusion list. They now are.
- Windowed BLAST seeding reported that the query's lineage saturates NCBI when every hit
  had been filtered out (the query's own record, or `--exclude` accessions), and then ran
  the negative-lineage search and the NCBI Virus download on that false explanation. It
  now needs a kept hit, as whole-query seeding already did.
- A taxon with no RefSeq genome stopped the run with "NCBI Virus returned no genomes"
  instead of broadening to complete genomes, which it did for a taxon with one or two
  RefSeq genomes. It now broadens; the error remains for a taxon with no genome at all.
- Without `--deep-typing`, the panel was typed by mining headers again, ignoring the exact
  NCBI Virus lineages recorded during recruitment in `ncbi_lineages.tsv`; names without a
  digit-bearing token of four or more characters (`B.1`, `A`, `XBB`) were lost and those
  genomes skipped lineage-aware handling. Both typing paths now read the sidecar.
- When `--curate` was requested (always, for `detect` and `build-panel`) but never ran --
  a seeded panel with no coverage gap downloads nothing for it to curate -- the run wrote
  no panel table and no log line. The log and the report's "Reference panel" section now
  state that curation did not run.
- Pool selection kept every genome without comment when all of them were siblings of the
  query, and returned a normal-looking panel. It now logs a warning, as BLAST seeding
  already did for the same situation.
- Regional pool selection treated a query between one and two selection windows long
  (1.5-3 kb at the default) as a single window, which is whole-genome selection. Such a
  query is now cut into two halves, with the whole-genome ranking as a fallback if neither
  half matches. Longer queries are unaffected.
- `--fetch-limit` was part of the NCBI Virus cache key although it does not change what
  is fetched (it only sets when a notice is logged), so changing it re-downloaded an
  identical genome set. It no longer is. A cache written by an earlier version is keyed
  differently and is fetched once more on first use.
- Nearest-reference typing took a Nextclade tip labelled `?` for a clade named `?`; it
  now uses the shared set of "no clade" markers. The per-clade consensus step in the
  Nextclade pool builder skipped only `NA`, so tips labelled `?` or `example` each got a
  consensus genome for a clade that does not exist; it uses the same set.
- The candidate-pool reader accepted four file extensions, so a `.fas` or `.fsa` genome
  was left out of the pool without notice. It now reads a pool as every other command
  reads a collection, and rejects a file that is not FASTA.

Callers (audit 4). Caller-statistics fixes from the fourth audit
(`docs/superpowers/specs/2026-10-12-audit-swot-gap.md`, items A1-A14, B2, B4). The
lineage-cluster changes (A1, A2), the descent floor (A3, A6), the HMM segment padding
(B4, A10), the last-window fix and the Bootscan null can change which regions are reported
on panels with near-duplicate references, at the genome ends, or wherever a padded segment
moves across the sign test's significance threshold; the others change reported values
only. Measured
effects on the offline harnesses are listed after the bullets.

- **Near-duplicate references could silence the HMM.** A lineage cluster's per-window
  evidence was summed over its members, so three copies of the backbone outweighed a
  single-genome donor and the HMM lost a 1 kb tract (HMM among the callers in 3 of 10
  simulated runs). The evidence is now the member mean (6 of 10). (A1)
- **A redundant donor lineage was split into two regions and lost Bootscan.** Clusters are
  now computed once and shared: the merge treats two members of one cluster as one donor,
  and Bootscan counts support per cluster. One row in 10 of 10 runs against 7 of 10, and
  Bootscan among the callers in 10 of 10 against 1 of 10. (A2)
- **The HMM sign test returned p = 0 above 1000 discordant sites.** It is now exact at
  every size, in log space, and never below the smallest normal double. (A9)
- **3SEQ discarded short tracts before testing them, and both triplet callers pruned the
  Benjamini-Hochberg family.** 3SEQ now tests every donor with at least 20 discriminating
  sites, using an exact dynamic program; MaxChi counts shallow donors in its family with
  p = 1 rather than dropping them. (A3, A6)
- **A merged region reported one member's identity over the union span**, and its
  `n_windows` was a maximum over different units. Both are restated over the merged span;
  `n_windows` counts base-pair windows for every caller. (A11, A8)
- **HMM `n_windows` and `posterior_support` described the untrimmed segment**; they now
  describe the reported span. Segment ends are also padded by the median window step
  rather than by the gap between the first two windows, which is not typical under
  informative-site windowing. That padding moves the segment, and so can change whether
  the HMM's sign test over it reaches significance -- it changes which HMM regions are
  called, not only the values reported for them. Whether a segment widens or narrows
  depends on which of the two gaps is larger, so neither direction is systematic. The
  padded width is also what `--min-region` is compared against, so a segment near that
  floor can cross it either way; on the curated datasets the shift was 1460 columns
  against a default floor of 500, and no segment crossed. (B4, A10)
- **Bootscan's permutation null dropped wins** when the window count was not a multiple of
  the block, giving anti-conservative p-values. (A5)
- **The last `(width - window) % step` columns were in no window**, so a tract or a
  coverage gap there was invisible to Bootscan, the barcode caller and the coverage
  diagnostic. A final window flush with the end of the alignment is added.
- **A tract at the first discriminating site had a zero-width breakpoint interval**
  (3SEQ, MaxChi, GENECONV). (A14)
- **Memory:** the PHI gamete matrices, the Hudson-Kaplan pair list and GENECONV's
  permutations no longer exceed the scaling caps (136 MB to 25 MB, 381 MB to 0.25 MB,
  680 MB to 250 MB on the test inputs); MaxChi's null no longer overflows int64 at
  bacterial scale. (A13, A7)

What this measurably changed, on the simulated and curated harnesses (10 replicates):

- False regions on the clonal specificity panels fell at both agreement settings once the
  shallow donors entered the Benjamini-Hochberg family: 17 of 40 runs with a false region
  (20 regions) to 16 of 40 (19 regions) at `--min-methods 1`, and 2 of 40 to 1 of 40 at
  `--min-methods 2`. No new 3SEQ or MaxChi false region appeared.
- On the `near_identical` clonal control the median-step segment padding added one HMM
  false region at `--min-methods 1`: the same 1 run in 10 carries a false region before
  and after, but that run now reports two rather than one. At `--min-methods 2` there is
  none either way.
- `n_windows` moves in 13 of the 21 curated regions and in both directions, and
  `mean_sim_minor` in 4 of them -- two in the third decimal (0.8982 to 0.8996, 0.8347 to
  0.8363) and two in the fourth -- both by design.
- One curated region (`sarscov2_xbb`) loses its HMM vote and 406 bp of extent. The cause
  is the segment padding above, not the restated statistics: reverting only the
  median-step pad restores the HMM to the row's callers (`3seq,maxchi` becomes
  `hmm,3seq,maxchi`) and the extent with it. The lost vote is a correction rather than a
  regression. Over the segment's own window union (22000-29734) the sign test gives 9
  against 7 discordant sites, p = 0.40, under either pad; the old pad reached p = 0.03
  only by clipping five major-favouring sites out of that union, so the call it granted
  was not supported by the segment's own evidence. The region is still called by three
  callers and the dataset still passes.
- The padding affects informative-site windowing only. Under base-pair windowing the
  first inter-window gap, the median and the step are the same number, so all three of
  its paths are no-ops: measured on this panel, 291 windows with a first gap of 100 and a
  median of 100, against 12 windows with a first gap of 1577 and a median of 847 under
  informative-site windowing, which is the 1460-column shift above. This is why the
  simulated panels are byte-identical and only curated informative-site datasets moved.
- No curated verdict changes, and `hcv_clonal_1b` fails before and after (a pre-existing
  12 bp MaxChi region on a clonal control).

### Changed

- **`--min-region`, when given, now bounds every reported region**, whichever caller found
  it. It was read only by the HMM, and on the untrimmed segment, so a 175 bp 3SEQ region
  and a 49 bp MaxChi region were reported under `--min-region 5000`. Left unset the
  behaviour is unchanged: the HMM keeps its segment floor at the window size, the
  heuristic caller keeps its own, and the site callers are not bounded, because finding
  tracts shorter than a window is what they are for. The help text said otherwise.
  (B2, A4)
- **An HMM left a single state by lineage clustering no longer counts toward
  `--min-methods`**, as an unrunnable barcode caller does not. When clustering pools every
  reference into one lineage the HMM has no second state and cannot call a region whatever
  the data, so a gate it could not help meet only suppressed what the other callers found:
  `--method hmm,3seq --min-methods 2` on five references 0.25 % apart with a real 6 kb
  tract reported nothing and read as a negative. The clamp, its warning and the provenance
  line now treat it as the barcode caller is treated. With the default four callers and
  `--min-methods 2` nothing changes. The record also says "could test" rather than "ran",
  since a one-state HMM runs and cannot test. (B7)
- **`--merge-gap` now also sets how far apart two coverage gaps may be and still be
  merged**, and so shapes the donor-absent regions built from them. It was read only by
  the heuristic caller, and the coverage diagnostic always merged within one window size
  whatever the option said. The default is unchanged, so only a run that sets it is
  affected. (B6)
- **Some seed options that cannot act are refused.** `detect` and `build-panel` used
  `--candidate-pool` and ignored `--nextclade`/`--nextclade-dataset` when both were given;
  `--pool-consensus` did nothing without a Nextclade pool; and `fill-references` ignored
  `--candidate-pool` without `--seed-source local`. Each of these three combinations is now
  a usage error. Seed options are still inert when `--collection` is given, because
  seeding then does not run; they are not refused in that case.
- **A candidate pool directory must hold only FASTA files.** A samtools `.fai` index, a
  README or a checksum file in the pool directory used to be skipped because the reader
  filtered by extension; it now stops the run with an error naming the file. Move such a
  file out of the pool directory.
- **One new reference is downloaded per coverage gap**, as documented. Every sub-interval
  search contributed all of its hits, so at the defaults one round could add up to 20
  genomes for a 1000 bp gap and 145 for a 6 kb gap, all of which went into the next
  alignment. The best candidate is the one with the highest identity over the searched
  sub-interval. `find-references` writes the full candidate list to
  `candidate_references.tsv`; the fill loop (`fill-references`, `detect` and `build-panel`)
  does not record the candidates it did not take. This applies to
  `find-references --download`, `fill-references`, `detect` and `build-panel`. With
  `--curate`, a gap whose best candidate is a sibling of the query gets its next
  candidate in the following round rather than in the same one; `find-references`, which
  runs once, leaves such a gap without a download.
- **The fill loop stops on the coverage shortfall.** A round must reduce the under-covered
  length, weighted by how far each gap falls below the threshold, by at least
  `--min-improvement` of the previous round's value. **`--min-improvement` is now a
  fraction (0-1, default unchanged at 0.01)**, where it was an absolute gain in the worst
  gap's best similarity. `fill_summary.tsv` gains a last column, `coverage_shortfall`.

- **HMM region coordinates are the span the donor explains, not the HMM segment.** A
  segment ran on through columns where donor and major are identical (a 600 bp difference
  was reported as 11.7 kb when a sibling recombinant was in the panel; now 1.0 kb) and its
  ends sat on the window grid. The reported span now runs from the first to the last site
  in the segment where the query matches the donor and not the major, and reaches the
  alignment end for a terminal tract. **`msa_start`, `msa_end`, `query_start`,
  `query_end`, `length_bp`, `mean_sim_minor` and `mean_sim_major` change for regions the
  HMM called**; which regions are called, and their p-values, do not. On the curated
  validation set six HMM rows moved by 10-776 bp; one end (`hcv_2k1b`, 3127 to 3113) moved
  14 bp further from the published breakpoint.
- **On near-identical panels the HMM can contribute again.** Lineage clustering pooled an
  entire panel less than about 1.5 % divergent into one lineage, leaving the HMM a single
  state. Under informative-site windowing the pairwise agreement is now measured on the
  informative columns. On simulated panels of this kind the HMM joined the call of the
  tract in 10 of 10 runs (0 of 10 before), and on clonal ones the HMM alone reported a
  false region in 1 of 20 runs at `--min-methods 1` (0 of 20 before, and 0 of 20 at
  `--min-methods 2`). On one hybrid-harness case (SARS-CoV-2, 0.4 % divergent) the HMM no
  longer joins the 3SEQ call, which is still reported.
- **The PHI window is capped at a tenth of the informative sites.** `--phi-window` is an
  upper bound; the report states the window used. With the fixed 100-rank window the test
  could not reject on alignments with few informative sites (p = 1 on the shipped
  `divergent` example, now 0.001). **The PHI p-value, and with it the
  `parent_free_support` flag, can change on any alignment with fewer than 1000 informative
  sites**; the per-site PHI profile is smoothed over the same, narrower, window.
  Alignments with 1000 or more informative sites are unaffected. On the Jaya 2023
  benchmark 13 of 185 testable recombining alignments are significant (8 of 73 before:
  more alignments are testable, and the share significant is 7 % against 11 %), with
  specificity 36/36 (16/16 before); 2 of 40 clonal simulations have p < 0.05.
- `RecombinationSignal.phi_p` is `float | None`.
- `seaborn` is no longer a dependency; nothing imported it.

### Added

- `validation/run_regimes.py`: a self-contained harness for near-identical panels,
  sibling tracts, and panels carrying a redundant backbone or donor lineage -- the regimes
  the specificity harness does not sample.
- `validation/characterise_false_regions.py`, and a decision record on the false-positive
  rate at the default agreement gate
  (`docs/superpowers/specs/2026-10-01-min-methods-default-decision.md`).

### Documentation

- 3SEQ and MaxChi are no longer described as independent votes. MaxChi scores the tract
  3SEQ's maximum descent finds, so the two are one triplet signal tested twice: on 300
  simulated clonal walks their p-values have a Spearman correlation of 0.95 (0.88 over
  the walks MaxChi permutes), and at a nominal 5 % both reject in 3.0 % of walks where
  independence would give 0.12 %. Each remains calibrated on its own. The open
  `--min-methods` decision record gains these numbers and the options for the barcode
  caller, which counts toward the agreement gate without a significance test. (A12, A14)
- `docs/detection-methods.md` describes lineage clustering (including its limit on panels
  below about 1.5 % divergence under base-pair windowing), the `breakpoint` coverage kind,
  the reported span of an HMM region and the PHI window cap.
- `validation/README.md` states that the specificity harness defaults to `--min-methods 2`
  while the CLI defaults to 1, gives the measured rate at each, and records that
  `hcv_clonal_1b` currently fails.

## [1.2.0] - 2026-10-01

### Fixed

- **On near-identical panels, values reported as "similarity" were not similarity** (#67).
  Below ~8 % divergence the HMM windows over polymorphic columns only, where a reference
  98-99 % identical to the query typically scores 0.4-0.8. Those values were written to
  `similarity_windows.tsv`, `similarity_stats.tsv`, `window_winners.tsv` and both similarity
  plots, so a panel of near-identical genomes was shown at 50-70 % "similarity". Everything
  named similarity is now identity over all comparable columns, whichever windowing ran. The
  windows the HMM segmented are kept, under their own name: a new
  `informative_site_windows.tsv`, an `informative_sites_top{N}` plot, and a separate, labelled
  track in the report. **In this mode the contents of the three similarity tables change**;
  on panels above the threshold (base-pair windowing) every output is byte-identical.
- **A near-identical donor could be flagged `donor_undercovered`.** In the same mode the HMM
  reported a region's `mean_sim_minor` / `mean_sim_major` as the mean of its informative-site
  windows, and the coverage check compared that against a base-pair threshold -- so a donor
  99.8 % identical to the query over the region was marked as a poor match. Those columns are
  now identity over the region's columns, as the site-based callers already reported them.
  Which regions are called, and their coordinates and p-values, are unchanged.
- A donor-absent row's `mean_sim_major` and `margin` are now read from the same base-pair
  scan that called the coverage gap, rather than mixing it with informative-site values.
- **A collection inside the output directory was deleted.** `curate-panel`, `fill-references`
  and `build-panel` rebuild `<output>/collection` on every run by clearing it and copying the
  input in. Given `--collection <output>/collection` -- the natural way to continue from an
  earlier run -- the clear removed the input before the copy, and the run then failed with the
  references gone. Such a collection (the working copy itself, a parent of it, or a directory
  inside it) is now refused before anything is touched.
- **Every file in a collection directory was staged as a genome**, including hidden ones. A
  `.DS_Store` left by Finder could be picked as the backbone and then crash the aligner.
  Hidden files are now skipped, and any other file that does not start with `>` is rejected
  by name before alignment.
- **`--deep-typing --cache-dir X` built the Nextclade pool in the cache root itself** rather
  than in a per-dataset directory, so a second dataset was answered with the first one's tips
  (references typed against the wrong dataset), or the build failed once the root held
  another cache.
- **A cache directory without a manifest blocked its own replacement.** The fresh fetch could
  not be installed over it and was discarded, so every run re-downloaded the set and then
  reported an empty candidate pool. A manifest-less directory is now replaced; a complete
  cache written by a concurrent run is still kept.
- A failed `efetch` no longer leaves an empty file in the collection, where the next round
  would stage it as a genome.
- The `panel.msa.fasta` published by `detect` / `fill-references` now carries its provenance
  sidecar, so `run_provenance.json` names the aligner and version for those commands too.
- `reassort`: a segment name taken from a FASTA header is sanitised before it is used as a
  file name. A header containing `/` failed with "Unexpected error", and one containing `../`
  wrote outside the temporary directory.
- `sibeliaz`: a sequence ID shared by two genomes is rejected. SibeliaZ names rows by sequence
  ID alone, so the genomes were merged under one label and the other vanished from the MSA.
- A gzipped query is now readable everywhere a query is read (it was accepted by staging but
  failed in the readers that look inside it), and Windows line endings are tolerated.
- A reference already in the collection under a versioned name (`NC_045512.2`) is recognised
  when a search returns the bare accession, instead of being downloaded again.
- Seeding no longer reports "only siblings found" when every BLAST search failed or returned
  nothing; that outcome switched the seed source under a false explanation. It now ends with
  the "could not seed" error, after the per-search failure warnings.
- A missing external binary is reported as such (`MissingBinaryError`) from every tool call,
  and a missing `--query` in `detect`, `build-panel` and `fill-references` as a missing file,
  rather than as `Unexpected error: [Errno 2]`.

### Added

- **Dependency vulnerability scanning and Dependabot.** A `Security` workflow runs `pip-audit`
  over the installed dependency tree on every pull request, on `main`, and weekly on a
  schedule -- most advisories arrive while nobody is touching the repository, so a
  pull-request-only scan would miss them. It is blocking: a security check whose failures can
  be ignored is decoration. `.github/dependabot.yml` proposes monthly updates for both pip
  dependencies and GitHub Actions, excluding major bumps, which stay a deliberate decision.
- **Every GitHub Action is pinned to a commit SHA** rather than a floating major tag, with the
  tag kept as a trailing comment so the pins stay readable and Dependabot can update both
  together. A tag can be moved; a commit cannot, so this is what stops a supply-chain change
  from silently altering what runs in CI.

### Changed

- **The type check is blocking in CI.** The 17-error backlog is cleared (annotation-only
  changes; no defect among them), so `mypy` no longer runs with `continue-on-error`. Making
  it blocking showed that on Python 3.12 and 3.13 it had not been checking anything: the
  configuration pinned `python_version = "3.11"`, under which mypy rejects the numpy 2.5
  stubs and stops. The pin is removed, so each CI job checks under its own interpreter.
- **Dependabot no longer raises the lower bound of a range requirement.** With ranges and no
  lockfile its default is to move the floor to the newest release (`typer>=0.12` became
  `>=0.27.1`), which narrows what an existing environment may hold and changes nothing for a
  fresh install. `versioning-strategy: increase-if-necessary` limits pip updates to releases
  that fall outside the declared range.
- **The `pytest` cap is widened to `<10.0`.** pytest 8.x carries PYSEC-2026-1845, whose fix is
  9.0.3, so the previous `<9.0` cap blocked a security fix -- worse than the surprise major it
  was guarding against. Verified: the full suite passes on pytest 9.1.1.

- **Remote BLAST results are cached on disk.** BLAST is the slowest step in recruitment --
  minutes per search, and now deliberately paced at one submission per 10 seconds to respect
  NCBI's guidance -- and every re-run, resumed run, or later round asking about the same
  window re-issued all of them. Results are cached under `$TESSERA_CACHE/blast`, keyed by the
  query sequence and every parameter that changes what comes back, so a hit is only ever
  reused for an identical request. Empty results are cached too: "nothing matched here" took
  just as long to establish. Unlike the panel caches this one expires (14 days, configurable
  with `TESSERA_BLAST_CACHE_DAYS`, `0` to disable) because `nt` grows continuously and a hit
  list reused indefinitely would quietly stop reflecting the database being searched. Reuse
  is logged with the entry's age. `--cache-dir` now reaches this cache as well, so isolating
  the cache directory isolates all of it.

## [1.1.0] - 2026-08-28

### Added

- **`CITATION.cff` and `CONTRIBUTING.md`.** Tessera had no machine-readable citation despite
  benchmarking against published datasets, and no written contribution guide -- in particular
  no statement of the rule that a change touching a caller or a default must be validated on
  `validation/run_hybrids.py` rather than on the unit suite alone.
- **`py.typed`.** About 91 % of public functions are fully annotated, but without the marker
  none of that was visible to a consumer -- including the third-party aligner adapters the
  `tessera.aligners` entry-point group invites. mypy is configured and runs in CI, permissive
  and non-blocking for now, since a never-type-checked codebase has an existing backlog.
- **A macOS CI runner and a packaging job.** Every check previously ran on Ubuntu only, while
  the project is developed on Darwin/arm64 and `docs/aligners.md` records platform-specific
  aligner problems there. A separate job builds the wheel, runs `twine check`, and installs it
  into a clean environment to confirm the console script works, so packaging breakage appears
  on a pull request rather than at release time.

### Changed

- **Three computations that grew unbounded with the panel are now bounded.** The PHI
  incompatibility matrix is quadratic in informative columns and had no guard at all; above
  25 M cells the columns are thinned to a sorted, seeded random subsample, which costs
  resolution rather than validity (`n_informative` reports the count actually used). Reference
  clustering, an all-pairs comparison gated only by a *lower* bound, gained an upper one.
  MaxChi's permutation null allocated several `(20000, sites)` int64 arrays per candidate
  donor; the permutations are now evaluated in batches, which bounds peak memory **without**
  reducing the null -- capping the permutation count instead was tried and reverted, because
  it engaged above 1000 discriminating sites (an ordinary panel) and dropped the specificity
  harness's positive control from 3/3 to 1/3.
- **The version is single-sourced** from `src/tessera/__init__.py`; it was duplicated in
  `pyproject.toml` with nothing enforcing agreement. The deprecated PEP 639 `license` form is
  replaced, and classifiers, keywords and `[project.urls]` are added.

### Fixed

- **Donor re-attribution no longer crashes when no major parent is resolved.** `run_recomb`
  passed `major_parent` (which is `None` when no reference wins a plurality of windows) into
  `lineage_of`, which raises `AttributeError` on `None`. Reachable with `--reattribute-donors`
  on a typed panel. Found by adding mypy.
- **`.gitignore` no longer lists tracked files.** `README.md` and `example_data/` were both
  ignored *and* tracked; it worked, but any new file under `example_data/` would have been
  silently dropped from `git add`, and `tests/conftest.py` depends on that directory.

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

- **NCBI requests are now paced, identified and retried.** Recruitment issued remote BLAST
  submissions back to back across gaps and sub-tiles, against a shared public service whose
  guidance is not to be contacted more than once every 10 seconds, with no retries and no way
  to supply an API key. A tool installed on many machines that ignores those limits is
  indistinguishable from abuse, and the practical result is a blocked address for whoever is
  running it. Requests to each service are now spaced (10 s between BLAST submissions,
  configurable with `TESSERA_BLAST_INTERVAL`; NCBI's 3-or-10 per second for E-utilities
  depending on whether a key is set), `NCBI_API_KEY` and a contact address are passed to
  EDirect through the environment -- which is how `efetch` and `datasets` accept them -- and
  transient failures are retried with exponential backoff instead of losing the coverage gap
  or seed window they were serving. A failure that persists is still raised: a search that
  could not run and a search that found nothing are different outcomes.
  Note that `--email` never reached the BLAST search and still cannot: Biopython's `qblast`
  accepts neither an e-mail nor an API key, so the BLAST URL API is unauthenticated and
  pacing is the only lever there. This is now stated in `docs/reference-panels.md` rather
  than implied otherwise.
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

[Unreleased]: https://github.com/FOI-Bioinformatics/tessera/compare/v1.2.0...HEAD
[1.2.0]: https://github.com/FOI-Bioinformatics/tessera/compare/v1.1.0...v1.2.0
[1.1.0]: https://github.com/FOI-Bioinformatics/tessera/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/FOI-Bioinformatics/tessera/releases/tag/v1.0.0
