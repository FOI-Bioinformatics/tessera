# Fourth audit: findings, SWOT and gap analysis

**Date:** 2026-10-12. **Audited commit:** `4b35e39` (branch `fix-audit-caller-behaviour`, the
plan A/B/C fixes applied, not yet merged). **Baseline measured at that commit:** 757 tests
pass (`pytest -m "not requires_binary"`), `mypy` clean over 79 source files, 89 % line
coverage, `ruff` clean over tracked files.

This is the fourth audit of Tessera. The first three covered inputs and caches (v1.1.0), the
output contract and the informative-site identity defect (v1.2.0), and the post-1.2.0 pass
whose three plans are the branch under audit here. This one re-read the source in five
parallel passes -- callers, orchestration and reporting, panel discovery, core/aligners/
converters, CLI and reassortment -- with instructions to report only defects the earlier
plans do not already cover.

Nothing in the repository was changed. Every finding below is labelled CONFIRMED (reproduced,
or traced end to end in the source) or SUSPECTED (mechanism established, not driven from a
real input). Reproduction scripts were written to a session scratchpad outside the repository.
The caller modules were read twice, by two independent passes. Both reported; the two agree on
the cluster-size emission defect (A1) and the Bootscan padding (A5), each reproduced
separately, and the rest of their findings are disjoint and merged below.

**What this audit could not reach.** No external aligner binary is installed, so the
progressiveMauve, sibeliaz, cactus and minimap2 findings rest on the code plus synthetic MAF,
XMFA, SAM and FASTA inputs, not on a run of the real tool. The network paths of `detect`,
`fill-references` and `build-panel` were not exercised (one BLAST taxon probe went out
inadvertently through `reassort`). Live Nextclade tree contents, large-input memory behaviour
and bacterial-scale alignments were not tested. The two validation harnesses that need an
aligner were not run, so nothing here is gated on `run_hybrids.py`.

## 1. Findings

### 1.1 Callers (`src/tessera/recomb/`)

**A1. The HMM emission is weighted by lineage-cluster size, so adding redundant references
can silence a true call. CONFIRMED, high.** `clusters.py:185-188` sums each cluster's
per-window numerators and denominators, and `hmm.py:126-134` then pads every state up to
`n_max = den.max(axis=0)` with `log(0.5)` per missing site. Cluster membership, not data
coverage, sets `n_max`, so a pooled cluster of k genomes carries k times the evidence while
each singleton state is charged a flat penalty for the shortfall. Reproduced: a 12 kb query
with a 1 kb tract from donor D against a panel of A, D, E, F reports the region; adding three
near-duplicates of the backbone A, which carry no new information, takes the HMM from one
region to none, and `--no-cluster-lineages` finds it again. In the window at the tract centre
the pooled A state (size 4, identity 0.944) emits -1177.3 against D (size 1, identity 0.997)
at -2100.7, because D is charged 3000 padded sites. At defaults the region survives on 3SEQ
and MaxChi, so the method list drops from `hmm,3seq,maxchi` to `3seq,maxchi` and the HMM
breakpoint interval is lost; a tract only the HMM can localise would be lost outright. The
docstring at `clusters.py:7-9` states that the pooled view leaves behaviour unchanged, which
holds only when every cluster is a singleton. The exposure is base-pair windowing; the
identity contrast under informative-site windowing was large enough to resist the flip.

**A2. Redundancy in the donor lineage splits one event into two regions and defeats
Bootscan. CONFIRMED, high on an untyped recruited panel.** Lineage clustering is applied only
inside `_call_regions_hmm` (`regions.py:280-283`), so the HMM names the cluster
representative while 3SEQ and MaxChi name whichever genome their own ranking picks, and
`ensemble.py:32-41` groups by genome label unless a lineage map exists. Separately
`bootscan.py:143` requires a single genome to reach 0.70 bootstrap support, which two
near-duplicates split. Reproduced with one true tract and the donor lineage represented by
one, then two, then three genomes 0.4 % apart: one genome gives a single region called by
four methods; two or three give two rows with different donors, one `hmm` and one
`3seq,maxchi`, and Bootscan alone drops from one region to none. Under `--min-methods 2` the
HMM row is suppressed and the surviving row names the other genome. The cluster membership
needed to group these is already computed in `clustered_view`; only its size reaches the
ensemble.

**A3. The `_MIN_DESCENT` floor blinds 3SEQ and MaxChi to short but significant tracts.
CONFIRMED, medium.** `threeseq.py:157` sets the floor at 4 and `threeseq.py:189-190` and
`maxchi.py:142-143` apply it before any p-value is computed, although the comment at
`threeseq.py:154-155` states that the p-value gate does the real work above the floors.
Reproduced: a query matching the donor at exactly three discriminating sites in a 189-site
walk has an exact p of 0.000169 and is discarded. At 1000 sites a depth of three gives
6e-06. This contradicts the module docstring's near-identical-panel power claim.

**A4. `--min-region` does not bound the length of a reported HMM region. CONFIRMED, medium
(contract, not statistics).** `regions.py:300` filters on the untrimmed segment and
`regions.py:327` then replaces the span with the trimmed tract from `_tract_span`. Reproduced
at the default floor of 1000: a segment of 11 200 bp passes the filter and the reported
region is 132 bp. Keeping the sign test on the segment is the deliberate choice recorded as
item C2; the consequence for the option's meaning is undocumented, and the help text in
`cli/cmd_recomb.py:70-73` is now false.

**A5. Bootscan's block permutation does not hold the win count fixed. CONFIRMED, low,
anti-conservative.** `bootscan.py:69-73` pads the window vector to whole blocks, shuffles,
then truncates to `n_win`, discarding wins that land past the end whenever `n_win` is not a
multiple of the block size, although the docstring at line 49 states the count is fixed.
Fewer wins give shorter null runs and a smaller p-value. Measured with one run of 12 wins at
block 10: `n_win` 100 gives a null mean of 12.00 and p 0.0930; `n_win` 103 gives 11.39 and p
0.0855.

**A6. The depth filter prunes the Benjamini-Hochberg family on a p-dependent statistic.
Arithmetic CONFIRMED, reachability SUSPECTED, low.** `threeseq.py:191-198` and
`maxchi.py:143-150` build the candidate list after the depth filter, so the family size
counts survivors only. Depth is monotone in p at fixed counts, so the least significant
candidates are removed and the rest inflated: a raw p of 0.004 gives q 0.2 among 50 tested
donors and q 0.012 among three survivors. How often a donor clearing 20 discriminating sites
has depth below 4 was not measured.

**A7. The MaxChi chi-square denominator can overflow int64. SUSPECTED, low.**
`maxchi.py:86` computes `r1 * r2 * tot_major * tot_minor`, which wraps above roughly 110 000
discriminating sites with balanced counts; a negative denominator fails the positivity guard,
forces that permutation's statistic to zero and drives p to the permutation floor. The
observed statistic uses Python integers, so the two sides disagree. Overflow confirmed in
isolation at 120 000 sites; not reached through `_exceedances`, because the null's own tract
stays short. Bacterial-scale alignments only.

**A8. The merged `n_windows` is a maximum over incommensurable units. CONFIRMED, low.**
`ensemble.py:111` takes `max(n_windows)` across callers, but the HMM and Bootscan count
windows while 3SEQ and MaxChi count discriminating sites.

**A9. The HMM sign-test p-value underflows to exactly zero above 1000 discordant sites.
CONFIRMED, high.** `stats.py:76` uses the exact binomial tail only up to n = 1000 and the
continuity-corrected normal approximation above it, and that returns 0.0 beyond a z of about
8.3. The boundary is a cliff, reproduced directly: 900 against 100 gives 6.7e-162, while 909
against 101, 1500 against 300 and 1827 against 200 all give 0.0. On a 25 % divergent panel a
called region with 1827 discordant sites therefore reports a p-value and q-value of zero
beside a weaker region at 2.1e-235 in the same run, so the strongest evidence carries the
least informative number. Zero q-values also tie, and `ensemble._best` then orders group
members arbitrarily and can carry the wrong member's test name and similarities onto the
merged row. This is the same defect class the changelog records as fixed for the 3SEQ exact
p-value (item C4), left in place for the sign test. The remedy is the same: compute the tail
in log space, with lgamma for the exact branch and a Mills-ratio asymptotic for the normal
branch.

**A10. The HMM segment half-width is the window step, not half the window width. CONFIRMED,
medium.** `hmm.py:182` sets `half` from the gap between the first two window positions. Under
base-pair windowing at 1000/100 a segment over windows 10 to 14 is reported as columns 1400
to 2000 rather than 1000 to 2400, and a single-window segment gets a 200-column span for a
window covering 1000. Under informative-site windowing it reads one unrepresentative gap and
applies it genome-wide: on a panel with a hypervariable 5' region the first gap is 33 columns
where half the median window width is 237, and a five-window segment is reported as 1137
columns against a true coverage of 2330. The project already wrote `median_window_step` for
exactly this flaw, and `tests/unit/test_clusters.py:156` names this expression as what the old
code used; the HMM was never migrated. The span feeds the `--min-region` gate and bounds the
trimming search, so marginal regions are dropped.

**A11. A merged region reports one member's identity over the union span. CONFIRMED,
medium.** `ensemble.py:107` takes the span as the union of the group while
`ensemble.py:112-114` keeps the best member's minor and major similarities and margin. An HMM
member over columns 1000 to 2000 at 0.99 merged with a 3SEQ member over 1900 to 6000 at 0.60
yields a region spanning 1000 to 6000 that reports a minor identity of 0.99, an identity over
a fifth of the span. This is the same class as the region-length defect already fixed in plan
B; the remedy is to recompute the identity over the merged span.

**A12. 3SEQ and MaxChi are not independent votes. CONFIRMED, medium, scientific rather than a
code defect.** Both score the same maximum-drawdown tract from the same triplet walk. Over 300
simulated clonal walks their p-values have a Spearman correlation of 0.97, and at alpha 0.05
the joint rejection rate of 0.030 equals 3SEQ's own, where independence would give 0.0015.
Both are individually well calibrated, at 0.03 and 0.05 against a nominal 0.05.
`docs/detection-methods.md` describes them as genuinely independent votes and treats a region
significant under both as corroborated twice, so an agreement gate of two satisfied by this
pair adds nothing. This bears directly on the open `--min-methods` decision, where the
remaining false regions at the agreement gate were exactly a 3SEQ and MaxChi pair.

**A13. Three memory sites sit outside the audited scaling caps. CONFIRMED, low.** At the PHI
subsampling cap the four float64 gamete matrices in the incompatibility routine are live
together, peaking at 892 MB against the roughly 225 MB its own comment budgets
(`diagnostics.py:131`); the Hudson-Kaplan bound materialises every incompatible pair as a
Python tuple, adding 320 MB at 30 % density (`diagnostics.py:213`); and the GENECONV
permutation allocates its full 40-million-cell matrix unbatched at 730 MB per donor
(`geneconv.py:52`), where MaxChi deliberately batches at 20 million. The scaling-guard tests
cover none of the three. This matters on an 18 GB development machine where the MAFFT
out-of-memory incident already occurred.

**A14. The barcode caller votes toward the agreement gate with no significance test.
CONFIRMED, low.** `barcode.py:155` sets no p-value, q-value or test name, yet the caller
counts toward `--min-methods`, which the documentation presents as a precision control. A
related small defect: a terminal tract gets a zero-width breakpoint interval at
`threeseq.py:218` and `geneconv.py:121`, because the lower bound collapses to the first column
when the peak is at index zero.

Two further items from this pass belong with the test gaps in section 3:
`similarity.py:192-193` leaves the final `(width - window) mod step` columns in no window
(499 of 5999 in one configuration), which the HMM mitigates through `_tract_span` but the
site callers and the coverage diagnostic do not; and `diagnostics.py:213` materialises every
incompatible site pair as a Python tuple, bounded only indirectly by the 5000-column matrix
cap.

### 1.2 Orchestration and reporting (`src/tessera/recomb/`)

**B1. `--no-phi` makes the report state that the alignment is too poor to test. CONFIRMED,
high.** `run.py:304-309` leaves `signal` as `None` both when the test was skipped and when it
returned nothing, and `report_html.py:501-506` renders a single message for that case: "Too
few informative sites in the alignment for a parent-free recombination test." Reproduced on
`example_data/divergent.msa.fasta`, which has 159 informative columns and a PHI p of about
0.001 when the test runs. Two knock-on effects: `corroborating_intervals(None)` returns an
empty list (`diagnostics.py:307`), so every region carries `parent_free_support = no` and a
test that never ran is rendered as a negative; and `run.py:529` adds the PHI provenance keys
only when a signal exists, so `run_provenance.json` records nothing about the test being off.
This is the pattern item B6 forbade for the barcode caller, recurring elsewhere.

**B2. `--min-region` is honoured by two of the seven callers. CONFIRMED, high.** The value is
read at `regions.py:300` and `regions.py:426`, `clusters.py:122` and `siblings.py:131,154`.
3SEQ, MaxChi, Bootscan, GENECONV and barcode never read it. Reproduced: with
`--min-region 5000`, a 175 bp 3SEQ region and a 49 bp MaxChi region are reported. Scoped
options elsewhere carry an explicit `(hmm)` marker in their help; this one does not.

**B3. A `low_information` coverage gap is promoted to the headline as a possible missing
reference. CONFIRMED, moderate.** `report_context.py:57-65` excludes only the breakpoint kind
from the caveat set, so a gap that `coverage.py:16-17` defines as "few comparable bases, so
the gap is uncertain rather than informative" reaches `_caveat_html`. Reproduced with a query
whose 2.5 kb stretch is 97 % N: the header advises running `find-references` to search NCBI,
when the real cause is that the query is unsequenced there. Recruitment is unaffected, so
this is report-only.

**B4. After HMM span trimming, `n_windows` and `posterior` still describe the untrimmed
segment. CONFIRMED, moderate.** `regions.py:338,341` keep the segment values while
`regions.py:327` has replaced the span, and `ensemble.py:111` propagates them. A row's window
count cannot be reconciled with its own length. `docs/detection-methods.md:75-83` documents
the same choice for the sign test and p-value but is silent on these two columns.

**B5. The provenance states base-pair geometry under informative-site windowing. CONFIRMED,
moderate.** `run.py:504` writes the window and step from the parameters unconditionally, so
the record holds `1000 / 100` beside `windowing = informative-site (40 sites/window, step
5)`. The glossary entry at `report_assets.py:107-109` describes only base-pair scans.

**B6. `--margin` and `--merge-gap` are inert outside `--method heuristic`. CONFIRMED,
moderate.** `params.margin` is read only at `regions.py:398` and `params.merge_gap` only at
`regions.py:416`; `run.py:443-447` builds the coverage parameters without `merge_gap`, so gap
merging always uses the window size. Neither help string carries a scope marker, so
`--margin 0.2` on a default run changes nothing and says nothing.

**B7. A single-state HMM is reported as a clean negative. CONFIRMED, moderate.**
`run.py:359` allows only the barcode caller to be recorded as "not run". When clustering pools
every reference into one lineage the HMM segments one state and returns a non-`None` major,
so the run reads as a negative rather than as a caller that could not run. Reproduced with
five references 0.25 % apart and a real 8 kb tract under `--no-informative-sites`: no regions
and a verdict of no recombination, beside a PHI p of 0.000999 in the same report. The default
path is protected by commit `a53ea53`.

**B8. The pair plot colours a dataset as a different dataset. CONFIRMED, low.**
`report_plots.py:153` builds its palette from the two plotted labels while the HTML and the
top-N plot use the ranked list, and colours are positional. Measured: one donor is green in
the report and orange in the pair plot, where orange denotes a different reference elsewhere.

**B9. The headline backbone is not the reconciled major parent. SUSPECTED, low.**
`report_html.py:61-66` uses the first region's major parent, or the first dataset, while the
reconciled major from `run.py:380` never reaches the report, so the verdict and the provenance
can name different backbones.

**Smaller reporting items, confirmed by reading.** Donor-absent regions are appended after
the method breakdown is final (`run.py:485`), so `recombination_methods.tsv` has fewer rows
than `recombination_regions.tsv` and a positional join breaks.
`report_text.py:252-255` points at a provenance key for coordinate conventions that the
provenance dict does not contain. `report_html.py:402-403` states that every caller applies
Benjamini-Hochberg, which the barcode caller does not. `report_assets.py:121-123` attributes
the breakpoint interval to the HMM posterior even for a site-caller-only region.
`report_assets.py:99-100` defines the major parent as the reference winning the most windows,
while the reported one is the HMM's majority state after sibling exclusion and pooling.
`report_text.py:157-158` still describes the superseded gap-wide parent rule, and
`example_data/README.md` still advises `--phi-window 5` for a case the PHI cap now handles
automatically.

**Checked and clean.** Every label reaching the HTML passes `html.escape` with quoting,
including inside the plotly payload, so a label containing a script tag cannot break out.
The agreement filter and the consensus region list stay parallel. HMM similarities are
computed on base-pair rows, so the informative-site scale does not leak into the reported
means.

### 1.3 Core, aligners and converters

**C1. The progressiveMauve backbone row holds only the columns the first query aligned to.
CONFIRMED, high.** Each per-query projection carries a reference row gapped wherever that
query did not align, and `progressivemauve.py:114-115` writes the reference row from the
first projection only. Reproduced: with a 20 bp backbone and two queries covering opposite
halves, half the backbone is lost, and which half depends on list order. The row stays full
width, so the ragged-row guard at `similarity.py:172` never fires; downstream the backbone is
both the coordinate anchor and a donor candidate. The repository's own rectangularity test
exhibits this and asserts only widths. The adapter already reads the reference sequence at
line 56.

**C2. Duplicate sequence IDs inside one genome silently misplace blocks. CONFIRMED,
medium-high.** Guards exist for an ID shared between genomes (`sibeliaz.py:171-181`) and for
colliding file labels (`io.py:312`), but both contig layout tables are dicts keyed on the ID.
A backbone with two records named `c1` yields a reference length of 8 with a single offset, so
the backbone row comes back half gapped and a donor aligned to the first contig is written
onto the second contig's columns. The minimap2 adapter (`minimap2.py:51-59`) builds its
offsets the same way, where the reference row is correct and every query row shifts. Both
produce a full-width, well-formed, wrong alignment with no warning.

**C3. The minimap2 adapter multiplies its thread and memory budget by the thread count.
CONFIRMED, medium.** `minimap2.py:62,69,77` passes `-t <threads>` to each process and also
runs the mapping pool with the same thread count, so `--threads 16` starts 16 processes of 16
threads, each building its own reference index. The MAFFT adapter solved exactly this in
commit `44e6b3f`; minimap2 was not changed, and it is the backend used for the 200 kb
orthopox and mpox cases.

**C4. The Nextclade pool cache is installed with a bare rename. CONFIRMED, medium-low.**
`nextclade.py:352` uses `os.replace` instead of the `atomic_cache_dir` helper that exists for
this and handles the collision explicitly (`cache.py:163-182`). Replacing a non-empty
directory raises an OSError that surfaces as "Unexpected error". Two triggers: two concurrent
runs on one dataset, and a cache directory holding only files the reader does not recognise,
which reads as absent so every run rebuilds and then fails. The pool cache also carries no
completeness manifest, unlike the NCBI Virus cache.

**C5. The aligner documentation overstates the unplaced-genome warning. CONFIRMED,
low-medium.** `docs/aligners.md` promises that a genome the aligner cannot place is named in a
warning. That warning exists only in `maf_to_fasta.py:141-147`, so for sibeliaz and cactus.
Under minimap2 and MAFFT the genome becomes an all-gap row with nothing logged, surfacing
later at `run.py:289` as a message about its alphabet.

**C6. A MAF row cut before its text column bypasses the truncation guard. CONFIRMED, low.**
`maf_to_fasta.py:215` keeps a row only when it has seven fields, so a truncated row is dropped
during parsing and the width guard written for exactly this condition never sees it. The
genome becomes an all-gap row reported as sharing no block with the backbone, which reads as
divergence rather than as a truncated file.

**C7. A plain file at the working collection path leaks a raw error. CONFIRMED, low.**
`io.py:160` calls `rmtree` without checking that the destination is a directory, giving
`NotADirectoryError` where the surrounding code is built to raise a user-input error.

**Suspected and minor.** `cactus.py:95-97` maps dots to underscores in sample names, so two
accessions differing only in that character collapse to one row. Child processes are started
without redirecting standard input (`binaries.py:47`, `process.py:96,116`), and alignment has
no default timeout, so a tool that waits on input would hang rather than fail; no current
call site triggers it. The `flank > 0` path in `xmfa_to_fasta.py:162-186` is unreachable and
internally wrong. `collection_genomes` skips genomes in subdirectories. The logger replaces a
file handler without closing it.

**Two carried-over items settled.** The progressiveMauve leaf-name defect is already fixed:
rows are identified by position and labelled from the staged file. The sorted-contig-order
defect in the MAF converter is real but reachable only through the cactus backend, which does
not supply the contig order, and both the previous audit spec and `docs/aligners.md` record
that choice.

**Caller internals verified correct, so a later pass need not re-walk them.** The exact 3SEQ
descent p-value matches brute-force enumeration over all arrangements for every m and n up to
7 and depth up to 5; the last-maximum peak rule agrees between the 3SEQ walk and the MaxChi
permutation; the Benjamini-Hochberg routine is a correct step-up capped at 1.0; the base-pair
window bounds are reconstructed exactly, odd widths included; four-gamete missing-data
handling and the Hudson-Kaplan greedy disjointness test are right; and both triplet callers
are individually calibrated under the clonal null.

**Core internals verified correct.** SAM CIGAR consumption including
hard clips and padding, minus-strand MAF block placement, CRLF handling in the FASTA and XMFA
readers, the minimap2 argument order, and the BLAST cache key.

### 1.4 Panel discovery (`src/tessera/discover/`)

**D1. The fill loop downloads every new accession from every sub-tile, not one per gap.
CONFIRMED, high.** `run.py:272` fetches all deduplicated accessions although its own docstring
and `iterate.py:1027` both state one reference per gap. `_tiles` (`run.py:125`) splits a gap
into half-overlapping sub-intervals, each contributing up to the per-search hit cap, so at the
defaults a 1000 bp gap can add 60 genomes and a 6 kb gap 435, times the number of searched
gaps. Reproduced: 40 candidates from one gap produced 40 downloads. This feeds the next
alignment, which is the MAFFT memory path just bounded, and makes the added-references column
of `fill_summary.tsv` a list hundreds long.

**D2. The stop test watches the worst gap while the search targets the longest, so the loop
exits after round 2. CONFIRMED, high.** `iterate.py:434` compares the minimum best-similarity
over all gaps against the previous round, while `iterate.py:445` searches only the longest
few. An unsearched short gap pins the metric, so improvement always stalls. Reproduced with
`--max-rounds 5`: three long gaps improved from about 0.70 to about 0.90 and the run stopped
at round 2 reporting that coverage had stopped improving at an unchanged 0.550. The
already-computed undercovered base count, or the mean over the searched gaps, would be
faithful.

**D3. The default typing path discards the lineage sidecar it just wrote. CONFIRMED,
medium-high.** The non-deep branch at `iterate.py:546` builds the lineage map without the
structured rows, so lineages are re-mined from FASTA headers by a rule that needs a token of
four or more characters containing a digit, while `<output>/ncbi_lineages.tsv`, written during
recruitment at `iterate.py:485`, holds the exact values and is read only by the deep branch
(`iterate.py:531`). Measured on seven realistic headers, five lost their lineage on the
default path and all seven survived when the sidecar was supplied. Untyped genomes skip
lineage-aware selection and recombinant exclusion. This is the same loss as item C6 with a
cheaper remedy.

**D4. Windowed seeding reports saturation when every hit was filtered out. CONFIRMED,
medium.** `iterate.py:1000` computes the flag from hits taken before filtering
(`iterate.py:963`) while scanning the kept hits, so an empty kept list makes the test
vacuously true. Reproduced: a pool where every hit is the query's own record, or the only hit
is excluded by the user, reports saturation; no hits at all does not, and the whole-query
seeder gets it right on the same input. The consequence is a full negative-lineage sweep and
an entire taxon download, logged with a false explanation.

**D5. Zero RefSeq genomes abort the run instead of broadening. CONFIRMED, medium.** The fetch
at `iterate.py:776` requests RefSeq only and `pool.py:171` raises when the download holds no
genomes, so the broadening on the next line is unreachable for the zero case. Reproduced: a
taxon with no RefSeq and 50 complete genomes aborts; the same taxon with two RefSeq genomes
broadens and returns all 50.

**D6. `--curate` can be requested, never run, and never reported. CONFIRMED, medium.** The
pending flag starts false for a freshly seeded collection by design and is set only once a
round downloads something (`iterate.py:389,461`), while a seeded sibling-only panel has no
coverage gap and breaks out of round 1 (`iterate.py:431`). Reproduced: curation returned no
panel rows, and because the report gates its panel section on that dict (`iterate.py:576`)
there is no panel table, no lineage file and no log line. The masking panel is the case
curation exists for.

**D7. Pool selection keeps every sibling, silently, when all of them are siblings.
CONFIRMED, medium.** `pool.py:328` guards the drop with a condition that fails when the
sibling set is the whole pool. Reproduced with four genomes at 99 % identity over 99 %
coverage: none dropped, no warning, a normal-looking panel. The analogous seeding path logs
clearly.

**D8. `--fetch-limit` is in the cache key but does not change the fetch. CONFIRMED, low.**
`iterate.py:762` puts it in the cache scope and the fetch never passes it, so its only effects
are a log line and a cache slot: changing the flag forces a full re-download of an identical
genome set.

**D9. Regional selection degenerates to whole-genome selection for a short query. CONFIRMED,
low.** `pool.py:338` divides the query length by the window with integer division, so a 2.9 kb
genome at the default 1500 bp window gets one window spanning the whole query. This silently
removes the feature for small genomes such as hepatitis B and hepatitis delta.

**D10. The non-clade marker set is not used by the consumer it was written for. CONFIRMED,
low.** `nextclade.py:50` defines four markers and states that consumers should treat all four
as "no clade"; `lineage_assign.py:121` hardcodes two of them.

**Suspected.** The negative-lineage seeder BLASTs the first 1500 bases while its docstring
claims a whole-query search (`iterate.py:872`), so for a recombinant query the dominant token
can be one parent's lineage and the exclusion then removes a lineage the search wants.
Reference typing resolves its Nextclade dataset from whichever untyped genome sorts first
(`lineage_assign.py:153`), so a divergent outgroup there types the collection against another
organism; the identity floor rejects the result, so the cost is a wasted call and an untyped
panel.

**Recorded items confirmed as described.** The RSV-B and HIV-2 alias fallbacks behave as
documented (substring match, first insertion-order hit wins), and curated-away siblings are
re-downloaded each round because the exclusion list is never extended.

### 1.5 CLI and reassortment

**E1. Four commands crash with a raw traceback when the output path names an existing file.
CONFIRMED, medium.** `main.py:131-133` creates the output directory through the logger, and
`cmd_detect.py:122`, `cmd_build_panel.py:93`, `cmd_type_lineages.py:52` and
`cmd_reassort.py:56` call it outside their error handler; only `recomb` checks the path first.
Observed: a traceback on four commands, a clean message on `recomb`, and an "Unexpected
error" on `fill-references`. The previous audit's item B11 covered `recomb` alone.

**E2. `reassort --margin` has no upper bound, and a wide margin inverts the verdict.
CONFIRMED, medium.** `cmd_reassort.py:62` checks only the lower bound, and the near-best set
(`constellation.py:31-36`) grows with the margin, so once it covers the ANI spread every pair
shares a strain and every query reads clonal. Measured on one dataset: 0.5 gives a
reassortant verdict with two groups; 40, 100 and 1e6 all give clonal with one. A margin of
zero does the reverse, making a 0.01-point tie-break read as reassortant.

**E3. A reassortant verdict can be reported with a single unnamed parent group. CONFIRMED,
medium.** `constellation.py:62-68` builds groups from concordant edges only and then lets any
discordant pair set the verdict, so a middle segment concordant with both ends collapses all
three into one component while the ends disagree. The user sees a reassortant verdict beside
one group with no parent strains, and every row of the output names group 0. The documented
explanation for an empty parent group (`docs/reference-panels.md:347-351`) covers only the
variant that stays clonal.

**E4. The evidence for a reassortant or clonal call is computed and discarded. CONFIRMED,
medium.** `cmd_reassort.py:110-111` logs the pairwise notes only for an undetermined verdict,
so a reassortant call records nowhere which pair disagreed or on which strains.

**E5 to E13, lower severity, all confirmed.** The `--seed-source` help omits the `nextclade`
value the code accepts and the documentation describes. `--candidate-pool` silently overrides
`--nextclade`, and `--pool-consensus` is inert under any other seed source.
`build-panel` cannot set three options that shape the panel it exists to build
(`--keep-recombinant-lineages`, `--deep-typing`, `--pool-consensus`). The two coverage
fractions are range-checked on `recomb` only, so `--coverage-floor 50` is a clean error on one
command and makes every window a gap on two others. The candidate-pool reader accepts four
extensions where every other collection reader accepts any non-hidden file, so a `.fas` genome
is dropped silently. `reassort --ani-floor 0.8` is accepted where 80 is meant, and `--margin`
is a 0-1 similarity fraction on `recomb` but ANI percentage points on `reassort`.
`find-references -q` is an alignment label declared as a path. `msa --reference` silently
overrides `--query-as-backbone`. Six count options are unbounded on the discover commands, and
`--top-gaps 0` searches nothing and exits cleanly. Minor: `msa` overwrites an existing
alignment without notice; `reassort` has no thread option, so segment scans align
single-threaded; `--metric` has one legal value; a computed return value is discarded; and the
per-segment clade is tie-broken by file name with no indication that several clades sat inside
the margin.

**E14. The exported constellation function can return a reassortant verdict from an empty
candidate list. SUSPECTED, low.** `constellation.py:51-58` sets the discordant flag when one
segment has no near-best set, which is unreachable through the CLI but unguarded in the
public API.

### 1.6 Repository, harnesses and documentation (own checks)

**F1. The specificity harness hand-parses its arguments, so `--help` runs the whole harness.
CONFIRMED, low.** `validation/run_specificity.py:221-230` scans `sys.argv` for two flags with
no parser, so `--help` starts a ten-replicate run, an unknown flag is ignored silently, and a
missing or non-numeric value raises. The documented invocation in `CLAUDE.md` and
`CONTRIBUTING.md` is this script.

**F2. The documented `PATH` setup for the aligner environment breaks the harness. CONFIRMED,
low.** `CLAUDE.md` and `validation/README.md` both prepend the aligner environment to `PATH`,
which also replaces the Python interpreter. Appending is the form that works, and the session
record of a previous run notes this.

**F3. `CONTRIBUTING.md` states that mypy is not blocking in CI and that a backlog exists.
CONFIRMED, low.** CI runs it as a blocking step and it is clean over 79 files.

**F4. The README development section omits the validation path from the lint command.
CONFIRMED, cosmetic.** CI and `CONTRIBUTING.md` both lint `src tests validation`.

**F5. The measured false-region rate at the shipped default is not user-facing. CONFIRMED,
low.** The open decision record states that the README already describes the agreement-gate
trade, but the README does not mention `--min-methods` at all. The CLI help describes the
trade qualitatively and gives no rate, while the harness measured a false region in 17 of 40
clonal runs at that default.

**F6. Permutation seeds are fixed per call and not settable.** Runs are reproducible, which is
the right default, but a borderline p-value cannot be checked for seed sensitivity.

**F7. Untracked scratch scripts under the validation data directory fail the documented lint
command locally.** They are git-ignored, so CI is unaffected.

## 2. SWOT

### Strengths

- **Measured validation, reported honestly.** Specificity and the matching positive control
  are reported together with Wilson intervals, the harness default is documented as differing
  from the CLI default, and a currently failing real clonal control is recorded in the
  validation README rather than hidden.
- **The dependency-free design holds.** The PHI/Rmin, 3SEQ, MaxChi, Bootscan and GENECONV
  callers are numpy reimplementations, there is no `shell=True` anywhere, external tools go
  through one choke point, and the runtime dependency list is six capped packages.
- **Aligner backends are genuinely pluggable** through an entry-point registry that core never
  bypasses, and the registry is exercised by tests.
- **Engineering baseline is strong for a research tool**: 757 tests against 13 700 lines of
  source, 89 % line coverage with a CI floor, mypy clean and blocking, a three-version Linux
  matrix plus macOS, SHA-pinned actions, a container build in CI, and provenance written per
  run.
- **An audit culture exists.** Three previous audits left written specs, per-item designs,
  recorded revisions after verification, and an open decision record that states what is still
  the maintainer's call.

### Weaknesses

- **Panel redundancy is handled inconsistently across the statistical core.** The HMM emission
  is weighted by cluster size (A1), clustering is applied in one caller only so the ensemble
  splits a redundant donor (A2), Bootscan's support threshold is split by near-duplicates
  (A2), and the pool selector keeps every sibling when all of them are siblings (D7). The
  tool's own recruitment story is that more references help; in these paths they hurt.
- **"Could not test" is still rendered as "tested, found nothing"** in paths the previous
  audit did not reach: a disabled parent-free test (B1), a single-state HMM (B7), an
  uninformative coverage gap promoted to a missing-reference caveat (B3), a saturation flag
  computed from pre-filter hits (D4), and curation that never ran and is never mentioned (D6).
- **Options are accepted that do nothing, or that do less than their help says.**
  `--min-region` on five callers (B2), `--margin` and `--merge-gap` outside one caller (B6),
  `--pool-consensus` outside one seed source, `--fetch-limit` entirely (D8), and
  `--candidate-pool` silently overriding `--nextclade`.
- **The converter and adapter layer has no backbone invariant.** The progressiveMauve backbone
  row can lose half its columns (C1), a duplicate contig ID shifts donor rows onto the wrong
  contig (C2), and a truncated MAF row becomes an all-gap row reported as divergence (C6). All
  three produce a well-formed alignment that the user cannot falsify.
- **Resource budgets are per-adapter rather than shared.** The MAFFT memory fix was not
  carried to minimap2 (C3), and the fill loop's download fan-out is unbounded (D1).
- **The recruitment loop stops early and over-fetches at the same time** (D1, D2), and the
  default typing path throws away structured lineages it already wrote to disk (D3).
- **Reported significance is not always usable.** The sign test returns exactly zero above a
  thousand discordant sites (A9), a merged region reports one member's identity over the union
  span (A11), the merged window count mixes two units (A8), and the two callers the
  documentation offers as independent corroboration score the same tract (A12).
- **The HMM's significance gate is not aware of the scan**, which is the mechanism behind the
  false-region rate at the shipped default; the decision record is open and the rate is not in
  the README (F5).
- **Two shipped aligner backends have never been run against their real binary** in CI or in
  the harnesses, and the only end-to-end binary coverage is MAFFT.
- **No API or CLI reference documentation**, open since the August audit.

### Opportunities

- **The two high-severity caller defects share one fix surface.** Cluster membership is
  already computed in `clustered_view`; passing it to the ensemble groups a redundant donor's
  rows, and averaging rather than summing within a cluster removes the emission weighting.
- **The infrastructure to gate a scan-aware HMM test already exists** (`run_specificity.py`,
  `run_regimes.py`, `characterise_false_regions.py`, the alpha sweep), so option B of the open
  decision can be measured without new harness work.
- **Helpers for three findings are already in the tree**: the MAFFT shared thread budget for
  minimap2, `atomic_cache_dir` for the Nextclade pool, and the lineage sidecar for the default
  typing path.
- **Several invariant tests are cheap and would have caught high-severity defects**: the
  concatenated backbone row equals the reference for every backend; every option changes some
  output; a reported region respects `--min-region`; a test that never ran is not a negative.
- **Several findings have a remedy the project has already written once**: the 3SEQ log-space
  p-value for the sign test, `median_window_step` for the HMM segment span, the region-length
  recomputation for merged similarities, and the MAFFT shared budget for minimap2.
- **The validation harnesses are a publishable asset.** A specificity rate with intervals
  beside a matched positive control is more than most comparable tools report.

### Threats

- **A real clonal control currently fails** (`hcv_clonal_1b`, a 12 bp single-caller region)
  while the shipped default reports a false region in roughly four of ten simulated clonal
  runs, and neither number is user-facing. This is the kind of result a reviewer finds first.
- **Silent wrong alignments are unfalsifiable by the user.** C1, C2 and C6 all yield a
  full-width alignment with no warning, and the scan then reports divergence that is an
  artefact.
- **Adding references can remove a detection** (A1), which contradicts the documented
  recruitment workflow and is the opposite of what a user expects from a larger panel.
- **Scope has outgrown the verification budget.** Five aligner backends, nine commands, a
  reassortment module and twelve harnesses rest on one maintainer, one end-to-end binary in
  CI, and an audit backlog now four passes deep with the third pass unmerged.
- **The unmerged branch carries real behaviour changes.** Its final gate was taken from
  simulated harnesses only, with `run_hybrids.py` recorded as an executor gate that has not
  been run on the current tip.

## 3. Gap analysis

Ordered by the product of impact and cost to close. "Cheapest next step" is deliberately
smaller than a fix: in several cases the next step is a measurement, because unit-level
reasoning about caller behaviour has been contradicted by the harnesses twice.

### G1. Panel redundancy distorts the callers (highest priority)

**Evidence:** A1, A2, D7, and the already-recorded C1 collapse. **Impact:** a true event can
be lost (A1) or split into two contradictory rows (A2); both get worse as a panel grows, and
both are invisible in the output. **Why now:** these are the only findings in this audit that
can remove a detection rather than add or mislabel one. **Cheapest next step:** reproduce A1
in a unit test with clusters of unequal size feeding the emission matrix, then measure
averaging-within-cluster on `run_specificity.py` at both gates and on `run_regimes.py` before
touching defaults. Pass cluster membership to the ensemble and assert that one event on a
redundant donor lineage yields one region.

### G2. Reported p-values, and what corroboration means

**Evidence:** A9, A11, A12, A8, and the already-fixed 3SEQ underflow. **Impact:** the HMM's
strongest regions report a p-value and q-value of exactly zero, which also ties them in the
ensemble's ordering; a merged region reports an identity measured over a fifth of its span;
and the two callers the documentation offers as independent corroboration score the same
tract, so an agreement gate satisfied by that pair adds no information. The second point
bears on the open default-gate decision, whose remaining false regions at the agreement gate
were exactly that pair. **Cheapest next step:** compute the sign-test tail in log space and
add a test above a thousand trials; recompute the merged region's identity over its own span;
and correct the independence claim in `docs/detection-methods.md`, stating which caller pairs
are genuinely distinct.

### G3. Backbone and alignment integrity has no invariant

**Evidence:** C1, C2, C6, plus the recorded cactus contig-order item. **Impact:** a wrong
alignment is a wrong scientific result, and all three failures look well-formed. **Cheapest
next step:** one test per backend asserting that the concatenated backbone row equals the
reference sequence, and a duplicate-sequence-ID rejection at staging. Both are cheap and the
first reuses an existing fixture that already reproduces C1.

### G4. "Not testable" must not be reported as "negative"

**Evidence:** B1, B3, B7, D4, D6. **Impact:** the failure mode the project has already
committed to avoiding, now found in five further places; a user reads a skipped test as
evidence of absence. **Cheapest next step:** make the distinction explicit once, in the
reporting context rather than per call site, so that a caller or test that did not run has a
state of its own, and add the five assertions named in the test gaps.

### G5. The recruitment loop over-fetches and stops early

**Evidence:** D1, D2, D3, D5, D6. **Impact:** hundreds of genomes downloaded per round feed the
alignment step whose memory was just bounded; the loop then stops while coverage is still
improving; and the panel it produces is partly untyped. **Cheapest next step:** cap the
download to the documented one reference per gap, switch the stop test to the undercovered
base count already computed, and pass the lineage sidecar on the default typing branch. All
three are local and unit-testable.

### G6. Options that do not do what they say

**Evidence:** B2, B6, A4, D8, E5-E13. **Impact:** a user who sets a threshold and sees no
change cannot tell a no-op from a null result; the `--min-region` help is now false.
**Cheapest next step:** a convention plus a test shape, not a sweep: every option either
changes an output in a test or carries an explicit scope marker in its help. Then fix the
help text for the scoped ones and remove or implement `--fetch-limit` and `--metric`.

### G7. Resource budgets and failure handling in the adapters

**Evidence:** C3, C4, C7, E1, D1. **Impact:** an avoidable out-of-memory on the 200 kb cases
on a development machine with 18 GB, an aborted run on a cache race, and tracebacks where
clean messages exist elsewhere. **Cheapest next step:** carry the MAFFT shared-budget helper
to minimap2, switch the pool cache to `atomic_cache_dir`, and move the four `get_logger`
calls inside the error handler.

### G8. The open default-gate decision, and what the user is told

**Evidence:** F5, the open decision record, the failing `hcv_clonal_1b` control. **Impact:**
the headline false-positive rate of the shipped default is documented only inside the
validation directory and an internal decision record. **Cheapest next step:** independent of
which option is chosen, state the measured rate and the recommended setting for a redundant
panel in the README, and align the harness and CLI defaults or say plainly that they differ.
The audit's own recommendation remains option C now and option B as a design task gated on
`run_hybrids.py`.

### G9. Verification breadth

**Evidence:** no binary coverage for progressiveMauve, cactus or sibeliaz; `run_hybrids.py`
not run on the current branch tip; the reassortment module's intermediate case untested
(E3); the harness entry points hand-parse arguments (F1); no API or CLI reference docs.
**Impact:** two backends are shipped as supported on the strength of code reading alone, and
the branch about to merge has not met its own recorded gate. **Cheapest next step:** run
`run_hybrids.py` on the branch tip before merging, add the transitive-discordant constellation
test, and give the harness entry points a real argument parser.

### Test gaps named by the audit, consolidated

Caller side: the sign test above a thousand trials, which no test covers, which is why A9
survived three audits; the HMM segment span under a window step smaller than the window
width, which `test_hmm.py` cannot see because its fixture sets step equal to window; a merged
region's identity against its own span; the PHI window a test names against the one
`effective_phi_window` leaves it with (`test_diagnostics.py:105` asks for 8 and gets 2); the
Hudson-Kaplan and GENECONV memory paths, which the scaling guards omit;
`clustered_view` with clusters of unequal size feeding the emission matrix; two
near-duplicates of the donor in an untyped panel yielding one region; what the
pre-Benjamini-Hochberg floors exclude, pinned by an exact p-value; the permuted win count
when the window count is not a multiple of the block; any relation between `--min-region` and
a reported HMM region length; the trailing unscanned columns; a dense incompatibility matrix;
the exact-versus-permutation switch.

Reporting side: the wording of the signal section under `--no-phi` and that the PHI column is
not a negative; a low-information gap reaching the headline; `n_windows` and `posterior`
against a trimmed span; `--min-region` against a site caller; the provenance window line under
informative-site windowing; one dataset's colour across the three figures; a single-state HMM
recorded as "not run"; the HTML verdict's backbone equal to the provenance major parent;
`--margin` and `--merge-gap` being inert outside the heuristic caller.

Core side: the backbone row equal to the reference for every backend; one genome with two
records sharing a sequence ID; the minimap2 worker count; the short-MAF-row branch; a source
size disagreeing with the declared contig length; the Nextclade pool install under
concurrency and with stray files; a non-directory at the working collection path.

Discover side: how many genomes one round may add; the stop test with more gaps than the
search covers; saturation with zero surviving hits; zero RefSeq genomes; what a run reports
when curation never runs; an all-sibling pool; the window count for a short query; a pool
build over an existing cache; the default and deep typing branches compared on one panel.

CLI side: an existing file as the output path on the five uncovered commands; the margin upper
bound and zero; the transitive-discordant constellation shape; the pairwise evidence reaching
the user; option interactions (`--candidate-pool` with `--nextclade`, `msa --reference` with
`--query-as-backbone`); range checks on the discover commands; the two collection readers
accepting the same file set; a one-clade consensus panel and an empty regions file in the
segment scan.

## Revisions after the plans were written

The four implementation plans (`docs/superpowers/plans/2026-10-12-audit-4-index.md`) were
written against `4b35e39`, and every planned test was run there. That checked the findings
above against the code a second time and corrected several. The plans and the index hold
where they disagree with this document; the main corrections are:

- **The aligner statement in the introduction is wrong.** The `recomfi-aln` environment
  holds mafft, minimap2 and sibeliaz 1.2.7; only progressiveMauve and cactus are absent. The
  audit passes ran without that environment on `PATH`.
- **New: the default sibeliaz backend writes an empty alignment on macOS.** Its wrapper uses
  two GNU-only commands, and the adapter attributes the empty output to memory (plan D Task
  3). Whether this explains the recorded orthopox out-of-memory note was not tested.
- **D1:** at the default `--max-hits 5` the bounds are 20 and 145 genomes per gap per round,
  not 60 and 435. **D9:** hepatitis B already gets two windows; the affected range is 1.5 to
  3 kb. **E1:** seven commands are affected, not five. **C7:** the raw error needs a previous
  run's marker or a symlink at the path.
- **B1** also covers PHI runs with too few informative sites or not testable at their
  window; **B9** reproduces in the no-region case and is confirmed.
- **A3's remedy** as stated costs about ten times the runtime; **A5's example numbers** do
  not reproduce, though the defect does; **A14** also applies at `maxchi.py:170`; **A10**
  read literally changes the HMM's significance, not only its coordinates, and is referred
  to the open `--min-methods` decision.
- **Found while planning:** a donor-absent region's parent-free support is never computed;
  a donor copy 0.03 % away is not pooled under informative-site windowing; GENECONV shares
  the 3SEQ/MaxChi floor and family pruning.
