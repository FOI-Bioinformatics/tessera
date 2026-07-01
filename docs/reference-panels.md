# Building and curating the reference panel

Detection is only as good as the donor panel. For a recombinant query the most useful
panel holds the *distinct parental lineages* and none of the query's own siblings
(whole-genome relatives that match it everywhere and mask the event). This guide
covers recruiting a panel, finding missing donors, and curating for detection.

The one-shot `tessera detect` command does all of this automatically; reach for the
pieces below when you want to inspect or control the panel.

## Build the panel separately from detection

To recruit the donor panel without running detection -- so you can inspect the panel,
or re-run detection with different window/HMM settings without re-fetching or
re-aligning -- use `build-panel`, then `recomb`:

```
tessera build-panel --query CRF01_AE.fasta --output panel/ --email you@example.org
tessera recomb --msa panel/panel.msa.fasta --query CRF01_AE --output panel/
```

`build-panel` uses the same parent-recruiting, sibling-dropping recruitment as
`detect` but stops at the panel, writing `panel/collection/` (the donor genomes), a
stable `panel/panel.msa.fasta` alignment, and `panel/panel_lineages.tsv`. The exact
`tessera recomb` follow-up command is logged on completion. The panel is recruited for
this query (donors are its regional best matches), so it is not a reusable
organism-wide reference set. Pass `--collection` to grow an existing collection
instead of recruiting from scratch. The full-featured `fill-references` command takes
`--no-report` for the same effect.

## Typed lineage names in the report

When `detect`, `build-panel`, or `fill-references` recruit a panel, each reference is
typed with a genotype mined from its genome header -- the NCBI lineage note
(`pangolinClassification` / `isolate.lineage`) or a designation token in the GenBank
title (e.g. `GII.P16-GII.4`, `CRF01_AE`, `B.1.1.7`). The mapping is written to
`<output>/lineages.tsv`, and the report then names parents by lineage, e.g.
"GII.P16-GII.4 (MK573073)" instead of the bare accession. A standalone `tessera
recomb` picks up a `lineages.tsv` sitting beside the output or the MSA, so the
two-step workflow carries the typed names through. To override or supply names (for
schemes that titles do not carry well, such as HIV pure subtypes or a curated CRF
reference set), pass `--lineage-map accession_to_genotype.tsv` (an
`accession<TAB>genotype` table) to any of these commands; user-supplied names take
precedence over mined ones. With no source, references fall back to bare accessions
and the report is unchanged.

The query is also typed from its own header; when its lineage differs from the file
name the verdict states it ("The query is typed as GII.P16-GII.4"). For a SARS-CoV-2
query whose lineage is a Pango recombinant (an `X` lineage such as `XBB.1.5`), Tessera
looks up the designated parents in the Pango `alias_key.json` (fetched once into the
cache) and adds a cross-check block, so the recruited parents can be compared with the
designated ones. The lookup is best-effort: a network failure simply skips it.

## Is a reference missing?

Tessera always reports the *closest* reference, even when every reference is far from
the query -- so a recombinant whose true donor is not in the collection is still
assigned to the least-bad reference. To catch this, the scan reports **reference
coverage**: stretches where even the closest reference is below a best-similarity
threshold (adaptive by default, or set with `--coverage-floor`). These are written to
`coverage_gaps.tsv`, shaded on the report, and flagged with a caveat banner; a region
labelled `divergent` (ample comparable bases, yet a poor match) is the signature of an
absent reference, as opposed to `low_information` (too few comparable bases to judge).

### Find one missing reference

When a gap is found, `tessera find-references` searches NCBI for the missing genome:
it BLASTs the under-covered query subsequence against `nt`, reports candidate
references (accession, identity, whether already in your collection), and -- with
`--download` -- fetches the best new one into a collection directory to re-run with.

```
tessera find-references --msa msa.fasta --query my_query \
    --collection collection/ --email you@example.org --download collection/
```

This contacts NCBI over the network; `--download` needs Entrez Direct (`conda install
-c bioconda entrez-direct`). BLAST often returns the query's own GenBank record (the
MSA labels the query by name, not accession), so a near-identical, near-full-length
hit is auto-skipped; use `--keep-self-hits` to keep it, or `--exclude <accession>` to
drop specific records.

A gap longer than `--subtile` (default 400 bp) is searched in overlapping
sub-intervals, not as one subsequence. A short divergent tract inside a longer gap is
otherwise diluted by the gap's high-matching flanks: one whole-gap BLAST returns the
flanks' donor and hides the tract's. (On a norovirus capsid gap this is the difference
between recruiting only the flank reference and also recruiting the divergent core's
distinct donor.) Resolution has two levers: `--window-size` controls whether a short
region is *detected* as a gap at all (a wide window averages it away), and `--subtile`
controls whether its donor is *found* once detected. Set `--subtile 0` to disable;
very short intervals (< ~150 bp) lose BLAST specificity.

### Fill references repeatedly

To do this repeatedly until coverage stops improving, `tessera fill-references` runs
the whole cycle -- build MSA, scan, BLAST, download -- for several rounds, growing a
copy of the collection each time:

```
tessera fill-references --query my_query.fasta \
    --collection collection/ --output filled/ --aligner mafft --email you@example.org
```

It stops when the gaps close, when no new reference can be found, or when a round no
longer improves the worst gap (so a genuinely hypervariable region is reported, not
chased forever). Each round is recorded in `filled/fill_summary.tsv`, and the query's
own record is auto-excluded from its FASTA header. This needs an aligner and Entrez
Direct, and rebuilds the alignment every round.

Omit `--collection` to **start fresh** with no suggested references: the first round
seeds the collection from an NCBI search, then the loop fills the remaining gaps as
usual.

```
tessera fill-references --query CRF01_AE.fasta --output filled/ --aligner mafft \
    --seed-mode parents --curate --email you@example.org
```

## Seed modes: where the starting references come from

How the seed is chosen (`--seed-mode`) matters for a recombinant query. A whole-query
search returns the query's closest *whole-genome* relatives -- which, for a
recombinant, are its own-lineage **siblings**; they cover the whole query and mask the
recombination. The parents that donated each region are *regional* best matches, so by
default seeding searches per region:

- `whole` -- one whole-query search; the closest whole-genome relatives. Best for a
  non-recombinant query.
- `windowed` (default) -- search each window (`--seed-window`, default 1500 bp) and
  seed the per-window best hits. Safe for any query; surfaces regional matches.
- `parents` -- like windowed, but drop near-identical full-coverage hits (the
  siblings) so each region contributes its best *divergent* source. This recruits the
  parental lineages of a recombinant query; pair it with `--curate`. If a region has
  only siblings in NCBI, it falls back to the best hit so the panel is never empty.

`parents` mode helps when the parents are discoverable by BLAST -- the common case of
a newly observed recombinant whose own form is still rare. When the query's own
lineage saturates `nt` (for example HIV-1 CRF01_AE, with thousands of genomes), every
region's top hits are siblings and the parents never surface; the run reports this
("no distinct parental lineage found in NCBI") rather than masking it. BLAST ranks by
similarity, so it cannot escape an over-represented lineage.

## Seed source: seeding from a finite pool

For an over-represented lineage, seed from a **finite pool** instead
(`--seed-source`):

- `local` -- a directory of genomes you control (`--candidate-pool dir/`: a RefSeq
  dump, a subtype reference set, or a lab collection). Because the pool is finite it is
  dereplicated and the query's siblings (near-identical genome-wide) are dropped, so
  the divergent parents are no longer out-ranked. The query is then matched region by
  region and the best-covering genomes are kept -- a focused parental panel, offline
  and reproducible.
- `ncbi-virus` -- a taxon-scoped set fetched with the `datasets` CLI (`--taxon
  "HIV-1"`, auto-detected if omitted). By default it pulls the RefSeq representative
  set (one genome per lineage, diverse by construction); add `--source-complete` for
  all complete genomes, which are then dereplicated.
- `nextclade` -- a pool reconstructed from a [Nextclade
  dataset](https://docs.nextstrain.org/projects/nextclade/en/stable/user/datasets.html).
  Tessera auto-detects the dataset from the query (`nextclade sort` when the CLI is
  installed, otherwise BLAST taxon detection mapped to a dataset), or you pass it
  explicitly with `--nextclade-dataset <path>` (e.g. `nextstrain/sars-cov-2/XBB`,
  `community/neherlab/hiv-1/hxb2`). Every reference-tree tip is reconstructed from the
  dataset reference plus its mutations and labelled by clade, so the report names
  parents by clade. Fetched pools are cached per dataset version. Add `--pool-consensus`
  to compete **one denoised consensus genome per clade** (the majority base per column
  across the clade's tips) instead of every tip -- a stable per-lineage reference (the
  jpHMM principle) that is harder to lose to a near-identical adjacent-clade genome at low
  divergence. The consensus pool is cached separately. Use it for low-divergence,
  near-identical panels rather than as an always-on default: the harness attribution
  measurement (`validation/attribution-results.md`) found that one-genome-per-clade
  competition can *degrade* backbone (major-parent) attribution at moderate divergence
  (~7-9%), so it is opt-in.

```
tessera fill-references --query CRF01_AE.fasta --output filled/ --aligner mafft \
    --seed-source local --candidate-pool subtype_refs/ --curate
```

```
tessera detect --query CRF01_AE.fasta --output out/ --nextclade
tessera fill-references --query q.fasta --output out/ --seed-source nextclade \
    --nextclade-dataset nextstrain/sars-cov-2/XBB
```

A Nextclade dataset is a clade-typing reference tree, so the pool spans clade
representatives rather than clean recombinant-parental endpoints; it suits a broad
screen. Single-segment datasets (most flu segments) yield segment-length pool genomes,
which is logged. The `nextclade` CLI is optional (it only sharpens dataset
auto-detection).

The sibling drop uses an absolute identity cutoff, which suits divergent-lineage
recombinants (HIV subtypes ~12 % apart). For closely related parents that have no
masking sibling (e.g. SARS-CoV-2 sublineages), add `--seed-keep-siblings`. NCBI Virus
/ `datasets` needs `conda install -c conda-forge ncbi-datasets-cli`.

## Auto-diversify

By default (`--auto-diversify`), BLAST seeding will **switch to the `ncbi-virus`
diversity path automatically** when it finds only siblings -- i.e. when the query's
lineage saturates `nt` and no parental lineage can be recruited by similarity. A broad
fetch is capped (`--fetch-limit`, default 2000) and dereplicated; for a heavily
sequenced taxon the capped sample may miss lineages, so a curated `--candidate-pool`
is recommended (and the run says so). Whether the diversity panel actually contains the
parents depends on the taxon: it works when they are genotype/lineage representatives,
less so for fine genotype-specific recombinants. Disable with `--no-auto-diversify`.
This complements the caller-side defence: `tessera recomb` excludes whole-genome
siblings from the competition (`--exclude-siblings`, on by default), so a sibling that
still slips into the panel cannot win every window and mask the event.

## Curating a panel for detection

Auto-filling maximises *coverage* by recruiting the query's closest genomes -- but for
a recombinant query those closest genomes are its own-lineage relatives, which match it
everywhere and **mask** the recombination. (Filling an HIV-1 CRF01_AE query recruited
CRF01_AE relatives that hid the underlying subtype-A/E mosaic.) Detection needs
*diversity*: distinct parental lineages, and none of the query's siblings.

`tessera curate-panel` curates a collection for detection. It uses
[skani](https://github.com/bluenote-1577/skani) to measure each reference's
genome-wide identity (ANI) and how much of the query it covers, then drops the query's
siblings and dereplicates near-duplicates with
[skDER](https://github.com/raufs/skDER):

```
tessera curate-panel --query CRF01_AE.fasta --collection collection/ --output curated/
```

A *sibling* is a reference whose ANI to the query exceeds the backbone's by a margin
**and** that covers most of the query -- a whole-genome relative. A regional parent
(for example an env-only reference) has high ANI but low coverage and is kept. The rule
is relative to the backbone (the query's closest whole-genome match), so it needs no
per-organism identity cutoff: it flags HIV subtype relatives (~12 % apart) and
SARS-CoV-2 sublineages (<1 % apart) alike. The curated `curated/collection/` and a
`panel_lineages.tsv` (each reference's role and ANI/coverage) are written; rebuild with
`tessera msa` then `tessera recomb`.

The same curation runs inside the fill loop with `fill-references --curate`, which
keeps the growing panel diverse and sibling-free each round and adds a "Reference
panel" section to the report. Both need skani (and skDER for dereplication):
`conda install -c bioconda skani skder`.

## Make a collection lineage-ready (`type-lineages`)

`tessera type-lineages -c genomes/ -o out/` assigns a lineage to every genome and writes
`out/lineages.tsv`, so an NCBI or local set is typed like a Nextclade pool. Each genome takes
the best label available: a header/metadata designation, else the clade of its nearest
Nextclade-dataset reference tip (skani ANI, above `--ref-ani-floor`), else a de-novo lineage
from ANI clustering (`--cluster-ani`). Downstream `recomb`/`detect`/`fill-references` read the
`lineages.tsv` beside the MSA or output, so the typed names carry through. Needs skani; the
nearest-neighbour step fetches the Nextclade dataset on first use.

### Type the recruited panel automatically (`--deep-typing`)

`detect` and `fill-references` type their recruited panel at the end of a run. By default this
is header/metadata mining only. Pass `--deep-typing` to run the full `type-lineages` ladder on
the recruited references instead, so genomes that carry no clean header lineage still get a
Nextclade nearest-neighbour or de-novo label. It writes the same `lineages.tsv` and feeds the
lineage-aware panel selection. It is off by default because it adds a Nextclade/skani pass at
the end of the run; it needs skani (the command errors early if skani is absent). References
fetched from NCBI Virus are additionally tagged with the source `ncbi-datasets` -- their
structured datasets lineage, captured during download in an `ncbi_lineages.tsv` sidecar,
rather than a token mined back out of the header.
