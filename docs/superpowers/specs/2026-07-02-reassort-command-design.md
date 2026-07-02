# `tessera reassort`: per-segment reassortment detection -- design

## Context

Tessera detects **intragenic** recombination by anchoring a query to a single backbone genome and
scanning sliding windows. Segmented viruses (influenza, bunyaviruses, reoviruses) also undergo
**reassortment** -- a whole-segment swap between parents -- which the single-backbone model cannot
represent: the harder-harness `reassort_flu` probe (a concatenated `HA ++ NA` query) produced an
empty regions table because one backbone cannot anchor two segments (the NA half fell off). There
is no per-segment or reassortment concept anywhere in `src/`.

Reassortment is a **per-segment** phenomenon: for each segment, which parent did the query get it
from? A reassortant is a query whose segments trace to **different parents**. This is orthogonal to
(and simpler than) the intragenic scan -- each segment is analysed independently, sidestepping the
single-backbone problem. This spec adds a dedicated `tessera reassort` command that does exactly
that: **per-segment nearest-reference typing + a strain-consistency reassortment call**.

Decisions (from brainstorming): a **new `reassort` command**, **assignment-only** v1 (per-segment
intragenic recombination is a later cycle); references come from **Nextclade per-segment datasets**
(auto-typed); the reassortment call is by **nearest-reference-strain consistency**.

## Feasibility (confirmed)

Nextclade flu trees label tips by **strain name** (`A/Darwin/6/2021`), and the same strain is typed
in multiple segment datasets -- for H3N2 HA vs NA, **477 strains overlap by name**. So a query's
per-segment nearest tips can be compared by strain across datasets: segments tracing to the same
strain/lineage are clonal; segments tracing to disjoint strains/lineages are reassortment.

## Command

`tessera reassort` -- assign each segment of a segmented query to its nearest reference lineage and
report the per-segment genotype + a reassortment verdict.

Inputs (Typer, following `cmd_detect.py`): `-q/--query` (a **multi-FASTA**, one record per segment),
`-o/--output`, `--email` (NCBI, for dataset auto-detection), optional `--dataset seg=path` overrides
(map a segment record name to an explicit Nextclade dataset), `--ani-floor` (a segment below this
ANI to every tip is "unassigned"). Needs skani.

## Architecture

New module `src/tessera/reassort/` (keeps the segment logic out of the single-backbone `recomb`
path). Core = `assign_segments(query_records, ...) -> ReassortmentResult`:

For each query segment record:
1. **Resolve its dataset** -- an explicit `--dataset` override, else `resolve_dataset` on that
   segment's sequence (nextclade-sort / BLAST-alias), so HA resolves to the HA dataset, NA to NA.
   A segment whose dataset does not resolve is recorded `unassigned` (non-fatal).
2. **Reconstruct the clade + strain-labelled tips** -- `build_pool(dataset)` for the tip genome
   files, plus an `accession -> (strain, clade)` map read from the tree (`node["name"]` is the
   strain; `_accession_of` / `node_clade` already exist). Cache per `path@tag` as elsewhere.
3. **Nearest tip by ANI** -- `skani_query_ani(segment_fasta, tips)`; the highest-ANI tip that clears
   `--ani-floor` gives the segment's `(strain, clade, ANI)`. Below the floor -> `unassigned`.

Then the **reassortment call** (`_call_reassortment`):
- Collect each segment's nearest `(strain, clade)`.
- The query is **clonal** (not reassortant) when a single reference strain is the nearest, or a
  near-nearest (within a small ANI margin of the best), for **every** assigned segment -- i.e. one
  parent strain explains all segments. To test this, for each segment keep the top-k nearest strains
  (not just the top 1); the query is clonal iff the intersection of the per-segment top-k strain
  sets is non-empty.
- Otherwise it is a **reassortant**: no single strain explains all segments. Report each segment's
  nearest strain/clade as its independent origin, and group segments by donor strain/lineage.
- Segments left `unassigned` (no dataset / below the floor) are reported as such and excluded from
  the clonal/reassortant decision (with a note), never silently dropped.

## Output
`ReassortmentResult`: per-segment rows `(segment, dataset, nearest_strain, clade, ani, status)` and a
verdict (`clonal` / `reassortant` / `undetermined` when too few segments assign). Written as:
- `reassortment.tsv` -- the per-segment table.
- a text summary -- the genotype mosaic (e.g. `HA:3C.2a1b (A/X) | NA:N2.2 (A/Y) -> reassortant`) and
  the verdict, via the CLI.

## Data flow
```
multi-FASTA query -> per segment: resolve dataset -> build_pool tips + strain map
                                -> nearest tip by skani ANI -> (strain, clade, ANI)
   -> clonal iff one strain is near across all assigned segments; else reassortant
   -> reassortment.tsv + text mosaic
```

## Reuse (no reimplementation)
`resolve_dataset`, `build_pool`, `_accession_of`, `node_clade` (`discover/nextclade.py`);
`nextclade_cache` (`core/cache.py`); `skani_query_ani`, `skani_available` (`discover/panel.py`);
`read_fasta`, `write_fasta_record`, `strip_sequence_extension` (`core/io.py`); the `cmd_*.py` +
`main.py` registration pattern; `UserInputError`/`stage_errors`.

## Edge cases
- **Single-record query** (one "segment") -> `undetermined` (reassortment needs >= 2 segments).
- **A segment's dataset does not resolve** -> that segment is `unassigned`; the rest still assign.
- **A segment below `--ani-floor` to every tip** -> `unassigned` (a divergent segment with no close
  reference).
- **skani missing** -> a clear `UserInputError` up front (as `--curate` does).
- **All segments clonal to one strain** -> `clonal`, the honest common case (not every segmented
  query is a reassortant).

## Testing / verification
- **CI unit tests (no binaries):** `_call_reassortment` on synthetic per-segment
  `(strain, clade, ani)` inputs -- clonal when a strain is shared across all segments' top-k;
  reassortant when the top-k strain sets are disjoint; `undetermined` for < 2 assigned segments;
  `unassigned` segments excluded from the decision. The dataset-resolution / ANI steps are stubbed.
- **Opt-in end-to-end (needs skani/Nextclade, not CI):** a synthetic reassortant multi-FASTA
  (H3N2 HA from one clade + NA from another) -> `reassortant` with the two segments attributed to
  different strains/clades; a clonal control (HA + NA of the same strain) -> `clonal`.
- **Harness:** flip the `reassort_flu` frontier probe to call `tessera reassort` on the segmented
  query and assert `reassortant` (replacing the concatenated-junction approach that XFAILed);
  record in `attribution-results.md` that the reassortment known-limitation is closed by the new
  command.
- `ruff check src tests validation` clean; CI coverage floor (80%) holds.

## Scope / non-goals
- v1 = per-segment nearest-reference assignment + strain-consistency reassortment call, on a
  multi-FASTA query, via Nextclade per-segment datasets. No per-segment **intragenic** recombination
  scan (a later cycle). No change to the existing `detect`/`recomb` single-backbone path.
- Reassortment among **internal** genes beyond what the resolved datasets cover is out of scope
  (only segments whose datasets resolve are assigned).
- No new runtime dependency (skani/Nextclade tips already optional; native tip reconstruction).
