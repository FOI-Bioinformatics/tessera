# Auto-run the lineage-assignment ladder in recruitment (+ ncbi-datasets source tag) -- design

## Context

`type-lineages` (`discover/lineage_assign.py`) assigns lineages through a three-tier ladder:
tier-1 metadata/header mining, tier-2 Nextclade nearest-neighbour, tier-3 de-novo ANI
clustering. But that ladder is reachable only from the standalone `type-lineages` command.
Inside `fill-references` / `detect` recruitment, `_type_panel` (`discover/iterate.py:361`) types
the recruited references with **tier-1 only** (`build_lineage_map`: user TSV + mined title
token). So freshly downloaded NCBI genomes that a Nextclade dataset could place, or that share
de-novo clusters, are left with only whatever a title token yields.

This design lets recruitment optionally run the **full ladder** once at the end, enriching
`lineages.tsv`, the report, the Pango cross-check, and the lineage-aware panel selection that a
later `tessera recomb` consumes. A second, smaller component fixes a provenance gap: the
structured NCBI-datasets lineage is currently mined back out of the FASTA header note and
tagged `title`, never `ncbi-datasets`.

Two follow-up items from prior cycles, taken together:
- **A (#3):** auto-run the assignment ladder in recruitment.
- **B (#4):** emit the `ncbi-datasets` lineage source as its own tag.

## Decisions (settled in brainstorming)

- **Placement:** run the ladder **once at the end** of recruitment (replace `_type_panel`'s
  tier-1-only call), not during the rounds. Cheapest and simplest; the reduction inside the
  rounds keeps its header-mined labels.
- **Activation:** an **opt-in flag** `--deep-typing/--no-deep-typing`, **default off on both**
  `detect` and `fill-references`. Off reproduces today's tier-1 typing exactly.
- **Part B scope:** both pieces -- the `assign_lineages` pass-through param *and* the recruitment
  capture -- so recruited NCBI-Virus genomes are correctly tagged `ncbi-datasets` end to end.

## Part A -- the ladder in recruitment

### `FillParams` and CLI
- `FillParams.deep_typing: bool = False`, forwarded through `FillParams.for_detection`.
- `--deep-typing/--no-deep-typing` (default off) on `cli/cmd_detect.py` and
  `cli/cmd_fill_references.py`, threaded into `FillParams`. Help: "Type the recruited panel with
  the full lineage ladder (Nextclade nearest-neighbour + de-novo ANI clustering) instead of just
  header mining. Needs skani; off by default."

### `fill_references` guard
Mirror the existing `--curate` check: if `params.deep_typing and not skani_available()`, raise a
`UserInputError` up front (the ladder's tier-2/tier-3 need skani; failing early is clearer than a
mid-run crash).

### `_type_panel` branch (`discover/iterate.py:361`)
When `deep_typing`, replace the tier-1 `build_lineage_map` call with the full ladder:

```
rows = assign_lineages(
    coll_files,
    user_lineage_map=params.lineage_map,
    taxon=params.taxon,
    nextclade_dataset=params.nextclade_dataset,
    datasets_rows=_read_ncbi_lineages(params.output),   # Part B; None when absent
    email=params.email,
    cache_dir=params.cache_dir,
    logger=logger,
)
```

The query row (`genotype_from_title` + `("query", ...)`) and `write_lineage_map` are appended
exactly as today, so the `lineages.tsv` format and every downstream consumer are unchanged --
only the labels get richer. When `deep_typing` is off, `_type_panel` is byte-for-byte today's
tier-1 path.

### Data flow / error handling
Unchanged pipeline; `_type_panel` swaps its one typing call for the ladder when the flag is on.
The ladder's tier-2 already treats "no Nextclade dataset resolves" as a normal non-fatal skip
(falls through to de-novo), so the only new hard failure is the opt-in-triggered missing-skani
guard.

## Part B -- ncbi-datasets source tag

### Piece 1: `assign_lineages` pass-through
Add `datasets_rows: Iterable[tuple[str, str]] | None = None` to `assign_lineages`
(`discover/lineage_assign.py`), forwarded to `build_lineage_map(datasets_rows=...)`.
`build_lineage_map`'s existing priority (user > ncbi-datasets > title) means a structured row
overrides a mined title token for the same accession and is tagged `ncbi-datasets`.

### Piece 2: capture in recruitment
- **Write at fetch:** `fetch_ncbi_virus` (`discover/pool.py`) already computes
  `labels = _lineage_labels(data_report.jsonl)`. It writes those as a 2-column sidecar
  `ncbi_lineages.tsv` (`accession<TAB>lineage`) in its `dest` before the jsonl is discarded.
- **Consolidate:** a helper `_record_ncbi_lineages(output, dest)` merges `dest/ncbi_lineages.tsv`
  into one well-known file `<output>/ncbi_lineages.tsv` (dedup by accession), called from the
  recruitment paths that pull NCBI-Virus genomes into the collection (the seed branch in
  `_seed_collection` and the pool cache path `_cached_ncbi_panel`).
- **Read at typing:** `_read_ncbi_lineages(output)` reads `<output>/ncbi_lineages.tsv` into a list
  of `(accession, lineage)` rows for `assign_lineages(datasets_rows=...)`; returns `None` when the
  file is absent.

Genomes from BLAST or a local `--candidate-pool` contribute no rows and keep title/de-novo
typing. When `deep_typing` is off, none of Part B runs.

### Data flow
```
fetch_ncbi_virus -> dest/ncbi_lineages.tsv
   -> _record_ncbi_lineages merges to <output>/ncbi_lineages.tsv
      -> _type_panel (deep_typing) reads it -> assign_lineages(datasets_rows=...)
         -> build_lineage_map tags those accessions 'ncbi-datasets' in lineages.tsv
```

## Reuse
`assign_lineages`, `assign_by_nearest_reference`, `cluster_by_ani` (`discover/lineage_assign.py`);
`build_lineage_map`, `write_lineage_map`, `genotype_from_title`, `titles_from_collection`,
`LINEAGES_TSV` (`recomb/typing.py`); `skani_available` (`discover/panel.py`); `_lineage_labels`,
`fetch_ncbi_virus` (`discover/pool.py`); the `cmd_*.py` typer-option pattern and `FillParams`
plumbing used by `--curate` / `--reattribute-donors`.

## Edge cases
- **`deep_typing` off** -> `_type_panel` is exactly today's tier-1; neither part runs.
- **No NCBI-Virus fetch** (local pool / pure BLAST) -> no sidecar; Piece 1 inert, no error.
- **No Nextclade dataset for the taxon** -> tier-2 skips, tier-3 de-novo labels the rest.
- **skani missing with `--deep-typing`** -> hard `UserInputError` up front (like `--curate`).
- **Sidecar row for a genome later curated out** -> unused, harmless.
- **Stale `<output>/ncbi_lineages.tsv` from a prior run in the same output dir** -> recruitment
  clears/overwrites the collection each run; the sidecar is written under the same fresh `output`,
  so it is rebuilt per run (merge is additive within one run only).

## Testing / verification
- **CI unit tests (no binaries):**
  - `assign_lineages` with tool calls stubbed: a supplied `datasets_rows` entry wins over a mined
    title token for the same accession and carries source `ncbi-datasets`.
  - Sidecar round-trip: `fetch_ncbi_virus` (fetch stubbed) writes `dest/ncbi_lineages.tsv`;
    `_record_ncbi_lineages` merges/dedups into `<output>/ncbi_lineages.tsv`; `_read_ncbi_lineages`
    reads it back; absent-file returns `None`.
  - `_type_panel` on the `deep_typing` path (with `assign_lineages` stubbed) passes the read
    `datasets_rows` through and still writes a well-formed `lineages.tsv` including the query row;
    `deep_typing` off calls the tier-1 path unchanged.
  - CLI: `--deep-typing/--no-deep-typing` default off on both commands and threaded to
    `FillParams` (via `for_detection` for `detect`).
- **Opt-in end-to-end (needs skani/Entrez, not CI):** `tessera detect ... --deep-typing` on a
  small taxon yields a `lineages.tsv` with `nextclade-nn` and/or `ncbi-datasets` sources present.
- `ruff check src tests` clean; CI coverage floor (80%) holds.

## Scope / non-goals
- Typing runs **once at the end** of recruitment; the per-round `select_regional` reduction is
  unchanged (still header-mined labels).
- No change to detection, attribution, the report layout, or the panel-selection algorithm.
- No new runtime dependency (skani/Nextclade tips are already optional; native tip
  reconstruction, no `nextclade` binary).
- Retyping already-typed user collections or re-running the ladder mid-round is out of scope.
