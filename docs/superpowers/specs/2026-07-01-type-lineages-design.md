# `tessera type-lineages`: create lineage-ready datasets from NCBI/local data -- design

## Context

Tessera's detection and panel building work best on **lineage-labelled** references (a
Nextclade pool gives every genome a clade). But when references come from NCBI or a local
directory they usually carry **no clean lineage**, so the pipeline falls back to ANI
dereplication -- a proxy for diversity that (as the workflow audit and the `rsv_a`/`hiv1`
measurements showed) drops true parents and keeps masking recombinants. The fix is upstream:
a dedicated command that **assigns lineages to an arbitrary genome collection**, making it
"lineage-ready" the way a Nextclade dataset is, so everything downstream (panel selection,
attribution, the report) can use real labels instead of ANI guesses.

This is modelled on how Nextstrain does it -- Nextclade assigns a clade by placing each query
at its **nearest neighbour on a clade-labelled reference tree** -- but implemented natively
(no `nextclade` binary), and extended with a de-novo fallback for taxa that have no dataset.
Dereplication is part of it: it reduces the reference tips before matching, and clustering is
the de-novo mechanism.

Scope this cycle: **the command only** (produce `lineages.tsv`). The consumer change --
lineage-aware panel selection with a **toggle to keep/exclude recombinant lineages** (so a
user can ask for the best genome even when it is a known hybrid) -- is the immediate
follow-up cycle (it consumes this command's output).

## The command

`tessera type-lineages` -- assign a lineage/genotype to every genome in a collection and write
the shared `lineages.tsv` sidecar (`label<TAB>genotype<TAB>source`, the existing format).

Inputs (Typer, following `cmd_build_panel.py`): `-c/--collection` (genome dir), `-o/--output`,
optional `--taxon` / `--nextclade-dataset` (for the reference step; auto-detected otherwise),
`--email`, and thresholds (`--ref-ani-floor`, `--cluster-ani`).

**Assignment ladder** -- each genome takes the best label available, most-authoritative first:

1. **Metadata / header** -- reuse `build_lineage_map` (`recomb/typing.py:157`): a user TSV, then
   the NCBI-datasets lineage rows (`discover/pool.py:_lineage_labels`, which reads
   `isolate.lineage` / `pangolinClassification`), then a mined title token
   (`genotype_from_title`). Source `user` / `ncbi-datasets` / `title`. Cheapest, no tools.
2. **Reference nearest-neighbour** (Nextclade-style, native) -- for genomes still untyped, if a
   Nextclade dataset resolves (`resolve_dataset`), reconstruct its clade-labelled tips
   (`build_pool`), **dereplicate the tips per clade** (`dereplicate(..., clade_of=...)`, reuse
   the clade-aware helper), then assign each untyped genome the clade of its **nearest tip by
   `skani_query_ani`** when the top ANI clears `--ref-ani-floor`. Source `nextclade-nn`.
3. **De-novo clustering** -- for genomes still untyped (no dataset / below the floor): cluster
   the remaining collection all-vs-all and label each cluster a de-novo lineage
   (`denovo_1`, ...). Source `denovo`.

Recombinant lineages are **flagged** in the output (a shared `is_recombinant_lineage`) but not
excluded here -- exclusion is the consumer's toggle (next cycle).

## Components

### `src/tessera/discover/lineage_assign.py` (new -- the substance)
- `assign_by_nearest_reference(genomes, labeled_tips, *, ani_floor, logger) -> dict[str, str]`
  -- per genome, `skani_query_ani` vs the (deduped, clade-labelled) tips; return the nearest
  tip's clade when ANI >= `ani_floor`.
- `skani_triangle(genomes, logger) -> list[tuple[Path, Path, float]]` -- new thin wrapper over
  `skani triangle` (all-vs-all ANI; the current wrapper is query-vs-refs only). Reuses the
  `run_tool(SKANI, ...)` + staging pattern from `panel.py`.
- `cluster_by_ani(genomes, edges, *, threshold) -> dict[str, list[Path]]` -- **pure**:
  connected components of the graph of pairs with ANI >= `threshold`; each component is a
  de-novo lineage. Unit-testable without binaries.
- `assign_lineages(params, logger) -> list[tuple[str, str, str]]` -- orchestrates the ladder and
  returns lineage-map rows.

### `src/tessera/recomb/typing.py`
- Add `is_recombinant_lineage(name) -> bool` (move the regex verbatim from
  `validation/run_hybrids.py:187-191`; also treat Pango `X*` via the existing `expand_recombinant`
  sense). Shared single source of truth; the harness imports it.

### `src/tessera/cli/cmd_type_lineages.py` (new)
- `@app.command(name="type-lineages")` wrapping `assign_lineages` + `write_lineage_map`; register
  by importing it in `cli/main.py` (the side-effect import block, `:107-115`).

### Reuse (no reimplementation)
`build_lineage_map`, `titles_from_collection`, `genotype_from_title`, `write_lineage_map`,
`LINEAGES_TSV` (`recomb/typing.py`); `resolve_dataset`, `build_pool`, `collect_tips`
(`discover/nextclade.py`); `skani_query_ani`, `dereplicate`, `SKANI`/`run_tool`
(`discover/panel.py`); `_lineage_labels` (`discover/pool.py`); `expand_recombinant`,
`load_alias_key` (`recomb/pango.py`); the `cmd_*.py` + `main.py` registration pattern.

## Data flow

```
collection -> (1) metadata/header  \
              (2) nearest-ref clade  >-- per-genome label + source --> lineages.tsv
              (3) de-novo cluster   /         (is_recombinant flag)
```

Downstream (detect / recomb / fill-references) already read `lineages.tsv` via
`load_lineage_map`, so a typed collection is immediately usable -- no other change needed for
the command to deliver value.

## Testing / verification
- **CI unit tests** (no binaries): `cluster_by_ani` (connected components, threshold edges,
  singletons); `assign_by_nearest_reference` with a stubbed `skani_query_ani` (nearest tip above
  floor -> its clade; below floor -> unassigned); `is_recombinant_lineage` (CRF01_AE, URF, XBB,
  a plain clade, "recombinant"); the ladder precedence (metadata beats nearest-ref beats
  de-novo) with the tool calls stubbed.
- **Opt-in end-to-end** (needs `skani`; not CI): `tessera type-lineages -c genomes/ -o out/`
  produces a `lineages.tsv` whose labels match known lineages on a small curated set; then
  `tessera recomb`/`detect` consumes it and names parents by lineage.
- `ruff check src tests` clean; CI coverage floor (80%) holds.

## Follow-up cycles (not in this plan)
- **Lineage-aware panel selection** consuming this output: one representative per lineage, with
  a `--keep-recombinant-lineages/--no-...` **toggle** (default exclude for clean detection; keep
  when the user wants the best genome even if a known hybrid). This is what closes `rsv_a` and
  `hiv1` on the harness; it is a separate, smaller cycle.
- Tier-2 wiring so `fill-references`/`detect` optionally run `type-lineages` on freshly
  downloaded NCBI genomes automatically.

## Branch disposition
Started fresh from `main`. The superseded `clade-aware-derep-spec` branch (the ANI patch) is
abandoned; its `rsv_a`/`hiv1` measurement lives in that branch's history and the SDD ledger.

## Scope / non-goals
- The `type-lineages` command only (assign + write `lineages.tsv`).
- No change this cycle to panel selection, the recombinant toggle, detection, or the report.
- No new runtime dependency (skani/skDER already optional external tools); native reference
  placement (no `nextclade` binary required).
