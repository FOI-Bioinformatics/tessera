# Lineage-aware panel selection + recombinant toggle -- design

## Context

The workflow audit and the measured attribution runs (`validation/attribution-results.md`)
showed that ANI dereplication is the wrong panel-reduction primitive for **typed** references:
it is a proxy for diversity that (a) drops a true parent lineage when it is near-identical to a
neighbour (`rsv_a`: `A.D.1.8` collapsed away) and (b) keeps recombinant CRF lineages that carry
both parents' segments and mask the true parents (`hiv1`). The clade-aware ANI patch that only
addressed (a) produced a **swap, not a gain** (16/18: `rsv_a` in, `hiv1` out) and was abandoned.

The fix is to select the panel **by lineage** when lineage labels exist: keep one representative
per lineage and **exclude recombinant lineages** (which are never valid parental endpoints for a
clean detection). This is the consumer of the just-merged `type-lineages` command (which
produces the `lineages.tsv` labels) and the piece that closes **both** failures. Per the user,
the recombinant exclusion is a **toggle** (default exclude), because sometimes the best-matching
genome is wanted even if it is a known hybrid. Activation is **default-on when a lineage map is
present**, gated by the `--compare` harness showing a clean 17/18.

## Goal

At `select_regional`, when references carry lineage labels, build the panel by lineage instead of
by ANI: one query-closest representative per **non-recombinant** lineage; recombinant (CRF/URF/X)
lineages excluded unless the toggle keeps them; untyped genomes fall back to today's ANI
dereplication. Backward compatible: with no lineage map, behaviour is exactly today's.

## Components

### `src/tessera/discover/pool.py` -- `select_regional`
Two new keyword params: `lineage_of: dict[str, str] | None = None` (genome label -> lineage) and
`keep_recombinant: bool = False`. Reorder the whole-genome query ANI (`skani_query_ani`) to run
*before* the reduction step so it is available to pick representatives. Replace the
`dereplicate_pool` block with a new `_lineage_select` when `lineage_of` is given; otherwise the
existing global `dereplicate` runs unchanged.

New helper (same module):
```
def _lineage_select(genomes, lineage_of, query_ani, *, keep_recombinant, derep_ani, logger)
    -> list[Path]
```
- Partition into typed (`lineage_of.get(strip_sequence_extension(g.name))` truthy) and untyped.
- Drop typed genomes whose lineage `is_recombinant_lineage`, unless `keep_recombinant`.
- Group the rest by lineage; keep the member with the highest `query_ani[g][0]` per lineage.
- Untyped genomes: `dereplicate(untyped, ani=derep_ani, logger=logger)` (the ANI fallback).
- Return the union; the caller's sibling-drop / related-filter / per-window steps are unchanged.

### `src/tessera/discover/iterate.py` -- product wiring
`_select_from` builds `lineage_of` from the pool genomes' headers (`titles_from_collection` +
`genotype_from_title`, as it already does elsewhere) and passes it plus
`keep_recombinant=params.keep_recombinant` to `select_regional`. `FillParams` gains
`keep_recombinant: bool = False`, forwarded through `FillParams.for_detection`.

### CLI -- the toggle
`--keep-recombinant-lineages/--no-keep-recombinant-lineages` (default off = exclude) on
`cli/cmd_detect.py` and `cli/cmd_fill_references.py`, threaded into `FillParams`. Help: "Keep
recombinant (CRF/URF/X) lineages in the panel (default: exclude them, since they mask the true
parents). Turn on to get the best-matching genome even if it is a known hybrid."

### `validation/run_hybrids.py` -- measure it
In `_build_and_score`'s tip branch, build `lineage_of` from the tips (reuse `clade_of_label`) and
pass it to `select_regional` (with `keep_recombinant=False`). This closes `rsv_a` (its
non-recombinant `A.D.1.8` lineage is kept) and `hiv1` (CRF lineages excluded).

## Reuse
`is_recombinant_lineage` (`recomb/typing.py`, on main); `skani_query_ani` and `dereplicate`
(`discover/panel.py`, the latter now the untyped fallback); `titles_from_collection` /
`genotype_from_title` (`recomb/typing.py`); `clade_of_label` (harness).

## Data flow

```
pool ─▶ query ANI (skani) ─▶  lineage_of present?
                               ├─ yes: drop recombinant lineages (unless keep), keep one
                               │       query-closest rep per lineage; untyped ─▶ ANI derep
                               └─ no:  global ANI dereplicate (today)
      ─▶ sibling-drop ─▶ related-filter ─▶ per-window selection (unchanged)
```

## Edge cases
- **No lineage map** -> current global path (strict superset).
- **`keep_recombinant=True`** -> recombinant lineages compete like any other (one rep each).
- **A lineage with one member** -> that member is the rep.
- **All-untyped panel** -> `_lineage_select` routes everything to the ANI fallback (== today).
- **skani/skDER absent** -> unchanged: `select_regional` still raises its existing "needs skani"
  error before reduction.

## Testing / verification
- **CI unit tests** (`tests/unit/test_lineage_select.py`, stub `skani_query_ani`/`dereplicate`):
  `_lineage_select` keeps one query-closest rep per non-recombinant lineage; excludes CRF/URF/X
  by default and keeps them when `keep_recombinant=True`; routes untyped genomes to the ANI
  fallback; is the fallback path when `lineage_of` is empty.
- **Harness gate** (`recomfi-aln` env, opt-in): `python validation/run_hybrids.py --compare` ->
  confirm `rsv_a` and `hiv1` both PASS and the other 16 hold (clean **17/18**); record in
  `validation/attribution-results.md`. Ship default-on only if clean.
- `ruff check src tests validation` clean; CI coverage floor (80%) holds.

## Scope / non-goals
- Only `select_regional`'s reduction step, the product/harness wiring, and the CLI toggle change.
  No change to detection, attribution, sibling exclusion, the per-window selection, or the report.
- The BLAST-recruited `fill-references` rounds (which do not go through `select_regional`) are out
  of scope; auto-running `type-lineages` inside recruitment is a separate follow-up.
- No new runtime dependency.
