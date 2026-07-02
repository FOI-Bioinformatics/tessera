# `tessera reassort` v2: intragenic per-segment scan (`--scan-segments`) -- design

## Context

`tessera reassort` (v1.1) is assignment-only: it types each query segment against its own Nextclade
dataset and calls whole-segment reassortment (which parent each segment came from). It does not look
*inside* a segment. A segment can itself be a recombinant -- a crossover between two lineages within
that one segment -- which is the ordinary intragenic recombination the `detect`/`recomb` single-
backbone scan already finds. A full picture of a segmented virus needs both: reassortment *across*
segments and recombination *within* each segment.

This spec adds an opt-in `--scan-segments` flag to `reassort` that runs the existing `run_recomb`
pipeline on each assigned segment against that segment's own Nextclade panel. Decisions (from
brainstorming): a **flag on `reassort`** (not a separate command); **the full `run_recomb` report per
segment**; a **missing aligner is an up-front error** (not a silent no-op); the per-segment panel is
the **per-clade consensus** (as `detect` uses).

## Command

`tessera reassort -q segments.fasta -o out/ --scan-segments` -- assign each segment (as today) and,
for each *assigned* segment, additionally scan it for intragenic recombination against its Nextclade
panel. Existing options (`-q/-o`, `--email`, `--dataset`, `--ani-floor`, `--margin`) are unchanged.
New options: `--scan-segments` (bool, default off), `--aligner` (default `mafft`, validated against
the aligner registry, only used when scanning).

## Architecture

New module **`src/tessera/reassort/scan.py`** -- keeps the MSA/recomb machinery out of the
assignment path (`assign.py` stays skani-only when `--scan-segments` is off).

```
@dataclass
class SegmentScan:
    segment: str
    scanned: bool          # False if it could not be scanned (single-clade panel, MSA/recomb error)
    recombinant: bool      # True if run_recomb called at least one region
    n_regions: int
    note: str              # human-readable status, e.g. "2 regions", "single-clade panel", error text

def scan_segment(
    segment: str, seq: str, dataset, out_dir: Path, *,
    aligner: str, cache_dir: Path | None, logger: logging.Logger,
) -> SegmentScan
```

`scan_segment`:
1. Build the segment's panel: `build_pool(dataset, cache_dir=<consensus cache>, per_clade_consensus=True)`.
   Fewer than two clades -> `SegmentScan(scanned=False, note="single-clade panel")` (nothing to
   recombine between).
2. Write the segment query FASTA and each panel consensus into `out_dir/<segment>/` (sanitised
   segment name for the directory).
3. `build_msa(MsaParams(query, collection, output=<segment>/panel.msa.fasta, aligner=aligner))`.
4. `run_recomb(RecombParams(msa, output=out_dir/<segment>/, query=<segment>, organism=<segment>,
   window_size, window_step, methods=DEFAULT_METHODS))` -- the full pipeline, writing
   `recombination_regions.tsv`, the profile, the report, and the plot into the per-segment directory.
5. Parse `recombination_regions.tsv`; `n_regions` counts present (non-`donor_absent`) regions;
   `recombinant = n_regions > 0`.

Window sizing: port the harness's length-aware `window_params(genome_len)` into `scan.py` (segments
are ~1-2.5 kb, so the `detect` default 1000/100 window would give too few windows). It returns
`(window, step, selection_window)` scaled to the segment length.

**Integration in `assign.py`.** `assign_segments` gains `scan_segments: bool = False` and
`aligner: str = "mafft"`. When `scan_segments` is on:
- `_type_segment` additionally returns the resolved `dataset` (or `None`) so the scan reuses it --
  no second `resolve_dataset` (which for auto mode would re-run `nextclade sort`).
- After the assignment loop, for each assigned segment, call `scan_segment(...)` with the reused
  dataset, collecting `SegmentScan`s onto `ReassortmentResult.scans: list[SegmentScan]`.
- Unassigned segments are not scanned (no panel); they get `SegmentScan(scanned=False,
  note="unassigned")` so every input segment has a row.

`ReassortmentResult` gains `scans: list[SegmentScan] = field(default_factory=list)` (empty when
`--scan-segments` is off, so existing consumers and the harness are unaffected).

## New requirement and error handling

- **Aligner (new for this flag only).** The scan needs an aligner (mafft/minimap2) for the MSA;
  assignment still needs only skani. The CLI validates `--aligner` against
  `aligner_registry.names()` and calls the adapter's `preflight()` up front; if the aligner is
  unknown or its binary is absent, raise `UserInputError` **before** any typing (mirroring the
  `skani_available()` check). This is the "error up front" decision.
- **Per-segment scan failure is non-fatal.** Any exception inside `scan_segment` (MSA failure,
  recomb error, empty panel) is caught, logged, and recorded as `SegmentScan(scanned=False, note=<why>)`;
  the reassortment verdict and the other segments' scans are still produced.
- **Cache separation.** The consensus panel uses a cache directory distinct from the tip pool that
  assignment already cached under the same dataset key (append a `-consensus` suffix), so the two
  `build_pool` modes cannot collide.

## Output

- Unchanged: `reassortment.tsv`, `constellation.tsv`, the verdict + constellation mosaic.
- New per assigned segment: `out/<segment>/` with the full `run_recomb` outputs (regions TSV,
  profile, report, plot).
- New rollup `out/segment_scan.tsv`: columns `segment`, `intragenic_recombination` (`yes`/`no`/`n/a`),
  `n_regions`, `note`. Plus a summary log line, e.g.
  `HA: intragenic recombination (2 regions) | NA: none | PB2: not scanned (single-clade panel)`.

## Data flow

```
reassort --scan-segments:
  assign each segment (as today), keeping the resolved dataset for assigned ones
  per assigned segment:
    build_pool(per_clade_consensus=True) -> per-clade panel   (>=2 clades, else skip)
    build_msa(segment query + panel)     -> out/<segment>/panel.msa.fasta
    run_recomb(...)                      -> out/<segment>/{regions, profile, report, plot}
    parse regions -> SegmentScan
  write segment_scan.tsv + rollup summary; the reassortment verdict is unchanged
```

## Testing / verification

CI unit tests (no binaries):
- Rollup: a synthetic `recombination_regions.tsv` -> `SegmentScan` (regions counted, `recombinant`
  set, `donor_absent` rows excluded); an empty regions file -> `recombinant=False`.
- Stubbed `scan_segment` orchestration: monkeypatch `build_pool`, `build_msa`, `run_recomb`, and the
  regions parse to assert the per-segment directory layout, the single-clade-panel skip, and that a
  `build_msa`/`run_recomb` exception yields `scanned=False` (non-fatal), not a raised error.
- Stubbed `assign_segments(scan_segments=True)`: assigned segments get a `SegmentScan`, unassigned
  ones get `scanned=False, note="unassigned"`, and the reassortment verdict is unchanged from the
  no-scan path.
- CLI: `--scan-segments` and `--aligner` are exposed (parameter introspection); an unknown
  `--aligner` raises before typing.

Opt-in end-to-end (aligner env, not CI):
- A synthetic intragenic-recombinant flu HA segment (a within-HA mosaic of two clades) paired with a
  clonal NA, run through `reassort --scan-segments`: the HA `SegmentScan` shows recombination, NA
  shows none. Record the measured result in `validation/attribution-results.md`; state it faithfully
  (a scan that misses is a real finding, not retuned away).

`ruff check src tests validation` clean; CI coverage floor (80%) holds.

## Scope / non-goals

- v2 adds per-segment intragenic **detection** by reusing `run_recomb`. It does **not** change the
  reassortment call, and it does **not** attempt cross-segment donor attribution. A segment that is
  both reassorted and internally recombinant is reported as both, without a combined model.
- No change to `detect`/`recomb`. The scan reuses them unchanged.
- No new runtime dependency in the package: the aligner is the same external, optional backend
  `detect` already uses; it becomes required only when `--scan-segments` is passed.
- Deferred audit items #9 (example-tip clade label) and #10 (transient-network vs unassigned) remain
  documented limitations, untouched here.
