# Nextclade as a recomfi seed-source

Date: 2026-06-28
Status: approved design, ready for implementation
Branch: typed-lineages

## Goal

Make it possible to run recomfi recombination detection against a reference pool
derived from any Nextclade dataset, with the dataset auto-detected from the query.
This generalises the one-off HIV-1 pool experiment (`reports/hiv_nextclade_pool/`)
into a first-class `--seed-source nextclade` inside recomfi, usable across the
~79 official `nextstrain/*` datasets and the community datasets (SARS-CoV-2, flu
A/B segments, RSV, mpox, measles, dengue, mumps, rubella, ebola, hantavirus, WNV,
yellow fever, hMPV, VZV, HIV-1, …).

A Nextclade dataset is a clade-typing reference tree, not a set of clean
recombinant-parental endpoints; see `memory/hiv-nextclade-pool-finding.md` for the
documented limitation (a Nextclade-only pool could not reproduce the CRF01_AE A/E
mosaic). The feature makes the workflow reproducible; it does not claim the pool is
ideal for every recombination question.

## Non-goals

- Installing or wrapping the full Nextclade clade-assignment workflow. The
  `nextclade` CLI is an *optional* aid to dataset detection only.
- A dataset-listing/browse command (YAGNI; auto-detect plus an explicit override
  cover the need).
- Reconstructing insertions relative to the dataset reference (dropped; see
  Limitations).

## Architecture

The Nextclade pool reuses the exact path the existing `ncbi-virus` seed-source
uses: fetch genomes -> `select_regional` -> copy into the collection. The caller,
sibling-drop, dereplication, MSA and HMM stages are unchanged. Only the
"where do the candidate genomes come from" step is new.

```
query ──► resolve_dataset ──► build_pool ──► select_regional ──► collection/ ──► msa ──► recomb
          (nextclade.py)      (nextclade.py)  (pool.py, reused)    (iterate.py)
```

### New module: `src/recomfi/discover/nextclade.py`

Self-contained; `iterate.py` only calls its two public functions. No concrete
Nextclade logic leaks into the rest of the package.

Dataclass:

```python
@dataclass
class NextcladeDataset:
    path: str        # e.g. "community/neherlab/hiv-1/hxb2"
    tag: str         # newest version tag, e.g. "2026-04-14--11-55-23Z"
    files: dict      # file-role -> filename, from the index entry
```

Public functions:

- `resolve_dataset(query: Path, override: str | None, *, email, logger) -> NextcladeDataset`
  Hybrid detection, in priority order:
  1. `override` given (`--nextclade-dataset PATH`): validate the path exists in the
     index; resolve its newest tag. Raise `UserInputError` if unknown.
  2. else if the `nextclade` binary is on PATH: run `nextclade sort` on the query
     and parse the best-matching dataset path from its TSV output (handles
     multi-segment pathogens correctly). Resolve that path against the index.
  3. else: call the existing `pool.detect_taxon(query)` (BLAST top-hit organism)
     and map the organism string to a dataset path via `_ALIAS` (keyword -> path).
     Raise `UserInputError` with guidance if no alias matches (suggest passing
     `--nextclade-dataset`).

- `build_pool(dataset: NextcladeDataset, *, cache_dir, logger) -> list[Path]`
  1. If the cache dir for `path@tag` already holds genomes, return them (no
     network, no reconstruction).
  2. Otherwise download `reference.fasta`, `tree.json`, and the examples file
     (`sequences.fasta`/`example_sequences.fasta`, per the index `files` entry).
  3. Reconstruct every tree tip: walk root->tip accumulating
     `branch_attrs.mutations.nuc`, apply substitutions and deletions to the
     reference (1-based positions), strip gaps -> one ungapped genome.
  4. Add the example sequences (already full FASTA).
  5. Write each genome as `>{accession} {clade}` into the cache dir; return paths.
  6. Guard: skip a reconstructed sequence shorter than a small floor (truncated).

Helpers (module-private):

- `fetch_index() -> dict` — GET `https://data.clades.nextstrain.org/v3/index.json`
  via `urllib.request.urlopen` (same pattern as `recomb/pango.py`), parse JSON,
  flatten the nested `collections`/`datasets` tree to `path -> entry`. Newest tag
  is `entry["versions"][0]["tag"]`.
- `_dataset_file_url(dataset, role) -> str` —
  `https://data.clades.nextstrain.org/v3/{path}/{tag}/{filename}`.
- `_clade_of(node_attrs) -> str` — priority list
  `["Nextclade_pango", "clade_nextstrain", "clade_membership", "subclade",
  "LANL_subtype"]`, each value possibly `{"value": ...}` or a bare scalar; fall
  back to the universal `clade_membership`; else `"NA"`.
- `_accession_of(node) -> str` — `node_attrs.accession` (scalar or `{"value":}`)
  else `node["name"]`; sanitised to `[\w.]`.
- `_apply_mutations(ref_chars, nuc_muts)` — regex `^([A-Za-z])(\d+)([A-Za-z-])$`;
  out-of-range positions ignored.
- `_ALIAS: dict[str, str]` — keyword (lowercased organism substring) -> dataset
  path, covering the core pathogens. Examples:
  `"human immunodeficiency virus" -> "community/neherlab/hiv-1/hxb2"`,
  `"severe acute respiratory syndrome coronavirus 2" -> "nextstrain/sars-cov-2/XBB"`,
  `"monkeypox"/"mpox" -> "nextstrain/mpox/all-clades"`,
  `"measles" -> "nextstrain/measles/genome/WHO-2012"`,
  `"respiratory syncytial virus" -> "nextstrain/rsv/a/EPI_ISL_412866"`,
  `"dengue" -> "nextstrain/dengue/all"`, `"mumps" -> "nextstrain/mumps/genome"`,
  `"rubella" -> "nextstrain/rubella/genome"`, `"west nile" -> "nextstrain/wnv/all-lineages"`,
  `"yellow fever" -> "nextstrain/yellow-fever/prM-E"`,
  `"ebola"/"zaire ebolavirus" -> "nextstrain/orthoebolavirus/ebov"`,
  `"metapneumovirus" -> "nextstrain/hmpv/all-clades/NC_039199"`.

Tool capability for the optional binary:

```python
NEXTCLADE = ToolCapabilities(
    name="nextclade", conda=("bioconda::nextclade",),
    required_binaries=(BinarySpec("nextclade", version_args=("--version",)),),
)
```
`nextclade` is checked with `shutil.which`, never required; absence simply routes
detection to the BLAST fallback.

### Cache: `src/recomfi/core/cache.py`

Add, mirroring `ncbi_virus_cache`:

```python
def nextclade_cache(path: str, tag: str, *, override=None) -> Path:
    key = hashlib.sha1(f"{path}@{tag}".encode()).hexdigest()[:12]
    return cache_root(override) / "nextclade" / f"{_slug(path)}_{key}"
```

### Wiring: `src/recomfi/discover/iterate.py`

- `FillParams` gains `nextclade_dataset: str | None = None` (override; `None` =
  auto-detect).
- `_seed_from_pool` gains a `nextclade` branch (alongside `local` / NCBI):
  ```python
  if params.seed_source == "nextclade" and not force_ncbi:
      genomes = _fetch_nextclade(params, logger)
  ```
- New `_fetch_nextclade(params, logger) -> list[Path]`:
  ```python
  from .nextclade import build_pool, resolve_dataset
  ds = resolve_dataset(params.query, params.nextclade_dataset,
                       email=params.email, logger=logger)
  cache = nextclade_cache(ds.path, ds.tag, override=params.cache_dir)
  return build_pool(ds, cache_dir=cache, logger=logger)
  ```
  Selection then runs through the existing `_select_from` (`select_regional`) and
  `_copy_into`, exactly like the `ncbi-virus` path.

`_seed_collection` already routes `seed_source in ("local", "ncbi-virus")` to
`_seed_from_pool`; extend that membership test to include `"nextclade"`.

### CLI

- `cli/cmd_fill_references.py`:
  - add `nextclade` to the `--seed-source` `_require_choice` set.
  - add `--nextclade-dataset TEXT` option -> `FillParams.nextclade_dataset`.
- `cli/cmd_detect.py` and `cli/cmd_build_panel.py`:
  - add `--nextclade/--no-nextclade` (auto-detect) and `--nextclade-dataset TEXT`
    (override). When either is set, the preset uses `seed_source="nextclade"` and
    passes `nextclade_dataset` through; `candidate_pool` still wins if also given
    (document precedence: explicit local pool > nextclade > blast).

Headline usage:

```
recomfi detect --query CRF01_AE.fasta --output out/ --nextclade
recomfi detect --query q.fasta --output out/ --nextclade-dataset nextstrain/sars-cov-2/XBB
recomfi fill-references --query q.fasta --output out/ --seed-source nextclade
```

## Data flow

1. `detect`/`build-panel`/`fill-references` build `FillParams` with
   `seed_source="nextclade"` (+ optional `nextclade_dataset`).
2. `fill_references` -> `_seed_collection` -> `_seed_from_pool` (nextclade branch)
   -> `_fetch_nextclade`.
3. `resolve_dataset` picks the dataset (override / `nextclade sort` / BLAST+alias).
4. `build_pool` returns cached or freshly-reconstructed genome paths.
5. `select_regional` dereplicates, drops the query's siblings, and keeps the
   per-window best donors (unchanged).
6. The focused panel is copied into the collection; the normal MSA + HMM detection
   and typed-lineage report follow (clade names come from the genome headers).

## Error handling

- Unknown override path / no alias match / index unreachable -> `UserInputError`
  with an actionable message (suggest `--nextclade-dataset`, list a few examples).
- A network failure mid-download raises `UserInputError`; a partial cache dir is
  not treated as complete (write to a temp dir, move into place on success — same
  discipline as the NCBI Virus fetch).
- `nextclade sort` failure (non-zero exit / unparseable TSV) is non-fatal: log and
  fall through to the BLAST+alias path.
- A tree with no usable tips (e.g. a tree-less dataset) -> `UserInputError`
  pointing at the examples-only fallback or a different dataset.

## Testing (binary-free, mocked; matches the existing suite)

`tests/unit/test_nextclade.py`:
- Reconstruction: tiny `tree.json` + reference -> a tip's substitutions/deletions
  applied, gaps stripped, header `>{accession} {clade}`.
- Clade labeling: `Nextclade_pango` beats `clade_membership`; universal
  `clade_membership` fallback; `"NA"` when none.
- Dataset resolution: override validated against a mocked index; `nextclade sort`
  parse via mocked `run_tool`; BLAST->alias fallback via mocked `detect_taxon`;
  unknown organism raises `UserInputError`.
- Cache reuse: a populated cache dir short-circuits `build_pool` (no `urlopen`).
- All HTTP through a mocked `urllib.request.urlopen` (as in `test_pango.py`); the
  `nextclade` binary is never required.

`tests/unit/test_iterate.py`:
- `seed_source="nextclade"` routes through a mocked `_fetch_nextclade` into
  `select_regional` and seeds the collection.

## Dependencies & docs

- No new Python dependencies (urllib + stdlib JSON/regex). `nextclade` CLI
  optional. skani/skDER/mafft/datasets unchanged.
- README: add a "Nextclade-derived pool" subsection under the seed-source docs,
  including the segment-length caveat and the clade-naming benefit.

## Limitations

- Insertions relative to the dataset reference are dropped (reconstruction is in
  reference coordinates). Negligible for these panels.
- Single-gene/segment datasets (most flu segments, some others) yield
  segment-length pool genomes; recombination calling on one short segment is
  weaker than on a full genome. Auto-detect will pick them; the run logs a note
  rather than blocking.
- Pool composition (clade representatives, no recombinant-parental fragments) can
  prevent recovery of a specific recombinant's parental mosaic — see
  `memory/hiv-nextclade-pool-finding.md`.
