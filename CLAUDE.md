# CLAUDE.md

Guidance for agents working in this repository.

## What this is

Tessera detects recombination in a query genome against a reference panel. It builds a
reference-anchored pseudo-MSA (one backbone sequence), scans it in sliding windows, and
calls recombinant regions with three complementary methods. It is a dependency-free
Python package (numpy / pandas / matplotlib / plotly / biopython / typer); alignment is
delegated to external aligner binaries through a plugin registry.

For the user-facing overview see [README.md](README.md). For depth:

- [docs/detection-methods.md](docs/detection-methods.md) -- the HMM, 3SEQ, and PHI/Rmin
  callers, low-divergence informative-site windowing, and every output file.
- [docs/reference-panels.md](docs/reference-panels.md) -- recruiting, finding, filling,
  and curating the donor panel; seed modes and sources.
- [docs/aligners.md](docs/aligners.md) -- the aligner backends.

## Layout

- `src/tessera/cli/` -- Typer app; one `cmd_*.py` per command (`detect`, `build-panel`,
  `msa`, `recomb`, `find-references`, `fill-references`, `curate-panel`). Importing a
  `cmd_*` module registers its command as a side effect.
- `src/tessera/recomb/` -- the scan: `similarity` (windowed distances), `hmm`,
  `threeseq`, `diagnostics` (PHI + Rmin), `regions`, `coverage`, `typing`, `report`,
  `run` (orchestration).
- `src/tessera/discover/` -- panel recruitment (BLAST, NCBI Virus, Nextclade, curation).
- `src/tessera/aligners/` -- aligner adapters, discovered via the `tessera.aligners`
  entry-point group (`pyproject.toml`); core never imports a specific adapter.
- `src/tessera/msa/`, `converters/`, `core/` -- alignment build, format converters,
  shared IO/cache/logging/errors.
- `tests/unit/`, `tests/integration/` -- pytest suite. `validation/` -- opt-in harnesses
  that exercise the real aligner path (not part of CI).

## Build, test, lint

```
pip install -e ".[dev]"
ruff check src tests validation
pytest -m "not requires_binary"     # fast suite, no aligner binaries needed
pytest                              # full suite (needs aligners on PATH)
```

The `requires_binary` marker gates tests that need an external aligner. Ruff config is
in `pyproject.toml` (line length 100, `E/F/I/UP/B`).

## Conventions

- Dependency-free by design: no new runtime dependency without discussion; reimplement
  in numpy rather than add a binary. The PHI/Rmin and 3SEQ callers follow this.
- Aligner adapters are added via the entry-point registry, never imported directly by
  core.
- Use modest scientific language in code and documentation; avoid overstating results.
- Reported numbers must be faithful (state what passes, what fails, what was skipped).

## Validation harnesses

`validation/run_hybrids.py` builds synthetic recombinants from Nextclade datasets and
checks detection across ~24 pathogens; `validation/run_validation.py` runs curated
published-recombination datasets. Both need an aligner env on `PATH`:

```
export PATH="$HOME/miniforge3/envs/recomfi-aln/bin:$PATH"   # mafft/minimap2/skani/skDER
python validation/run_hybrids.py
```

The aligner conda env is named `recomfi-aln` (predates the rename; it is external
infrastructure, not part of the package). See `validation/README.md`.

## Notes

- The on-disk cache is `~/.cache/tessera` (override with `$TESSERA_CACHE`).
- The package was renamed from `recomfi` to `tessera`; the git repo directory is still
  `Code/recomfi`. Do not "fix" the directory name or the `recomfi-aln` env name.
- Commit messages end with the project's co-author trailer; branch before committing on
  the default branch.
