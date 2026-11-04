# Contributing to Tessera

## Setup

Tessera needs **Python 3.13 or newer** (`requires-python` in `pyproject.toml`), which
is the one version CI tests. The package uses 3.12 syntax, so 3.11 and older fail at
import with a `SyntaxError` rather than with a dependency message; 3.12 itself will
still import, but it is neither supported nor tested, and nothing keeps it working.

```
python3.13 -m venv .venv            # or any interpreter >= 3.13
.venv/bin/pip install -e ".[dev]"
.venv/bin/ruff check src tests validation
.venv/bin/pytest -m "not requires_binary"   # fast suite, no aligner binaries needed
```

`.venv` is git-ignored. If `python` on your PATH is already 3.13 or newer, the plain
forms work as before:

```
pip install -e ".[dev]"
ruff check src tests validation
pytest -m "not requires_binary"
```

The fast suite takes a few minutes (about four, last timed on an Apple-silicon laptop in a
single process) and needs no external tools, network, or downloaded data. For a quick loop while editing, run
the one test module you changed. The full suite additionally needs an aligner on `PATH`; see
[docs/aligners.md](docs/aligners.md).

For the aligner-dependent work, `environment.yml` creates a conda environment named
`tessera` with the default aligners (SibeliaZ, MAFFT, minimap2):

```
conda env create -f environment.yml
conda activate tessera
```

Some `validation/run_*.py` docstrings refer to an environment named `recomfi-aln`. That
is one developer's local environment name and predates the rename; any environment with
the aligners on `PATH` works.

## Before opening a pull request

- `ruff check src tests validation` is clean.
- `pytest -m "not requires_binary"` is green, and the test count has not gone down.
- `mypy` is clean. CI runs it as a blocking step over `src/tessera`.
- New behaviour has a test that fails without the change. A test that passes either way
  documents intent but does not protect it.
- A user-visible change has an entry under `## [Unreleased]` in `CHANGELOG.md`: what was
  wrong, what changed, and what was measured or left alone.

## Releasing

The version lives in `src/tessera/__init__.py`, and `CITATION.cff` must carry the same
version (a test fails when they differ, and the release workflow checks the tag against
both). To cut a release: move the `[Unreleased]` entries under a new `## [X.Y.Z] - date`
heading in `CHANGELOG.md` (the date the release is published), update `__version__` and
`CITATION.cff` (`version` and `date-released`), and tag `vX.Y.Z`.

## Changing detection behaviour

**Anything that touches a caller, a default, or region calling must be validated on the
real harness before it lands**, not only on the unit suite:

```
export PATH="$PATH:$HOME/miniforge3/envs/tessera/bin"   # append, never prepend
python validation/run_hybrids.py
```

This is not ceremony. A change that was unit-green and simulation-backed once shipped a
default that cost three pathogens their detections outright, and had to be reverted
(PR #47). Simulation is reliable for showing that a failure mode exists; it is not
reliable for deciding what a default should be, because it cannot reproduce how the
callers' applicability domains overlap at real divergence.

`validation/run_specificity.py` is the exception that needs nothing: it simulates its own
already-aligned clonal data, so it runs with no aligner, no network, and no downloaded
data. Every region it reports is a false positive by construction.

```
python validation/run_specificity.py --reps 3
```

The harness defaults to `--min-methods 2`, the CLI default is 1, and the two report
very different rates; pass `--min-methods 1` to measure the CLI default.

## Conventions

- **Few runtime dependencies by design.** No new one without discussion; reimplement
  in numpy rather than add a binary. The PHI/Rmin, 3SEQ, MaxChi and Bootscan
  callers all follow this.
- **Aligners are plugins.** New backends are added through the `tessera.aligners`
  entry-point group; core never imports a specific adapter.
- **External tools go through `core/process.run_tool`** — one choke point, no
  `shell=True`, arguments as lists. It also stops a timed-out tool's whole process group,
  so a wrapper script cannot leave its children running. The one exception is the
  `--version` probe in `core/binaries.py`, a bounded (30 s) call that follows the same
  no-shell, list-arguments rule.
- **Errors are `TesseraError` subclasses** raised from library code. Library code never
  calls `sys.exit`; `cli/main.py` turns exceptions into clean exits.
- **Modest scientific language** in code and documentation. State what passes, what
  fails, and what was skipped; do not overstate a result.
- **Commit messages** explain why the change is right, not just what changed: the
  problem, the evidence, and any trade-off taken. See the existing history for the shape.

## Reporting a problem

Include the Tessera version, the command, and `run_provenance.json` and `tessera.log` from
the output directory. That file records the parameters, the caller configuration, and — when the
alignment came from `tessera msa` — the aligner and its version, which is usually the
first thing worth knowing.
