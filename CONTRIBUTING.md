# Contributing to Tessera

## Setup

```
pip install -e ".[dev]"
ruff check src tests validation
pytest -m "not requires_binary"     # fast suite, no aligner binaries needed
```

The fast suite runs in a few seconds and needs no external tools, network, or downloaded
data. The full suite additionally needs an aligner on `PATH`; see
[docs/aligners.md](docs/aligners.md).

For the aligner-dependent work, `environment.yml` creates a conda environment named
`tessera` with every backend:

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
- `mypy` does not add new errors. It is not yet blocking in CI — the codebase has an
  existing backlog — but it should not grow.
- New behaviour has a test that fails without the change. A test that passes either way
  documents intent but does not protect it.

## Changing detection behaviour

**Anything that touches a caller, a default, or region calling must be validated on the
real harness before it lands**, not only on the unit suite:

```
export PATH="$HOME/miniforge3/envs/tessera/bin:$PATH"
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

## Conventions

- **Dependency-free by design.** No new runtime dependency without discussion;
  reimplement in numpy rather than add a binary. The PHI/Rmin, 3SEQ, MaxChi and Bootscan
  callers all follow this.
- **Aligners are plugins.** New backends are added through the `tessera.aligners`
  entry-point group; core never imports a specific adapter.
- **External tools go through `core/process.run_tool`** — one choke point, no
  `shell=True`, arguments as lists.
- **Errors are `TesseraError` subclasses** raised from library code. Library code never
  calls `sys.exit`; `cli/main.py` turns exceptions into clean exits.
- **Modest scientific language** in code and documentation. State what passes, what
  fails, and what was skipped; do not overstate a result.
- **Commit messages** explain why the change is right, not just what changed: the
  problem, the evidence, and any trade-off taken. See the existing history for the shape.

## Reporting a problem

Include the Tessera version, the command, and `run_provenance.json` from the output
directory. That file records the parameters, the caller configuration, and — when the
alignment came from `tessera msa` — the aligner and its version, which is usually the
first thing worth knowing.
