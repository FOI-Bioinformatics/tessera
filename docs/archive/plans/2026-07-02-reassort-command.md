# `tessera reassort` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A `tessera reassort` command that assigns each segment of a multi-FASTA query to its nearest Nextclade reference (skani ANI) and calls reassortment by nearest-reference-strain consistency.

**Architecture:** New `src/tessera/reassort/` module with a pure `_call_reassortment` (unit-tested) and an I/O `assign_segments` (per-segment: resolve dataset -> build_pool tips -> nearest tip by skani -> strain+clade). A `cli/cmd_reassort.py` wraps it. The existing `recomb` single-backbone path is untouched.

**Tech Stack:** Python 3, typer, skani + Nextclade (already-optional external tools), pytest, ruff. Design spec: `docs/superpowers/specs/2026-07-02-reassort-command-design.md`.

## Global Constraints

- v1 = per-segment nearest-reference assignment + strain-consistency reassortment call; **no** per-segment intragenic recombination and **no** change to `detect`/`recomb`.
- For flu Nextclade tips `_accession_of` returns the (sanitized) **strain name** and `build_pool` writes each tip file `<strain>.fasta` with header `>{strain} {clade}` -- so nearest-tip -> strain (filename via `strip_sequence_extension`) + clade (header's second token).
- Reuse `resolve_dataset`, `build_pool` (`discover/nextclade.py`), `skani_query_ani`/`skani_available` (`discover/panel.py`), `nextclade_cache` (`core/cache.py`), `read_fasta`/`write_fasta_record`/`strip_sequence_extension` (`core/io.py`), `first_header`/`genotype`-style header parsing, the `cmd_*.py`+`main.py` registration pattern, `UserInputError`/`stage_errors`.
- No new runtime dependency. Modest language; ruff clean (line length 100, E/F/I/UP/B); CI coverage floor (80%) holds.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. Branch `reassort-command` (checked out); do not touch main.

---

### Task 1: `reassort` core -- `_call_reassortment` + `assign_segments`

**Files:** Create `src/tessera/reassort/__init__.py`, `src/tessera/reassort/assign.py`. Test: `tests/unit/test_reassort.py`.

**Interfaces:**
- Produces: `SegmentAssignment` (dataclass: `segment: str`, `dataset: str`, `strain: str|None`, `clade: str|None`, `ani: float`, `status: str`); `ReassortmentResult` (dataclass: `segments: list[SegmentAssignment]`, `verdict: str`); `_call_reassortment(per_segment_topk, assigned) -> str`; `assign_segments(query, *, dataset_overrides, ani_floor, email, cache_dir, logger) -> ReassortmentResult`.
- Consumes: `resolve_dataset`, `build_pool`, `skani_query_ani`, `nextclade_cache`, `read_fasta`, `write_fasta_record`, `strip_sequence_extension`.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_reassort.py`:

```python
"""Unit tests for the reassortment call (per-segment inputs; no binaries)."""

from __future__ import annotations

from tessera.reassort.assign import _call_reassortment


def test_clonal_when_a_strain_is_near_across_all_segments():
    # HA and NA both have strain S1 in their top-k -> one parent explains both -> clonal.
    topk = {"HA": {"S1", "S2"}, "NA": {"S1", "S3"}}
    assert _call_reassortment(topk, assigned=["HA", "NA"]) == "clonal"


def test_reassortant_when_topk_strain_sets_are_disjoint():
    topk = {"HA": {"S1", "S2"}, "NA": {"S3", "S4"}}
    assert _call_reassortment(topk, assigned=["HA", "NA"]) == "reassortant"


def test_undetermined_with_fewer_than_two_assigned_segments():
    assert _call_reassortment({"HA": {"S1"}}, assigned=["HA"]) == "undetermined"


def test_unassigned_segments_excluded_from_the_decision():
    # NA unassigned -> decision on HA + PB2; they share S1 -> clonal.
    topk = {"HA": {"S1"}, "PB2": {"S1"}}
    assert _call_reassortment(topk, assigned=["HA", "PB2"]) == "clonal"
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/unit/test_reassort.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'tessera.reassort'`.

- [ ] **Step 3: Create the module + `_call_reassortment`**

`src/tessera/reassort/__init__.py`:

```python
"""Per-segment reassortment detection (the ``tessera reassort`` command)."""

from .assign import ReassortmentResult, SegmentAssignment, assign_segments

__all__ = ["ReassortmentResult", "SegmentAssignment", "assign_segments"]
```

`src/tessera/reassort/assign.py` (the pure caller first):

```python
"""Assign each segment of a multi-FASTA query to its nearest Nextclade reference and
call reassortment by nearest-reference-strain consistency."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

DEFAULT_ANI_FLOOR = 80.0  # a segment below this ANI to every tip is left unassigned
TOP_K = 5                 # nearest strains kept per segment for the consistency test


@dataclass
class SegmentAssignment:
    segment: str
    dataset: str
    strain: str | None
    clade: str | None
    ani: float
    status: str  # "assigned" | "unassigned"


@dataclass
class ReassortmentResult:
    segments: list[SegmentAssignment] = field(default_factory=list)
    verdict: str = "undetermined"


def _call_reassortment(per_segment_topk: dict[str, set[str]], assigned: list[str]) -> str:
    """clonal iff one reference strain is in the top-k nearest of EVERY assigned segment
    (one parent explains all segments); reassortant iff no strain is shared across them;
    undetermined for < 2 assigned segments."""
    considered = [s for s in assigned if s in per_segment_topk]
    if len(considered) < 2:
        return "undetermined"
    shared = set.intersection(*(per_segment_topk[s] for s in considered))
    return "clonal" if shared else "reassortant"
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/unit/test_reassort.py -v`
Expected: PASS (4 passed).

- [ ] **Step 5: Add `assign_segments` (the I/O path)**

Append `assign_segments` to `assign.py`. For each query record (segment):
- Write the segment to a temp FASTA; resolve its dataset (`dataset_overrides.get(record_name)` else `resolve_dataset(segment_fasta, None, email=email, logger=logger)`); on failure record `SegmentAssignment(status="unassigned", dataset="?", strain=None, clade=None, ani=0.0)` and continue.
- `tips = build_pool(dataset, cache_dir=nextclade_cache(dataset.path, dataset.tag, override=cache_dir), logger=logger)`.
- `ani = skani_query_ani(segment_fasta, tips, logger)`; rank tips by ANI desc; keep the top-`TOP_K` whose ANI >= `ani_floor` as the segment's strain set (`strip_sequence_extension(tip.name)`); the best is `(strain, clade, ani)` where `clade` is the second token of the tip's `first_header`.
- If no tip clears the floor -> `unassigned`.
Collect `per_segment_topk[segment] = {strains}` for assigned segments; build `SegmentAssignment`s; `verdict = _call_reassortment(per_segment_topk, [assigned segment names])`. Return `ReassortmentResult`. Raise `UserInputError` up front if `not skani_available()`.

- [ ] **Step 6: Run tests + ruff**

Run: `python -m pytest tests/unit/test_reassort.py -v && ruff check src/tessera/reassort tests/unit/test_reassort.py`
Expected: PASS, ruff clean.

- [ ] **Step 7: Commit** ("Add the reassort core: per-segment assignment + strain-consistency call"; trailer).

---

### Task 2: `tessera reassort` CLI command

**Files:** Create `src/tessera/cli/cmd_reassort.py`; modify `src/tessera/cli/main.py` (register the import). Test: `tests/unit/test_cli_reassort.py`.

**Interfaces:** Consumes `assign_segments`. Produces the `reassort` command writing `reassortment.tsv` + a text summary.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_cli_reassort.py`:

```python
"""The reassort command registers and renders help."""

from __future__ import annotations

from typer.testing import CliRunner

from tessera.cli.main import app


def test_reassort_in_help():
    res = CliRunner().invoke(app, ["--help"])
    assert res.exit_code == 0 and "reassort" in res.output


def test_reassort_help_renders():
    res = CliRunner().invoke(app, ["reassort", "--help"])
    assert res.exit_code == 0 and "--query" in res.output
```

- [ ] **Step 2: Run to verify failure** (`reassort` absent from help).

- [ ] **Step 3: Add `cmd_reassort.py`**

Model on `cmd_detect.py`: `@app.command(name="reassort")` with `query: Path (-q/--query)`, `output: Path (-o/--output)`, `email: str | None`, `dataset: list[str]` (repeatable `seg=path` overrides, parsed to a dict), `ani_floor: float = DEFAULT_ANI_FLOOR`. Body: `read_fasta(query)` -> records (one per segment); `assign_segments(...)`; write `output/reassortment.tsv` (columns `segment dataset nearest_strain clade ani status`) and print the mosaic + verdict (e.g. `HA:J.2.2 (A/X) | NA:B.4 (A/Y) -> reassortant`). Wrap in `stage_errors(logger)`.

- [ ] **Step 4: Register in `main.py`**

Add `from . import cmd_reassort  # noqa: F401` to the side-effect import block (with the other `cmd_*` imports), so the command registers on app import.

- [ ] **Step 5: Run the test + `--help` checks + ruff**

Run: `python -m pytest tests/unit/test_cli_reassort.py -v && python -m tessera.cli.main reassort --help | grep -- --query && ruff check src/tessera/cli/cmd_reassort.py`
Expected: PASS, help shows `--query`, ruff clean.

- [ ] **Step 6: Commit** ("Add the tessera reassort CLI command"; trailer).

---

### Task 3: Harness -- flip `reassort_flu` to the new command; end-to-end; record

**Files:** Modify `validation/run_hybrids.py` (the `reassortant` frontier branch), `validation/attribution-results.md`, `docs/reference-panels.md` (a short command note).

- [ ] **Step 1: Point the frontier reassortant probe at `tessera reassort`**

Change `_run_frontier_case`'s `reassortant` branch: instead of the concatenated-query MSA/junction approach, write the query as a **multi-FASTA** (one record per segment: `HA_consensus[X]`, `NA_consensus[Y]`), call `assign_segments(...)` (or shell `tessera reassort`), and set the verdict from the result: `XPASS` when `verdict == "reassortant"` (the true reassortant is detected), else `XFAIL`. Keep `_prepare_frontier_reassortant`'s consensus building; just emit a 2-record multi-FASTA and route through the new command. Remove the now-unused junction/`_score_frontier_reassortant` path (or leave it dead-code-free).

- [ ] **Step 2: End-to-end gate (opt-in, aligner env)**

```bash
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
python validation/run_hybrids.py --frontier reassort_flu
```
Expected: `reassort_flu` now returns **reassortant** (XPASS) -- the HA and NA segments trace to different strains/clades. Also try a clonal control (HA + NA of the same strain) manually to confirm it reports `clonal`. If the strain-consistency call is too strict/loose on real flu tips (the 477-shared-strain caveat), tune `TOP_K`/`--ani-floor` for validity, not to force a result.

- [ ] **Step 3: Record**

Update `validation/attribution-results.md`: the reassortment known-limitation is now closed by `tessera reassort` (per-segment nearest-strain consistency), replacing the concatenated-junction probe that XFAILed. Document `tessera reassort` briefly in `docs/reference-panels.md`. State the measured verdict faithfully.

- [ ] **Step 4: Commit** ("Detect reassortment via tessera reassort; close the frontier known-limitation"; trailer).

---

## Self-Review

**Spec coverage:** the `reassort` module with `_call_reassortment` (pure) + `assign_segments` (per-segment Nextclade typing) -> Task 1; the CLI command -> Task 2; the harness flip + e2e + record -> Task 3. Assignment-only; no `recomb` change.

**Placeholder scan:** complete code for `_call_reassortment`, the dataclasses, and the module/CLI structure; `assign_segments` and the CLI body are I/O-heavy and precisely described (exercised by the opt-in gate + `--help`). No TODO/TBD.

**Type consistency:** `_call_reassortment(dict[str, set[str]], list[str]) -> str`; `assign_segments(...) -> ReassortmentResult` whose `segments` are `SegmentAssignment`; the CLI reads those to write the TSV + mosaic. Nearest-tip -> strain via `strip_sequence_extension(tip.name)`, clade via the tip header's second token -- consistent with how `build_pool` writes tip files.

## Execution Handoff

Two options: **1. Subagent-Driven (recommended)** or **2. Inline Execution**. Which approach?
