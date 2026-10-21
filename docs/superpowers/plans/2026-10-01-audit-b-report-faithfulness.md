# Post-1.2.0 Audit, Plan B: Reporting Faithfulness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make every Tessera output state what the run did -- the callers that ran, the settings that shaped the result, what the coverage gaps mean, and what could not be tested -- without changing which regions are called.

**Architecture:** All changes sit in the reporting layer (`recomb/run.py` provenance and orchestration, `recomb/coverage.py`, the `recomb/report_*.py` writers) and in CLI input validation (`cli/`). No caller's region calls, coordinates or p-values change. One region flag changes: a coverage gap that a called breakpoint explains is relabelled `breakpoint` after region calling, so it no longer sets `donor_undercovered` or lowers the confidence wording; `call_coverage_gaps` itself is untouched, so reference recruitment behaves as before. "Could not test" becomes explicit in two places: a barcode caller on an untyped panel is `not run`, and a PHI test with too few informative sites is `not testable`.

**Tech Stack:** Python 3.11+, numpy, pandas, matplotlib, plotly, typer, biopython; pytest, ruff, mypy.

**Spec:** `docs/superpowers/specs/2026-10-01-post-1.2.0-audit-design.md`, items B1-B11. Read it first; the plan argues from it. Plans A and C cover the other items and touch other code; nothing here depends on them.

## Global Constraints

- No new runtime dependency. `seaborn` is removed, not replaced.
- Modest scientific language in code, docs and messages; reported numbers must be faithful (state what passes, what fails, what was skipped).
- A new behaviour needs a test that fails without the change.
- Do not change any caller's region calls, the PHI default window (`--phi-window 100`), or `call_coverage_gaps`. Those belong to plan C.
- "Could not test" must never be reported as "tested, found nothing".
- No test may reach the network. Tests that drive `detect`, `fill-references`, `build-panel`, `find-references`, `type-lineages` or `reassort` replace the pipeline entry point with a function that fails the test (`type-lineages` contacts NCBI before it validates anything else).
- Ruff: line length 100, rules `E/F/I/UP/B`. `mypy src` is blocking in CI.
- Run every command from the repository root; several tests use `example_data/...` relative paths.
- Append the aligner env to `PATH` (`PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"`), never prepend it: prepending swaps the Python interpreter.
- Work on branch `fix-audit-report-faithfulness`; commit messages end with the trailer `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`. Do not merge the PR: the maintainer does that.
- Baseline on `main` (commit `457bdfb`): `586 passed, 1 deselected`; ruff and mypy clean. After the last task: `671 passed, 1 deselected`.

## Review Focus

Inputs and conditions the spec implies but does not spell out, most likely first. Each has a test in the task that owns the code.

1. **A genuinely divergent stretch next to a breakpoint.** A missing third source adjacent to a called region must stay `divergent` and keep caveating the region; only a gap the two parents explain is relabelled. Task 3, `test_gap_the_two_parents_do_not_explain_stays_divergent`.
2. **A region whose parent label is not a row of the alignment** (the `n/a` placeholder when no backbone was resolved). Gap marking must skip it, not raise `KeyError`. Task 3, `test_region_naming_a_label_outside_the_alignment_is_skipped`.
3. **An output directory that already holds files from an earlier run** in another plot format. The footer must list what this run wrote, not what is on disk. Task 10, `test_footer_lists_the_plots_in_the_format_they_were_written`.
4. **A typed panel on which the barcode caller still cannot run** (every reference in one clade). It must be `not run` with the reason given, not treated as typed-and-negative. Task 7, `test_barcode_not_run_on_a_typed_panel_without_marked_clades`.
5. **Downstream readers of the PHI header.** `validation/run_hybrids.py` parses the first line of `recombination_profile.tsv`; it must still find the Rmin when the p-value is `NA`. Task 8, `test_parse_signal_reads_an_untestable_phi_header`.

## File Structure

| File | Responsibility after this plan | Tasks |
|---|---|---|
| `src/tessera/recomb/run.py` | orchestration; `caller_description`; provenance settings; breakpoint-gap marking call; barcode "not run" handling | 2, 3, 5, 7, 8 |
| `src/tessera/recomb/coverage.py` | gap calling (unchanged) plus `mark_breakpoint_gaps` | 3 |
| `src/tessera/recomb/report_context.py` | `caveat_gaps`, `methods_not_run`, `alpha` | 3, 7, 8 |
| `src/tessera/recomb/report_html.py` | union length, method table, PHI section, footer, methods text | 3, 4, 7, 8, 9, 10 |
| `src/tessera/recomb/report_text.py` | coverage console note, methods TSV, profile TSV header | 3, 7, 8 |
| `src/tessera/recomb/report_plots.py` | `region_labels`, donor-absent bands | 9 |
| `src/tessera/recomb/report.py` | `pair_datasets`, companion-file list | 7, 9, 10 |
| `src/tessera/recomb/report_assets.py` | glossary and references | 10 |
| `src/tessera/recomb/barcode.py` | returns no major parent when it cannot run | 7 |
| `src/tessera/recomb/diagnostics.py` | `phi_p` is `None` when not testable | 8 |
| `src/tessera/cli/main.py` and `cmd_*.py` | input validation | 6, 10, 12 |
| `src/tessera/core/binaries.py`, `aligners/sibeliaz.py` | version probe | 11 |
| `pyproject.toml` | drop `seaborn` | 11 |
| `docs/detection-methods.md`, `docs/reference-panels.md`, `example_data/README.md`, `validation/README.md`, `CHANGELOG.md` | documentation | 3, 7, 8, 9, 10, 11, 13 |
| `tests/unit/test_caller_provenance.py`, `test_breakpoint_gaps.py`, `test_report_summary.py`, `test_cli_input_checks.py`, `test_report_signal.py`, `test_report_plots.py`, `test_report_stale_text.py`, `test_dependencies.py` | new test files | 2-12 |

Edit steps are written as "replace this text with that text". Every "replace" block occurs exactly once in its file at the point the task is reached; the tasks must be done in order, because later tasks edit text earlier tasks introduced. When a step says "append to the end of" a test file, leave two blank lines between the existing content and the appended block.

Every task below was replayed on a fresh checkout of `457bdfb` before this plan was written: the tests of each task fail for the stated reason before its implementation step and pass after it, `ruff` and `mypy` are clean after each task, and the "Expected" lines are the recorded outputs.

---

### Task 1: Branch and baseline

**Files:** none changed.

**Interfaces:**
- Consumes: a clean checkout of `main` at or after `457bdfb`.
- Produces: branch `fix-audit-report-faithfulness`; baseline harness numbers saved outside the repository for Task 13.

- [ ] **Step 1: Create the branch**

```bash
git switch main && git pull --ff-only
git switch -c fix-audit-report-faithfulness
```

- [ ] **Step 2: Confirm the baseline**

Run: `ruff check src tests validation && pytest -m "not requires_binary" -q && mypy src`
Expected: `All checks passed!`, `586 passed, 1 deselected`, `Success: no issues found in 79 source files`.

- [ ] **Step 3: Record the "before" specificity numbers**

Task 3 changes the `donor_undercovered` flag, so the harness is run before and after. It needs no aligner, network or data; each command takes a few minutes.

```bash
mkdir -p /tmp/tessera-plan-b
python validation/run_specificity.py --reps 3 | tee /tmp/tessera-plan-b/before_gate2.txt
python validation/run_specificity.py --reps 3 --min-methods 1 | tee /tmp/tessera-plan-b/before_gate1.txt
```

Expected (the spec's "Measured context"): `TOTAL 0/12 (0%, CI 0%-24%) 0 --` at the default gate; `TOTAL 7/12 (58%, CI 32%-81%) 8 hmm=8` with `--min-methods 1`; positive control `detected 3/3 | correct donor 3/3 | median breakpoint error 55 bp` in both. If your numbers differ, stop and report them before changing anything: the plan's "after" comparison assumes this baseline.

- [ ] **Step 4: Record the "before" hybrids result, if the harness can run here**

`validation/run_hybrids.py` needs the Nextclade pools under `validation/data/hybrids/` (and the pool cache in `~/.cache/tessera`; without the cache it contacts the Nextclade dataset server) and mafft, skani and skDER on `PATH`. It takes a long time. If those are in place:

```bash
PY=$(which python)
PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin" $PY validation/run_hybrids.py \
    | tee /tmp/tessera-plan-b/hybrids_before.txt | tail -5
```

Expected: a final line of the form `N/M passed  (S skipped, E error)`. If the data or the tools are missing, skip this step and note that Task 13 Step 4 cannot be run either.

Nothing to commit.

---

### Task 2: Describe every caller under its own name (B1)

`run.py` describes every caller other than hmm and 3seq with the heuristic caller's text, so a default run records `hmm + 3seq + heuristic + heuristic` in `run_provenance.json` and the report. The record also omits four settings that change what is reported.

**Files:**
- Create: `tests/unit/test_caller_provenance.py`
- Modify: `src/tessera/recomb/run.py`

**Interfaces:**
- Consumes: the local `min_agree` already computed in `run_recomb`.
- Produces: `caller_description(method: str, params: RecombParams) -> str` (module level in `recomb/run.py`, replacing the nested `_caller_desc`); provenance keys `min methods (agreement gate)`, `sibling exclusion`, `lineage clustering`, `donor re-attribution`.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_caller_provenance.py`:

```python
"""The run provenance names the callers that ran and the settings that shape the result.

`run_provenance.json` and the report's "Run parameters" table are what a reader uses to
reproduce a run. A caller described under another caller's name, or a setting that
changes which regions are reported but is not recorded, makes the record misleading.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tessera.recomb.regions import CALLERS
from tessera.recomb.run import RecombParams, caller_description, run_recomb


def _params(example_data: Path, output: Path, **kwargs) -> RecombParams:
    return RecombParams(
        msa=example_data / "divergent.msa.fasta", output=output, query="query",
        window_size=300, window_step=30, plot_format="png", **kwargs,
    )


@pytest.mark.parametrize("method", CALLERS)
def test_each_caller_is_described_under_its_own_name(example_data, tmp_path, method) -> None:
    description = caller_description(method, _params(example_data, tmp_path))
    assert description.startswith(method)


@pytest.mark.parametrize("method", [m for m in CALLERS if m != "heuristic"])
def test_only_the_heuristic_caller_is_called_heuristic(example_data, tmp_path, method) -> None:
    # maxchi, bootscan, geneconv and barcode used to fall through to the heuristic text.
    assert "heuristic" not in caller_description(method, _params(example_data, tmp_path))


def test_default_run_names_its_four_callers(example_data, tmp_path, logger) -> None:
    out = tmp_path / "out"
    run_recomb(_params(example_data, out), logger)
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert "heuristic" not in run["caller"]
    for name in ("hmm", "3seq", "maxchi", "bootscan"):
        assert name in run["caller"]


def test_provenance_records_the_settings_that_change_what_is_reported(
    example_data, tmp_path, logger
) -> None:
    out = tmp_path / "default"
    run_recomb(_params(example_data, out), logger)
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["min methods (agreement gate)"] == "1"
    assert run["sibling exclusion"] == "on"
    assert run["lineage clustering"] == "on"
    assert run["donor re-attribution"] == "off"

    out = tmp_path / "changed"
    run_recomb(
        _params(example_data, out, min_methods=2, exclude_siblings=False,
                cluster_lineages=False, reattribute_donors=True, reattribute_margin=0.05),
        logger,
    )
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["min methods (agreement gate)"] == "2"
    assert run["sibling exclusion"] == "off"
    assert run["lineage clustering"] == "off"
    assert run["donor re-attribution"] == "on (margin 0.05)"
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_caller_provenance.py -q`
Expected: FAIL --

```text
ImportError: cannot import name 'caller_description' from 'tessera.recomb.run'
ERROR tests/unit/test_caller_provenance.py
1 error in 0.61s
```

- [ ] **Step 3: Implement**

In `src/tessera/recomb/run.py`, replace:

```python
def _select_windowing(bp_result, params: RecombParams, query_label: str, logger):
```

with:

```python
def caller_description(method: str, params: RecombParams) -> str:
    """One caller's name and the settings that govern it, for the run provenance.

    Every caller is described under its own name. The text is what a reader uses to
    tell which tests produced the regions, so a caller must never be described as a
    different one.
    """
    alpha = f"alpha {params.alpha:g}"
    if method == "hmm":
        return f"hmm (jump-rate {params.jump_rate:g}, {alpha})"
    if method == "3seq":
        return f"3seq (triplet max-descent test, {alpha})"
    if method == "maxchi":
        return f"maxchi (chi-square triplet test, scan-aware permutation, {alpha})"
    if method == "bootscan":
        return f"bootscan (bootstrap support, block-permutation run-length test, {alpha})"
    if method == "geneconv":
        return f"geneconv (longest donor-match run, permutation test, {alpha})"
    if method == "barcode":
        return "barcode (clade-marker attribution on a typed panel; no significance test)"
    min_region = params.min_region if params.min_region is not None else params.window_size
    merge_gap = params.merge_gap if params.merge_gap is not None else params.window_size
    return f"{method} (min {min_region} / margin {params.margin} / merge {merge_gap})"


def _select_windowing(bp_result, params: RecombParams, query_label: str, logger):
```

In `src/tessera/recomb/run.py`, replace:

```python
    def _caller_desc(method: str) -> str:
        if method == "hmm":
            return f"hmm (jump-rate {params.jump_rate:g}, alpha {params.alpha:g})"
        if method == "3seq":
            return f"3seq (triplet max-descent test, alpha {params.alpha:g})"
        min_region = params.min_region if params.min_region is not None else params.window_size
        merge_gap = params.merge_gap if params.merge_gap is not None else params.window_size
        return f"heuristic (min {min_region} / margin {params.margin} / merge {merge_gap})"

    caller_desc = " + ".join(_caller_desc(m) for m in params.methods)
```

with:

```python
    caller_desc = " + ".join(caller_description(m, params) for m in params.methods)
```

In `src/tessera/recomb/run.py`, replace:

```python
        "windowing": windowing,
        "major parent": major_parent or "n/a",
```

with:

```python
        "windowing": windowing,
        # Settings that change which regions are reported. Without them a run with
        # --min-methods 2 or --no-cluster-lineages is indistinguishable from a default
        # run in the record.
        "min methods (agreement gate)": str(min_agree),
        "sibling exclusion": "on" if params.exclude_siblings else "off",
        "lineage clustering": "on" if params.cluster_lineages else "off",
        "donor re-attribution": (
            f"on (margin {params.reattribute_margin:g})" if params.reattribute_donors else "off"
        ),
        "major parent": major_parent or "n/a",
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_caller_provenance.py tests/unit/test_provenance.py tests/integration/test_recomb_pipeline.py -q`
Expected: PASS -- `29 passed`

Run: `ruff check src tests validation && mypy src`
Expected: `All checks passed!` and `Success: no issues found in 79 source files`

- [ ] **Step 5: Commit**

```bash
git add tests/unit/test_caller_provenance.py src/tessera/recomb/run.py
git commit -m "Name each caller correctly in the run provenance" -m "Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 3: Label breakpoint-straddling coverage gaps (B2)

A window straddling a breakpoint between divergent parents has a best similarity below the coverage threshold, is called a `divergent` gap, and marks the region `donor_undercovered`. On `example_data/divergent.msa.fasta` the report says "low confidence. Possible missing reference" for a donor identical to the query. The gap is relabelled after region calling; `call_coverage_gaps` is not changed.

One refinement of the spec's rule, kept deliberately: besides lying within one window width of a called region boundary, the gap must be *explained* by the region's two parents (the query matches at least one of them at a fraction of comparable columns that reaches the coverage threshold). Distance alone would also relabel a genuinely divergent stretch that happens to sit next to a breakpoint.

**Files:**
- Create: `tests/unit/test_breakpoint_gaps.py`
- Modify: `src/tessera/recomb/coverage.py`
- Modify: `src/tessera/recomb/run.py`
- Modify: `src/tessera/recomb/report_context.py`
- Modify: `src/tessera/recomb/report_html.py`
- Modify: `src/tessera/recomb/report_text.py`
- Modify: `docs/detection-methods.md`
- Modify: `example_data/README.md`

**Interfaces:**
- Consumes: `CoverageGap`, `Region`, `similarity._canonical_mask`; `bp_result.rows`, `query_label`, `coverage_threshold` in `run_recomb`.
- Produces: `BREAKPOINT_KIND = "breakpoint"` and `mark_breakpoint_gaps(gaps: list[CoverageGap], regions: list[Region], rows: dict[str, np.ndarray], query: str, window_size: int, threshold: float) -> int` in `recomb/coverage.py`; `ReportContext.caveat_gaps` (property, `list[CoverageGap]`).

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_breakpoint_gaps.py`:

```python
"""A coverage gap produced by a window straddling a called breakpoint is not a
missing reference.

A window that spans a breakpoint between two divergent parents matches neither parent
well on its own, so its best similarity dips below the coverage threshold even though
both parents are in the panel. Such a gap is relabelled ``breakpoint``; it must not
caveat the region, become a donor-absent region, or reach the report headline.
"""

from __future__ import annotations

import csv
import re

import numpy as np

from tessera.recomb.coverage import (
    BREAKPOINT_KIND,
    CoverageGap,
    gaps_as_regions,
    mark_breakpoint_gaps,
)
from tessera.recomb.regions import Region
from tessera.recomb.run import RecombParams, run_recomb

WIDTH = 2000
WINDOW = 400


def _enc(seq: str) -> np.ndarray:
    return np.frombuffer(seq.encode("ascii"), dtype=np.uint8)


def _rows(query: str) -> dict[str, np.ndarray]:
    """Parents that differ at every tenth column; the query is given by the caller."""
    major = "A" * WIDTH
    minor = "".join("C" if i % 10 == 0 else "A" for i in range(WIDTH))
    return {"query": _enc(query), "major": _enc(major), "minor": _enc(minor)}


def _mosaic() -> str:
    """Major up to column 1000, minor after it: a clean breakpoint at 1000."""
    minor = "".join("C" if i % 10 == 0 else "A" for i in range(WIDTH))
    return "A" * 1000 + minor[1000:]


def _region() -> Region:
    return Region(
        minor_parent="minor", major_parent="major", msa_start=1000, msa_end=WIDTH,
        query_start=1000, query_end=WIDTH, n_windows=10,
        mean_sim_minor=1.0, mean_sim_major=0.9, margin=0.1,
    )


def _gap(start: int, end: int, kind: str = "divergent") -> CoverageGap:
    return CoverageGap(
        msa_start=start, msa_end=end, query_start=start, query_end=end,
        length_bp=end - start, n_windows=2, best_label="major", mean_best=0.94, kind=kind,
    )


def test_gap_straddling_a_called_breakpoint_is_relabelled() -> None:
    gaps = [_gap(800, 1200)]
    n = mark_breakpoint_gaps(gaps, [_region()], _rows(_mosaic()), "query", WINDOW, 0.95)
    assert n == 1
    assert gaps[0].kind == BREAKPOINT_KIND


def test_breakpoint_gap_is_not_bridged_to_a_donor_absent_region() -> None:
    gaps = [_gap(800, 1200)]
    mark_breakpoint_gaps(gaps, [_region()], _rows(_mosaic()), "query", WINDOW, 0.95)

    class _Result:
        positions: list[int] = []
        similarities: dict[str, list[float]] = {}

    assert gaps_as_regions(gaps, _Result(), "major") == []


def test_gap_the_two_parents_do_not_explain_stays_divergent() -> None:
    # Next to the breakpoint, but the query carries a base neither parent has at one
    # column in five: a third source, not a straddling window.
    query = list(_mosaic())
    for i in range(800, 1200, 5):
        query[i] = "G"
    gaps = [_gap(800, 1200)]
    n = mark_breakpoint_gaps(gaps, [_region()], _rows("".join(query)), "query", WINDOW, 0.95)
    assert n == 0
    assert gaps[0].kind == "divergent"


def test_gap_far_from_every_region_boundary_stays_divergent() -> None:
    gaps = [_gap(100, 500)]  # more than one window away from 1000 and from 2000
    n = mark_breakpoint_gaps(gaps, [_region()], _rows(_mosaic()), "query", WINDOW, 0.95)
    assert n == 0
    assert gaps[0].kind == "divergent"


def test_low_information_gap_is_left_alone() -> None:
    gaps = [_gap(800, 1200, kind="low_information")]
    mark_breakpoint_gaps(gaps, [_region()], _rows(_mosaic()), "query", WINDOW, 0.95)
    assert gaps[0].kind == "low_information"


def test_donor_absent_region_does_not_explain_a_gap() -> None:
    region = _region()
    region.donor_absent = True
    gaps = [_gap(800, 1200)]
    assert mark_breakpoint_gaps(gaps, [region], _rows(_mosaic()), "query", WINDOW, 0.95) == 0


def test_region_naming_a_label_outside_the_alignment_is_skipped() -> None:
    # A region can carry a placeholder parent ("n/a" when no backbone was resolved).
    region = _region()
    region.major_parent = "n/a"
    gaps = [_gap(800, 1200)]
    assert mark_breakpoint_gaps(gaps, [region], _rows(_mosaic()), "query", WINDOW, 0.95) == 0
    assert gaps[0].kind == "divergent"


# --- the shipped example, end to end --------------------------------------

def test_clean_recombinant_is_not_reported_as_a_missing_reference(
    example_data, tmp_path, logger
) -> None:
    """`example_data/divergent.msa.fasta`: the donor matches the query exactly over the
    insert and four callers agree. The only coverage gaps are the windows straddling
    the two breakpoints."""
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=example_data / "divergent.msa.fasta", output=out, query="query",
                     window_size=300, window_step=30, plot_format="png"),
        logger,
    )
    regions = list(csv.DictReader((out / "recombination_regions.tsv").open(), delimiter="\t"))
    assert len(regions) == 1
    assert regions[0]["donor_undercovered"] == "no"
    assert regions[0]["donor_absent"] == "no"

    lines = (out / "coverage_gaps.tsv").read_text().splitlines()
    gaps = list(csv.DictReader([ln for ln in lines if not ln.startswith("#")], delimiter="\t"))
    assert gaps and {g["kind"] for g in gaps} == {BREAKPOINT_KIND}

    report = (out / "report.html").read_text()
    verdict = re.sub(r"<[^>]+>", "", re.search(r'<p class="verdict">.*?</p>', report).group(0))
    assert "high confidence" in " ".join(verdict.split())
    assert "Possible missing reference" not in report
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_breakpoint_gaps.py -q`
Expected: FAIL --

```text
ImportError: cannot import name 'BREAKPOINT_KIND' from 'tessera.recomb.coverage'
ERROR tests/unit/test_breakpoint_gaps.py
1 error in 0.49s
```

- [ ] **Step 3: Implement**

In `src/tessera/recomb/coverage.py`, replace:

```python
A gap is labelled ``divergent`` (ample comparable bases, the query really is far
from all references -- a likely missing reference) or ``low_information`` (few
comparable bases, so the gap is uncertain rather than informative).
"""
```

with:

```python
A gap is labelled ``divergent`` (ample comparable bases, the query really is far
from all references -- a likely missing reference) or ``low_information`` (few
comparable bases, so the gap is uncertain rather than informative).

Once regions have been called, a ``divergent`` gap that sits on a called breakpoint
and is explained by that region's two parents together is relabelled ``breakpoint``
(:func:`mark_breakpoint_gaps`): a window straddling a switch between divergent
parents matches neither of them well on its own, which is a property of the window,
not a missing reference.
"""
```

In `src/tessera/recomb/coverage.py`, replace:

```python
from .regions import Region
from .similarity import WindowSimilarity
```

with:

```python
from .regions import Region
from .similarity import WindowSimilarity, _canonical_mask

# The kind given to a gap that a called breakpoint explains. See mark_breakpoint_gaps.
BREAKPOINT_KIND = "breakpoint"
```

In `src/tessera/recomb/coverage.py`, replace:

```python
    kind: str  # "divergent" | "low_information"
```

with:

```python
    kind: str  # "divergent" | "low_information" | "breakpoint"
```

In `src/tessera/recomb/coverage.py`, replace:

```python
def flag_undercovered_regions(regions: list[Region], threshold: float) -> None:
```

with:

```python
def _explained_fraction(
    rows: dict[str, np.ndarray], query: str, parents: tuple[str, str], start: int, end: int
) -> float:
    """Fraction of comparable columns in ``[start, end)`` where the query matches at
    least one of ``parents`` (``nan`` when nothing is comparable)."""
    q = rows[query][start:end]
    q_canon = _canonical_mask(q)
    comparable = np.zeros(q.size, dtype=bool)
    matched = np.zeros(q.size, dtype=bool)
    for label in parents:
        ref = rows[label][start:end]
        canon = q_canon & _canonical_mask(ref)
        comparable |= canon
        matched |= canon & (q == ref)
    n = int(np.count_nonzero(comparable))
    return float(np.count_nonzero(matched) / n) if n else float("nan")


def mark_breakpoint_gaps(
    gaps: list[CoverageGap],
    regions: list[Region],
    rows: dict[str, np.ndarray],
    query: str,
    window_size: int,
    threshold: float,
) -> int:
    """Relabel ``divergent`` gaps that a called breakpoint explains; return how many.

    A window straddling a breakpoint between two divergent parents is part one parent
    and part the other, so its similarity to either alone falls below the coverage
    threshold although both are in the panel. A gap is relabelled ``breakpoint`` when

    - it lies within one window width of a boundary of a called, donor-present region
      (the only place a straddling window can be), and
    - the region's donor and major parent together explain it: the fraction of
      comparable columns where the query matches at least one of them reaches
      ``threshold``.

    The second condition keeps a genuinely divergent stretch that happens to sit next
    to a breakpoint labelled ``divergent``. The gaps are mutated in place; call this
    after region calling and before :func:`gaps_as_regions`, which skips every kind
    other than ``divergent``.
    """
    relabelled = 0
    for gap in gaps:
        if gap.kind != "divergent":
            continue
        for region in regions:
            if region.donor_absent:
                continue
            parents = (region.minor_parent, region.major_parent)
            if any(label not in rows for label in parents):
                continue
            near = any(
                boundary - window_size <= gap.msa_start
                and gap.msa_end <= boundary + window_size
                for boundary in (region.msa_start, region.msa_end)
            )
            if not near:
                continue
            explained = _explained_fraction(rows, query, parents, gap.msa_start, gap.msa_end)
            if not isnan(explained) and explained >= threshold:
                gap.kind = BREAKPOINT_KIND
                relabelled += 1
                break
    return relabelled


def flag_undercovered_regions(regions: list[Region], threshold: float) -> None:
```

In `src/tessera/recomb/run.py`, replace:

```python
    gaps_as_regions,
    reconcile_gaps,
)
```

with:

```python
    gaps_as_regions,
    mark_breakpoint_gaps,
    reconcile_gaps,
)
```

In `src/tessera/recomb/run.py`, replace:

```python
    flag_undercovered_regions(regions, coverage_threshold)
    if coverage_gaps:
        logger.info(
            "Reference coverage: %d region(s) where the closest reference is below "
            "%.3f -- a better reference may be missing.",
            len(coverage_gaps), coverage_threshold,
        )
```

with:

```python
    flag_undercovered_regions(regions, coverage_threshold)
    # A window straddling a called breakpoint matches neither parent well on its own.
    # That is not a missing reference, so such gaps are relabelled before they can
    # caveat a region or be bridged to a donor-absent one. Recruitment
    # (fill-references / find-references) calls call_coverage_gaps directly and is
    # unaffected.
    n_breakpoint = mark_breakpoint_gaps(
        coverage_gaps, regions, bp_result.rows, query_label,
        params.window_size, coverage_threshold,
    )
    n_poor = len(coverage_gaps) - n_breakpoint
    if n_poor:
        logger.info(
            "Reference coverage: %d region(s) where the closest reference is below "
            "%.3f -- a better reference may be missing.",
            n_poor, coverage_threshold,
        )
    if n_breakpoint:
        logger.info(
            "Reference coverage: %d low-similarity stretch(es) sit on a called breakpoint "
            "and are explained by the two parents there; not treated as missing references.",
            n_breakpoint,
        )
```

In `src/tessera/recomb/report_context.py`, replace:

```python
from .coverage import CoverageGap
```

with:

```python
from .coverage import BREAKPOINT_KIND, CoverageGap
```

In `src/tessera/recomb/report_context.py`, replace:

```python
        return self.coverage_gaps or []
```

with:

```python
        return self.coverage_gaps or []

    @property
    def caveat_gaps(self) -> list[CoverageGap]:
        """The gaps that may mean a missing reference: every kind except ``breakpoint``.

        A breakpoint gap is a window straddling a called breakpoint. It is listed in the
        coverage table and ``coverage_gaps.tsv`` under its own kind, but it is not a
        poorly covered stretch, so the headline caveat, the mosaic and the plots leave
        it out.
        """
        return [g for g in self.gaps if g.kind != BREAKPOINT_KIND]
```

In `src/tessera/recomb/report_html.py`, replace:

```python
        f'bases to judge.</p>'
    )
```

with:

```python
        f'bases to judge; <strong>breakpoint</strong> = windows straddling a called '
        f'breakpoint, where the two parents together explain the query (not a missing '
        f'reference, and not counted in the caveat above).</p>'
    )
```

In `src/tessera/recomb/report_html.py`, replace:

```python
    fig = build_interactive_figure(
        ctx.site_result, datasets, regions, ctx.gaps,
```

with:

```python
    fig = build_interactive_figure(
        ctx.site_result, datasets, regions, ctx.caveat_gaps,
```

In `src/tessera/recomb/report_html.py`, replace:

```python
    gaps = ctx.gaps
    lineage_map = ctx.lineage_map
    threshold = ctx.coverage_threshold
    fig = build_interactive_figure(result, datasets, regions, gaps)
```

with:

```python
    gaps = ctx.gaps
    # Breakpoint gaps are tabulated but are not poorly covered stretches.
    caveat_gaps = ctx.caveat_gaps
    lineage_map = ctx.lineage_map
    threshold = ctx.coverage_threshold
    fig = build_interactive_figure(result, datasets, regions, caveat_gaps)
```

In `src/tessera/recomb/report_html.py`, replace:

```python
        f"{_caveat_html(gaps, threshold)}</header>"
```

with:

```python
        f"{_caveat_html(caveat_gaps, threshold)}</header>"
```

In `src/tessera/recomb/report_html.py`, replace:

```python
        f"{_mosaic_html(regions, colors, s, gaps, lineage_map)}</section>"
```

with:

```python
        f"{_mosaic_html(regions, colors, s, caveat_gaps, lineage_map)}</section>"
```

In `src/tessera/recomb/report_text.py`, replace:

```python
from .coverage import CoverageGap
```

with:

```python
from .coverage import BREAKPOINT_KIND, CoverageGap
```

In `src/tessera/recomb/report_text.py`, replace:

```python
    print_formatted_table(rows, header=COVERAGE_HEADER, echo=echo)
    echo("  ^ the closest reference here is poor; the true source may be missing. "
         "Run 'tessera find-references' to search NCBI.")
    echo("")
```

with:

```python
    print_formatted_table(rows, header=COVERAGE_HEADER, echo=echo)
    if any(g.kind != BREAKPOINT_KIND for g in gaps):
        echo("  ^ the closest reference here is poor; the true source may be missing. "
             "Run 'tessera find-references' to search NCBI.")
    if any(g.kind == BREAKPOINT_KIND for g in gaps):
        echo("  breakpoint = windows straddling a called breakpoint; the two parents "
             "together explain the query there (not a missing reference).")
    echo("")
```

In `docs/detection-methods.md`, replace:

```markdown
| `coverage_gaps.tsv` | Stretches where even the closest reference is a poor match -- possible missing references |
```

with:

```markdown
| `coverage_gaps.tsv` | Stretches where even the closest reference is below the best-similarity threshold, with a `kind`: `divergent` (the query is far from every reference -- a possible missing reference), `low_information` (too few comparable bases to judge), or `breakpoint` (windows straddling a called breakpoint, where the region's two parents together explain the query; not a missing reference, and it does not caveat the region) |
```

In `example_data/README.md`, replace:

```markdown
Both callers call `parent_B` over the insert (q-value ~1e-29) with a sharp breakpoint,
so the region is flagged as agreeing (high confidence); the similarity plot shows an
obvious crossover.
```

with:

```markdown
The four default callers all call `parent_B` over the insert with a sharp breakpoint,
so the region is flagged as agreeing (high confidence); the similarity plot shows an
obvious crossover. The two short stretches listed under reference coverage are the
windows straddling the breakpoints (kind `breakpoint`), not missing references.
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_breakpoint_gaps.py tests/unit/test_coverage.py tests/unit/test_coverage_reconcile.py tests/unit/test_output_contract.py tests/integration/test_recomb_pipeline.py -q`
Expected: PASS -- `41 passed`

Run: `ruff check src tests validation && mypy src`
Expected: `All checks passed!` and `Success: no issues found in 79 source files`

Then run the shipped example and read the result:

```bash
tessera recomb --msa example_data/divergent.msa.fasta --query query \
    --output /tmp/tessera-plan-b/divergent --window-size 300 --window-step 30 --plot-format png
cut -f1,2,5,6,21,22,23 /tmp/tessera-plan-b/divergent/recombination_regions.tsv
cat /tmp/tessera-plan-b/divergent/coverage_gaps.tsv
```

Expected: one region `parent_B  parent_A  960  2010  no  no  hmm,3seq,maxchi,bootscan`; two gaps (`810-1170`, `1860-2190`) both of kind `breakpoint`; the log line `2 low-similarity stretch(es) sit on a called breakpoint ...`; and the report verdict ends "high confidence." with no "Possible missing reference" box.

- [ ] **Step 5: Commit**

```bash
git add tests/unit/test_breakpoint_gaps.py src/tessera/recomb/coverage.py src/tessera/recomb/run.py src/tessera/recomb/report_context.py src/tessera/recomb/report_html.py src/tessera/recomb/report_text.py docs/detection-methods.md example_data/README.md
git commit -m "Do not report breakpoint-straddling windows as missing references" -m "Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 4: Report the union of overlapping regions (B3)

The verdict and the "Query recombinant" card add region lengths. The ensemble keeps overlapping regions that name different donors as separate rows, so a 2.2 kb union was reported as "3.3 kb (54.6 %)".

**Files:**
- Create: `tests/unit/test_report_summary.py`
- Modify: `src/tessera/recomb/report_html.py`

**Interfaces:**
- Consumes: nothing new.
- Produces: `_union_length(spans: list[tuple[int, int]]) -> int` in `recomb/report_html.py`.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_report_summary.py`:

```python
"""The report's headline numbers agree with the regions they summarise."""

from __future__ import annotations

import numpy as np

from tessera.recomb.regions import Region
from tessera.recomb.report_html import _summary, _verdict_html


class _Result:
    """The one attribute ``_summary`` reads: a 6000-base, gap-free query."""

    query = "query"
    query_cumulative = np.arange(6001)


def _region(minor: str, start: int, end: int, *, donor_absent: bool = False) -> Region:
    return Region(
        minor_parent=minor, major_parent="backbone", msa_start=start, msa_end=end,
        query_start=start, query_end=end, n_windows=5,
        mean_sim_minor=0.99, mean_sim_major=0.92, margin=0.07,
        qvalue=1e-9, support=0.95, methods=("hmm", "3seq"), donor_absent=donor_absent,
    )


def test_overlapping_regions_are_counted_once() -> None:
    """The ensemble keeps overlapping regions that name different donors separate.
    Adding their lengths reported a 2.2 kb union as 3.3 kb (54.6 % of the query)."""
    regions = [_region("donorB", 1900, 3100), _region("donorA", 2025, 4100)]
    s = _summary(_Result(), regions, ["backbone", "donorA", "donorB"])
    assert s["recomb_bp"] == 2200  # the union 1900-4100, not 1200 + 2075
    assert round(s["pct"], 1) == 36.7
    assert s["n_regions"] == 2  # both regions are still listed


def test_recombinant_fraction_never_exceeds_the_query() -> None:
    regions = [_region(f"donor{i}", 0, 6000) for i in range(3)]
    s = _summary(_Result(), regions, ["backbone"])
    assert s["recomb_bp"] == 6000
    assert s["pct"] == 100.0


def test_touching_and_nested_regions() -> None:
    touching = [_region("donorA", 1000, 2000), _region("donorB", 2000, 3000)]
    assert _summary(_Result(), touching, ["backbone"])["recomb_bp"] == 2000
    nested = [_region("donorA", 1000, 4000), _region("donorB", 2000, 2500)]
    assert _summary(_Result(), nested, ["backbone"])["recomb_bp"] == 3000


def test_disjoint_regions_still_add_up() -> None:
    regions = [_region("donorA", 500, 1500), _region("donorB", 3000, 3500)]
    assert _summary(_Result(), regions, ["backbone"])["recomb_bp"] == 1500


def test_donor_absent_regions_stay_out_of_the_recombinant_span() -> None:
    regions = [_region("donorA", 500, 1500), _region("x", 3000, 5000, donor_absent=True)]
    s = _summary(_Result(), regions, ["backbone"])
    assert s["recomb_bp"] == 1000
    assert s["n_absent"] == 1


def test_verdict_states_the_union() -> None:
    regions = [_region("donorB", 1900, 3100), _region("donorA", 2025, 4100)]
    s = _summary(_Result(), regions, ["backbone", "donorA", "donorB"])
    verdict = _verdict_html(s, "query", {})
    assert "2.2 kb" in verdict
    assert "36.7%" in verdict
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_report_summary.py -q`
Expected: FAIL --

```text
FAILED tests/unit/test_report_summary.py::test_overlapping_regions_are_counted_once
FAILED tests/unit/test_report_summary.py::test_recombinant_fraction_never_exceeds_the_query
FAILED tests/unit/test_report_summary.py::test_touching_and_nested_regions
FAILED tests/unit/test_report_summary.py::test_verdict_states_the_union
4 failed, 2 passed in 0.45s
```

- [ ] **Step 3: Implement**

In `src/tessera/recomb/report_html.py`, replace:

```python
def _summary(
    result: WindowSimilarity, regions: list[Region], datasets: list[str]
) -> dict:
```

with:

```python
def _union_length(spans: list[tuple[int, int]]) -> int:
    """Total length covered by ``spans`` (half-open intervals), overlaps counted once."""
    total = 0
    covered_to: int | None = None
    for start, end in sorted(spans):
        if end <= start:
            continue
        if covered_to is None or start > covered_to:
            total += end - start
            covered_to = end
        elif end > covered_to:
            total += end - covered_to
            covered_to = end
    return total


def _summary(
    result: WindowSimilarity, regions: list[Region], datasets: list[str]
) -> dict:
```

In `src/tessera/recomb/report_html.py`, replace:

```python
    recomb_bp = sum(max(0, r.query_end - r.query_start) for r in present)
```

with:

```python
    # The union, not the sum: overlapping regions that name different donors are kept
    # as separate rows, and adding their lengths counts the shared stretch twice.
    recomb_bp = _union_length([(r.query_start, r.query_end) for r in present])
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_report_summary.py tests/unit/test_report_typed.py -q`
Expected: PASS -- `9 passed`

Run: `ruff check src tests validation && mypy src`
Expected: `All checks passed!` and `Success: no issues found in 79 source files`

- [ ] **Step 5: Commit**

```bash
git add tests/unit/test_report_summary.py src/tessera/recomb/report_html.py
git commit -m "Report the union of overlapping regions in the verdict" -m "Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 5: Keep the method comparison in step with re-attribution (B4)

`--reattribute-donors` re-labels the region but not the per-method breakdown rows, so `recombination_methods.tsv` and the report's method table name the donor from before re-attribution.

The spec suggests returning an old-to-new mapping. That is not needed: `reattribute_donors` returns exactly one region per input region in the same order, and `regions` and `method_breakdown` are parallel lists (built together, sorted alike, filtered together), so the rows can be updated in step with a strict `zip`.

**Files:**
- Modify: `src/tessera/recomb/run.py`
- Test (extend): `tests/unit/test_reattribute.py`

**Interfaces:**
- Consumes: `regions` and `method_breakdown` in `run_recomb` being parallel lists.
- Produces: nothing new.

- [ ] **Step 1: Write the failing tests**

In `tests/unit/test_reattribute.py`, replace the imports:

```python
import numpy as np

from tessera.recomb.reattribute import reattribute_donors
from tessera.recomb.regions import Region
```

with:

```python
import csv
import random
from pathlib import Path

import numpy as np

from tessera.recomb.reattribute import reattribute_donors
from tessera.recomb.regions import Region
from tessera.recomb.run import RecombParams, run_recomb

from ..conftest import write_fasta
```

Append to the end of `tests/unit/test_reattribute.py`:

```python
# --- the run's other outputs follow the re-attribution ---------------------

def _reattribution_panel(tmp_path: Path) -> tuple[Path, dict[str, str]]:
    """A cowpox backbone with a variola insert at 2000-4000. The insert is copied from
    ``variolaA``, but ``variolaA`` shares clade V with two distant genomes, so clade V's
    consensus matches the insert worse than clade W's (``variolaB`` alone)."""
    rng = random.Random(11)

    def mutate(seq: str, frac: float) -> str:
        chars = list(seq)
        for i in range(len(chars)):
            if rng.random() < frac:
                chars[i] = rng.choice("ACGT")
        return "".join(chars)

    base = "".join(rng.choice("ACGT") for _ in range(6000))
    cowpox, variola = mutate(base, 0.03), mutate(base, 0.08)
    variola_b = mutate(variola, 0.01)
    query = list(cowpox)
    query[2000:4000] = list(variola[2000:4000])
    msa = write_fasta(tmp_path / "panel.fasta", {
        "query": "".join(query), "cowpox": cowpox, "variolaA": variola,
        "variolaB": variola_b, "junk1": mutate(base, 0.15), "junk2": mutate(base, 0.15),
    })
    lineages = {"cowpox": "CPX", "variolaA": "V", "junk1": "V", "junk2": "V", "variolaB": "W"}
    return msa, lineages


def test_method_comparison_names_the_reattributed_donor(tmp_path: Path, logger) -> None:
    """Re-attribution re-labelled the region but not the per-method breakdown, so
    `recombination_methods.tsv` and the report's method table kept the old donor."""
    msa, lineages = _reattribution_panel(tmp_path)
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=msa, output=out, query="query", plot_format="png",
                     lineage_map=lineages, reattribute_donors=True, cluster_lineages=False),
        logger,
    )
    regions = list(csv.DictReader((out / "recombination_regions.tsv").open(), delimiter="\t"))
    methods = list(csv.DictReader((out / "recombination_methods.tsv").open(), delimiter="\t"))
    assert [r["minor_parent"] for r in regions] == ["variolaB"]  # re-attributed from variolaA
    assert [m["minor_parent"] for m in methods] == ["variolaB"]
    assert [(m["query_start"], m["query_end"]) for m in methods] == [
        (r["query_start"], r["query_end"]) for r in regions
    ]
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_reattribute.py -q`
Expected: FAIL --

```text
FAILED tests/unit/test_reattribute.py::test_method_comparison_names_the_reattributed_donor
1 failed, 7 passed in 2.72s
```

- [ ] **Step 3: Implement**

In `src/tessera/recomb/run.py`, replace:

```python
            margin=params.reattribute_margin, logger=logger,
        )
```

with:

```python
            margin=params.reattribute_margin, logger=logger,
        )
        # The breakdown rows were built from the regions before re-attribution and are
        # written to recombination_methods.tsv and the report's method table. Keep the
        # donor they name in step with the region. reattribute_donors returns one region
        # per input region in the same order, so the two lists stay parallel.
        for region, row in zip(regions, method_breakdown, strict=True):
            row["minor_parent"] = region.minor_parent
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_reattribute.py tests/unit/test_ensemble.py -q`
Expected: PASS -- `24 passed`

Run: `ruff check src tests validation && mypy src`
Expected: `All checks passed!` and `Success: no issues found in 79 source files`

- [ ] **Step 5: Commit**

```bash
git add tests/unit/test_reattribute.py src/tessera/recomb/run.py
git commit -m "Keep the method comparison in step with donor re-attribution" -m "Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 6: Reject a --lineage-map that does not exist (B5)

`--lineage-map` pointing at a file that does not exist is ignored (exit 0, untyped report) in five commands, because the readers treat a missing lineage file as "no typed names". The barcode caller and donor re-attribution then do nothing, silently.

This task creates `tests/unit/test_cli_input_checks.py`; Task 12 appends to it. Its `no_pipeline` fixture is what keeps these tests off the network.

**Files:**
- Create: `tests/unit/test_cli_input_checks.py`
- Modify: `src/tessera/cli/main.py`
- Modify: `src/tessera/cli/cmd_recomb.py`
- Modify: `src/tessera/cli/cmd_type_lineages.py`
- Modify: `src/tessera/cli/cmd_detect.py`
- Modify: `src/tessera/cli/cmd_build_panel.py`
- Modify: `src/tessera/cli/cmd_fill_references.py`

**Interfaces:**
- Consumes: `_require_file` in `cli/main.py`.
- Produces: `_require_lineage_map(path: Path | None) -> None` in `cli/main.py`. Task 12 extends the import lists this task writes in `cmd_detect.py`, `cmd_build_panel.py`, `cmd_fill_references.py` and `cmd_recomb.py`.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_cli_input_checks.py`:

```python
"""Inputs the CLI must reject before any work starts.

Each case used to be accepted: the option was silently ignored, or the run failed
later under "Unexpected error". None of these tests may reach the network -- the
panel-building entry points are replaced with a function that fails the test if the
command gets that far.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

from tessera.cli.main import app

runner = CliRunner()

EXAMPLE = "example_data/divergent.msa.fasta"


def _must_not_run(*args, **kwargs):
    raise AssertionError("the command reached the pipeline; it should have been rejected")


@pytest.fixture
def no_pipeline(monkeypatch):
    """Fail the test if a command gets past validation into panel building or typing."""
    monkeypatch.setattr("tessera.discover.iterate.fill_references", _must_not_run)
    monkeypatch.setattr("tessera.discover.lineage_assign.assign_lineages", _must_not_run)
    monkeypatch.setattr("tessera.discover.run.find_references", _must_not_run)


def _query(tmp_path: Path) -> Path:
    path = tmp_path / "query.fasta"
    path.write_text(">query\nACGTACGTACGT\n")
    return path


def _collection(tmp_path: Path) -> Path:
    directory = tmp_path / "collection"
    directory.mkdir()
    (directory / "ref.fasta").write_text(">ref\nACGTACGTACGT\n")
    return directory


def _rejected(result, needle: str) -> None:
    assert result.exit_code == 1, result.output
    assert "Unexpected error" not in result.output
    assert needle in " ".join(result.output.split())


# --- --lineage-map must exist ---------------------------------------------

def test_recomb_rejects_a_missing_lineage_map(tmp_path: Path) -> None:
    """A mistyped path used to be ignored: the run exited 0 with an untyped report, and
    the barcode caller and donor re-attribution silently did nothing."""
    result = runner.invoke(app, [
        "recomb", "--msa", EXAMPLE, "--query", "query",
        "--output", str(tmp_path / "out"), "--lineage-map", str(tmp_path / "nope.tsv"),
    ])
    _rejected(result, "--lineage-map file not found")
    assert not (tmp_path / "out" / "report.html").exists()


def test_recomb_rejects_a_lineage_map_that_is_a_directory(tmp_path: Path) -> None:
    result = runner.invoke(app, [
        "recomb", "--msa", EXAMPLE, "--query", "query",
        "--output", str(tmp_path / "out"), "--lineage-map", str(tmp_path),
    ])
    _rejected(result, "is a directory, not a file")


def test_type_lineages_rejects_a_missing_lineage_map(tmp_path: Path, no_pipeline) -> None:
    result = runner.invoke(app, [
        "type-lineages", "--collection", str(_collection(tmp_path)),
        "--output", str(tmp_path / "out"), "--lineage-map", str(tmp_path / "nope.tsv"),
    ])
    _rejected(result, "--lineage-map file not found")


@pytest.mark.parametrize("command", ["detect", "fill-references", "build-panel"])
def test_panel_commands_reject_a_missing_lineage_map(
    tmp_path: Path, no_pipeline, command: str
) -> None:
    result = runner.invoke(app, [
        command, "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        "--lineage-map", str(tmp_path / "nope.tsv"),
    ])
    _rejected(result, "--lineage-map file not found")
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_cli_input_checks.py -q`
Expected: FAIL --

```text
FAILED tests/unit/test_cli_input_checks.py::test_recomb_rejects_a_missing_lineage_map
FAILED tests/unit/test_cli_input_checks.py::test_recomb_rejects_a_lineage_map_that_is_a_directory
FAILED tests/unit/test_cli_input_checks.py::test_type_lineages_rejects_a_missing_lineage_map
FAILED tests/unit/test_cli_input_checks.py::test_panel_commands_reject_a_missing_lineage_map[detect]
FAILED tests/unit/test_cli_input_checks.py::test_panel_commands_reject_a_missing_lineage_map[fill-references]
FAILED tests/unit/test_cli_input_checks.py::test_panel_commands_reject_a_missing_lineage_map[build-panel]
6 failed in 2.52s
```

- [ ] **Step 3: Implement**

In `src/tessera/cli/main.py`, replace:

```python
def _require_directory(path: Path, label: str) -> None:
```

with:

```python
def _require_lineage_map(path: Path | None) -> None:
    """Reject a ``--lineage-map`` that does not exist.

    The readers treat a missing lineage file as "no typed names", which is right for
    the automatically discovered ``lineages.tsv`` and wrong for a path the user typed:
    the run would finish with an untyped report, and the barcode caller and donor
    re-attribution would do nothing, without a word.
    """
    if path is not None:
        _require_file(path, "--lineage-map file")


def _require_directory(path: Path, label: str) -> None:
```

In `src/tessera/cli/cmd_recomb.py`, replace:

```python
    _require_file,
    _require_range,
```

with:

```python
    _require_file,
    _require_lineage_map,
    _require_range,
```

In `src/tessera/cli/cmd_recomb.py`, replace:

```python
        _require_file(msa, "MSA file")
```

with:

```python
        _require_file(msa, "MSA file")
        _require_lineage_map(lineage_map)
```

In `src/tessera/cli/cmd_type_lineages.py`, replace:

```python
from .main import _require_directory, app, get_logger, stage_errors
```

with:

```python
from .main import _require_directory, _require_lineage_map, app, get_logger, stage_errors
```

In `src/tessera/cli/cmd_type_lineages.py`, replace:

```python
        _require_directory(collection, "Collection directory")
```

with:

```python
        _require_directory(collection, "Collection directory")
        _require_lineage_map(lineage_map)
```

In `src/tessera/cli/cmd_detect.py`, replace:

```python
from .main import _require_choice, _require_file, app, get_logger, stage_errors
```

with:

```python
from .main import (
    _require_choice,
    _require_file,
    _require_lineage_map,
    app,
    get_logger,
    stage_errors,
)
```

In `src/tessera/cli/cmd_detect.py`, replace:

```python
        _require_file(query, "Query file")
```

with:

```python
        _require_file(query, "Query file")
        _require_lineage_map(lineage_map)
```

In `src/tessera/cli/cmd_build_panel.py`, replace:

```python
from .main import _require_choice, _require_file, app, get_logger, stage_errors
```

with:

```python
from .main import (
    _require_choice,
    _require_file,
    _require_lineage_map,
    app,
    get_logger,
    stage_errors,
)
```

In `src/tessera/cli/cmd_build_panel.py`, replace:

```python
        _require_file(query, "Query file")
```

with:

```python
        _require_file(query, "Query file")
        _require_lineage_map(lineage_map)
```

In `src/tessera/cli/cmd_fill_references.py`, replace:

```python
    _require_file,
    app,
```

with:

```python
    _require_file,
    _require_lineage_map,
    app,
```

In `src/tessera/cli/cmd_fill_references.py`, replace:

```python
        _require_file(query, "Query file")
```

with:

```python
        _require_file(query, "Query file")
        _require_lineage_map(lineage_map)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_cli_input_checks.py tests/unit/test_cli_commands.py tests/unit/test_cli_validation.py -q`
Expected: PASS -- `54 passed`

Run: `ruff check src tests validation && mypy src`
Expected: `All checks passed!` and `Success: no issues found in 79 source files`

- [ ] **Step 5: Commit**

```bash
git add tests/unit/test_cli_input_checks.py src/tessera/cli/main.py src/tessera/cli/cmd_recomb.py src/tessera/cli/cmd_type_lineages.py src/tessera/cli/cmd_detect.py src/tessera/cli/cmd_build_panel.py src/tessera/cli/cmd_fill_references.py
git commit -m "Reject a --lineage-map that does not exist" -m "Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 7: Report a barcode caller that could not run as not run (B6)

On an untyped panel the barcode caller returns no regions and names the first record of the alignment as major parent. A `--method barcode` run then reports that record as the backbone and a clean negative; in an ensemble the method table shows barcode as `no`.

**Files:**
- Modify: `src/tessera/recomb/barcode.py`
- Modify: `src/tessera/recomb/run.py`
- Modify: `src/tessera/recomb/report_context.py`
- Modify: `src/tessera/recomb/report_text.py`
- Modify: `src/tessera/recomb/report.py`
- Modify: `src/tessera/recomb/report_html.py`
- Modify: `docs/detection-methods.md`
- Test (extend): `tests/unit/test_barcode.py`

**Interfaces:**
- Consumes: `majors: dict[str, str | None]` in `run_recomb`; `recombinant_msa` and `write_fasta` from `tests/conftest.py`.
- Produces: `call_regions_barcode` returns `([], None, [])` when it cannot run; locals `not_run: tuple[str, ...]` and `n_ran: int` in `run_recomb`; `ReportContext.methods_not_run: tuple[str, ...] = ()`; `write_methods_tsv(..., methods_not_run=())`; `_method_comparison_html(..., methods_not_run=())` and `_method_section(..., methods_not_run=())`; provenance key `callers not run`. Task 8 edits the `methods_not_run=not_run,` line; Task 10 edits the `write_methods_tsv(...)` call block.

- [ ] **Step 1: Write the failing tests**

In `tests/unit/test_barcode.py`, replace the imports:

```python
from pathlib import Path

import numpy as np

from tessera.recomb.analyze import analyze
from tessera.recomb.barcode import clade_markers
from tessera.recomb.regions import RegionParams, call_regions
from tessera.recomb.similarity import compute_similarity

from ..conftest import write_fasta
```

with:

```python
import csv
import json
import logging
from pathlib import Path

import numpy as np
import pytest

from tessera.core.errors import UserInputError
from tessera.recomb.analyze import analyze
from tessera.recomb.barcode import clade_markers
from tessera.recomb.regions import RegionParams, call_regions
from tessera.recomb.run import RecombParams, run_recomb
from tessera.recomb.similarity import compute_similarity

from ..conftest import recombinant_msa, write_fasta
```

Append to the end of `tests/unit/test_barcode.py`:

```python
# --- "could not run" is not "found nothing" --------------------------------

def test_barcode_names_no_major_parent_when_it_cannot_run(tmp_path: Path) -> None:
    """It used to return the first record of the alignment as the major parent, which a
    barcode-only run then reported as the backbone."""
    result = compute_similarity(str(_panel(tmp_path)), "query", window_size=300, window_step=30)
    params = RegionParams.with_defaults(300, method="barcode")  # no lineage map
    regions, major, _ = call_regions(result, analyze(result), 300, params)
    assert regions == []
    assert major is None


def _run(tmp_path: Path, logger, **kwargs) -> Path:
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=recombinant_msa(tmp_path, recombinant=True), output=out,
                     query="query", plot_format="png", **kwargs),
        logger,
    )
    return out


def test_barcode_only_run_on_an_untyped_panel_is_refused(tmp_path: Path, logger) -> None:
    """Nothing was tested, so there is no result to report -- not a clean negative."""
    with pytest.raises(UserInputError, match="typed references"):
        _run(tmp_path, logger, methods=("barcode",))
    assert not (tmp_path / "out" / "report.html").exists()


def _collecting_logger() -> tuple[logging.Logger, list[str]]:
    """A logger outside the ``tessera`` hierarchy that keeps its warnings in a list.

    Not ``caplog``: once a CLI test has configured the ``tessera`` logger it stops
    propagating to the root logger, and whether that has happened depends on test order.
    """
    messages: list[str] = []

    class _Collect(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            messages.append(record.getMessage())

    log = logging.getLogger("test_barcode.collect")
    log.handlers = [_Collect(level=logging.WARNING)]
    log.propagate = False
    return log, messages


def test_ensemble_reports_barcode_as_not_run_on_an_untyped_panel(tmp_path: Path) -> None:
    log, warnings = _collecting_logger()
    out = _run(tmp_path, log, methods=("3seq", "barcode"))
    assert any("barcode" in message and "could not run" in message for message in warnings)

    rows = list(csv.DictReader((out / "recombination_methods.tsv").open(), delimiter="\t"))
    assert rows, "3seq should still call the insert"
    assert {row["3seq"] for row in rows} == {"yes"}
    assert {row["barcode"] for row in rows} == {"not run"}  # was "no"

    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert "barcode" in run["callers not run"]
    assert "not run" in (out / "report.html").read_text()


def test_agreement_gate_counts_only_the_callers_that_ran(tmp_path: Path, logger) -> None:
    """--min-methods is clamped to the callers actually run. A caller that could not run
    cannot agree, so it must not make the gate unreachable."""
    out = _run(tmp_path, logger, methods=("3seq", "barcode"), min_methods=2)
    rows = list(csv.DictReader((out / "recombination_regions.tsv").open(), delimiter="\t"))
    assert [row["methods"] for row in rows] == ["3seq"]


def test_barcode_not_run_on_a_typed_panel_without_marked_clades(tmp_path: Path) -> None:
    """Typed, but every reference is in one clade: there are no clade markers to
    compete, so the caller cannot run and the warning says why."""
    log, warnings = _collecting_logger()
    one_clade = {"A": "L1", "B": "L1", "other": "L1"}
    out = _run(tmp_path, log, methods=("3seq", "barcode"), lineage_map=one_clade)
    assert any("fewer than two typed clades" in message for message in warnings)
    rows = list(csv.DictReader((out / "recombination_methods.tsv").open(), delimiter="\t"))
    assert {row["barcode"] for row in rows} == {"not run"}
```

The warning test uses its own logger instead of `caplog`: once a CLI test has configured the `tessera` logger it stops propagating to the root logger, so `caplog` would make the test depend on test order.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_barcode.py -q`
Expected: FAIL --

```text
FAILED tests/unit/test_barcode.py::test_barcode_names_no_major_parent_when_it_cannot_run
FAILED tests/unit/test_barcode.py::test_barcode_only_run_on_an_untyped_panel_is_refused
FAILED tests/unit/test_barcode.py::test_ensemble_reports_barcode_as_not_run_on_an_untyped_panel
FAILED tests/unit/test_barcode.py::test_agreement_gate_counts_only_the_callers_that_ran
FAILED tests/unit/test_barcode.py::test_barcode_not_run_on_a_typed_panel_without_marked_clades
5 failed, 4 passed in 1.84s
```

- [ ] **Step 3: Implement**

In `src/tessera/recomb/barcode.py`, replace:

```python
    that donor clade. Returns ``(regions, major, [])`` to match ``call_regions``.
    """
```

with:

```python
    that donor clade. Returns ``(regions, major, [])`` to match ``call_regions``.

    When the caller cannot run -- the panel is untyped, or fewer than two clades carry
    enough markers -- it returns ``([], None, [])``. ``major is None`` is the signal that
    nothing was tested; the pipeline reports the caller as not run rather than as having
    found nothing.
    """
```

In `src/tessera/recomb/barcode.py`, replace:

```python
    labels = list(result.similarities)
    default_major = labels[0] if labels else None
    lineage_map = getattr(params, "lineage_map", None)
    if not lineage_map:
        return [], default_major, []
    cols, alleles, rep = clade_markers(result.rows, result.query, lineage_map)
    if not cols:
        return [], default_major, []
```

with:

```python
    lineage_map = getattr(params, "lineage_map", None)
    if not lineage_map:
        return [], None, []
    cols, alleles, rep = clade_markers(result.rows, result.query, lineage_map)
    if not cols:
        return [], None, []
```

In `src/tessera/recomb/run.py`, replace:

```python
    major_parent, per_major = reconcile_major(majors, window_wins=analysis_bp.winners_with_ties)
```

with:

```python
    # The barcode caller needs typed references and returns no major parent when it
    # cannot run. That is "could not test", which must not be reported as "tested and
    # found nothing": refuse a run that selected nothing else, and say so otherwise.
    not_run = tuple(m for m in params.methods if m == "barcode" and majors[m] is None)
    if not_run:
        reason = (
            "fewer than two typed clades carry enough characteristic markers"
            if lineage_map else
            "it needs typed references (a lineage map: --lineage-map, or a lineages.tsv "
            "beside the output or the MSA) and none was found"
        )
        if len(not_run) == len(params.methods):
            raise UserInputError(
                f"The barcode caller could not run: {reason}. No other caller was "
                "selected, so this scan could not test for recombination. Supply typed "
                "references or choose another --method."
            )
        logger.warning(
            "The barcode caller could not run (%s); it is reported as 'not run', not as "
            "a negative.", reason,
        )
    n_ran = len(params.methods) - len(not_run)

    major_parent, per_major = reconcile_major(majors, window_wins=analysis_bp.winners_with_ties)
```

In `src/tessera/recomb/run.py`, replace:

```python
    min_agree = max(1, min(params.min_methods, len(params.methods)))
```

with:

```python
    min_agree = max(1, min(params.min_methods, n_ran))
```

In `src/tessera/recomb/run.py`, replace:

```python
    if excluded_siblings:
        provenance["excluded siblings (query's own lineage)"] = ", ".join(
```

with:

```python
    if not_run:
        provenance["callers not run"] = ", ".join(not_run) + " (needs typed references)"
    if excluded_siblings:
        provenance["excluded siblings (query's own lineage)"] = ", ".join(
```

In `src/tessera/recomb/run.py`, replace:

```python
            methods_run=params.methods, method_breakdown=method_breakdown, per_major=per_major,
```

with:

```python
            methods_run=params.methods, method_breakdown=method_breakdown, per_major=per_major,
            methods_not_run=not_run,
```

In `src/tessera/recomb/report_context.py`, replace:

```python
    methods_run: tuple[str, ...] = ()
```

with:

```python
    methods_run: tuple[str, ...] = ()
    # Selected callers that could not run (barcode on an untyped panel). They stay in
    # ``methods_run`` so the method table keeps a column for them, marked "not run".
    methods_not_run: tuple[str, ...] = ()
```

In `src/tessera/recomb/report_text.py`, replace:

```python
def write_methods_tsv(
    breakdown: list[dict], methods_run: tuple[str, ...], output_dir: Path,
    logger: logging.Logger,
) -> None:
    """Write the per-region x per-method agreement matrix (the ensemble breakdown)."""
```

with:

```python
def write_methods_tsv(
    breakdown: list[dict], methods_run: tuple[str, ...], output_dir: Path,
    logger: logging.Logger, methods_not_run: tuple[str, ...] = (),
) -> None:
    """Write the per-region x per-method agreement matrix (the ensemble breakdown).

    A cell is ``yes`` / ``no`` for a caller that ran, and ``not run`` for one that was
    selected but could not run -- a ``no`` there would read as a caller that looked and
    found nothing.
    """
```

In `src/tessera/recomb/report_text.py`, replace:

```python
            cells = [("yes" if m in called else "no") for m in methods_run]
```

with:

```python
            cells = [
                "not run" if m in methods_not_run else ("yes" if m in called else "no")
                for m in methods_run
            ]
```

In `src/tessera/recomb/report.py`, replace:

```python
        write_methods_tsv(ctx.method_breakdown, ctx.methods_run, output_dir, logger)
```

with:

```python
        write_methods_tsv(
            ctx.method_breakdown, ctx.methods_run, output_dir, logger, ctx.methods_not_run
        )
```

In `src/tessera/recomb/report_html.py`, replace:

```python
def _method_comparison_html(
    breakdown: list[dict], methods_run: tuple[str, ...], per_major: dict[str, str],
    lineage_map: LineageMap | None = None,
) -> str:
    """A compact region x method agreement matrix; omitted for a single-method run."""
    if len(methods_run) < 2:
        return ""
```

with:

```python
def _method_comparison_html(
    breakdown: list[dict], methods_run: tuple[str, ...], per_major: dict[str, str],
    lineage_map: LineageMap | None = None, methods_not_run: tuple[str, ...] = (),
) -> str:
    """A compact region x method agreement matrix; omitted for a single-method run."""
    if len(methods_run) < 2:
        return ""
    not_run_note = ""
    if methods_not_run:
        names = ", ".join(html.escape(m) for m in methods_not_run)
        not_run_note = (
            f'<p class="cap"><strong>Not run:</strong> {names}. The caller needs typed '
            f'references and none were available, so it tested nothing; its column below '
            f'is not a negative result.</p>'
        )
```

In `src/tessera/recomb/report_html.py`, replace:

```python
        f'Backbone per method &mdash; {majors}.</p>'
    )
    if not breakdown:
```

with:

```python
        f'Backbone per method &mdash; {majors}.</p>{not_run_note}'
    )
    if not breakdown:
```

In `src/tessera/recomb/report_html.py`, replace:

```python
        cells = "".join(
            f'<td class="num">{"&check;" if m in b["per_method_support"] else "&middot;"}</td>'
            for m in methods_run
        )
```

with:

```python
        cells = "".join(
            '<td class="num">not run</td>' if m in methods_not_run else
            f'<td class="num">{"&check;" if m in b["per_method_support"] else "&middot;"}</td>'
            for m in methods_run
        )
```

In `src/tessera/recomb/report_html.py`, replace:

```python
def _method_section(
    method_breakdown: list[dict] | None, methods_run: tuple[str, ...],
    per_major: dict[str, str] | None, lineage_map: LineageMap | None,
) -> str:
    """Wrap the method-comparison table in a report section (empty for one method)."""
    body = _method_comparison_html(
        method_breakdown or [], methods_run, per_major or {}, lineage_map
    )
```

with:

```python
def _method_section(
    method_breakdown: list[dict] | None, methods_run: tuple[str, ...],
    per_major: dict[str, str] | None, lineage_map: LineageMap | None,
    methods_not_run: tuple[str, ...] = (),
) -> str:
    """Wrap the method-comparison table in a report section (empty for one method)."""
    body = _method_comparison_html(
        method_breakdown or [], methods_run, per_major or {}, lineage_map, methods_not_run
    )
```

In `src/tessera/recomb/report_html.py`, replace:

```python
        f"{_method_section(ctx.method_breakdown, ctx.methods_run, ctx.per_major, lineage_map)}"
```

with:

```python
        f"{method_section}"
```

In `src/tessera/recomb/report_html.py`, replace:

```python
    doc = (
        '<!DOCTYPE html>
```

with:

```python
    method_section = _method_section(
        ctx.method_breakdown, ctx.methods_run, ctx.per_major, lineage_map,
        ctx.methods_not_run,
    )

    doc = (
        '<!DOCTYPE html>
```

In `docs/detection-methods.md`, replace:

```markdown
the genome-level callers at low divergence. It is silent on untyped panels, and opt-in
(needs typed references). It composes with `--pool-consensus` but needs only a typed panel.
```

with:

```markdown
the genome-level callers at low divergence. It is opt-in and needs typed references. On an
untyped panel it cannot run, which is not the same as finding nothing: in an ensemble it is
logged and shown as `not run` in the method comparison, and a run that selected only
`barcode` is refused. It composes with `--pool-consensus` but needs only a typed panel.
```

In `docs/detection-methods.md`, replace:

```markdown
one row per region with a Y/n per method and the parent-free flag |
```

with:

```markdown
one row per region with `yes` / `no` per method (`not run` for a selected caller that could not run) and the parent-free flag |
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_barcode.py tests/unit/test_ensemble.py tests/unit/test_run.py tests/unit/test_output_contract.py tests/unit/test_report_typed.py -q`
Expected: PASS -- `53 passed`

Run: `ruff check src tests validation && mypy src`
Expected: `All checks passed!` and `Success: no issues found in 79 source files`

- [ ] **Step 5: Commit**

```bash
git add tests/unit/test_barcode.py src/tessera/recomb/barcode.py src/tessera/recomb/run.py src/tessera/recomb/report_context.py src/tessera/recomb/report_text.py src/tessera/recomb/report.py src/tessera/recomb/report_html.py docs/detection-methods.md
git commit -m "Report a barcode caller that could not run as not run" -m "Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 8: PHI: use the run's alpha and say when it is not testable (B7)

Two defects in the parent-free signal. (i) The report states and applies alpha 0.05 whatever `--alpha` is, while the per-region flag uses `--alpha`. (ii) With `z` informative sites and a window of at least `z - 1` ranks, every pair of sites is inside the window, the statistic cannot change under permutation, and p is 1 for any data; this is reported as "no significant signal".

The default window is not changed here (that is plan C, item C5).

**Files:**
- Create: `tests/unit/test_report_signal.py`
- Modify: `src/tessera/recomb/diagnostics.py`
- Modify: `src/tessera/recomb/run.py`
- Modify: `src/tessera/recomb/report_context.py`
- Modify: `src/tessera/recomb/report_text.py`
- Modify: `src/tessera/recomb/report_html.py`
- Modify: `docs/detection-methods.md`
- Modify: `example_data/README.md`
- Test (extend): `tests/unit/test_diagnostics.py`
- Test (extend): `tests/unit/test_harness_scoring.py`

**Interfaces:**
- Consumes: `ReportContext` construction in `run_recomb` (the line Task 7 wrote); `corroborating_intervals`, which already returns `[]` for `phi_p is None`.
- Produces: `RecombinationSignal.phi_p: float | None` (`None` = not testable); `ReportContext.alpha: float = 0.05`.

- [ ] **Step 1: Write the failing tests**

Append to the end of `tests/unit/test_diagnostics.py`:

```python
# --- "not testable" is not "no signal" --------------------------------------

def test_phi_is_not_testable_when_every_site_pair_is_inside_the_window() -> None:
    """With z informative sites and a window of at least z - 1 ranks, the statistic
    averages over every pair of sites, so reordering the sites cannot change it and the
    permutation p-value is 1 whatever the data. That was reported as "no signal"."""
    rows = _block_alignment(recombinant=True)  # 24 informative columns
    signal = recombination_signal(rows, "s0", lambda c: c, window=100, seed=1)
    assert signal is not None
    assert signal.n_informative == 24
    assert signal.phi_p is None  # was 1.0
    assert signal.rmin >= 1  # Rmin does not depend on the window
    assert corroborating_intervals(signal, alpha=0.05) == []


def test_phi_becomes_testable_one_rank_below_the_site_count() -> None:
    rows = _block_alignment(recombinant=True)  # z = 24, so z - 1 = 23
    at_limit = recombination_signal(rows, "s0", lambda c: c, window=23, seed=1)
    below = recombination_signal(rows, "s0", lambda c: c, window=22, seed=1)
    assert at_limit is not None and below is not None
    assert at_limit.phi_p is None
    assert below.phi_p is not None
```

Create `tests/unit/test_report_signal.py`:

```python
"""The parent-free signal section states the alpha the run used, and says "not
testable" when the PHI test could not have rejected."""

from __future__ import annotations

import json
import logging

from tessera.recomb.diagnostics import RecombinationSignal
from tessera.recomb.report_html import _signal_html
from tessera.recomb.report_text import write_profile_tsv
from tessera.recomb.run import RecombParams, run_recomb


def _signal(phi_p: float | None, n_informative: int = 500) -> RecombinationSignal:
    return RecombinationSignal(
        n_informative=n_informative, phi_p=phi_p, phi_observed=0.1, phi_window=100,
        rmin=2, rmin_intervals=[(100, 200), (900, 1000)], profile=[(150, 150, 0.2)],
    )


def test_signal_section_says_not_testable() -> None:
    html = _signal_html(_signal(None, n_informative=54), alpha=0.05)
    assert "not testable" in html
    assert "p =" not in html
    assert "no significant recombination signal" not in html
    assert "54 informative sites" in html
    assert "--phi-window" in html  # tells the reader what to change


def test_profile_tsv_header_marks_an_untestable_phi(tmp_path) -> None:
    write_profile_tsv(_signal(None, n_informative=54), tmp_path, logging.getLogger("tessera"))
    header = (tmp_path / "recombination_profile.tsv").read_text().splitlines()[0]
    assert header.split("\t")[1] == "NA"
    assert "not testable" in header


def test_run_reports_an_untestable_phi(example_data, tmp_path, logger) -> None:
    # cryptic_insert has 54 informative sites: fewer than the default window of 100.
    out = tmp_path / "cryptic"
    run_recomb(
        RecombParams(msa=example_data / "cryptic_insert.msa.fasta", output=out,
                     query="query", plot_format="png"),
        logger,
    )
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["recombination signal (PHI)"].startswith("not testable")
    assert "not testable" in (out / "report.html").read_text()


def test_report_states_the_alpha_the_run_used(example_data, tmp_path, logger) -> None:
    """The section always printed "alpha 0.05" and judged significance at 0.05, while
    the per-region PHI flag used --alpha."""
    out = tmp_path / "divergent"
    run_recomb(
        RecombParams(msa=example_data / "divergent.msa.fasta", output=out, query="query",
                     window_size=300, window_step=30, plot_format="png", alpha=0.2,
                     phi_window=20),
        logger,
    )
    report = (out / "report.html").read_text()
    assert "(alpha 0.2;" in report
    assert "(alpha 0.05;" not in report
```

Append to the end of `tests/unit/test_harness_scoring.py`:

```python
def test_parse_signal_reads_an_untestable_phi_header(tmp_path):
    """The profile header carries NA when the PHI test could not have rejected; the
    harness must still find the Rmin that follows it."""
    import logging

    from tessera.recomb.diagnostics import RecombinationSignal
    from tessera.recomb.report_text import write_profile_tsv

    signal = RecombinationSignal(
        n_informative=54, phi_p=None, phi_observed=0.1, phi_window=100, rmin=2,
    )
    write_profile_tsv(signal, tmp_path, logging.getLogger("tessera"))
    assert rh.parse_signal(tmp_path / "recombination_profile.tsv") == ("NA", "2")
```

`validation/run_benchmark.py` and `run_coalescent_benchmark.py` already treat a `None` p-value as "untestable" and count it separately, so their denominators change for small alignments; that is the intended reading.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_diagnostics.py tests/unit/test_report_signal.py tests/unit/test_harness_scoring.py -q`
Expected: FAIL --

```text
FAILED tests/unit/test_diagnostics.py::test_phi_is_not_testable_when_every_site_pair_is_inside_the_window
FAILED tests/unit/test_diagnostics.py::test_phi_becomes_testable_one_rank_below_the_site_count
FAILED tests/unit/test_report_signal.py::test_signal_section_says_not_testable
FAILED tests/unit/test_report_signal.py::test_profile_tsv_header_marks_an_untestable_phi
FAILED tests/unit/test_report_signal.py::test_run_reports_an_untestable_phi
FAILED tests/unit/test_report_signal.py::test_report_states_the_alpha_the_run_used
FAILED tests/unit/test_harness_scoring.py::test_parse_signal_reads_an_untestable_phi_header
7 failed, 26 passed in 3.02s
```

- [ ] **Step 3: Implement**

In `src/tessera/recomb/diagnostics.py`, replace:

```python
    phi_p: float  # PHI permutation p-value (one-sided; small = recombination)
```

with:

```python
    # PHI permutation p-value (one-sided; small = recombination). None when the test
    # could not have rejected: see recombination_signal.
    phi_p: float | None
```

In `src/tessera/recomb/diagnostics.py`, replace:

```python
    informative columns (too little variation to test). The PHI p-value is reported
    as-is; the significance threshold is a reporting concern, applied downstream.
    ``query_label`` is accepted for interface symmetry; the statistics use every
    sequence in ``rows``.
    """
```

with:

```python
    informative columns (too little variation to test). The PHI p-value is reported
    as-is; the significance threshold is a reporting concern, applied downstream.
    ``query_label`` is accepted for interface symmetry; the statistics use every
    sequence in ``rows``.

    ``phi_p`` is ``None`` when the test is not testable at this window: with ``z``
    informative columns and ``window >= z - 1`` every pair of columns falls inside the
    window, the statistic is the mean over all pairs, and no reordering of the columns
    can change it -- the permutation p-value would be 1 whatever the data. Rmin, the
    intervals and the profile do not depend on the permutation and are still returned.
    """
```

In `src/tessera/recomb/diagnostics.py`, replace:

```python
    p, observed = phi_pvalue(incompatible, window, seed=seed)
```

with:

```python
    p: float | None
    if z - 1 > window:
        p, observed = phi_pvalue(incompatible, window, seed=seed)
    else:
        p, observed = None, phi(incompatible, window)
```

In `src/tessera/recomb/run.py`, replace:

```python
        if signal is not None:
            logger.info(
                "Recombination signal (parent-free): PHI p=%.4g, Rmin=%d (%d informative "
                "sites).", signal.phi_p, signal.rmin, signal.n_informative,
            )
```

with:

```python
        if signal is not None and signal.phi_p is None:
            logger.info(
                "Recombination signal (parent-free): PHI not testable (%d informative "
                "site(s) do not exceed the window of %d ranks; lower --phi-window), Rmin=%d.",
                signal.n_informative, signal.phi_window, signal.rmin,
            )
        elif signal is not None:
            logger.info(
                "Recombination signal (parent-free): PHI p=%.4g, Rmin=%d (%d informative "
                "sites).", signal.phi_p, signal.rmin, signal.n_informative,
            )
```

In `src/tessera/recomb/run.py`, replace:

```python
        provenance["recombination signal (PHI)"] = (
            f"p={signal.phi_p:.4g} ({signal.n_informative} informative sites, "
            f"window {signal.phi_window})"
        )
```

with:

```python
        phi_text = "not testable" if signal.phi_p is None else f"p={signal.phi_p:.4g}"
        provenance["recombination signal (PHI)"] = (
            f"{phi_text} ({signal.n_informative} informative sites, "
            f"window {signal.phi_window})"
        )
```

In `src/tessera/recomb/run.py`, replace:

```python
            methods_not_run=not_run,
```

with:

```python
            methods_not_run=not_run, alpha=params.alpha,
```

In `src/tessera/recomb/report_context.py`, replace:

```python
    signal: RecombinationSignal | None = None
```

with:

```python
    signal: RecombinationSignal | None = None
    # The run's significance level, so the report judges the PHI p-value at the same
    # alpha the per-region corroboration used.
    alpha: float = 0.05
```

In `src/tessera/recomb/report_text.py`, replace:

```python
        fo.write(
            f"# PHI p-value\t{signal.phi_p:.4g}\t(window {signal.phi_window} "
            f"informative sites, {signal.n_informative} sites, Rmin {signal.rmin})\n"
        )
```

with:

```python
        # "NA" when the PHI test could not have rejected (too few informative sites
        # for the window). The Rmin stays last in the note: the harness reads it there.
        phi_p = "NA" if signal.phi_p is None else f"{signal.phi_p:.4g}"
        note = "" if signal.phi_p is not None else "not testable at this window; "
        fo.write(
            f"# PHI p-value\t{phi_p}\t({note}window {signal.phi_window} "
            f"informative sites, {signal.n_informative} sites, Rmin {signal.rmin})\n"
        )
```

In `src/tessera/recomb/report_html.py`, replace:

```python
    significant = signal.phi_p < alpha
    verdict = (
        '<strong>significant recombination signal</strong>' if significant
        else 'no significant recombination signal'
    )
```

with:

```python
    if signal.phi_p is None:
        # Every pair of informative sites is inside one window, so the permutation test
        # cannot reject whatever the data. Say so instead of printing p = 1.
        p_cell = "not testable"
        verdict = (
            f'{signal.n_informative} informative sites do not exceed the window of '
            f'{signal.phi_window} ranks, so the permutation test cannot reject here; '
            f'this is not evidence against recombination. Lower '
            f'<span class="mono">--phi-window</span> to test'
        )
    else:
        p_cell = f"p = {signal.phi_p:.4g}"
        verdict = (
            '<strong>significant recombination signal</strong>' if signal.phi_p < alpha
            else 'no significant recombination signal'
        ) + f' (alpha {alpha:g}; {signal.n_informative} informative sites)'
```

In `src/tessera/recomb/report_html.py`, replace:

```python
        f'<tr><td class="lbl">PHI test</td><td class="num strong">p = {signal.phi_p:.4g}</td>'
        f'<td class="lbl">{verdict} (alpha {alpha:g}; {signal.n_informative} informative '
        f'sites)</td></tr>'
```

with:

```python
        f'<tr><td class="lbl">PHI test</td><td class="num strong">{p_cell}</td>'
        f'<td class="lbl">{verdict}</td></tr>'
```

In `src/tessera/recomb/report_html.py`, replace:

```python
        f"{_signal_html(ctx.signal)}</section>"
```

with:

```python
        f"{_signal_html(ctx.signal, ctx.alpha)}</section>"
```

In `docs/detection-methods.md`, replace:

```markdown
parent-free outputs. The diagnostic runs for every `--method`; disable with
`--no-phi`, or widen its window with `--phi-window`.
```

with:

```markdown
parent-free outputs. The diagnostic runs for every `--method`; disable with
`--no-phi`, or widen its window with `--phi-window`.

The PHI p-value is judged at the run's `--alpha`, in the report and for the per-region
flag alike. The test is **not testable** when the alignment has too few informative sites
for the window (`--phi-window`, default 100 site ranks; at least window + 2 sites are
needed): every pair of sites then falls inside one window, so reordering the sites cannot
change the statistic and the permutation p-value would be 1 whatever the data. Tessera
reports this as `not testable` (`NA` in the `recombination_profile.tsv` header) rather
than as a non-significant result, and no region is flagged `parent_free_support`; lower
`--phi-window` to test such an alignment. Just above that limit the test is defined but
has little power: the shipped `divergent` example (159 informative sites) gives p = 1 at
the default window and p = 0.001 at `--phi-window 20`.
```

In `docs/detection-methods.md`, replace:

```markdown
| `recombination_profile.tsv` | Parent-free signal: header with the PHI p-value and Rmin, then
```

with:

```markdown
| `recombination_profile.tsv` | Parent-free signal: header with the PHI p-value (`NA` when not testable) and Rmin, then
```

In `example_data/README.md`, replace:

```markdown
Both runs also report the parent-free PHI / Rmin signal in `recombination_profile.tsv`.
```

with:

```markdown
Both runs also write the parent-free PHI / Rmin signal to `recombination_profile.tsv`.
The cryptic example has 54 informative sites, fewer than the default PHI window of 100, so
its PHI test is reported as `not testable`; add `--phi-window 5` to test it.
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_diagnostics.py tests/unit/test_report_signal.py tests/unit/test_harness_scoring.py tests/unit/test_benchmark_scoring.py tests/unit/test_coalescent_benchmark.py tests/integration -q`
Expected: PASS -- `59 passed, 1 skipped`

Run: `ruff check src tests validation && mypy src`
Expected: `All checks passed!` and `Success: no issues found in 79 source files`

- [ ] **Step 5: Commit**

```bash
git add tests/unit/test_diagnostics.py tests/unit/test_harness_scoring.py tests/unit/test_report_signal.py src/tessera/recomb/diagnostics.py src/tessera/recomb/run.py src/tessera/recomb/report_context.py src/tessera/recomb/report_text.py src/tessera/recomb/report_html.py docs/detection-methods.md example_data/README.md
git commit -m "Judge PHI at the run's alpha and report an untestable test as such" -m "Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 9: Plots: donor-absent bands, the pair plot, and colours (B8)

The plots label a donor-absent region "recombinant: <backbone>" in the backbone's colour; `similarity_pair` shows the two leading window winners, which are the backbone and a near-duplicate of it when the panel holds one; and with `--top-n 1` the donor is outside the colour map and drawn grey.

**Files:**
- Create: `tests/unit/test_report_plots.py`
- Modify: `src/tessera/recomb/report_plots.py`
- Modify: `src/tessera/recomb/report.py`
- Modify: `src/tessera/recomb/report_html.py`
- Modify: `docs/detection-methods.md`

**Interfaces:**
- Consumes: `Region.donor_absent`, `Region.length_bp`.
- Produces: `region_labels(datasets: list[str], regions: list[Region]) -> list[str]` and `ABSENT_LABEL = "donor absent"` in `recomb/report_plots.py`; `pair_datasets(regions: list[Region], ranked: list[str]) -> list[str]` in `recomb/report.py`. Task 10 edits the block containing `pair = pair_datasets(...)`.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_report_plots.py`:

```python
"""The plots label what the regions table says, in the same colours."""

from __future__ import annotations

from pathlib import Path

from matplotlib.figure import Figure

from tessera.recomb.regions import Region
from tessera.recomb.report import pair_datasets
from tessera.recomb.report_plots import (
    GREY,
    _color_map,
    _palette,
    _shade_regions,
    build_interactive_figure,
    region_labels,
)
from tessera.recomb.similarity import compute_similarity

from ..conftest import recombinant_msa


def _region(minor: str, major: str = "A", *, start: int = 2000, end: int = 4000,
            donor_absent: bool = False) -> Region:
    return Region(
        minor_parent=minor, major_parent=major, msa_start=start, msa_end=end,
        query_start=start, query_end=end, n_windows=10,
        mean_sim_minor=0.99, mean_sim_major=0.95, margin=0.04, donor_absent=donor_absent,
    )


# --- colours cover every label a region names ------------------------------

def test_region_labels_appends_region_parents_to_the_plotted_datasets() -> None:
    labels = region_labels(["A"], [_region("B"), _region("other", donor_absent=True)])
    assert labels == ["A", "B"]  # the donor is added; a donor-absent stand-in is not
    assert region_labels(["A", "B"], [_region("B")]) == ["A", "B"]  # no duplicates


def test_donor_outside_the_top_n_keeps_a_colour(tmp_path: Path) -> None:
    """With --top-n 1 only the backbone was in the colour map, so the donor band and
    its label were drawn grey."""
    result = compute_similarity(str(recombinant_msa(tmp_path, recombinant=True)), "query")
    fig = build_interactive_figure(result, ["A"], [_region("B")])
    donor_colour = _color_map(["A", "B"])["B"]
    band = next(a for a in fig.layout.annotations if a.text == "B")
    assert band.font.color == donor_colour
    assert donor_colour != GREY


# --- a donor-absent region is not a recombinant from the backbone ----------

def test_static_plot_does_not_call_a_donor_absent_region_recombinant() -> None:
    """The stand-in label of a donor-absent region is the closest reference, often the
    backbone itself: the legend read "recombinant: <backbone>" in the backbone's colour."""
    ax = Figure().subplots()  # no pyplot state, no display backend needed
    _shade_regions(ax, [_region("A", donor_absent=True)], _palette(["A", "B"]))
    labels = ax.get_legend_handles_labels()[1]
    assert labels == ["donor absent (no close reference)"]


def test_static_plot_still_labels_a_called_donor() -> None:
    ax = Figure().subplots()
    _shade_regions(ax, [_region("B")], _palette(["A", "B"]))
    labels = ax.get_legend_handles_labels()[1]
    assert labels == ["recombinant: B"]


def test_interactive_plot_marks_a_donor_absent_region(tmp_path: Path) -> None:
    result = compute_similarity(str(recombinant_msa(tmp_path, recombinant=True)), "query")
    fig = build_interactive_figure(result, ["A", "B"], [_region("A", donor_absent=True)])
    texts = [a.text for a in fig.layout.annotations]
    assert texts == ["donor absent"]
    assert fig.layout.annotations[0].font.color == GREY


# --- the pair plot is major versus leading minor ---------------------------

def test_pair_plot_shows_the_backbone_and_the_leading_donor() -> None:
    """It showed the top two window winners, which are the backbone and a
    near-duplicate of it when the panel holds one."""
    regions = [_region("variola", "cowpox", start=2000, end=2600),
               _region("camelpox", "cowpox", start=4000, end=4100)]
    assert pair_datasets(regions, ["cowpox", "cowpox2"]) == ["cowpox", "variola"]


def test_pair_plot_falls_back_to_the_window_ranking_without_a_called_donor() -> None:
    assert pair_datasets([], ["cowpox", "cowpox2"]) == ["cowpox", "cowpox2"]
    absent = [_region("cowpox", "cowpox", donor_absent=True)]
    assert pair_datasets(absent, ["cowpox", "cowpox2"]) == ["cowpox", "cowpox2"]
```

The plot file names do not change (`similarity_top{N}` still counts the plotted datasets): colours come from the extended label list, the plotted lines do not.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_report_plots.py -q`
Expected: FAIL --

```text
ImportError: cannot import name 'pair_datasets' from 'tessera.recomb.report'
ERROR tests/unit/test_report_plots.py
1 error in 0.68s
```

- [ ] **Step 3: Implement**

In `src/tessera/recomb/report_plots.py`, replace:

```python
def _palette(datasets: list[str]):
```

with:

```python
def region_labels(datasets: list[str], regions: list[Region]) -> list[str]:
    """``datasets`` followed by any parent of a called region that is not among them.

    The plots draw the top-N datasets, but a region can name a parent outside that
    list (``--top-n 1`` leaves out every donor). Colours are assigned from this
    extended list so a region is never drawn in the fallback grey; the plotted
    datasets come first, so their colours do not depend on which regions were called.
    Donor-absent regions are skipped: their label is a stand-in, not a called parent.
    """
    labels = list(datasets)
    for region in regions:
        if region.donor_absent:
            continue
        for label in (region.major_parent, region.minor_parent):
            if label not in labels:
                labels.append(label)
    return labels


def _palette(datasets: list[str]):
```

In `src/tessera/recomb/report_plots.py`, replace:

```python
def _shade_regions(ax, regions: list[Region], colors: dict) -> None:
    seen: set[str] = set()
    for region in regions:
        color = colors.get(region.minor_parent, "grey")
        label = f"recombinant: {region.minor_parent}" if region.minor_parent not in seen else None
        seen.add(region.minor_parent)
        ax.axvspan(region.msa_start, region.msa_end, color=color, alpha=0.12, label=label)
```

with:

```python
ABSENT_LABEL = "donor absent"


def _shade_regions(ax, regions: list[Region], colors: dict) -> None:
    seen: set[str] = set()
    for region in regions:
        if region.donor_absent:
            # The stand-in label is the closest reference -- often the backbone itself --
            # so naming it here would read as "recombinant from the backbone".
            key, color = ABSENT_LABEL, GREY
            text = f"{ABSENT_LABEL} (no close reference)"
        else:
            key, color = region.minor_parent, colors.get(region.minor_parent, GREY)
            text = f"recombinant: {region.minor_parent}"
        label = text if key not in seen else None
        seen.add(key)
        ax.axvspan(region.msa_start, region.msa_end, color=color, alpha=0.12, label=label)
```

In `src/tessera/recomb/report_plots.py`, replace:

```python
    subset = df.loc[available]
    colors = _palette(available)
```

with:

```python
    subset = df.loc[available]
    colors = _palette(region_labels(available, regions))
```

In `src/tessera/recomb/report_plots.py`, replace:

```python
    colors = _palette([seq1, seq2])
```

with:

```python
    colors = _palette(region_labels([seq1, seq2], regions))
```

In `src/tessera/recomb/report_plots.py`, replace:

```python
    df = result.to_dataframe()
    colors = _color_map(datasets)
    fig = go.Figure()
```

with:

```python
    df = result.to_dataframe()
    colors = _color_map(region_labels(datasets, regions))
    fig = go.Figure()
```

In `src/tessera/recomb/report_plots.py`, replace:

```python
    for region in regions:
        color = colors.get(region.minor_parent, GREY)
        fig.add_vrect(
            x0=region.msa_start, x1=region.msa_end,
            fillcolor=color, opacity=0.12, line_width=0, layer="below",
            annotation_text=("" if region.minor_parent in seen else region.minor_parent),
            annotation_position="top left",
            annotation_font_size=11, annotation_font_color=color,
        )
        seen.add(region.minor_parent)
```

with:

```python
    for region in regions:
        if region.donor_absent:
            key, color = ABSENT_LABEL, GREY
        else:
            key, color = region.minor_parent, colors.get(region.minor_parent, GREY)
        fig.add_vrect(
            x0=region.msa_start, x1=region.msa_end,
            fillcolor=color, opacity=0.12, line_width=0, layer="below",
            annotation_text=("" if key in seen else key),
            annotation_position="top left",
            annotation_font_size=11, annotation_font_color=color,
        )
        seen.add(key)
```

In `src/tessera/recomb/report.py`, replace:

```python
- ``similarity_pair.{fmt}``      static major-vs-minor pairwise plot
```

with:

```python
- ``similarity_pair.{fmt}``      static pairwise plot: the major parent against the
                                 donor of the longest called region (the two leading
                                 window winners when no donor was called)
```

In `src/tessera/recomb/report.py`, replace:

```python
def write_reports(
```

with:

```python
def pair_datasets(regions: list[Region], ranked: list[str]) -> list[str]:
    """The two datasets of the pairwise plot: major parent, then the leading donor.

    The leading donor is the donor of the longest called, donor-present region. Without
    one there is no minor parent to show and the plot falls back to ``ranked`` (the two
    leading window winners) -- which is not "major versus minor" when the panel holds a
    near-duplicate of the backbone, hence the region-based choice whenever possible.
    """
    present = [r for r in regions if not r.donor_absent]
    if not present:
        return ranked
    longest = max(present, key=lambda r: r.length_bp)
    return [longest.major_parent, longest.minor_parent]


def write_reports(
```

In `src/tessera/recomb/report.py`, replace:

```python
    pair = rank_datasets(ranking, 2)
```

with:

```python
    pair = pair_datasets(regions, rank_datasets(ranking, 2))
```

In `src/tessera/recomb/report.py`, replace:

```python
    "write_html_report", "write_reports",
]
```

with:

```python
    "write_html_report", "write_reports", "pair_datasets",
]
```

In `src/tessera/recomb/report_html.py`, replace:

```python
from .report_plots import GREY, _color_map, build_interactive_figure
```

with:

```python
from .report_plots import GREY, _color_map, build_interactive_figure, region_labels
```

In `src/tessera/recomb/report_html.py`, replace:

```python
    colors = _color_map(datasets)
    s = _summary(result, regions, datasets)
```

with:

```python
    # Colour every label a region names, not only the top-N: with --top-n 1 the donor
    # is outside the plotted datasets and its swatch and mosaic segment were grey.
    colors = _color_map(region_labels(datasets, regions))
    s = _summary(result, regions, datasets)
```

In `docs/detection-methods.md`, replace:

```markdown
| `similarity_pair.{fmt}` | Static plot of the major vs leading minor parent, region shaded |
```

with:

```markdown
| `similarity_pair.{fmt}` | Static plot of the major parent against the donor of the longest called region, regions shaded (the two leading window winners when no donor was called) |
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_report_plots.py tests/unit/test_report_typed.py tests/unit/test_output_contract.py tests/integration -q`
Expected: PASS -- `47 passed, 1 skipped`

Run: `ruff check src tests validation && mypy src`
Expected: `All checks passed!` and `Success: no issues found in 79 source files`

- [ ] **Step 5: Commit**

```bash
git add tests/unit/test_report_plots.py src/tessera/recomb/report_plots.py src/tessera/recomb/report.py src/tessera/recomb/report_html.py docs/detection-methods.md
git commit -m "Label donor-absent regions and the pair plot for what they show" -m "Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 10: Replace stale report, help and documentation text (B9)

Text that describes an earlier version of the tool: the report footer hard-codes `.pdf` plots and lists files that were not written; the methods paragraph and references describe a two-caller ensemble; `--method` help says "all but the legacy heuristic" and omits `geneconv`; `docs/reference-panels.md` describes a fetch cap that no longer exists; `docs/detection-methods.md` does not describe lineage clustering.

The clustering section documents a limitation found by the audit (spec item C1: on a panel below about 1.5 % divergence all references pool and the HMM cannot call). It is described as it behaves today; plan C revises it.

**Files:**
- Create: `tests/unit/test_report_stale_text.py`
- Modify: `src/tessera/recomb/report.py`
- Modify: `src/tessera/recomb/report_html.py`
- Modify: `src/tessera/recomb/report_assets.py`
- Modify: `src/tessera/cli/cmd_recomb.py`
- Modify: `src/tessera/cli/cmd_detect.py`
- Modify: `src/tessera/cli/cmd_fill_references.py`
- Modify: `docs/detection-methods.md`
- Modify: `docs/reference-panels.md`
- Modify: `example_data/README.md`

**Interfaces:**
- Consumes: `pair_datasets` (Task 9) and the `write_methods_tsv` call (Task 7) inside `write_reports`; the plot functions returning `Path | None`.
- Produces: `write_html_report(..., ctx, companion_files: list[str] | None = None)` and `_footer_html(provenance, companion_files=None)`.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_report_stale_text.py`:

```python
"""Report and help text that describes what the run did, not an earlier version of it."""

from __future__ import annotations

import re
from pathlib import Path

from typer.testing import CliRunner

from tessera.cli.main import app
from tessera.recomb.report_assets import _REFERENCES
from tessera.recomb.run import RecombParams, run_recomb


def _footer(out: Path) -> str:
    report = (out / "report.html").read_text()
    return re.search(r"<footer>.*?</footer>", report, flags=re.S).group(0)


def _run(example_data: Path, out: Path, logger, **kwargs) -> Path:
    run_recomb(
        RecombParams(msa=example_data / "divergent.msa.fasta", output=out, query="query",
                     window_size=300, window_step=30, **kwargs),
        logger,
    )
    return out


def test_footer_lists_the_plots_in_the_format_they_were_written(
    example_data, tmp_path, logger
) -> None:
    out = tmp_path / "out"
    out.mkdir()
    (out / "similarity_pair.pdf").write_text("left over from an earlier run")
    footer = _footer(_run(example_data, out, logger, plot_format="png"))
    assert "similarity_pair.png" in footer
    assert "similarity_top3.png" in footer
    # Was hard-coded; and a file this run did not write is not listed even if present.
    assert ".pdf" not in footer


def test_footer_lists_only_files_that_exist(example_data, tmp_path, logger) -> None:
    out = _run(example_data, tmp_path / "out", logger, plot_format="png",
               methods=("hmm",), phi=False)
    footer = _footer(out)
    listed = re.findall(r"<code>(.*?)</code>", footer)
    assert listed, "the footer should list the companion files"
    for name in listed:
        assert (out / name).exists(), f"{name} is listed but was not written"
    # A single-caller run writes no method comparison; --no-phi writes no profile.
    assert "recombination_methods.tsv" not in listed
    assert "recombination_profile.tsv" not in listed
    assert "run_provenance.json" in listed


def test_methods_paragraph_describes_the_default_ensemble(
    example_data, tmp_path, logger
) -> None:
    report = (_run(example_data, tmp_path / "out", logger, plot_format="png")
              / "report.html").read_text()
    methods = re.search(r'<details class="methods">.*?</details>', report, flags=re.S).group(0)
    assert "(HMM and the 3SEQ triplet test)" not in methods  # the two-caller description
    for caller in ("3SEQ", "MaxChi", "Bootscan"):
        assert caller in methods


def test_references_cite_every_default_caller() -> None:
    titles = " ".join(title for title, _ in _REFERENCES)
    assert "MaxChi" in titles
    assert "Bootscan" in titles


def test_method_help_lists_every_caller_and_the_real_default() -> None:
    result = CliRunner().invoke(app, ["recomb", "--help"], env={"COLUMNS": "200"})
    assert result.exit_code == 0
    text = " ".join(result.output.replace("│", " ").split())
    assert "geneconv" in text  # run by `all`, and was missing from the list
    assert "all but the legacy heuristic" not in text  # geneconv and barcode are opt-in too
```

Not changed here: the comment `# cap a broad NCBI Virus fetch` at `src/tessera/discover/iterate.py:89` is also stale, but plan A rewrites that file; leave it to plan A to avoid a conflict.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_report_stale_text.py -q`
Expected: FAIL --

```text
FAILED tests/unit/test_report_stale_text.py::test_footer_lists_the_plots_in_the_format_they_were_written
FAILED tests/unit/test_report_stale_text.py::test_footer_lists_only_files_that_exist
FAILED tests/unit/test_report_stale_text.py::test_methods_paragraph_describes_the_default_ensemble
FAILED tests/unit/test_report_stale_text.py::test_references_cite_every_default_caller
FAILED tests/unit/test_report_stale_text.py::test_method_help_lists_every_caller_and_the_real_default
5 failed in 3.87s
```

- [ ] **Step 3: Implement**

In `src/tessera/recomb/report.py`, replace:

```python
    write_windows_tsv(result, per_window_winners, output_dir, logger)
    write_stats_tsv(analysis, output_dir, logger)
    write_winners_tsv(analysis, output_dir, logger)
    write_regions_tsv(regions, output_dir, logger)
    write_coverage_tsv(gaps, ctx.coverage_threshold, output_dir, logger)
    if ctx.signal is not None:
        write_profile_tsv(ctx.signal, output_dir, logger)
    if len(ctx.methods_run) > 1 and ctx.method_breakdown is not None:
        write_methods_tsv(
            ctx.method_breakdown, ctx.methods_run, output_dir, logger, ctx.methods_not_run
        )
```

with:

```python
    # The files written beside the report, in the order the footer lists them. Built as
    # they are written, so the footer never names a file that does not exist or a plot
    # in a format that was not asked for.
    companions: list[str] = []

    write_regions_tsv(regions, output_dir, logger)
    companions.append("recombination_regions.tsv")
    if len(ctx.methods_run) > 1 and ctx.method_breakdown is not None:
        write_methods_tsv(
            ctx.method_breakdown, ctx.methods_run, output_dir, logger, ctx.methods_not_run
        )
        companions.append("recombination_methods.tsv")
    write_coverage_tsv(gaps, ctx.coverage_threshold, output_dir, logger)
    companions.append("coverage_gaps.tsv")
    if ctx.signal is not None:
        write_profile_tsv(ctx.signal, output_dir, logger)
        companions.append("recombination_profile.tsv")
    write_winners_tsv(analysis, output_dir, logger)
    write_stats_tsv(analysis, output_dir, logger)
    write_windows_tsv(result, per_window_winners, output_dir, logger)
    companions += ["window_winners.tsv", "similarity_stats.tsv", "similarity_windows.tsv"]
```

In `src/tessera/recomb/report.py`, replace:

```python
    plot_top_n(result, top_datasets, regions, output_dir, plot_format, logger)
    if ctx.site_result is not None:
        write_site_windows_tsv(
            ctx.site_result, winners_per_window(ctx.site_result), output_dir, logger
        )
        plot_top_n(
            ctx.site_result, top_datasets, regions, output_dir, plot_format, logger,
            ylabel="Identity to query at informative sites",
            title=f"Identity at informative sites, query {result.query}",
            stem="informative_sites_top",
        )

    pair = pair_datasets(regions, rank_datasets(ranking, 2))
    plot_pairwise(result, pair, regions, output_dir, plot_format, logger)

    write_html_report(
        result, analysis, regions, top_datasets, provenance, output_dir, logger, ctx
    )
```

with:

```python
    plots = [plot_top_n(result, top_datasets, regions, output_dir, plot_format, logger)]
    if ctx.site_result is not None:
        write_site_windows_tsv(
            ctx.site_result, winners_per_window(ctx.site_result), output_dir, logger
        )
        companions.append("informative_site_windows.tsv")
        plots.append(plot_top_n(
            ctx.site_result, top_datasets, regions, output_dir, plot_format, logger,
            ylabel="Identity to query at informative sites",
            title=f"Identity at informative sites, query {result.query}",
            stem="informative_sites_top",
        ))

    pair = pair_datasets(regions, rank_datasets(ranking, 2))
    plots.append(plot_pairwise(result, pair, regions, output_dir, plot_format, logger))
    # A plot function returns None when it had nothing to draw and wrote no file.
    companions += [path.name for path in plots if path is not None]
    if (output_dir / "run_provenance.json").exists():
        companions.append("run_provenance.json")

    write_html_report(
        result, analysis, regions, top_datasets, provenance, output_dir, logger, ctx,
        companion_files=companions,
    )
```

In `src/tessera/recomb/report_html.py`, replace:

```python
def _footer_html(provenance: dict[str, str]) -> str:
    files = [
        "recombination_regions.tsv", "recombination_methods.tsv", "coverage_gaps.tsv",
        "recombination_profile.tsv", "window_winners.tsv", "similarity_stats.tsv",
        "similarity_windows.tsv", "similarity_top*.pdf", "similarity_pair.pdf",
    ]
    flist = ", ".join(f"<code>{f}</code>" for f in files)
    ver = html.escape(provenance.get("tessera version", ""))
    date = html.escape(provenance.get("date (UTC)", ""))
    return (
        f'<footer><div>Generated by Tessera <span class="mono">{ver}</span> &middot; '
        f'<span class="mono">{date}</span> UTC</div>'
        f'<div>Companion files in this folder: {flist}.</div></footer>'
    )
```

with:

```python
def _footer_html(provenance: dict[str, str], companion_files: list[str] | None = None) -> str:
    """The version line and the companion files that were written with this report.

    ``companion_files`` comes from the writer that produced them. With none given the
    sentence is left out: a list is only worth printing if it is the list of this run.
    """
    ver = html.escape(provenance.get("tessera version", ""))
    date = html.escape(provenance.get("date (UTC)", ""))
    companions = ""
    if companion_files:
        flist = ", ".join(f"<code>{html.escape(f)}</code>" for f in companion_files)
        companions = f"<div>Companion files in this folder: {flist}.</div>"
    return (
        f'<footer><div>Generated by Tessera <span class="mono">{ver}</span> &middot; '
        f'<span class="mono">{date}</span> UTC</div>{companions}</footer>'
    )
```

In `src/tessera/recomb/report_html.py`, replace:

```python
        '<p class="cap">Tessera segments the query against the reference panel with an HMM '
        '(jpHMM-style) and reports a region only when its donor beats the major parent on the '
        'sites that distinguish them (a sign test on discordant sites, immune to window '
        'overlap; Benjamini-Hochberg FDR across segments), with a posterior breakpoint '
        'interval. By default it runs an ensemble of callers (HMM and the 3SEQ triplet '
        'test) and merges their regions into one consensus, so a region found by more than '
        'one method is flagged as agreeing and treated as higher confidence (see the '
        'caller line under Run parameters). It remains an indicative screen, not a full '
        'phylogenetic test (e.g. GARD) -- confirm strong candidates.</p>'
```

with:

```python
        '<p class="cap">Tessera runs one or more region callers on the same alignment and '
        'merges their regions into one consensus; a region found by more than one caller '
        'is flagged as agreeing and treated as higher confidence. The callers that ran '
        'for this report are named in the caller line under Run parameters. The default '
        'ensemble is four callers: an HMM segmentation of the query against the reference '
        'panel (jpHMM-style), which reports a segment only when its donor beats the major '
        'parent on the sites that distinguish them (a one-sided sign test on discordant '
        'sites, with a posterior breakpoint interval); the 3SEQ triplet test; MaxChi; and '
        'Bootscan. Each caller applies Benjamini-Hochberg correction within its own '
        'candidates; nothing is corrected across callers or across the genome. This is an '
        'indicative screen, not a full phylogenetic test (e.g. GARD) -- confirm strong '
        'candidates.</p>'
```

In `src/tessera/recomb/report_html.py`, replace:

```python
    ctx: ReportContext,
) -> Path:
    """Write a single self-contained ``report.html``."""
```

with:

```python
    ctx: ReportContext,
    companion_files: list[str] | None = None,
) -> Path:
    """Write a single self-contained ``report.html``.

    ``companion_files`` are the names of the files written beside it, listed in the
    footer; ``write_reports`` supplies them.
    """
```

In `src/tessera/recomb/report_html.py`, replace:

```python
        f"{_footer_html(provenance)}"
```

with:

```python
        f"{_footer_html(provenance, companion_files)}"
```

In `src/tessera/recomb/report_assets.py`, replace:

```python
    ("Support",
     "The share of distinguishing (discordant) sites -- where the query matches one "
     "candidate parent but not the other -- that favour the donor. 0.5 = no "
     "preference, 1.0 = every distinguishing site favours the donor."),
    ("q-value",
     "The sign-test p-value after Benjamini-Hochberg correction across all candidate "
     "segments (false-discovery-rate control). A region is reported when q <= alpha."),
```

with:

```python
    ("Support",
     "The supporting statistic of the test behind the region; its meaning depends on "
     "the caller, and the statistic column of recombination_regions.tsv names it. For "
     "the HMM it is the share of distinguishing (discordant) sites -- where the query "
     "matches one candidate parent but not the other -- that favour the donor "
     "(0.5 = no preference, 1.0 = every distinguishing site favours the donor)."),
    ("q-value",
     "The p-value of the test behind the region after Benjamini-Hochberg correction "
     "across that caller's own candidates (false-discovery-rate control within one "
     "caller). A region is reported when q <= alpha. For a region several callers "
     "found, the most significant caller's value is shown."),
```

In `src/tessera/recomb/report_assets.py`, replace:

```python
    ("3SEQ triplet test",
     "Boni MF, Posada D, Feldman MW (2007). An exact nonparametric method for inferring "
     "mosaic structure in sequence triplets. Genetics 176(2):1035-1047."),
```

with:

```python
    ("3SEQ triplet test",
     "Boni MF, Posada D, Feldman MW (2007). An exact nonparametric method for inferring "
     "mosaic structure in sequence triplets. Genetics 176(2):1035-1047."),
    ("MaxChi",
     "Maynard Smith J (1992). Analyzing the mosaic structure of genes. Journal of "
     "Molecular Evolution 34(2):126-129."),
    ("Bootscan",
     "Salminen MO, Carr JK, Burke DS, McCutchan FE (1995). Identification of breakpoints "
     "in intergenotypic recombinants of HIV type 1 by bootscanning. AIDS Research and "
     "Human Retroviruses 11(11):1423-1425."),
```

In `src/tessera/cli/cmd_recomb.py`, replace:

```python
        "the default is hmm,3seq,maxchi,bootscan (all but the legacy heuristic). Callers: "
        "hmm (HMM segmentation + a discordant-site "
        "significance test), 3seq (scan-aware triplet max-drawdown test; strong at low "
        "divergence), maxchi (chi-square triplet test, complementary to 3seq), bootscan "
        "(distance + bootstrap support for the closest parent), barcode (clade-marker "
        "lineage attribution; needs typed references), heuristic (legacy margin/merge). "
        "Pass a single name (e.g. --method hmm) for one caller.",
```

with:

```python
        "the default is hmm,3seq,maxchi,bootscan; geneconv, barcode and heuristic are "
        "opt-in, and 'all' runs every one. Callers: "
        "hmm (HMM segmentation + a discordant-site "
        "significance test), 3seq (scan-aware triplet max-drawdown test; strong at low "
        "divergence), maxchi (chi-square triplet test, complementary to 3seq), bootscan "
        "(distance + bootstrap support for the closest parent), geneconv (longest "
        "uninterrupted donor-match run), barcode (clade-marker "
        "lineage attribution; needs typed references), heuristic (legacy margin/merge). "
        "Pass a single name (e.g. --method hmm) for one caller.",
```

In `src/tessera/cli/cmd_detect.py`, replace:

```python
        help="Region caller(s): a comma-separated list of hmm/3seq/maxchi/bootscan/"
        "heuristic, or 'all'. Several run as an ensemble and their regions are merged "
        "(default hmm,3seq,maxchi,bootscan).",
```

with:

```python
        help="Region caller(s): a comma-separated list of hmm/3seq/maxchi/bootscan/"
        "geneconv/barcode/heuristic, or 'all'. Several run as an ensemble and their "
        "regions are merged (default hmm,3seq,maxchi,bootscan).",
```

In `src/tessera/cli/cmd_fill_references.py`, replace:

```python
        help="Region caller(s) for the detection step: a comma-separated list of "
        "hmm/3seq/maxchi/bootscan/heuristic, or 'all'. Several run as an ensemble "
        "(default hmm,3seq,maxchi,bootscan).",
```

with:

```python
        help="Region caller(s) for the detection step: a comma-separated list of "
        "hmm/3seq/maxchi/bootscan/geneconv/barcode/heuristic, or 'all'. Several run as "
        "an ensemble (default hmm,3seq,maxchi,bootscan).",
```

In `docs/detection-methods.md`, replace:

```markdown
### Low-divergence panels (intra-species sets, DNA viruses)
```

with:

```markdown
### Near-duplicate references (lineage clustering)

A recruited panel often holds several near-identical genomes of one lineage. Competed
individually they tie in every window and fragment the call, so the HMM caller first
pools them (`--cluster-lineages`, on by default; hmm only). Two references are pooled
when their identity stays at or above 98.5 % in every window, with no region-sized run
below it. The pooled lineage competes as one state under the label of its best-covering
member, and a region names that member. Clustering is skipped for panels of fewer than 4
or more than 200 references.

One limitation follows from the absolute threshold. On a panel in which *every* pair of
references is at least 98.5 % identical in every window -- mpox, VZV and other sets below
roughly 1.5 % divergence -- all references pool into a single lineage and the HMM has
nothing left to compete, so it calls no region whatever the data. The site-based callers
(3SEQ, MaxChi) are unaffected. On such a panel, run with `--no-cluster-lineages` to keep
the HMM's vote. `run_provenance.json` records whether clustering was on.

### Low-divergence panels (intra-species sets, DNA viruses)
```

In `docs/reference-panels.md`, replace:

```markdown
lineage saturates `nt` and no parental lineage can be recruited by similarity. A broad
fetch is capped (`--fetch-limit`, default 2000) and dereplicated; for a heavily
sequenced taxon the capped sample may miss lineages, so a curated `--candidate-pool`
is recommended (and the run says so). Whether the diversity panel actually contains the
```

with:

```markdown
lineage saturates `nt` and no parental lineage can be recruited by similarity. A broad
fetch is not truncated: the whole set is downloaded and dereplicated locally, and the run
logs a notice when it exceeds `--fetch-limit` (default 2000) because that step can take
a few minutes. For a heavily sequenced taxon a curated `--candidate-pool` is the faster
route. Whether the diversity panel actually contains the
```

In `example_data/README.md`, replace:

```markdown
The default ensemble also runs 3SEQ, which pools the discriminating sites into an exact
triplet test and recovers the event (q-value ~1e-12, `methods` = 3seq):
```

with:

```markdown
The default ensemble also runs the site-based callers, which pool the discriminating
sites into triplet tests and recover the event (q-value ~1e-12, `methods` = 3seq,maxchi):
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_report_stale_text.py tests/unit/test_cli_commands.py tests/integration -q`
Expected: PASS -- `58 passed, 1 skipped`

Run: `ruff check src tests validation && mypy src`
Expected: `All checks passed!` and `Success: no issues found in 79 source files`

- [ ] **Step 5: Commit**

```bash
git add tests/unit/test_report_stale_text.py src/tessera/recomb/report.py src/tessera/recomb/report_html.py src/tessera/recomb/report_assets.py src/tessera/cli/cmd_recomb.py src/tessera/cli/cmd_detect.py src/tessera/cli/cmd_fill_references.py docs/detection-methods.md docs/reference-panels.md example_data/README.md
git commit -m "Bring report, help and documentation text up to date" -m "Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 11: Housekeeping: unused dependency, version probe, validation README (B10)

`seaborn` is a declared runtime dependency that nothing imports. `sibeliaz -v` is rejected by the wrapper (`illegal option -- v`, exit 1) and that line is recorded as the aligner version. `validation/README.md` says the clonal control reports 0 regions and "7 PASS, 0 FAIL", and its "default agreement gate" is the harness default, not the CLI default.

**Files:**
- Create: `tests/unit/test_dependencies.py`
- Modify: `src/tessera/core/binaries.py`
- Modify: `src/tessera/aligners/sibeliaz.py`
- Modify: `pyproject.toml`
- Modify: `validation/README.md`
- Test (extend): `tests/unit/test_binaries.py`

**Interfaces:**
- Consumes: nothing new.
- Produces: `BinarySpec.version_args: tuple[str, ...] | None` (`None` = the tool has no version option; it is not executed and its version is `unknown`).

- [ ] **Step 1: Write the failing tests**

Append to the end of `tests/unit/test_binaries.py`:

```python
# --- a failed probe is not a version --------------------------------------

def failing_binary(directory: Path, name: str, output: str) -> Path:
    """An executable that prints ``output`` to stderr and exits 1, like a tool given an
    option it does not have."""
    path = directory / name
    path.write_text(f'#!/bin/sh\nprintf %s "{output}" >&2\nexit 1\n')
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


def test_error_output_is_not_recorded_as_a_version(on_path) -> None:
    """`sibeliaz -v` prints "illegal option -- v" and exits 1; that line went into the
    provenance sidecar as the aligner version."""
    failing_binary(on_path, "sibeliaz", "sibeliaz: illegal option -- v")
    versions = check_binaries((BinarySpec("sibeliaz", version_args=("-v",)),))
    assert versions["sibeliaz"] == "unknown"


def test_failed_probe_that_still_prints_a_version_is_kept(on_path) -> None:
    # Some tools print their version in a usage message and exit non-zero.
    failing_binary(on_path, "usagey", "usagey 1.4.2 -- usage: usagey [options]")
    assert check_binaries((BinarySpec("usagey"),))["usagey"] == "1.4.2"


def test_tool_without_a_version_option_is_not_probed(on_path, tmp_path) -> None:
    marker = tmp_path / "was_run"
    path = on_path / "noversion"
    path.write_text(f'#!/bin/sh\ntouch "{marker}"\n')
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    versions = check_binaries((BinarySpec("noversion", version_args=None),))
    assert versions == {"noversion": "unknown"}
    assert not marker.exists()  # declared as having no version option: never executed


def test_sibeliaz_declares_no_version_probe() -> None:
    from tessera.aligners.sibeliaz import SibeliazAligner

    (spec,) = SibeliazAligner.capabilities.required_binaries
    assert spec.version_args is None
```

Create `tests/unit/test_dependencies.py`:

```python
"""Every declared runtime dependency is one the package imports.

Tessera is dependency-light by design, so a declared dependency nothing imports is a
cost with no benefit: it is installed for every user and audited in CI.
"""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

# Distribution name -> import name, where they differ.
IMPORT_NAME = {"biopython": "Bio"}


def test_every_runtime_dependency_is_imported() -> None:
    project = tomllib.loads((REPO / "pyproject.toml").read_text())["project"]
    names = [re.split(r"[<>=!~\[ ;]", dep, maxsplit=1)[0] for dep in project["dependencies"]]
    source = "\n".join(p.read_text() for p in (REPO / "src" / "tessera").rglob("*.py"))
    unused = [
        name for name in names
        if not re.search(rf"^\s*(?:import|from)\s+{IMPORT_NAME.get(name, name)}\b",
                         source, flags=re.M)
    ]
    assert unused == []
```

`test_failed_probe_that_still_prints_a_version_is_kept` passes before and after: it pins what must not change (a tool that prints its version in a usage message and exits non-zero).

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_binaries.py tests/unit/test_dependencies.py -q`
Expected: FAIL --

```text
FAILED tests/unit/test_binaries.py::test_error_output_is_not_recorded_as_a_version
FAILED tests/unit/test_binaries.py::test_tool_without_a_version_option_is_not_probed
FAILED tests/unit/test_binaries.py::test_sibeliaz_declares_no_version_probe
FAILED tests/unit/test_dependencies.py::test_every_runtime_dependency_is_imported
4 failed, 20 passed in 1.76s
```

- [ ] **Step 3: Implement**

In `src/tessera/core/binaries.py`, replace:

```python
    ``version_args`` is the argument vector that prints a version (e.g.
    ``("--version",)``). ``min_version`` is an optional dotted-string requirement.
    """

    name: str
    version_args: tuple[str, ...] = field(default=("--version",))
```

with:

```python
    ``version_args`` is the argument vector that prints a version (e.g.
    ``("--version",)``), or ``None`` for a tool that has no version option: it is then
    not executed at all and its version is recorded as ``unknown``. ``min_version`` is
    an optional dotted-string requirement.
    """

    name: str
    version_args: tuple[str, ...] | None = field(default=("--version",))
```

In `src/tessera/core/binaries.py`, replace:

```python
    blob = (proc.stdout or "") + (proc.stderr or "")
    parsed = _parse_version(blob)
    if parsed:
        return ".".join(map(str, parsed))
    return blob.strip().splitlines()[0] if blob.strip() else None
```

with:

```python
    blob = (proc.stdout or "") + (proc.stderr or "")
    parsed = _parse_version(blob)
    if parsed:
        return ".".join(map(str, parsed))
    # No dotted version in the output. A tool that exited cleanly said something about
    # itself (a build date, say) and that is worth keeping; one that exited non-zero
    # printed an error about the option, which is not a version.
    if proc.returncode != 0:
        return None
    return blob.strip().splitlines()[0] if blob.strip() else None
```

In `src/tessera/core/binaries.py`, replace:

```python
        reported = _query_version(spec.name, spec.version_args)
        versions[spec.name] = reported or "unknown"
```

with:

```python
        reported = (
            None if spec.version_args is None
            else _query_version(spec.name, spec.version_args)
        )
        versions[spec.name] = reported or "unknown"
```

In `src/tessera/aligners/sibeliaz.py`, replace:

```python
        required_binaries=(BinarySpec("sibeliaz", version_args=("-v",)),),
```

with:

```python
        # The sibeliaz wrapper has no version option ("-v" is rejected as an illegal
        # option), so there is nothing to probe; its version is recorded as unknown.
        required_binaries=(BinarySpec("sibeliaz", version_args=None),),
```

In `pyproject.toml`, delete this line:

```toml
    "seaborn>=0.13,<1.0",
```

In `validation/README.md`, replace:

```markdown
It is deliberately sensitive to the failure mode it was built for. With the default
agreement gate the scan is clean; dropping the gate with `--min-methods 1` surfaces the
single-caller regions again (measured at 4 replicates: 9/16 runs, 10 false regions,
almost all from one caller). Treat a non-zero total as a regression to explain.
```

with:

```markdown
It is deliberately sensitive to the failure mode it was built for. **The harness's own
default is `--min-methods 2`; the `tessera` CLI default is `--min-methods 1`**, so a clean
default harness run does not describe the shipped default. Measured at 3 replicates on
2026-10-01 (commit `457bdfb`):

| gate | runs with a false region | false regions | source |
|---|---|---|---|
| `--min-methods 2` (harness default) | 0/12 (CI 0-24 %) | 0 | -- |
| `--min-methods 1` (CLI default) | 7/12 (58 %, CI 32-81 %) | 8 | hmm = 8 |

The positive control was detected 3/3 with the correct donor at both gates (median
breakpoint error 55 bp). Treat a non-zero total at `--min-methods 2` as a regression to
explain, and quote the `--min-methods 1` row when describing what a default run reports.
```

In `validation/README.md`, replace:

```markdown
| `hcv_clonal_1b` | HCV ~9.4 kb | pure genotype-1b (non-recombinant control) | mafft | resolves to 1b throughout; **0 regions** (real-data specificity) |
```

with:

```markdown
| `hcv_clonal_1b` | HCV ~9.4 kb | pure genotype-1b (non-recombinant control) | mafft | resolves to 1b throughout, but **currently FAILS**: one 12 bp MaxChi-only region (GT2a donor, q = 0.045) is reported where none is expected (real-data specificity) |
```

In `validation/README.md`, replace:

```markdown
Each reproduces its published event (or, for the clonal control, its *absence* of
recombination) end-to-end; the current run is **7 PASS, 0 FAIL**
(`orthopox_example` SKIPs until its 7-genome collection is built).
```

with:

```markdown
The six recombinant datasets that ran reproduce their published events end-to-end
(`orthopox_example` SKIPs until its 7-genome collection is built). The clonal control does
not currently pass: the run on 2026-10-01 was **6 PASS, 1 FAIL** (`hcv_clonal_1b`), **1
SKIP**. The failing region rests on a single caller at the CLI default `--min-methods 1`;
it is recorded here as an open specificity item rather than hidden.
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_binaries.py tests/unit/test_dependencies.py tests/unit/test_plugins.py tests/unit/test_aligner_params.py -q`
Expected: PASS -- `32 passed`

Run: `ruff check src tests validation && mypy src`
Expected: `All checks passed!` and `Success: no issues found in 79 source files`

**The numbers in `validation/README.md` were measured on 2026-10-01** with these changes applied: `python validation/run_validation.py` gave PASS for `sarscov2_xbb`, `hiv1_crf`, `norovirus_gii`, `enterovirus_e11`, `hiv_crf02ag`, `hcv_2k1b`; FAIL for `hcv_clonal_1b` (one 12 bp MaxChi-only region, donor `GT2a_JFH1`, q = 0.0451); SKIP for `orthopox_example`. If the fetched datasets are present under `validation/data/` and mafft and minimap2 are available, re-run it and correct the README where your result differs:

```bash
PY=$(which python)
PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin" $PY validation/run_validation.py | tail -12
```

If the data are not present, keep the text as written: it states the date it was measured.

- [ ] **Step 5: Commit**

```bash
git add tests/unit/test_binaries.py tests/unit/test_dependencies.py src/tessera/core/binaries.py src/tessera/aligners/sibeliaz.py pyproject.toml validation/README.md
git commit -m "Drop the unused seaborn dependency and stop recording a failed version probe" -m "Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 12: Close the remaining CLI validation gaps (B11)

Inputs that are accepted and fail later, or are silently ignored: `find-references --msa <missing>` and `recomb -o <existing file>` give "Unexpected error"; `--max-rounds 0` exits 0 having built nothing; `reassort` accepts `--ani-floor 500`, `--margin -3` and a `--dataset` key that matches no segment; `type-lineages` rejects `.fasta.gz` / `.fas` collections the other commands accept; a multi-line failure note breaks `segment_scan.tsv`.

**Files:**
- Modify: `src/tessera/cli/main.py`
- Modify: `src/tessera/cli/cmd_recomb.py`
- Modify: `src/tessera/cli/cmd_find_references.py`
- Modify: `src/tessera/cli/cmd_detect.py`
- Modify: `src/tessera/cli/cmd_build_panel.py`
- Modify: `src/tessera/cli/cmd_fill_references.py`
- Modify: `src/tessera/cli/cmd_reassort.py`
- Modify: `src/tessera/cli/cmd_type_lineages.py`
- Test (extend): `tests/unit/test_cli_input_checks.py`

**Interfaces:**
- Consumes: `_require_lineage_map` and the import lists from Task 6; `core.io.read_fasta`, `core.io.collection_genomes`, `core.io._require_fasta`.
- Produces: `_require_output_directory(path: Path) -> None` and `_require_scan_windows(window_size: int, window_step: int) -> None` in `cli/main.py`.

- [ ] **Step 1: Write the failing tests**

Append to the end of `tests/unit/test_cli_input_checks.py`:

```python
# --- paths ------------------------------------------------------------------

def test_find_references_rejects_a_missing_msa(tmp_path: Path, no_pipeline) -> None:
    result = runner.invoke(app, [
        "find-references", "--msa", str(tmp_path / "nope.fasta"), "--query", "query",
        "--output", str(tmp_path / "out"),
    ])
    _rejected(result, "MSA file not found")


def test_recomb_rejects_an_output_path_that_is_a_file(tmp_path: Path) -> None:
    occupied = tmp_path / "out"
    occupied.write_text("not a directory\n")
    result = runner.invoke(app, [
        "recomb", "--msa", EXAMPLE, "--query", "query", "--output", str(occupied),
    ])
    _rejected(result, "is not a directory")
    assert occupied.read_text() == "not a directory\n"  # left untouched


# --- numeric options on the panel-building commands ------------------------

@pytest.mark.parametrize("command", ["detect", "fill-references", "build-panel"])
@pytest.mark.parametrize(
    ("option", "value"),
    [("--max-rounds", "0"), ("--window-size", "0"), ("--window-step", "0")],
)
def test_panel_commands_reject_out_of_range_options(
    tmp_path: Path, no_pipeline, command: str, option: str, value: str
) -> None:
    """`fill-references --max-rounds 0` exited 0 having built no alignment and no
    report; a bad window reached the scan only after seeding and alignment."""
    result = runner.invoke(app, [
        command, "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        option, value,
    ])
    _rejected(result, f"Invalid {option}")


@pytest.mark.parametrize(("option", "value"), [("--window-size", "0"), ("--window-step", "0")])
def test_find_references_rejects_out_of_range_windows(
    tmp_path: Path, no_pipeline, option: str, value: str
) -> None:
    result = runner.invoke(app, [
        "find-references", "--msa", EXAMPLE, "--query", "query",
        "--output", str(tmp_path / "out"), option, value,
    ])
    _rejected(result, f"Invalid {option}")


# --- reassort ---------------------------------------------------------------

def _segments(tmp_path: Path) -> Path:
    path = tmp_path / "segments.fasta"
    path.write_text(">HA\nACGTACGT\n>NA\nACGTACGT\n")
    return path


@pytest.fixture
def no_reassort(monkeypatch):
    monkeypatch.setattr("tessera.cli.cmd_reassort.assign_segments", _must_not_run)


@pytest.mark.parametrize(
    ("option", "value", "needle"),
    [("--ani-floor", "500", "Invalid --ani-floor"),
     ("--ani-floor", "-1", "Invalid --ani-floor"),
     ("--margin", "-3", "Invalid --margin")],
)
def test_reassort_rejects_out_of_range_options(
    tmp_path: Path, no_reassort, option: str, value: str, needle: str
) -> None:
    """A negative margin makes every near-best set empty (a clonal pair then reads as
    undetermined); an ANI floor above 100 leaves every segment unassigned. Both exited 0."""
    result = runner.invoke(app, [
        "reassort", "--query", str(_segments(tmp_path)), "--output", str(tmp_path / "out"),
        option, value,
    ])
    _rejected(result, needle)


def test_reassort_rejects_a_dataset_override_for_an_unknown_segment(
    tmp_path: Path, no_reassort
) -> None:
    """A typo in the segment name silently fell back to dataset auto-detection."""
    result = runner.invoke(app, [
        "reassort", "--query", str(_segments(tmp_path)), "--output", str(tmp_path / "out"),
        "--dataset", "HAA=nextstrain/flu/h3n2/ha",
    ])
    _rejected(result, "HAA")
    assert "HA, NA" in " ".join(result.output.split())  # names the segments that do exist


def test_segment_scan_tsv_keeps_one_row_per_segment(tmp_path: Path, monkeypatch) -> None:
    """An aligner error message spans several lines and may hold tabs; written verbatim
    into the note column it broke the table."""
    from tessera.reassort.assign import ReassortmentResult, SegmentAssignment
    from tessera.reassort.scan import SegmentScan

    def fake_assign(query, **kwargs):
        return ReassortmentResult(
            segments=[SegmentAssignment("HA", "dataset", None, None, 0.0, "unassigned")],
            verdict="undetermined", groups=[], pair_notes=[],
            scans=[SegmentScan("HA", False, False, 0,
                               "scan failed: mafft failed:\nline1\tx\nline2")],
        )

    monkeypatch.setattr("tessera.cli.cmd_reassort.assign_segments", fake_assign)
    out = tmp_path / "out"
    result = runner.invoke(app, [
        "reassort", "--query", str(_segments(tmp_path)), "--output", str(out),
    ])
    assert result.exit_code == 0, result.output
    lines = (out / "segment_scan.tsv").read_text().splitlines()
    assert len(lines) == 2  # header + one segment
    assert all(len(line.split("\t")) == 4 for line in lines)
    assert lines[1].split("\t")[3] == "scan failed: mafft failed: line1 x line2"


# --- type-lineages reads the same collection the other commands do ---------

def test_type_lineages_accepts_what_a_collection_may_hold(tmp_path: Path, monkeypatch) -> None:
    """Every other command treats each file in the collection as a genome, whatever its
    extension and gzip-compressed or not; type-lineages filtered on three suffixes and
    reported "No FASTA genomes found"."""
    import gzip

    collection = tmp_path / "collection"
    collection.mkdir()
    with gzip.open(collection / "a.fasta.gz", "wt") as fo:
        fo.write(">a\nACGTACGT\n")
    (collection / "b.fas").write_text(">b\nACGTACGT\n")
    (collection / ".DS_Store").write_text("not a genome")
    seen: list[str] = []

    def fake_assign(genomes, **kwargs):
        seen.extend(p.name for p in genomes)
        return [("a", "L1", "denovo"), ("b.fas", "L1", "denovo")]

    monkeypatch.setattr("tessera.discover.lineage_assign.assign_lineages", fake_assign)
    out = tmp_path / "out"
    result = runner.invoke(app, [
        "type-lineages", "--collection", str(collection), "--output", str(out),
    ])
    assert result.exit_code == 0, result.output
    assert seen == ["a.fasta.gz", "b.fas"]  # hidden files are not genomes
    assert (out / "lineages.tsv").exists()


def test_type_lineages_rejects_a_file_that_is_not_fasta(tmp_path: Path, no_pipeline) -> None:
    collection = _collection(tmp_path)
    (collection / "notes.txt").write_text("these are my notes\n")
    result = runner.invoke(app, [
        "type-lineages", "--collection", str(collection), "--output", str(tmp_path / "out"),
    ])
    _rejected(result, "does not look like a FASTA file")
```

`type-lineages` now reads a collection the way every other command does: each non-hidden file is a genome, and a file that is not FASTA is an error rather than something to skip. A collection holding a stray `notes.txt` was silently accepted before and is now rejected by name.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_cli_input_checks.py -q`
Expected: FAIL --

```text
FAILED tests/unit/test_cli_input_checks.py::test_find_references_rejects_a_missing_msa
FAILED tests/unit/test_cli_input_checks.py::test_recomb_rejects_an_output_path_that_is_a_file
FAILED tests/unit/test_cli_input_checks.py::test_panel_commands_reject_out_of_range_options[--max-rounds-0-detect]
FAILED tests/unit/test_cli_input_checks.py::test_panel_commands_reject_out_of_range_options[--max-rounds-0-fill-references]
FAILED tests/unit/test_cli_input_checks.py::test_panel_commands_reject_out_of_range_options[--max-rounds-0-build-panel]
FAILED tests/unit/test_cli_input_checks.py::test_panel_commands_reject_out_of_range_options[--window-size-0-detect]
FAILED tests/unit/test_cli_input_checks.py::test_panel_commands_reject_out_of_range_options[--window-size-0-fill-references]
FAILED tests/unit/test_cli_input_checks.py::test_panel_commands_reject_out_of_range_options[--window-size-0-build-panel]
FAILED tests/unit/test_cli_input_checks.py::test_panel_commands_reject_out_of_range_options[--window-step-0-detect]
FAILED tests/unit/test_cli_input_checks.py::test_panel_commands_reject_out_of_range_options[--window-step-0-fill-references]
FAILED tests/unit/test_cli_input_checks.py::test_panel_commands_reject_out_of_range_options[--window-step-0-build-panel]
FAILED tests/unit/test_cli_input_checks.py::test_find_references_rejects_out_of_range_windows[--window-size-0]
FAILED tests/unit/test_cli_input_checks.py::test_find_references_rejects_out_of_range_windows[--window-step-0]
FAILED tests/unit/test_cli_input_checks.py::test_reassort_rejects_out_of_range_options[--ani-floor-500-Invalid --ani-floor]
FAILED tests/unit/test_cli_input_checks.py::test_reassort_rejects_out_of_range_options[--ani-floor--1-Invalid --ani-floor]
FAILED tests/unit/test_cli_input_checks.py::test_reassort_rejects_out_of_range_options[--margin--3-Invalid --margin]
FAILED tests/unit/test_cli_input_checks.py::test_reassort_rejects_a_dataset_override_for_an_unknown_segment
FAILED tests/unit/test_cli_input_checks.py::test_segment_scan_tsv_keeps_one_row_per_segment
FAILED tests/unit/test_cli_input_checks.py::test_type_lineages_accepts_what_a_collection_may_hold
FAILED tests/unit/test_cli_input_checks.py::test_type_lineages_rejects_a_file_that_is_not_fasta
20 failed, 6 passed in 1.97s
```

- [ ] **Step 3: Implement**

In `src/tessera/cli/main.py`, replace:

```python
def _parse_key_values(items: list[str], label: str) -> dict[str, str]:
```

with:

```python
def _require_output_directory(path: Path) -> None:
    """Reject an output path that exists and is not a directory.

    The writers create the directory if it is missing; given an existing file they fail
    with ``[Errno 17] File exists`` from wherever the first output is written.
    """
    if Path(path).exists() and not Path(path).is_dir():
        raise UserInputError(f"Output path exists and is not a directory: {path}")


def _require_scan_windows(window_size: int, window_step: int) -> None:
    """The window options every scanning command shares; see :func:`_require_range`."""
    _require_range(window_size, "--window-size", lo=1)
    _require_range(window_step, "--window-step", lo=1)


def _parse_key_values(items: list[str], label: str) -> dict[str, str]:
```

In `src/tessera/cli/cmd_recomb.py`, replace:

```python
    _require_lineage_map,
    _require_range,
```

with:

```python
    _require_lineage_map,
    _require_output_directory,
    _require_range,
```

In `src/tessera/cli/cmd_recomb.py`, replace:

```python
        _require_lineage_map(lineage_map)
```

with:

```python
        _require_lineage_map(lineage_map)
        _require_output_directory(output)
```

In `src/tessera/cli/cmd_find_references.py`, replace:

```python
from .main import app, get_logger, stage_errors
```

with:

```python
from .main import _require_file, _require_scan_windows, app, get_logger, stage_errors
```

In `src/tessera/cli/cmd_find_references.py`, replace:

```python
    with stage_errors(logger):
        params = FindRefParams(
```

with:

```python
    with stage_errors(logger):
        _require_file(msa, "MSA file")
        _require_scan_windows(window_size, window_step)
        params = FindRefParams(
```

In `src/tessera/cli/cmd_detect.py`, replace:

```python
    _require_lineage_map,
    app,
```

with:

```python
    _require_lineage_map,
    _require_range,
    _require_scan_windows,
    app,
```

In `src/tessera/cli/cmd_detect.py`, replace:

```python
        _require_lineage_map(lineage_map)
```

with:

```python
        _require_lineage_map(lineage_map)
        # Checked here, before any network or aligner work: a round count of zero
        # builds nothing and still exits 0, and a bad window otherwise surfaces only
        # after the panel has been recruited and aligned.
        _require_range(max_rounds, "--max-rounds", lo=1)
        _require_scan_windows(window_size, window_step)
```

In `src/tessera/cli/cmd_build_panel.py`, replace:

```python
    _require_lineage_map,
    app,
```

with:

```python
    _require_lineage_map,
    _require_range,
    _require_scan_windows,
    app,
```

In `src/tessera/cli/cmd_build_panel.py`, replace:

```python
        _require_lineage_map(lineage_map)
```

with:

```python
        _require_lineage_map(lineage_map)
        # Checked here, before any network or aligner work: a round count of zero
        # builds nothing and still exits 0, and a bad window otherwise surfaces only
        # after the panel has been recruited and aligned.
        _require_range(max_rounds, "--max-rounds", lo=1)
        _require_scan_windows(window_size, window_step)
```

In `src/tessera/cli/cmd_fill_references.py`, replace:

```python
    _require_lineage_map,
    app,
```

with:

```python
    _require_lineage_map,
    _require_range,
    _require_scan_windows,
    app,
```

In `src/tessera/cli/cmd_fill_references.py`, replace:

```python
        _require_lineage_map(lineage_map)
```

with:

```python
        _require_lineage_map(lineage_map)
        # Checked here, before any network or aligner work: a round count of zero
        # builds nothing and still exits 0, and a bad window otherwise surfaces only
        # after the panel has been recruited and aligned.
        _require_range(max_rounds, "--max-rounds", lo=1)
        _require_scan_windows(window_size, window_step)
```

In `src/tessera/cli/cmd_reassort.py`, replace:

```python
from ..core.errors import UserInputError
from ..reassort import assign_segments
```

with:

```python
from ..core.errors import UserInputError
from ..core.io import read_fasta
from ..reassort import assign_segments
```

In `src/tessera/cli/cmd_reassort.py`, replace:

```python
from .main import _require_file, app, get_logger, stage_errors
```

with:

```python
from .main import _require_file, _require_range, app, get_logger, stage_errors
```

In `src/tessera/cli/cmd_reassort.py`, replace:

```python
        _require_file(query, "Query file")
        overrides: dict[str, str] = {}
        for item in dataset or []:
            if "=" not in item:
                raise UserInputError(f"--dataset must be SEGMENT=path, got {item!r}")
            seg, path = item.split("=", 1)
            overrides[seg.strip()] = path.strip()
```

with:

```python
        _require_file(query, "Query file")
        # ANI is a percentage. Above 100 nothing can be assigned; a negative margin
        # leaves every near-best set empty, so a clonal pair reads as undetermined.
        _require_range(ani_floor, "--ani-floor", lo=0.0, hi=100.0)
        _require_range(margin, "--margin", lo=0.0)
        overrides: dict[str, str] = {}
        for item in dataset or []:
            if "=" not in item:
                raise UserInputError(f"--dataset must be SEGMENT=path, got {item!r}")
            seg, path = item.split("=", 1)
            overrides[seg.strip()] = path.strip()
        if overrides:
            # An override is looked up by segment name; one that matches no record would
            # be ignored and that segment's dataset auto-detected instead.
            segments = [name for name, _seq in read_fasta(query) if name]
            unknown = sorted(set(overrides) - set(segments))
            if unknown:
                raise UserInputError(
                    f"--dataset names segment(s) not in the query: {', '.join(unknown)}. "
                    f"Segments in {query.name}: {', '.join(segments) or '(none)'}."
                )
```

In `src/tessera/cli/cmd_reassort.py`, replace:

```python
                    fo.write(f"{sc.segment}\t{flag}\t{sc.n_regions}\t{sc.note}\n")
```

with:

```python
                    # A failure note quotes the aligner's own message, which can span
                    # lines and hold tabs; keep it to one cell.
                    note = " ".join(sc.note.split())
                    fo.write(f"{sc.segment}\t{flag}\t{sc.n_regions}\t{note}\n")
```

In `src/tessera/cli/cmd_type_lineages.py`, replace:

```python
from ..core.errors import UserInputError
```

with:

```python
from ..core.errors import UserInputError
from ..core.io import _require_fasta, collection_genomes
```

In `src/tessera/cli/cmd_type_lineages.py`, replace:

```python
        genomes = sorted(
            p for p in collection.iterdir()
            if p.is_file() and p.suffix.lower() in (".fasta", ".fa", ".fna")
        )
        if not genomes:
            raise UserInputError(f"No FASTA genomes found in {collection}")
```

with:

```python
        # The same reading of a collection as every other command: each non-hidden
        # file is a genome, whatever its extension and gzip-compressed or not, and a
        # file that is not FASTA is an error rather than something to skip.
        genomes = collection_genomes(collection)
        if not genomes:
            raise UserInputError(f"No genome files found in {collection}")
        for genome in genomes:
            _require_fasta(genome)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_cli_input_checks.py tests/unit/test_cli_commands.py tests/unit/test_cli_validation.py tests/unit/test_cli_reassort.py -q`
Expected: PASS -- `79 passed`

Run: `ruff check src tests validation && mypy src`
Expected: `All checks passed!` and `Success: no issues found in 79 source files`

- [ ] **Step 5: Commit**

```bash
git add tests/unit/test_cli_input_checks.py src/tessera/cli/main.py src/tessera/cli/cmd_recomb.py src/tessera/cli/cmd_find_references.py src/tessera/cli/cmd_detect.py src/tessera/cli/cmd_build_panel.py src/tessera/cli/cmd_fill_references.py src/tessera/cli/cmd_reassort.py src/tessera/cli/cmd_type_lineages.py
git commit -m "Validate the remaining CLI inputs before any work starts" -m "Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 13: Changelog, full gate and harness comparison

**Files:**
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes: the "before" harness outputs saved in Task 1 (`/tmp/tessera-plan-b/before_gate2.txt`, `before_gate1.txt`).
- Produces: the pull request.

- [ ] **Step 1: Add the changelog entry**

In `CHANGELOG.md`, replace:

```markdown
## [Unreleased]
```

with:

```markdown
## [Unreleased]

Reporting fixes from the post-1.2.0 audit
(`docs/superpowers/specs/2026-10-01-post-1.2.0-audit-design.md`, items B1-B11). No region
call changes: which regions are found, and their coordinates and p-values, are as before.

### Fixed

- **The run provenance named the wrong callers.** MaxChi, Bootscan, GENECONV and the barcode
  caller were each described as `heuristic (min ... / margin ... / merge ...)` in
  `run_provenance.json` and the report, so a default run recorded
  `hmm + 3seq + heuristic + heuristic`. Each caller is now described under its own name, and
  the record gains the settings that change what is reported: the agreement gate
  (`--min-methods`), sibling exclusion, lineage clustering and donor re-attribution.
- **A clean recombinant between divergent parents was reported as a possible missing
  reference.** A window straddling a breakpoint matches neither parent well on its own, so
  its best similarity fell below the coverage threshold, the stretch was called a `divergent`
  coverage gap, and the region was marked `donor_undercovered` -- on the shipped
  `divergent` example the headline read "low confidence" for a donor identical to the query.
  A gap within one window of a called region boundary that the region's two parents together
  explain is now labelled `breakpoint` in `coverage_gaps.tsv` and the report; it does not
  caveat the region, is not turned into a donor-absent region, and is left out of the
  headline. **`donor_undercovered` and the confidence wording change for such regions.**
  Reference recruitment (`fill-references`, `find-references`) is unaffected.
- The report's "covering N kb (P %) of the query" added the lengths of overlapping regions
  that name different donors; it now reports their union.
- After `--reattribute-donors`, `recombination_methods.tsv` and the report's method table
  kept the donor from before re-attribution.
- `--lineage-map` pointing at a file that does not exist was ignored (exit 0, untyped
  report) by `recomb`, `type-lineages`, `detect`, `fill-references` and `build-panel`. It
  is now an error.
- **The barcode caller on an untyped panel read as a negative.** It named the first record
  of the alignment as major parent and its column in the method table said `no`. It now
  reports no major parent; in an ensemble it is logged and shown as `not run`, the
  agreement gate counts only callers that ran, and a run that selected only `barcode` is
  refused.
- The report judged the PHI p-value at alpha 0.05 whatever `--alpha` was.
- **A PHI test that could not reject was reported as "no signal".** With no more
  informative sites than the window can hold (`--phi-window`, default 100) the permutation
  p-value is 1 for any data. This is now reported as `not testable` (`NA` in the
  `recombination_profile.tsv` header, `phi_p = None` in the API).
- Plots labelled a donor-absent region "recombinant: <backbone>"; `similarity_pair` showed
  the two leading window winners rather than the major parent and the leading donor; with
  `--top-n 1` the donor was drawn grey.
- The report footer listed `.pdf` plots under `--plot-format png` and files that were not
  written; the methods text and references described a two-caller ensemble; `--method`
  help omitted `geneconv`.
- `sibeliaz` was probed with `-v`, which it rejects, and the error line was recorded as the
  aligner version. A failed probe is no longer recorded as a version.
- Input checks: `find-references --msa <missing>` and `recomb -o <existing file>` failed
  with "Unexpected error"; `--max-rounds 0` exited 0 having built nothing; `reassort`
  accepted `--ani-floor 500`, `--margin -3` and a `--dataset` key matching no segment;
  `type-lineages` rejected `.fasta.gz` / `.fas` collections; a multi-line aligner error
  broke `segment_scan.tsv`.

### Changed

- `RecombinationSignal.phi_p` is `float | None`.
- `seaborn` is no longer a dependency; nothing imported it.

### Documentation

- `docs/detection-methods.md` describes lineage clustering (including its limit on panels
  below about 1.5 % divergence), the `breakpoint` coverage kind and the untestable-PHI case.
- `validation/README.md` states that the specificity harness defaults to `--min-methods 2`
  while the CLI defaults to 1, gives the measured rate at each, and records that
  `hcv_clonal_1b` currently fails.
```

- [ ] **Step 2: Run the full gate**

Run: `ruff check src tests validation`
Expected: `All checks passed!`

Run: `pytest -m "not requires_binary" -q`
Expected: `671 passed, 1 deselected`

Run: `mypy src`
Expected: `Success: no issues found in 79 source files`

Then repeat the suite the way CI sees it -- a narrow terminal and no aligner or NCBI tools on `PATH` -- so a test that only passes because `skani` or `efetch` is installed locally is caught here:

```bash
PY=$(which python)
PATH="$(dirname "$PY"):/usr/bin:/bin" COLUMNS=80 "$PY" -m pytest -m "not requires_binary" -q
```

Expected: `671 passed, 1 deselected`

- [ ] **Step 3: Run the specificity harness again and compare**

```bash
python validation/run_specificity.py --reps 3 | tee /tmp/tessera-plan-b/after_gate2.txt
python validation/run_specificity.py --reps 3 --min-methods 1 | tee /tmp/tessera-plan-b/after_gate1.txt
diff /tmp/tessera-plan-b/before_gate2.txt /tmp/tessera-plan-b/after_gate2.txt
diff /tmp/tessera-plan-b/before_gate1.txt /tmp/tessera-plan-b/after_gate1.txt
```

Expected: both `diff` commands print nothing. This plan changes no region call, so the false-region counts and the positive control must be identical to the baseline (`0/12` and `7/12, hmm=8`; `detected 3/3 | correct donor 3/3 | median breakpoint error 55 bp`). Any difference means a caller's behaviour changed: stop and find out why before going further.

- [ ] **Step 4: Run the hybrids harness if its data and aligners are present**

Only if Task 1 Step 4 produced `/tmp/tessera-plan-b/hybrids_before.txt`. The harness scores detection, donor and backbone from `recombination_regions.tsv`; this plan changes the `donor_undercovered` flag only, so no case's verdict should move.

```bash
PY=$(which python)
PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin" $PY validation/run_hybrids.py \
    | tee /tmp/tessera-plan-b/hybrids_after.txt | tail -5
grep -E "passed  \(" /tmp/tessera-plan-b/hybrids_before.txt /tmp/tessera-plan-b/hybrids_after.txt
```

Expected: the same `N/M passed  (S skipped, E error)` line in both files. Read the per-case lines above it in the two files as well (run times differ between runs, so compare the verdicts, not the whole file) and list any case whose PASS/FAIL/SKIP changed. If the harness could not be run, say so in the pull request -- do not report it as run.

- [ ] **Step 5: Commit and open the pull request**

```bash
git add CHANGELOG.md
git commit -m "Record the reporting fixes in the changelog" -m "Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
git push -u origin fix-audit-report-faithfulness
gh pr create --base main --title "Reporting faithfulness: audit items B1-B11" --body "$(cat <<'BODY'
Implements plan B of the post-1.2.0 audit (docs/superpowers/specs/2026-10-01-post-1.2.0-audit-design.md, items B1-B11): outputs that misstated what a run did. No region call changes.

Verification
- ruff, mypy clean; pytest -m "not requires_binary": 671 passed (586 before).
- run_specificity.py --reps 3: <paste the TOTAL and positive-control lines, before and after, for both gates>
- run_hybrids.py: <paste the before/after comparison, or state that it was not run and why>

Behaviour a user will notice
- A region whose only coverage gaps are the windows straddling its breakpoints is no longer flagged donor_undercovered / "low confidence".
- --method barcode on an untyped panel is refused; in an ensemble it is shown as "not run".
- The PHI test is reported as "not testable" when there are too few informative sites for the window.
- A missing --lineage-map file is an error.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
BODY
)"
```

Replace the two `<paste ...>` placeholders in the body with the measured lines before running the command. Do not merge: hand the pull request to the maintainer.
