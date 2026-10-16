# Audit Plan F: Reporting Faithfulness, CLI and Reassortment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the report, the TSVs, `run_provenance.json`, the CLI and the reassortment call say what the run actually did -- in particular, report a test or caller that did not run as "not run" with its reason -- without changing which regions are called or their coordinates.

**Architecture:** One mapping, `not_run` (name -> reason), built in `run_recomb` and carried by `ReportContext`, replaces plan B's barcode-only `methods_not_run` / `not_run_reason` pair; the PHI test and a one-state HMM are added to it, and every writer reads it. The remaining reporting items are local text, ordering and colour fixes in the report writers. The CLI items add input checks and help text at the command layer only, the reassortment items change the pure constellation call and what the command writes, and the harness items replace hand-parsed `sys.argv` with `argparse`.

**Tech Stack:** Python 3.11+, numpy, pytest, ruff, mypy. No new dependency.

**Spec:** `docs/superpowers/specs/2026-10-12-audit-swot-gap.md` (items B1, B3, B5, B7, B8, B9, the "Smaller reporting items" paragraph; the documentation half of B4; help text and documentation only for B2, B6 and A4; E1-E4, E5 except the three plan-E sub-items, E14; F1-F5; the matching "Reporting side" and "CLI side" test gaps of section 3).

## Global Constraints

- No new runtime dependency.
- Modest scientific language in code, docs and messages; reported numbers must be faithful (state what passes, what fails, what was skipped).
- "Could not test" must never be reported as "tested, found nothing".
- A new behaviour needs a test that fails without the change; if a test is a guard that passes before and after, the plan says so.
- Branch: create `fix-audit4-reporting` from `fix-audit-caller-behaviour` (or from `main` once that branch is merged). Commit messages explain the problem, the evidence and the trade-off, and end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Append the aligner env to PATH, never prepend (prepending swaps the Python interpreter): `export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"`.
- The dev Mac has 18 GB RAM: run one harness at a time.
- Each task ends with `ruff check src tests validation`, `mypy`, and `pytest -m "not requires_binary" -q` green before its commit.
- No Unicode in code or docs beyond what files already contain; ASCII punctuation.
- **No task in this plan may change which regions are called or their coordinates.** The check: `python validation/run_specificity.py --reps 3 --min-methods 1` is run once in Task 1 (before any change) and once in Task 18 (after the last), never concurrently with anything else, and the two outputs must be identical apart from the two header lines. Measured while writing this plan, on `4b35e39` and on the full set of candidate changes below, both runs gave:

  ```
  TOTAL                7/12 (58%, CI 32%-81%)              8  hmm=8
  Positive control -- clade-B tract spliced at query 4000-7000:
    detected 3/3 | correct donor 3/3 | median breakpoint error 13 bp
  ```

  If the after-run differs, stop and find the task responsible (bisect the branch with the same command) before opening the pull request.
- Two things this plan deliberately does **not** do, because each would change which regions are reported; both belong to plan G:
  1. A one-state HMM (Task 3) is reported as not run, but it does **not** lower the `--min-methods` agreement gate the way an unrunnable barcode caller does.
  2. `--merge-gap` is **not** passed to the coverage diagnostic (B6). The coverage gaps feed the donor-absent regions, so wiring it through would change the reported regions whenever `--merge-gap` differs from the window size. Task 5 documents the option as heuristic-only instead.

## Review Focus

Inputs the spec implies that no task's tests would otherwise exercise, most likely to bite a user first. Each has a test in the task that owns the code.

1. **`--no-phi` on a single-caller run.** The method comparison is omitted for one caller, so the regions TSV is the only place that says "not run". *Task 2, `test_no_phi_on_a_single_caller_run_still_says_not_run`.*
2. **The barcode caller and a one-state HMM not run in the same run.** Both must appear in the provenance and in the method table, each with its own reason, and only the barcode may lower the gate. *Task 3, `test_barcode_and_one_state_hmm_are_both_reported_not_run`.*
3. **A query with both a divergent and a low-information stretch.** The missing-reference caveat must count only the divergent one and the low-information note must still appear. *Task 4, `test_divergent_and_low_information_gaps_get_separate_notes`.*
4. **A donor-absent region under `--min-methods 2`.** The agreement filter removes rows before the donor-absent rows are added; the two files must still pair row for row. *Task 8, `test_methods_rows_match_regions_rows_under_the_agreement_gate`.*
5. **A new, nested output directory.** The output-path check in `get_logger` must not stop a command from creating an output directory that does not exist yet. Guard, passes before and after. *Task 10, `test_a_new_nested_output_directory_is_still_created`.*

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `src/tessera/recomb/report_context.py` | modify | `not_run` mapping, `PHI` key, `methods_not_run` / `phi_not_run` / `missing_reference_gaps` / `low_information_gaps` properties, `major_parent` |
| `src/tessera/recomb/run.py` | modify | build `not_run`; PHI and HMM not-run reasons; provenance lines; donor-absent method rows; log split for coverage kinds |
| `src/tessera/recomb/regions.py` | modify | `call_regions(..., not_run=)`; the HMM records a one-state decode |
| `src/tessera/recomb/diagnostics.py` | modify | public `MIN_INFORMATIVE` |
| `src/tessera/recomb/coverage.py` | modify | `DIVERGENT_KIND`, `LOW_INFORMATION_KIND` constants |
| `src/tessera/recomb/report.py` | modify | pass `not_run` to the TSV writers; the pair plot gets the report's colour order |
| `src/tessera/recomb/report_text.py` | modify | `parent_free_support = not run`; console coverage notes by kind; stale comment |
| `src/tessera/recomb/report_html.py` | modify | signal section, not-run caveat, method table, low-information note, headline backbone, methods paragraph |
| `src/tessera/recomb/report_assets.py` | modify | glossary: major parent, window / step, breakpoint |
| `src/tessera/recomb/report_plots.py` | modify | `plot_pairwise(..., color_labels=)` |
| `src/tessera/cli/main.py` | modify | `get_logger` checks the output path; `_require_coverage_fractions`, `_require_search_counts` |
| `src/tessera/cli/cmd_*.py` | modify | output-path check, range checks, scope markers, `msa` and `build-panel` options, `reassort` bounds and evidence |
| `src/tessera/reassort/constellation.py` | modify | `PairCall`, `MAX_MARGIN`, within-group discordance, empty candidate lists |
| `src/tessera/reassort/assign.py`, `scan.py`, `__init__.py` | modify | carry `pairs`; `threads` for segment scans; `clades_within_margin` |
| `validation/run_specificity.py`, `run_regimes.py`, `characterise_false_regions.py` | modify | `argparse` |
| `tests/unit/test_not_run_reporting.py`, `test_coverage_caveat.py`, `test_cli_option_scope.py`, `test_report_text_faithful.py`, `test_report_backbone.py`, `test_methods_table_rows.py`, `test_cli_output_path.py`, `test_harness_arguments.py`, `test_docs_commands.py` | create | tests for the tasks below |
| `tests/unit/test_report_plots.py`, `test_constellation.py`, `test_cli_input_checks.py`, `test_reassort_assign.py` | modify | new tests; three existing `fake_scan` stubs gain a `threads` parameter (Task 15) |
| `docs/detection-methods.md`, `docs/reference-panels.md`, `README.md`, `CONTRIBUTING.md`, `CLAUDE.md`, `validation/README.md`, `example_data/README.md`, `CHANGELOG.md` | modify | each in the task whose code it describes |

## Tasks

### Task 1: Branch and baseline  (audit items: none; the region-call check)

Nothing is committed in this task.

**Files:**
- Create (untracked, git-ignored): `validation/data/audit-f/spec_before.txt`

**Interfaces:**
- Produces: `validation/data/audit-f/spec_before.txt`, compared in Task 18.

- [ ] **Step 1: Branch**

```bash
git switch fix-audit-caller-behaviour      # or main, once that branch is merged
git switch -c fix-audit4-reporting
pip install -e ".[dev]"
```

- [ ] **Step 2: Confirm the starting point is green**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: `All checks passed!`, `Success: no issues found in 79 source files`, `757 passed, 5 deselected`. If the count differs, note the number; later steps quote counts relative to it.

- [ ] **Step 3: Record the region-call baseline (nothing else running)**

```bash
mkdir -p validation/data/audit-f
python validation/run_specificity.py --reps 3 --min-methods 1 > validation/data/audit-f/spec_before.txt 2>&1
tail -4 validation/data/audit-f/spec_before.txt
```

Expected (measured on `4b35e39`, about 3 minutes): `TOTAL 7/12 (58%, CI 32%-81%) 8 hmm=8` and `detected 3/3 | correct donor 3/3 | median breakpoint error 13 bp`.

---

### Task 2: Report a PHI test that did not run as "not run"  (audit items: B1; the `not_run` mechanism)

`run.py:304-309` leaves `signal` as `None` both when `--no-phi` turned the test off and when there were too few informative sites, and `report_html.py:503-507` renders one message for both: "Too few informative sites in the alignment for a parent-free recombination test." On `example_data/divergent.msa.fasta` (159 informative columns) that is false. Every region then reads `parent_free_support = no` (`corroborating_intervals(None)` is empty, `diagnostics.py:307`), and `run.py:529-535` writes no PHI line to the provenance. The same "no" is written when PHI ran but was not testable at its window.

This task replaces plan B's barcode-only pair (`ReportContext.methods_not_run` tuple, `not_run_reason` string) with one mapping `not_run: dict[str, str]` (name -> reason). The barcode caller and the PHI test (key `PHI = "phi"`) are its first entries; Task 3 adds the HMM. `methods_not_run` becomes a property, so `report.py` and the barcode tests keep working.

**Files:**
- Modify: `src/tessera/recomb/report_context.py:13` (import), `:19` (constant), `:36-41` (fields), `:51` (properties)
- Modify: `src/tessera/recomb/diagnostics.py:43-44`
- Modify: `src/tessera/recomb/run.py:24`, `:29-35`, `:301-321`, `:356-378`, `:398`, `:523-535`, `:558-566`
- Modify: `src/tessera/recomb/report.py:39-41`, `:57-62`, `:101-107`
- Modify: `src/tessera/recomb/report_text.py:240-277`, `:280-304`
- Modify: `src/tessera/recomb/report_html.py:20`, `:279-345`, `:501-507`, `:615-618`, `:646-647`
- Modify: `docs/detection-methods.md` (PHI section, output-file table)
- Create: `tests/unit/test_not_run_reporting.py`

**Interfaces:**
- Produces: `report_context.PHI: str = "phi"`; `ReportContext.not_run: dict[str, str]` (replaces `methods_not_run` and `not_run_reason` as fields); properties `ReportContext.methods_not_run -> tuple[str, ...]` and `ReportContext.phi_not_run -> str | None`; `diagnostics.MIN_INFORMATIVE: int`; `write_regions_tsv(regions, output_dir, logger, phi_not_run: bool = False)`; `write_methods_tsv(breakdown, methods_run, output_dir, logger, methods_not_run=(), phi_not_run: bool = False)`; `_method_comparison_html(breakdown, methods_run, per_major, lineage_map=None, not_run: dict[str, str] | None = None)`; `_method_section(method_breakdown, methods_run, per_major, lineage_map, not_run=None)`; `_signal_html(signal, alpha=0.05, not_run: str | None = None)`. `report.py` re-exports `PHI`.
- Consumes: nothing new.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_not_run_reporting.py`:

```python
"""A test or caller that gave no result is reported as "not run", with its reason.

Audit 4, items B1 and B7: ``--no-phi`` made the report say the alignment had too few
informative sites, every region read ``parent_free_support = no``, and the provenance said
nothing; a single-state HMM read as a clean negative.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from tessera.recomb.run import RecombParams, run_recomb


def _rows(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(), delimiter="\t"))


def _divergent(example_data: Path, tmp_path: Path, logger, **kwargs) -> Path:
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=example_data / "divergent.msa.fasta", output=out, query="query",
                     window_size=300, window_step=30, plot_format="png", **kwargs),
        logger,
    )
    return out


def test_no_phi_is_reported_as_not_run(example_data, tmp_path, logger) -> None:
    out = _divergent(example_data, tmp_path, logger, phi=False)
    report = (out / "report.html").read_text()
    # divergent.msa.fasta has 159 informative columns; the old text blamed the alignment.
    assert "Too few informative sites" not in report
    assert "disabled with --no-phi" in report
    assert "not evidence against recombination" in report

    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["recombination signal (PHI)"] == "not run (disabled with --no-phi)"

    regions = _rows(out / "recombination_regions.tsv")
    assert regions, "the divergent example has a called region"
    assert {r["parent_free_support"] for r in regions} == {"not run"}  # was "no"
    methods = _rows(out / "recombination_methods.tsv")
    assert {r["parent_free_support"] for r in methods} == {"not run"}


def test_no_phi_on_a_single_caller_run_still_says_not_run(
    example_data, tmp_path, logger
) -> None:
    """Review focus 1: with one caller there is no method table, so the regions TSV is
    the only machine-readable place that says the test did not run."""
    out = _divergent(example_data, tmp_path, logger, phi=False, methods=("hmm",))
    assert not (out / "recombination_methods.tsv").exists()
    regions = _rows(out / "recombination_regions.tsv")
    assert regions and {r["parent_free_support"] for r in regions} == {"not run"}


def test_phi_that_ran_still_reports_yes_or_no(example_data, tmp_path, logger) -> None:
    """Guard (passes before and after): with the test on, the columns keep yes/no."""
    out = _divergent(example_data, tmp_path, logger, phi_window=20)
    regions = _rows(out / "recombination_regions.tsv")
    assert {r["parent_free_support"] for r in regions} <= {"yes", "no"}
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["recombination signal (PHI)"].startswith("p=")


def test_too_few_informative_sites_is_recorded_in_the_provenance(
    tmp_path: Path, logger
) -> None:
    """Below ten informative columns the test cannot run; the provenance said nothing."""
    import random

    rng = random.Random(1)
    base = [rng.choice("ACGT") for _ in range(3000)]
    seqs = {name: list(base) for name in ("query", "A", "B", "C")}
    for col in (500, 1500, 2500):  # three informative columns, split 2/2
        alt = "A" if base[col] != "A" else "C"
        seqs["query"][col] = seqs["A"][col] = alt
    msa = tmp_path / "few.fasta"
    msa.write_text("".join(f">{k}\n{''.join(v)}\n" for k, v in seqs.items()))
    out = tmp_path / "out"
    run_recomb(RecombParams(msa=msa, output=out, query="query", plot_format="png"), logger)
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["recombination signal (PHI)"].startswith("not run (fewer than 10 informative")
    assert "not run" in (out / "report.html").read_text()
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_not_run_reporting.py -q`
Expected: 3 failed, 1 passed. Observed on `4b35e39`:
- `test_no_phi_is_reported_as_not_run`: `assert 'Too few informative sites' not in '<!DOCTYPE h...'` -- the report contains "Too few informative sites in the alignment for a parent-free recombination test."
- `test_no_phi_on_a_single_caller_run_still_says_not_run`: `assert {'no'} == {'not run'}`.
- `test_too_few_informative_sites_is_recorded_in_the_provenance`: `KeyError: 'recombination signal (PHI)'`.
- `test_phi_that_ran_still_reports_yes_or_no` passes (guard).

- [ ] **Step 3: Implement**

`src/tessera/recomb/diagnostics.py:43-44`, after `_MIN_INFORMATIVE = 10`:

```python
# Below this many informative columns the statistics are not meaningful.
_MIN_INFORMATIVE = 10
# Public name for the reporting layer, which states the floor when the test cannot run.
MIN_INFORMATIVE = _MIN_INFORMATIVE
```

`src/tessera/recomb/report_context.py`: change the import to `from dataclasses import dataclass, field`; after `from .typing import LineageMap` add

```python

# The key under which the parent-free PHI test is named in ``ReportContext.not_run``.
PHI = "phi"
```

replace the two fields at lines 37-41

```python
    # Selected callers that could not run (barcode on an untyped panel). They stay in
    # ``methods_run`` so the method table keeps a column for them, marked "not run".
    methods_not_run: tuple[str, ...] = ()
    # Why they could not run, in the run's own words (the same text the log gives).
    not_run_reason: str = ""
```

with

```python
    # Callers and tests that produced no result, by name, with the reason in the run's
    # own words (the text the log gives). A caller named here stays in ``methods_run``,
    # so the method table keeps a column for it, marked "not run"; ``PHI`` is the
    # parent-free test. Whatever is named here is reported as "not run", never as a
    # negative.
    not_run: dict[str, str] = field(default_factory=dict)
```

and add before the `gaps` property:

```python
    @property
    def methods_not_run(self) -> tuple[str, ...]:
        """The selected callers that could not run, in ``methods_run`` order."""
        return tuple(m for m in self.methods_run if m in self.not_run)

    @property
    def phi_not_run(self) -> str | None:
        """Why the parent-free test gave no result, or ``None`` when it did."""
        return self.not_run.get(PHI)
```

`src/tessera/recomb/report.py`: import `from .report_context import PHI, ReportContext`, add `"PHI",` to `__all__`, and change the two TSV calls:

```python
    write_regions_tsv(regions, output_dir, logger, phi_not_run=ctx.phi_not_run is not None)
    companions.append("recombination_regions.tsv")
    if len(ctx.methods_run) > 1 and ctx.method_breakdown is not None:
        write_methods_tsv(
            ctx.method_breakdown, ctx.methods_run, output_dir, logger, ctx.methods_not_run,
            phi_not_run=ctx.phi_not_run is not None,
        )
```

`src/tessera/recomb/report_text.py`: add before `write_regions_tsv`

```python
def _parent_free_cell(flag: bool, phi_not_run: bool) -> str:
    """``parent_free_support`` for one row: ``not run`` when the PHI test gave no result.

    Without a PHI result no region can be corroborated, so ``no`` would read as a test
    that looked and found nothing.
    """
    if phi_not_run:
        return "not run"
    return "yes" if flag else "no"
```

give `write_regions_tsv` the parameter `phi_not_run: bool = False` (signature `def write_regions_tsv(regions: list[Region], output_dir: Path, logger: logging.Logger, phi_not_run: bool = False) -> None:`) and replace its last cell `"yes" if r.parent_free_support else "no",` with `_parent_free_cell(r.parent_free_support, phi_not_run),`. In `write_methods_tsv` add `phi_not_run: bool = False` after `methods_not_run`, add to the docstring "``parent_free_support`` is ``not run`` for the same reason when the PHI test gave no result.", and replace `"yes" if b["parent_free_support"] else "no",` with `_parent_free_cell(b["parent_free_support"], phi_not_run),`.

`src/tessera/recomb/report_html.py`: import `from .report_context import PHI, ReportContext`. Replace `_method_comparison_html` lines 279-294 (signature, docstring and the `not_run_note` block) with

```python
def _method_comparison_html(
    breakdown: list[dict], methods_run: tuple[str, ...], per_major: dict[str, str],
    lineage_map: LineageMap | None = None, not_run: dict[str, str] | None = None,
) -> str:
    """A compact region x method agreement matrix; omitted for a single-method run.

    ``not_run`` maps a caller (or ``PHI``) that gave no result to the reason; its
    column reads "not run", never a negative.
    """
    if len(methods_run) < 2:
        return ""
    not_run = not_run or {}
    methods_not_run = tuple(m for m in methods_run if m in not_run)
    phi_not_run = PHI in not_run
    not_run_note = "".join(
        f'<p class="cap"><strong>Not run:</strong> {html.escape(m)} '
        f'({html.escape(not_run[m])}). It tested nothing; its column below is not a '
        f'negative result.</p>'
        for m in methods_not_run
    )
```

and in its row loop replace `phi = "&check;" if b["parent_free_support"] else "&middot;"` with

```python
        phi = (
            "not run" if phi_not_run
            else "&check;" if b["parent_free_support"] else "&middot;"
        )
```

Replace `_method_section` (lines 330-345) signature and call with

```python
def _method_section(
    method_breakdown: list[dict] | None, methods_run: tuple[str, ...],
    per_major: dict[str, str] | None, lineage_map: LineageMap | None,
    not_run: dict[str, str] | None = None,
) -> str:
    """Wrap the method-comparison table in a report section (empty for one method)."""
    body = _method_comparison_html(
        method_breakdown or [], methods_run, per_major or {}, lineage_map, not_run,
    )
```

(the rest of the function is unchanged). Replace the head of `_signal_html` (lines 501-507) with

```python
def _signal_html(
    signal: RecombinationSignal | None, alpha: float = 0.05, not_run: str | None = None
) -> str:
    """The parent-free recombination-signal section: PHI p-value, Rmin, profile.

    ``not_run`` is the run's reason when the test gave no result (``--no-phi``, or too
    few informative sites). That is stated as such: a test that did not run is not a
    negative, and the reason must be the real one -- ``--no-phi`` is not a poor alignment.
    """
    if signal is None:
        reason = not_run or "no result was recorded"
        return (
            f'<p class="cap">The parent-free recombination test was <strong>not run</strong>: '
            f'{html.escape(reason)}. No region is corroborated or contradicted by it; the '
            f'PHI/Rmin column of the method comparison and <span class="mono">'
            f'parent_free_support</span> in the TSVs read <em>not run</em>. This is not '
            f'evidence against recombination.</p>'
        )
```

In `write_html_report` change the two call sites:

```python
    method_section = _method_section(
        ctx.method_breakdown, ctx.methods_run, ctx.per_major, lineage_map, ctx.not_run,
    )
```

```python
        f"{_signal_html(ctx.signal, ctx.alpha, ctx.phi_not_run)}</section>"
```

`src/tessera/recomb/run.py`: imports become `from .diagnostics import MIN_INFORMATIVE, corroborating_intervals, recombination_signal` and `from .report import (PHI, ReportContext, ...)`. Replace lines 304-317 (from `signal = None` to the `elif signal is not None:` line) with

```python
    # A test or caller that gives no result is recorded here, with the reason, and
    # reported as "not run" -- never as a negative. The report, the TSVs and the
    # provenance all read this one mapping.
    not_run: dict[str, str] = {}
    signal = None
    if not params.phi:
        not_run[PHI] = "disabled with --no-phi"
    else:
        signal = recombination_signal(
            bp_result.rows, query_label, bp_result.column_to_query,
            window=params.phi_window,
        )
        if signal is None:
            not_run[PHI] = (
                f"fewer than {MIN_INFORMATIVE} informative sites in the alignment, too few "
                "for the test"
            )
            logger.info("Recombination signal (parent-free): not run (%s).", not_run[PHI])
        elif signal.phi_p is None:
            not_run[PHI] = (
                f"not testable: every pair of the {signal.n_informative} informative sites "
                f"lies within the window of {signal.phi_window} ranks"
            )
            logger.info(
                "Recombination signal (parent-free): PHI not testable (every pair of the "
                "%d informative site(s) lies within the window of %d ranks; lower "
                "--phi-window), Rmin=%d.",
                signal.n_informative, signal.phi_window, signal.rmin,
            )
        else:
```

(the `logger.info("... PHI p=%.4g ...")` body under it is unchanged). Replace the barcode block at lines 359-378 with

```python
    callers_not_run = tuple(
        m for m in params.methods if m == "barcode" and majors[m] is None
    )
    if callers_not_run:
        reason = (
            "fewer than two typed clades carry enough characteristic markers"
            if lineage_map else
            "it needs typed references (a lineage map: --lineage-map, or a lineages.tsv "
            "beside the output or the MSA) and none was found"
        )
        for m in callers_not_run:
            not_run[m] = reason
        if len(callers_not_run) == len(params.methods):
            raise UserInputError(
                f"The barcode caller could not run: {reason}. No other caller was "
                "selected, so this scan could not test for recombination. Supply typed "
                "references or choose another --method."
            )
        logger.warning(
            "The barcode caller could not run (%s); it is reported as 'not run', not as "
            "a negative.", reason,
        )
    n_ran = len(params.methods) - len(callers_not_run)
```

At line 398 `gate_lowered = bool(not_run) and ...` becomes `gate_lowered = bool(callers_not_run) and min_agree < params.min_methods`. Replace lines 523-524 with

```python
    reported_not_run = [m for m in params.methods if m in not_run]
    if reported_not_run:
        provenance["callers not run"] = "; ".join(
            f"{m} ({not_run[m]})" for m in reported_not_run
        )
```

(for the barcode alone this gives the same `barcode (<reason>)` text as before). After the `if signal is not None:` block (line 535) add

```python
    else:
        provenance["recombination signal (PHI)"] = f"not run ({not_run[PHI]})"
```

and in the `ReportContext(...)` call replace `methods_not_run=not_run, not_run_reason=not_run_reason, alpha=params.alpha,` with `not_run=not_run, alpha=params.alpha,`.

`docs/detection-methods.md`, PHI section: replace the sentence "Should the test still be undefined, Tessera reports `not testable` (`NA` in the `recombination_profile.tsv` header) rather than a non-significant result, and no region is flagged `parent_free_support`." with

```markdown
Should the test still be undefined, Tessera reports `not testable` (`NA` in the
`recombination_profile.tsv` header) rather than a non-significant result. When the test
gives no result -- not testable, fewer than ten informative sites, or turned off with
`--no-phi` -- the report says so with the reason, `run_provenance.json` records it under
`recombination signal (PHI)` as `not run (<reason>)`, and `parent_free_support` reads
`not run` in both TSVs and in the method comparison. That is not a negative: no region
was corroborated or contradicted.
```

and in the output-file table add to the `recombination_regions.tsv` row "(`parent_free_support` is `not run` when the PHI test gave no result)" and to the `recombination_methods.tsv` row "; the parent-free flag is `not run` likewise".

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_not_run_reporting.py tests/unit/test_barcode.py tests/unit/test_report_signal.py tests/unit/test_output_contract.py -q`
Expected: all pass. The barcode tests read `run["callers not run"]` and the methods TSV, both unchanged for the barcode-only case.

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: clean, clean, baseline + 4 passed.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/report_context.py src/tessera/recomb/diagnostics.py \
        src/tessera/recomb/run.py src/tessera/recomb/report.py \
        src/tessera/recomb/report_text.py src/tessera/recomb/report_html.py \
        docs/detection-methods.md tests/unit/test_not_run_reporting.py
git commit -F - <<'EOF'
Report a PHI test that did not run as "not run", with its reason

With --no-phi the report said "Too few informative sites in the alignment for a
parent-free recombination test" -- on example_data/divergent.msa.fasta, which has 159
informative columns and a PHI p of about 0.001 when the test runs. Every region then
read parent_free_support = no, as if the test had looked and found nothing, and
run_provenance.json recorded nothing about the test at all. The same "no" was written
when the test ran but was not testable at its window, and when there were fewer than
ten informative sites.

The run now records every test or caller that gave no result in one mapping, name to
reason, which replaces plan B's barcode-only pair of fields. The report states the
reason, the provenance records "not run (<reason>)", and parent_free_support reads
"not run" in both TSVs and the method table. Only text changes: which regions are
called, their coordinates and p-values are untouched.

Audit 4, item B1.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 3: Report a one-state HMM as "not run"  (audit items: B7)

When lineage clustering pools every reference into one lineage (`regions.py:279-283`), the HMM decodes a single state, cannot switch, and returns a non-`None` major; `run.py:359` only knows the barcode caller as "not run", so the run reads as a clean negative. Reproduced while writing this plan: five references 0.25 % apart, a 20 kb query carrying a real 6 kb tract from one of them, `--method hmm --no-informative-sites`: no region, verdict "No recombination among the present references", beside a PHI p of 0.000999 in the same report. In the default ensemble the site callers find the tract and the HMM column reads `no`.

The fix records the one-state decode in the Task 2 `not_run` mapping; the decode, its major and its (empty) regions are unchanged. Per the Global Constraints, the HMM is **not** removed from the agreement-gate count.

**Files:**
- Modify: `src/tessera/recomb/regions.py:223-256` (`call_regions`), `:259-283` (`_call_regions_hmm`)
- Modify: `src/tessera/recomb/run.py:350` (the caller loop), after the barcode block (`n_ran`)
- Modify: `src/tessera/recomb/report_html.py:55` (new `_not_run_caveat` before `_summary`), `:126-135` (`_verdict_html`), `:606` and `:629-630` (`write_html_report`)
- Modify: `docs/detection-methods.md` (the lineage-clustering limitation paragraph)
- Modify: `tests/unit/test_not_run_reporting.py`

**Interfaces:**
- Consumes: `ReportContext.not_run`, `ReportContext.methods_not_run` (Task 2).
- Produces: `call_regions(result, analysis, window_size, params, *, not_run: dict[str, str] | None = None)`; `regions.HMM_ONE_STATE: str`; `_not_run_caveat(methods_run: tuple[str, ...], not_run: dict[str, str]) -> str`; `_verdict_html(..., untested: bool = False)`.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_not_run_reporting.py`:

```python


def _one_lineage_panel(tmp_path: Path) -> Path:
    """Five references 0.25 % apart and a query carrying a real 6 kb tract from one.

    Under base-pair windowing every pair of references stays above the 98.5 % pooling
    floor in every window, so lineage clustering pools them into a single HMM state.
    """
    import random

    rng = random.Random(3)
    anc = [rng.choice("ACGT") for _ in range(20000)]

    def mut(seq: list[str], rate: float) -> list[str]:
        return [rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c
                for c in seq]

    refs = {f"ref_{i}": mut(anc, 0.0025) for i in range(5)}
    query = list(refs["ref_0"])
    query[8000:14000] = refs["ref_3"][8000:14000]
    msa = tmp_path / "pool.fasta"
    msa.write_text(">query\n" + "".join(query) + "\n"
                   + "".join(f">{k}\n{''.join(v)}\n" for k, v in refs.items()))
    return msa


def test_single_state_hmm_is_reported_as_not_run(tmp_path: Path, logger) -> None:
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=_one_lineage_panel(tmp_path), output=out, query="query",
                     plot_format="png", methods=("hmm",), informative_sites=False),
        logger,
    )
    assert _rows(out / "recombination_regions.tsv") == []  # unchanged: nothing called
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["callers not run"].startswith("hmm (lineage clustering pooled every reference")
    report = (out / "report.html").read_text()
    assert "No recombination among the present references" not in report
    assert "not a negative result" in report
    assert "--no-cluster-lineages" in report


def test_single_state_hmm_is_not_run_in_the_method_table(tmp_path: Path, logger) -> None:
    """In an ensemble the other callers still call the tract; the HMM column must not
    read as a caller that looked and found nothing."""
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=_one_lineage_panel(tmp_path), output=out, query="query",
                     plot_format="png", informative_sites=False),
        logger,
    )
    rows = _rows(out / "recombination_methods.tsv")
    assert rows, "the site callers find the 6 kb tract"
    assert {r["hmm"] for r in rows} == {"not run"}  # was "no"
    regions = _rows(out / "recombination_regions.tsv")
    assert all("hmm" not in r["methods"] for r in regions)


def test_barcode_and_one_state_hmm_are_both_reported_not_run(tmp_path: Path, logger) -> None:
    """Review focus 2: two callers not run in one run, each with its own reason. Only
    the barcode lowers the agreement gate (a one-state HMM still counts; see plan G)."""
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=_one_lineage_panel(tmp_path), output=out, query="query",
                     plot_format="png", informative_sites=False,
                     methods=("hmm", "3seq", "barcode"), min_methods=2),
        logger,
    )
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["callers not run"].startswith("hmm (lineage clustering")
    assert "; barcode (it needs typed references" in run["callers not run"]
    assert run["min methods (agreement gate)"] == "2"
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_not_run_reporting.py -q -k "single_state or one_state"`
Expected: 3 failed. Observed on `4b35e39` (the Task 2 change does not affect these): `KeyError: 'callers not run'` (single-caller run), `assert {'no'} == {'not run'}` (ensemble), and for the third test `AssertionError: assert False` from `startswith('hmm (lineage clustering')`, because `callers not run` lists only `barcode (it needs typed references ...)`.

- [ ] **Step 3: Implement**

`src/tessera/recomb/regions.py`: `call_regions` gains a keyword-only parameter and passes it on:

```python
def call_regions(
    result: WindowSimilarity,
    analysis: AnalysisResult,
    window_size: int,
    params: RegionParams,
    *,
    not_run: dict[str, str] | None = None,
) -> tuple[list[Region], str | None, list[SiblingEvidence]]:
```

Append to its docstring: "A caller that ran but could not test (the HMM left with one state) records why in ``not_run`` under its own name, when a mapping is given; its result is unchanged." and change the first branch to `regions, major, siblings = _call_regions_hmm(result, window_size, params, not_run)`.

Before `def _call_regions_hmm` add

```python
# Why the HMM could not test when lineage clustering left it a single state.
HMM_ONE_STATE = (
    "lineage clustering pooled every reference into one lineage, so the HMM had no "
    "second state to switch to; rerun with --no-cluster-lineages"
)
```

give `_call_regions_hmm` the parameter `not_run: dict[str, str] | None = None` after `params`, and directly after the clustering block (after `work, cluster_size = clustered_view(result, clusters)`, at the indentation of the `if params.cluster_lineages ...` statement) add

```python
    if len(work.similarities) < 2 and not_run is not None:
        # One state cannot switch, so the decode below cannot call a region whatever
        # the data. The run continues unchanged; the caller is reported as not run.
        not_run["hmm"] = HMM_ONE_STATE
```

`src/tessera/recomb/run.py`: the caller loop passes the mapping

```python
        regs, major_m, sibs = call_regions(
            res, ana, params.window_size, _region_params(method), not_run=not_run
        )
```

and directly before `n_ran = len(params.methods) - len(callers_not_run)` add

```python
    if "hmm" in not_run:
        logger.warning(
            "The hmm caller could not test (%s); it is reported as 'not run', not as a "
            "negative.", not_run["hmm"],
        )
    # Only the barcode caller lowers the agreement gate. A one-state HMM is reported as
    # not run but still counted here: lowering the gate for it would change which
    # regions are reported, which this reporting change does not do.
```

The provenance block from Task 2 (`reported_not_run`) already lists every selected caller in `not_run`, in `--method` order.

`src/tessera/recomb/report_html.py`: add before `_summary`

```python
def _not_run_caveat(methods_run: tuple[str, ...], not_run: dict[str, str]) -> str:
    """A header note naming every selected caller that could not test, with the reason.

    Shown whatever the number of callers: the method comparison, which also marks them,
    is omitted for a single-caller run, and a caller that tested nothing must not leave
    the verdict reading as a negative.
    """
    callers = [m for m in methods_run if m in not_run]
    if not callers:
        return ""
    items = "; ".join(
        f'<strong>{html.escape(m)}</strong> ({html.escape(not_run[m])})' for m in callers
    )
    lead = (
        "No selected caller could test this alignment" if len(callers) == len(methods_run)
        else "Not every selected caller could test this alignment"
    )
    return (
        f'<div class="caveat"><span class="ic">&#9888;</span><div><strong>{lead}.</strong> '
        f'Not run: {items}. A caller that did not run tested nothing; a region it does not '
        f'list is not a negative result from it.</div></div>'
    )
```

Change `_verdict_html`'s signature and first branch:

```python
def _verdict_html(
    s: dict, query: str, colors: dict[str, str], lineage_map: LineageMap | None = None,
    query_lineage: str | None = None, untested: bool = False,
) -> str:
    major = html.escape(typed(s["major"], lineage_map))
    q_clause = _query_lineage_clause(query, query_lineage)
    if s["n_regions"] == 0 and untested:
        lead = ('No recombinant region was called, but no selected caller could test this '
                'alignment, so this is <strong>not a negative result</strong> (see the note '
                'below).')
        return f'<p class="verdict">{lead}{q_clause}{_absent_clause(s)}</p>'
    if s["n_regions"] == 0:
```

In `write_html_report`, after `s = _summary(...)` add

```python
    untested = bool(ctx.methods_run) and len(ctx.methods_not_run) == len(ctx.methods_run)
```

and in the header replace the verdict line with

```python
        f"{_verdict_html(s, result.query, colors, lineage_map, ctx.query_lineage, untested)}"
        f"{_not_run_caveat(ctx.methods_run, ctx.not_run)}"
```

`docs/detection-methods.md`, the lineage-clustering limitation paragraph: after "run with `--no-cluster-lineages` to keep the HMM's vote." add

```markdown
When clustering does leave the HMM a single state, the run says so: the HMM is reported
as `not run` with this reason in the report header, the method comparison and
`run_provenance.json` (`callers not run`), not as a caller that found nothing. It still
counts toward `--min-methods`.
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_not_run_reporting.py tests/unit/test_barcode.py tests/unit/test_clusters.py tests/unit/test_regions.py -q`
Expected: all pass.

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/regions.py src/tessera/recomb/run.py \
        src/tessera/recomb/report_html.py docs/detection-methods.md \
        tests/unit/test_not_run_reporting.py
git commit -F - <<'EOF'
Report a one-state HMM as not run instead of as a clean negative

When lineage clustering pools every reference into one lineage, the HMM decodes a
single state and cannot call a region whatever the data. The run reported that as a
negative. Reproduced with five references 0.25 % apart and a real 6 kb tract under
--method hmm --no-informative-sites: no region, "No recombination among the present
references", beside a PHI p of 0.000999 in the same report. In the default ensemble the
site callers find the tract and the HMM column read "no".

The HMM now records the one-state decode in the run's not-run mapping. The report
header names it with the reason and the remedy (--no-cluster-lineages), the method
table and run_provenance.json mark it "not run", and when no selected caller could test
the verdict says that this is not a negative result. The decode itself, its major
parent and its (empty) regions are unchanged.

Trade-off: unlike an unrunnable barcode caller, a one-state HMM still counts toward
--min-methods, because lowering the gate would change which regions are reported.
That decision is left to plan G.

Audit 4, item B7.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 4: A low-information coverage gap is not a possible missing reference  (audit items: B3; the stale breakpoint wording at `report_text.py:157-158`)

`report_context.py:57-65` leaves only the `breakpoint` kind out of the caveat set, so a `low_information` gap -- "few comparable bases, so the gap is uncertain rather than informative" (`coverage.py:16-17`) -- reaches `_caveat_html` and the headline advises running `find-references` to search NCBI. The console table (`report_text.py:153-155`) and the log (`run.py:460-466`) say the same. Reproduced: a query whose 2.5 kb stretch is 97 % `N` gives one `low_information` gap and the headline "Possible missing reference ... run tessera find-references to search NCBI". Recruitment calls `call_coverage_gaps` directly and is unaffected. The mosaic and the plots keep hatching both kinds (`caveat_gaps` is unchanged), since a stretch with no comparable bases is worth seeing.

The breakpoint wording in the console (`report_text.py:157-158`) and the coverage intro (`report_html.py:457-459`) still described the superseded gap-wide rule ("the two parents together explain the query there"); `mark_breakpoint_gaps` now requires a single-switch mosaic to reach the threshold in every under-covered window. Both are reworded here, beside the code they explain.

**Files:**
- Modify: `src/tessera/recomb/coverage.py:37-38`
- Modify: `src/tessera/recomb/report_context.py:16`, `:56-65`
- Modify: `src/tessera/recomb/report_html.py:451-460` (new `_low_information_html` before `_coverage_html`; intro text), `:630` (header)
- Modify: `src/tessera/recomb/report_text.py:15`, `:153-158`
- Modify: `src/tessera/recomb/run.py:16-23`, `:460-466`
- Modify: `docs/detection-methods.md` (output-file table, `coverage_gaps.tsv` row)
- Create: `tests/unit/test_coverage_caveat.py`

**Interfaces:**
- Produces: `coverage.DIVERGENT_KIND = "divergent"`, `coverage.LOW_INFORMATION_KIND = "low_information"`; `ReportContext.missing_reference_gaps -> list[CoverageGap]`, `ReportContext.low_information_gaps -> list[CoverageGap]`; `report_html._low_information_html(gaps: list[CoverageGap]) -> str`.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_coverage_caveat.py`:

```python
"""A low-information coverage gap is not reported as a possible missing reference.

Audit 4, item B3: a stretch of the query that is mostly ``N`` has too few comparable
bases to judge, which ``coverage.py`` labels ``low_information``. The report headline
still called it a possible missing reference and advised searching NCBI.
"""

from __future__ import annotations

import csv
import random
from pathlib import Path

from tessera.recomb.run import RecombParams, run_recomb


def _query_with_an_unsequenced_stretch(tmp_path: Path) -> Path:
    """2.5 kb of very poor sequencing: 97 % N, the surviving bases random miscalls."""
    rng = random.Random(11)
    anc = [rng.choice("ACGT") for _ in range(8000)]

    def mut(seq: list[str], rate: float) -> list[str]:
        return [rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c
                for c in seq]

    pa, pb, og = mut(anc, 0.05), mut(anc, 0.05), mut(anc, 0.12)
    query = list(pa)
    for i in range(3000, 5500):
        query[i] = "N" if rng.random() < 0.97 else rng.choice("ACGT")
    msa = tmp_path / "lowinfo.fasta"
    msa.write_text("".join(
        f">{name}\n{''.join(seq)}\n"
        for name, seq in (("query", query), ("parent_A", pa), ("parent_B", pb),
                          ("outgroup", og))
    ))
    return msa


def test_low_information_gap_is_not_a_missing_reference(tmp_path: Path, logger) -> None:
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=_query_with_an_unsequenced_stretch(tmp_path), output=out,
                     query="query", plot_format="png"),
        logger,
    )
    lines = (out / "coverage_gaps.tsv").read_text().splitlines()[1:]  # skip the preamble
    kinds = {row["kind"] for row in csv.DictReader(lines, delimiter="\t")}
    assert kinds == {"low_information"}, kinds  # the fixture holds no divergent stretch
    report = (out / "report.html").read_text()
    assert "Possible missing reference" not in report
    assert "find-references" not in report.split("</header>")[0]
    assert "too few comparable bases" in report.split("</header>")[0]


def test_console_table_does_not_advise_a_search_for_a_low_information_gap() -> None:
    from tessera.recomb.coverage import CoverageGap
    from tessera.recomb.report_text import print_coverage

    lines: list[str] = []
    gap = CoverageGap(msa_start=0, msa_end=1000, query_start=0, query_end=1000,
                      length_bp=1000, n_windows=5, best_label="ref", mean_best=0.4,
                      kind="low_information")
    print_coverage([gap], 0.95, echo=lines.append)
    text = "\n".join(lines)
    assert "find-references" not in text
    assert "too few comparable bases" in text


def test_divergent_and_low_information_gaps_get_separate_notes() -> None:
    """Review focus 3: with both kinds present the missing-reference caveat counts only
    the divergent gap, and the low-information note still appears."""
    from tessera.recomb.coverage import CoverageGap
    from tessera.recomb.report_context import ReportContext
    from tessera.recomb.report_html import _caveat_html, _low_information_html

    divergent = CoverageGap(msa_start=0, msa_end=1500, query_start=0, query_end=1500,
                            length_bp=1500, n_windows=6, best_label="ref", mean_best=0.70,
                            kind="divergent")
    low = CoverageGap(msa_start=3000, msa_end=5500, query_start=3000, query_end=5500,
                      length_bp=2500, n_windows=20, best_label="ref", mean_best=0.40,
                      kind="low_information")
    ctx = ReportContext(coverage_gaps=[divergent, low], coverage_threshold=0.95)
    caveat = _caveat_html(ctx.missing_reference_gaps, 0.95)
    assert "1 region" in caveat and "1.5 kb" in caveat
    assert "Too little sequence to judge" in _low_information_html(ctx.low_information_gaps)
    assert ctx.caveat_gaps == [divergent, low]  # the mosaic and plots still hatch both
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_coverage_caveat.py -q`
Expected: 3 failed. Observed on `4b35e39`: `assert 'Possible missing reference' not in ...` (the header reads "Possible missing reference. 1 region (2.9 kb of the query) is poorly covered ... run tessera find-references to search NCBI."); the console test fails on "find-references" in "^ the closest reference here is poor ... Run 'tessera find-references' to search NCBI."; the third test fails with `ImportError: cannot import name '_low_information_html' from 'tessera.recomb.report_html'`.

- [ ] **Step 3: Implement**

`src/tessera/recomb/coverage.py`, after `BREAKPOINT_KIND = "breakpoint"`:

```python
# The two kinds call_coverage_gaps assigns. Only a divergent gap suggests a missing
# reference; a low-information gap has too few comparable bases to say anything.
DIVERGENT_KIND = "divergent"
LOW_INFORMATION_KIND = "low_information"
```

`src/tessera/recomb/report_context.py`: import `from .coverage import BREAKPOINT_KIND, DIVERGENT_KIND, LOW_INFORMATION_KIND, CoverageGap`; after the `caveat_gaps` property add

```python
    @property
    def missing_reference_gaps(self) -> list[CoverageGap]:
        """The gaps that may mean a missing reference in the headline: ``divergent`` only.

        A ``low_information`` gap has too few comparable bases to judge -- an unsequenced
        or masked stretch of the query, say -- so it says nothing about the panel and is
        given its own note instead of the advice to search for a reference.
        """
        return [g for g in self.gaps if g.kind == DIVERGENT_KIND]

    @property
    def low_information_gaps(self) -> list[CoverageGap]:
        """The gaps with too few comparable bases to judge."""
        return [g for g in self.gaps if g.kind == LOW_INFORMATION_KIND]
```

`src/tessera/recomb/report_html.py`: add before `_coverage_html`

```python
def _low_information_html(gaps: list[CoverageGap]) -> str:
    """A note for stretches with too few comparable bases: uncertain, not a gap in the panel."""
    if not gaps:
        return ""
    n = len(gaps)
    word = "stretch" if n == 1 else "stretches"
    verb = "has" if n == 1 else "have"
    total = sum(max(0, g.query_end - g.query_start) for g in gaps)
    return (
        f'<div class="caveat"><span class="ic">&#9888;</span><div>'
        f'<strong>Too little sequence to judge.</strong> {n} {word} of the query '
        f'(<span class="mono">{_fmt_kb(total)}</span>) {verb} too few comparable bases '
        f'to compare with the references, for example an unsequenced or masked stretch. '
        f'Recombination there could be neither called nor excluded; this is not evidence '
        f'of a missing reference. See <em>Reference coverage</em> below.</div></div>'
    )
```

In `_coverage_html` replace the last two sentences of the intro

```python
        f'bases to judge; <strong>breakpoint</strong> = windows straddling a called '
        f'breakpoint, where the two parents together explain the query (not a missing '
        f'reference, and not counted in the caveat above).</p>'
```

with

```python
        f'bases to judge (not counted as a missing reference); <strong>breakpoint</strong> = '
        f'windows beside a called breakpoint in which a single switch between that '
        f'region\'s two parents reaches the threshold in every under-covered window (not a '
        f'missing reference, and not counted in the caveat above).</p>'
```

and in `write_html_report` replace `f"{_caveat_html(caveat_gaps, threshold)}</header>"` with

```python
        f"{_caveat_html(ctx.missing_reference_gaps, threshold)}"
        f"{_low_information_html(ctx.low_information_gaps)}</header>"
```

`src/tessera/recomb/report_text.py`: import `from .coverage import BREAKPOINT_KIND, DIVERGENT_KIND, LOW_INFORMATION_KIND, CoverageGap`; replace lines 153-158 with

```python
    if any(g.kind == DIVERGENT_KIND for g in gaps):
        echo("  divergent = the closest reference here is poor; the true source may be "
             "missing. Run 'tessera find-references' to search NCBI.")
    if any(g.kind == LOW_INFORMATION_KIND for g in gaps):
        echo("  low_information = too few comparable bases to judge (e.g. an unsequenced "
             "stretch); not evidence of a missing reference.")
    if any(g.kind == BREAKPOINT_KIND for g in gaps):
        echo("  breakpoint = windows beside a called breakpoint that a single switch "
             "between the region's two parents explains (not a missing reference).")
```

`src/tessera/recomb/run.py`: add `DIVERGENT_KIND, LOW_INFORMATION_KIND,` to the `from .coverage import (...)` list, and replace lines 460-466 with

```python
    n_poor = sum(1 for g in coverage_gaps if g.kind == DIVERGENT_KIND)
    n_low_info = sum(1 for g in coverage_gaps if g.kind == LOW_INFORMATION_KIND)
    if n_poor:
        logger.info(
            "Reference coverage: %d region(s) where the closest reference is below "
            "%.3f -- a better reference may be missing.",
            n_poor, coverage_threshold,
        )
    if n_low_info:
        logger.info(
            "Reference coverage: %d stretch(es) with too few comparable bases to judge "
            "(e.g. unsequenced); not treated as a missing reference.", n_low_info,
        )
```

`docs/detection-methods.md`, output-file table, `coverage_gaps.tsv` row: after "`low_information` (too few comparable bases to judge)" add "-- shown in the report as its own note, not as a possible missing reference --".

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_coverage_caveat.py tests/unit/test_breakpoint_gaps.py tests/unit/test_coverage.py -q`

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/coverage.py src/tessera/recomb/report_context.py \
        src/tessera/recomb/report_html.py src/tessera/recomb/report_text.py \
        src/tessera/recomb/run.py docs/detection-methods.md tests/unit/test_coverage_caveat.py
git commit -F - <<'EOF'
Do not report a low-information coverage gap as a missing reference

coverage.py labels a gap low_information when it holds too few comparable bases to
judge, which it describes as uncertain rather than informative. The report excluded only
breakpoint gaps from its headline caveat, so such a gap was announced as a possible
missing reference with the advice to search NCBI. Reproduced with a query whose 2.5 kb
stretch is 97 % N: the real cause is that the query is unsequenced there.

The headline caveat, the console hint and the log now count divergent gaps only, and a
low-information gap gets its own note saying recombination there could be neither
called nor excluded. The mosaic and the plots still hatch both kinds. The breakpoint
wording in the console and the coverage table, which still described the superseded
gap-wide rule, now describes the per-window single-switch rule mark_breakpoint_gaps
applies. Report text only; recruitment, which reads the gaps directly, is unaffected.

Audit 4, item B3 and the stale breakpoint text at report_text.py:157-158.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 5: Say which callers read `--min-region`, `--margin` and `--merge-gap`  (audit items: B2, B6 and A4 help text and documentation; B4 documentation)

Help text and documentation only; no option changes what it does. `--min-region` is read at `regions.py:316` (HMM segment, before trimming), `regions.py:442` (heuristic), `clusters.py:122` and `siblings.py:131,154`; 3SEQ, MaxChi, Bootscan, GENECONV and barcode never read it, and a reported HMM region can be shorter than it (A4: a segment of 11 200 bp passed at the default floor and the reported region was 132 bp). `--margin` is read only at `regions.py:414` and `--merge-gap` only at `regions.py:432`, both in the heuristic caller. The scoped options elsewhere carry a "(hmm)" marker; these three did not, so `--margin 0.2` on a default run changed nothing and said nothing.

`run.py:443-447` builds the coverage parameters without `merge_gap`. Passing it through would change the coverage gaps, which feed the donor-absent regions, so it is left for plan G (see Global Constraints); the help now says the coverage diagnostic ignores the option. For B4 the code is plan G's; here the documentation states that `n_windows` and `posterior` of an HMM region describe the segment.

**Overlap with plan G.** Plan G edits the same help strings, `RecombParams` comments and documentation paragraph, and changes what they must say: its Task 4 reports `n_windows` and `posterior` over the trimmed span, its Task 13 makes `--min-region` bound every reported region when given, and its Task 16 passes `--merge-gap` to the coverage diagnostic. The intended landing order is F before G (`docs/superpowers/plans/2026-10-12-audit-4-index.md`), and plan G's Tasks 4, 13 and 16 each carry an "If plan F Task 5 has landed" note with the replacement text. If plan G has landed first, do not reintroduce the statements it superseded: write this task's help text and documentation to match plan G's behaviour (`--min-region` bounds every region when given and defaults to the HMM-segment rule; `--merge-gap` also sets the coverage-gap merge distance; `n_windows` and `posterior` describe the reported span), keep the scope marker for `--margin`, and adjust the assertions in `tests/unit/test_cli_option_scope.py` to the same wording.

**Files:**
- Modify: `src/tessera/cli/cmd_recomb.py:70-83`
- Modify: `src/tessera/recomb/run.py:90-93` (comment only)
- Modify: `docs/detection-methods.md` (HMM section; output-file table)
- Create: `tests/unit/test_cli_option_scope.py`

**Interfaces:** none (help strings and comments).

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_cli_option_scope.py`:

```python
"""Options read by only some callers say so in their help.

Audit 4, items B2, B6 and A4: ``--min-region`` is read by the HMM (on the segment, before
the span is trimmed) and the heuristic caller only, and ``--margin`` and ``--merge-gap``
by the heuristic caller only, but their help did not say so. A user who set one and saw
no change could not tell a no-op from a null result. Introspects the parameters rather
than rendered ``--help``, which wraps differently across typer/rich/click versions.
"""

from __future__ import annotations

import typer

from tessera.cli.main import app


def _help(command: str, option: str) -> str:
    cmd = typer.main.get_command(app).commands[command]
    param = next(p for p in cmd.params if option in getattr(p, "opts", []))
    return " ".join((param.help or "").split())


def test_min_region_names_its_callers_and_what_it_bounds() -> None:
    text = _help("recomb", "--min-region")
    assert "(hmm, heuristic;" in text
    assert "A reported HMM region can be shorter" in text
    assert "3seq, maxchi, bootscan, geneconv and barcode do not apply it" in text


def test_margin_and_merge_gap_are_marked_heuristic_only() -> None:
    assert "(heuristic only;" in _help("recomb", "--margin")
    assert "a fraction, 0-1" in _help("recomb", "--margin")
    assert "(heuristic only;" in _help("recomb", "--merge-gap")
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_cli_option_scope.py -q`
Expected: 2 failed (`assert '(hmm, heuristic;' in 'Minimum recombinant region length in bp (default: the window size).'` and the `(heuristic only;` assertion).

- [ ] **Step 3: Implement**

`src/tessera/cli/cmd_recomb.py`, replace the three options at lines 70-83 with

```python
    min_region: int | None = typer.Option(
        None, "--min-region",
        help="Minimum length in alignment columns of an HMM segment before its reported "
        "span is trimmed to the sites that support the donor, and of a heuristic region "
        "(hmm, heuristic; also the run length used by lineage clustering and sibling "
        "exclusion). A reported HMM region can be shorter; 3seq, maxchi, bootscan, "
        "geneconv and barcode do not apply it. Default: the window size.",
    ),
    margin: float = typer.Option(
        0.0, "--margin",
        help="Minimum similarity (a fraction, 0-1) by which a minor parent must exceed "
        "the major parent to mark a window recombinant (heuristic only; the other "
        "callers ignore it).",
    ),
    merge_gap: int | None = typer.Option(
        None, "--merge-gap",
        help="Merge same-parent regions separated by at most this many bp (heuristic "
        "only; the other callers and the coverage diagnostic ignore it). Default: the "
        "window size.",
    ),
```

`src/tessera/recomb/run.py:90-93`, replace the comment `# Heuristic-only thresholds (None -> derive from the window size).` and the three fields with

```python
    # Region-length floor (None -> the window size). Read by the HMM, on the segment
    # before its span is trimmed, by its lineage clustering and sibling exclusion, and by
    # the heuristic caller; the site callers do not apply it.
    min_region: int | None = None
    # Heuristic-only thresholds. merge_gap is not passed to the coverage diagnostic,
    # which merges gaps within one window size.
    margin: float = 0.0
    merge_gap: int | None = None
```

`docs/detection-methods.md`, HMM section: after "The legacy `--method heuristic` (margin / merge-gap / min-region) is kept for comparison." add

```markdown
**Which options apply to which caller.** `--min-region` (default: the window size) is
read by the HMM and the heuristic caller, and as the run length of lineage clustering and
sibling exclusion. For the HMM it bounds the length of the segment before the reported
span is trimmed (next paragraph), so a reported HMM region can be shorter than
`--min-region`. 3SEQ, MaxChi, Bootscan, GENECONV and barcode do not apply it. `--margin`
(a similarity fraction, 0-1) and `--merge-gap` are read by the heuristic caller only; the
coverage diagnostic merges gaps within one window size whatever `--merge-gap` says. The
`--help` text of each scoped option names the callers that read it.
```

and at the end of the paragraph "**The reported span is not the segment.**" add

```markdown
`n_windows` and `posterior` describe the segment too, not the reported span: a trimmed
region can list more windows than its length holds, and its posterior is the mean over
the whole segment.
```

In the output-file table, `recombination_regions.tsv` row, after "`length_bp` / `length_msa`," add "`n_windows` and `posterior` (for an HMM region both describe the HMM segment, not the trimmed span),".

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_cli_option_scope.py tests/unit/test_cli_validation.py -q`

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`

- [ ] **Step 6: Commit**

```bash
git add src/tessera/cli/cmd_recomb.py src/tessera/recomb/run.py docs/detection-methods.md \
        tests/unit/test_cli_option_scope.py
git commit -F - <<'EOF'
Say which callers read --min-region, --margin and --merge-gap

--min-region is read by the HMM and the heuristic caller (and as a run length by
lineage clustering and sibling exclusion); 3SEQ, MaxChi, Bootscan, GENECONV and barcode
never read it, so with --min-region 5000 a 175 bp 3SEQ region is still reported. For the
HMM it bounds the segment before the span is trimmed, so a reported HMM region can be
much shorter (a 11 200 bp segment was reported as 132 bp at the default floor).
--margin and --merge-gap are read by the heuristic caller only. None of the three said
so, unlike the options that already carry a "(hmm)" marker, so a user who set one and
saw no change could not tell a no-op from a null result.

Help text, a code comment and documentation only; no option changes what it does.
Enforcing --min-region across callers, and passing --merge-gap to the coverage
diagnostic (which would change the donor-absent regions), are left to plan G. The
documentation also states that n_windows and posterior of an HMM region describe the
segment rather than the trimmed span (the documentation half of B4).

Audit 4, items B2, B6, A4 (help and documentation) and B4 (documentation).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 6: Provenance and glossary text that matches the run  (audit items: B5; smaller reporting items: coordinate key, Benjamini-Hochberg sentence, breakpoint and major-parent glossary, `example_data/README.md`)

- **B5.** `run.py:504` writes `window / step = 1000 / 100` beside `windowing = informative-site (40 sites/window, step 5)`. The base-pair scan still drives the similarity track, the coverage diagnostic and the site callers, so the line is kept and says which scan it describes.
- `report_text.py:252-255` says the coordinate conventions are in `run_provenance.json`; no key held them. A `coordinates` key now does (the convention documented in `docs/detection-methods.md`: 0-based, end-exclusive).
- `report_html.py:402-403` says every caller applies Benjamini-Hochberg; the barcode caller has no test.
- `report_assets.py:121-123` attributes every breakpoint interval to the HMM posterior; for a region only a site caller found it is that caller's own (`threeseq.py:217-218`, `maxchi.py:170-171`, `geneconv.py:121-122`: the informative sites either side of the tract start; `bootscan.py:192`, `barcode.py:164`: the region's own ends).
- `report_assets.py:99-100` (and `README.md`, "How detection works") define the major parent as the reference winning the most windows; the reported one is the HMM's majority state after sibling exclusion and pooling, reconciled across callers (`ensemble.reconcile_major`).
- `example_data/README.md:44-45` advises `--phi-window 5` for the cryptic example; the PHI window is now capped at a tenth of the informative sites, so the test runs without it (`tests/unit/test_report_signal.py` pins "54 informative sites, window 5").

**Files:**
- Modify: `src/tessera/recomb/run.py:504` (provenance dict)
- Modify: `src/tessera/recomb/report_text.py:251-255` (comment)
- Modify: `src/tessera/recomb/report_html.py:402-404`
- Modify: `src/tessera/recomb/report_assets.py:99-100`, `:107-109`, `:121-123`
- Modify: `README.md:139`, `example_data/README.md:44-45`, `docs/detection-methods.md` (Provenance section)
- Create: `tests/unit/test_report_text_faithful.py`

**Interfaces:**
- Produces: provenance keys `window / step` (text changes under informative-site windowing only) and `coordinates` (new).

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_report_text_faithful.py`:

```python
"""Report, glossary and provenance text that describes what the run did.

Audit 4, item B5 and the smaller reporting items: the provenance gave base-pair window
geometry under informative-site windowing, pointed at a coordinate convention it did not
record, and the methods text and glossary described the barcode caller, the breakpoint
interval and the major parent as something they are not.
"""

from __future__ import annotations

import json
from pathlib import Path

from tessera.recomb.report_assets import _GLOSSARY
from tessera.recomb.run import RecombParams, run_recomb


def _run(example_data: Path, tmp_path: Path, logger, **kwargs) -> tuple[dict, str]:
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=example_data / "divergent.msa.fasta", output=out, query="query",
                     window_size=300, window_step=30, plot_format="png", **kwargs),
        logger,
    )
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    return run, (out / "report.html").read_text()


def test_window_line_names_the_base_pair_scan_under_informative_windowing(
    example_data, tmp_path, logger
) -> None:
    run, _ = _run(example_data, tmp_path, logger, informative_sites=True)
    assert run["windowing"].startswith("informative-site")
    assert run["window / step"].startswith("300 / 30 (base-pair scan")


def test_window_line_is_unchanged_under_base_pair_windowing(
    example_data, tmp_path, logger
) -> None:
    """Guard (passes before and after)."""
    run, _ = _run(example_data, tmp_path, logger, informative_sites=False)
    assert run["window / step"] == "300 / 30"


def test_provenance_records_the_coordinate_convention(example_data, tmp_path, logger) -> None:
    run, _ = _run(example_data, tmp_path, logger)
    assert run["coordinates"] == "0-based, end-exclusive; MSA columns and query bases"


def test_methods_text_does_not_say_every_caller_applies_bh(
    example_data, tmp_path, logger
) -> None:
    _, report = _run(example_data, tmp_path, logger)
    assert "Each caller applies Benjamini-Hochberg" not in report
    assert "barcode caller has no significance test" in report


def test_glossary_describes_the_breakpoint_and_major_parent_faithfully() -> None:
    glossary = dict(_GLOSSARY)
    assert "site caller" in glossary["Breakpoint"]
    major = glossary["Major parent (backbone)"]
    assert "most windows overall" not in major
    assert "siblings" in major and "pooled" in major
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_report_text_faithful.py -q`
Expected: 4 failed, 1 passed (the base-pair guard). Observed on `4b35e39`: `AssertionError: assert False` where `False = '300 / 30'.startswith('300 / 30 (base-pair scan')`; `KeyError: 'coordinates'`; `'Each caller applies Benjamini-Hochberg'` found in the report; `assert 'site caller' in 'The query position where the source switches, with a posterior-derived uncertainty interval from the HMM.'`.

- [ ] **Step 3: Implement**

`src/tessera/recomb/run.py`, provenance dict: replace `"window / step": f"{params.window_size} / {params.window_step}",` with

```python
        # Under informative-site windowing the HMM did not use this geometry: the
        # base-pair scan still drives the similarity track, the coverage diagnostic and
        # the site callers, so the line says which scan it describes.
        "window / step": (
            f"{params.window_size} / {params.window_step}" if result is bp_result else
            f"{params.window_size} / {params.window_step} (base-pair scan: similarity "
            "track, coverage and the site callers; the HMM used informative-site windows)"
        ),
        "coordinates": "0-based, end-exclusive; MSA columns and query bases",
```

`src/tessera/recomb/report_text.py:252-255`, the comment ends: "The FDR method ("multiple testing") and the coordinate convention ("coordinates") are in run_provenance.json and the documentation instead."

`src/tessera/recomb/report_html.py:402-404`, replace

```python
        'Bootscan. Each caller applies Benjamini-Hochberg correction within its own '
        'candidates; nothing is corrected across callers or across the genome. This is an '
```

with

```python
        'Bootscan. Each caller with a significance test applies Benjamini-Hochberg '
        'correction within its own candidates; the opt-in barcode caller has no '
        'significance test. Nothing is corrected across callers or across the genome. '
        'This is an '
```

`src/tessera/recomb/report_assets.py`, three glossary entries:

```python
    ("Major parent (backbone)",
     "The reference the query follows over most of its length. When the HMM ran, it "
     "is the HMM's majority state after the query's own-lineage siblings are set aside "
     "and near-duplicate references are pooled into one lineage (named by one member); "
     "otherwise the backbone most callers report. It is usually, but not always, the "
     "reference that wins the most windows; run_provenance.json records it."),
```

```python
    ("Window / step",
     "The scan slides a fixed-width window along the alignment in fixed steps; each "
     "window is scored independently. On a near-identical panel the HMM instead uses "
     "windows holding a fixed number of informative (polymorphic) columns; the "
     "windowing line under Run parameters names the one used, and the window / step "
     "line then describes the base-pair scan behind the similarity track, the coverage "
     "diagnostic and the site callers."),
```

```python
    ("Breakpoint",
     "An interval for the query position where the source switches. For a region the "
     "HMM found it is derived from the HMM posterior. For a region only a site caller "
     "found it is that caller's own estimate: the informative sites either side of the "
     "tract start (3SEQ, MaxChi, GENECONV) or the region's own ends (Bootscan, "
     "barcode)."),
```

`README.md:139`, replace "The reference winning the most windows is the **major parent** (the backbone donor)." with "The **major parent** (the backbone) is the reference the query follows over most of its length: for the HMM, its majority state after the query's own-lineage siblings are set aside and near-duplicate references pooled."

`example_data/README.md:44-45`, replace "The cryptic example has 54 informative sites, fewer than the default PHI window of 100, so its PHI test is reported as `not testable`; add `--phi-window 5` to test it." with "The cryptic example has 54 informative sites, fewer than the default PHI window of 100. The window used is capped at a tenth of the informative sites (5 here), so the test runs without changing `--phi-window`, and the report states the window used."

`docs/detection-methods.md`, Provenance section, after "Its `run` block holds what the report's methods section shows;" add "it states the coordinate convention under `coordinates`, and under informative-site windowing `window / step` names the base-pair scan it describes (the HMM's windows are given under `windowing`);".

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_report_text_faithful.py tests/unit/test_provenance.py tests/unit/test_report_stale_text.py -q`

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/run.py src/tessera/recomb/report_text.py \
        src/tessera/recomb/report_html.py src/tessera/recomb/report_assets.py \
        README.md example_data/README.md docs/detection-methods.md \
        tests/unit/test_report_text_faithful.py
git commit -F - <<'EOF'
Make provenance and glossary text describe what the run did

Under informative-site windowing run_provenance.json gave "window / step 1000 / 100"
beside "windowing informative-site (40 sites/window, step 5)", with nothing to say the
first line described the base-pair scan behind the similarity track, the coverage
diagnostic and the site callers; it now says so. A comment in report_text.py pointed at
a coordinate convention in the provenance that was not there; a "coordinates" key now
records it (0-based, end-exclusive, as documented).

The methods text said every caller applies Benjamini-Hochberg; the barcode caller has
no test. The glossary attributed every breakpoint interval to the HMM posterior, which
is wrong for a region only a site caller found, and defined the major parent as the
reference winning the most windows, which the HMM's major after sibling exclusion and
pooling need not be (the README said the same). The example README still advised
--phi-window 5 for a case the PHI window cap now handles. Text only.

Audit 4, item B5 and the smaller reporting items.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 7: The headline names the run's major parent  (audit items: B9, suspected in the spec, confirmed here)

`report_html.py:61-66` takes the headline backbone from the first region's major or, with no region, from the first of the ranked datasets; the reconciled major from `run.py:380` never reaches the report. Confirmed while writing this plan: a clonal query, its identical reference with a 3 kb gapped stretch (`A_close`), and a 99.6 % relative with full coverage (`A_cov`) that therefore represents their pooled lineage in the HMM. `A_close` wins 90 windows and `A_cov` 20, so the guard in `reconcile_major` (which replaces a zero-win HMM major) does not fire: `run_provenance.json` says `major parent: A_cov` and the headline says "the query is most similar to A_close throughout". With regions present the relabelling in `ensemble._merge_group` normally makes the two agree, so the case is the no-region verdict, the cards, the mosaic legend and the highlighted statistics row, which all read `s["major"]`.

The headline wording "most similar to X throughout" becomes "follows X (the major parent) throughout", because the reconciled major is not always the most similar reference (in the reproduction it is not).

**Files:**
- Modify: `src/tessera/recomb/report_context.py:42-43` (after `per_major`)
- Modify: `src/tessera/recomb/report_html.py:55-66` (`_summary`), `:133-134` (`_verdict_html` no-region lead), `:606` (`write_html_report`)
- Modify: `src/tessera/recomb/run.py:564` (`ReportContext(...)`)
- Create: `tests/unit/test_report_backbone.py`

**Interfaces:**
- Produces: `ReportContext.major_parent: str | None = None`; `_summary(result, regions, datasets, major_parent: str | None = None) -> dict`.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_report_backbone.py`:

```python
"""The report headline names the run's reconciled major parent.

Audit 4, item B9 (suspected, confirmed while planning): with no region called, the headline
named the dataset that wins the most windows while run_provenance.json named the HMM's
major -- under lineage pooling a cluster's best-covering member, a different genome.
"""

from __future__ import annotations

import json
import random
import re
from pathlib import Path

from tessera.recomb.report_html import _summary
from tessera.recomb.run import RecombParams, run_recomb


def _panel(tmp_path: Path) -> Path:
    """A clonal query; its identical reference has a 3 kb gapped stretch, so a 99.6 %
    relative covers more of it and represents their pooled lineage in the HMM."""
    rng = random.Random(5)
    anc = [rng.choice("ACGT") for _ in range(12000)]

    def mut(seq: list[str], rate: float) -> list[str]:
        return [rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c
                for c in seq]

    a = mut(anc, 0.05)
    a_close = list(a)
    a_close[3000:6000] = ["-"] * 3000
    refs = {"A_close": a_close, "A_cov": mut(a, 0.004), "A_x": mut(a, 0.008),
            "A_y": mut(a, 0.008), "B": mut(anc, 0.05), "C": mut(anc, 0.05)}
    msa = tmp_path / "b9.fasta"
    msa.write_text(">query\n" + "".join(a) + "\n"
                   + "".join(f">{k}\n{''.join(v)}\n" for k, v in refs.items()))
    return msa


def test_headline_backbone_is_the_provenance_major_parent(tmp_path: Path, logger) -> None:
    out = tmp_path / "out"
    run_recomb(RecombParams(msa=_panel(tmp_path), output=out, query="query",
                            plot_format="png"), logger)
    major = json.loads((out / "run_provenance.json").read_text())["run"]["major parent"]
    verdict = re.search(r'class="verdict">(.*?)</p>', (out / "report.html").read_text())
    assert verdict is not None
    named = re.findall(r"<strong>(.*?)</strong>", verdict.group(1))
    assert named and named[0] == major, (named, major)


def test_summary_prefers_the_reconciled_major() -> None:
    class _Result:
        query_cumulative = [1000]

    assert _summary(_Result(), [], ["top_window_winner"], "reconciled")["major"] == "reconciled"
    assert _summary(_Result(), [], ["top_window_winner"])["major"] == "top_window_winner"
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_report_backbone.py -q`
Expected: 2 failed. Observed on `4b35e39`: `AssertionError: (['A_close'], 'A_cov')` and `TypeError: _summary() takes 3 positional arguments but 4 were given`.

- [ ] **Step 3: Implement**

`src/tessera/recomb/report_context.py`, after `per_major`:

```python
    # The run's reconciled major parent -- the one the provenance records and the
    # region rows are relabelled to. The headline names it; ``None`` when no caller
    # reported one, and the report then falls back to the leading dataset.
    major_parent: str | None = None
```

`src/tessera/recomb/report_html.py`, `_summary`:

```python
def _summary(
    result: WindowSimilarity, regions: list[Region], datasets: list[str],
    major_parent: str | None = None,
) -> dict:
    """Headline numbers: query length, backbone, donor list, recombinant fraction.

    The backbone is the run's reconciled ``major_parent`` when given, so the verdict and
    run_provenance.json name the same genome. The first region's major or the leading
    dataset is only a fallback: with lineage pooling the HMM's major is a cluster's
    representative, which need not be the reference that wins the most windows.
    """
    qcum = result.query_cumulative
    query_len = int(qcum[-1]) if len(qcum) else 0
    if major_parent:
        major = major_parent
    elif regions:
        major = regions[0].major_parent
    elif datasets:
        major = datasets[0]
    else:
        major = "n/a"
```

(the rest unchanged). In `_verdict_html`, the no-region lead becomes

```python
        lead = (f'No recombination among the present references &mdash; the query follows '
                f'<strong>{major}</strong> (the major parent) throughout.')
```

and in `write_html_report`: `s = _summary(result, regions, datasets, ctx.major_parent)`.

`src/tessera/recomb/run.py`, in the `ReportContext(...)` call, add `major_parent=major_parent,` after `per_major=per_major,`.

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_report_backbone.py tests/unit/test_report_summary.py tests/unit/test_report_typed.py -q`

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/report_context.py src/tessera/recomb/report_html.py \
        src/tessera/recomb/run.py tests/unit/test_report_backbone.py
git commit -F - <<'EOF'
Name the run's major parent in the report headline

The headline took its backbone from the first region's major parent or, with no
region, from the top-ranked dataset, while the reconciled major that the provenance
records and the region rows are relabelled to never reached the report. Confirmed on a
clonal query whose identical reference has a gapped stretch: lineage pooling makes a
fully covering 99.6 % relative the HMM's major, so run_provenance.json said A_cov while
the headline said "most similar to A_close throughout".

The report now takes the reconciled major from the run and uses the old rule only as a
fallback. Since that genome need not be the most similar reference, the no-region
verdict says the query "follows X (the major parent)" instead of "is most similar to
X". Report text only.

Audit 4, item B9 (suspected in the spec; reproduced here).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 8: One method-table row per region, donor-absent regions included  (audit items: smaller reporting item, `run.py:485`)

Donor-absent regions are appended to `regions` after the method breakdown is final (`run.py:485`), so `recombination_methods.tsv` has fewer rows than `recombination_regions.tsv` and a positional join pairs the wrong rows. Reproduced while writing this plan: a `parent_B` tract and a 2.5 kb tract from a donor absent from the panel give two region rows and one method row. The two lists are now extended and sorted together, so they stay parallel; a donor-absent row reads `no` for every caller that ran (no caller called it) and `not run` for one that did not. The region list and its order are unchanged: the pairs are sorted on the same key the regions were (`msa_start`), with a stable sort.

**Files:**
- Modify: `src/tessera/recomb/run.py:28` (import `Region`), `:224` (new `_absent_breakdown_row`), `:485`
- Modify: `src/tessera/recomb/report_html.py:313-326` (method-table donor cell)
- Modify: `docs/detection-methods.md` (output-file table, `recombination_methods.tsv` row)
- Create: `tests/unit/test_methods_table_rows.py`

**Interfaces:**
- Produces: `run._absent_breakdown_row(region: Region) -> dict` (keys `minor_parent`, `query_start`, `query_end`, `methods`, `per_method_support`, `parent_free_support`, `donor_absent`). A breakdown row may now carry `"donor_absent": True`; `write_methods_tsv` ignores the key, the HTML method table shows "missing (donor absent)".

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_methods_table_rows.py`:

```python
"""recombination_methods.tsv keeps one row per row of recombination_regions.tsv.

Audit 4, smaller reporting item: donor-absent regions were appended after the method
breakdown was final, so the methods file had fewer rows than the regions file and a
positional join paired the wrong rows.
"""

from __future__ import annotations

import csv
import random
from pathlib import Path

from tessera.recomb.run import RecombParams, run_recomb


def _panel_with_an_absent_donor(tmp_path: Path) -> Path:
    """A parent_B tract at 2000-4000 and a 2.5 kb tract from a donor not in the panel."""
    rng = random.Random(31)
    anc = [rng.choice("ACGT") for _ in range(12000)]

    def mut(seq: list[str], rate: float) -> list[str]:
        return [rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c
                for c in seq]

    pa, pb, og = mut(anc, 0.04), mut(anc, 0.04), mut(anc, 0.10)
    ghost = mut(anc, 0.30)
    query = list(pa)
    query[2000:4000] = pb[2000:4000]
    query[7000:9500] = ghost[7000:9500]
    msa = tmp_path / "absent.fasta"
    msa.write_text("".join(
        f">{n}\n{''.join(s)}\n"
        for n, s in (("query", query), ("parent_A", pa), ("parent_B", pb), ("outgroup", og))
    ))
    return msa


def _rows(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(), delimiter="\t"))


def test_methods_rows_match_regions_rows_with_a_donor_absent_region(
    tmp_path: Path, logger
) -> None:
    out = tmp_path / "out"
    run_recomb(RecombParams(msa=_panel_with_an_absent_donor(tmp_path), output=out,
                            query="query", plot_format="png"), logger)
    regions = _rows(out / "recombination_regions.tsv")
    methods = _rows(out / "recombination_methods.tsv")
    assert any(r["donor_absent"] == "yes" for r in regions), "fixture needs an absent row"
    assert len(methods) == len(regions)
    for region, row in zip(regions, methods, strict=True):
        assert (row["query_start"], row["query_end"]) == (
            region["query_start"], region["query_end"])
        if region["donor_absent"] == "yes":
            assert {row[m] for m in ("hmm", "3seq", "maxchi", "bootscan")} == {"no"}
    assert "donor absent" in (out / "report.html").read_text()


def test_methods_rows_match_regions_rows_under_the_agreement_gate(
    tmp_path: Path, logger
) -> None:
    """Review focus 4: --min-methods 2 removes rows before the donor-absent rows are
    added; the two files must still pair row for row."""
    out = tmp_path / "out"
    run_recomb(RecombParams(msa=_panel_with_an_absent_donor(tmp_path), output=out,
                            query="query", plot_format="png", min_methods=2), logger)
    regions = _rows(out / "recombination_regions.tsv")
    methods = _rows(out / "recombination_methods.tsv")
    assert any(r["donor_absent"] == "yes" for r in regions)
    assert [(r["query_start"], r["query_end"]) for r in regions] == [
        (m["query_start"], m["query_end"]) for m in methods]
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_methods_table_rows.py -q`
Expected: 2 failed. Observed on `4b35e39`: `assert 1 == 2` (one method row, two region rows) and `assert [('2024', '3989'), ('6200', '10300')] == [('2024', '3989')]`.

- [ ] **Step 3: Implement**

`src/tessera/recomb/run.py`: import `Region` (`from .regions import DEFAULT_METHODS, Region, RegionParams, call_regions`); add before `_write_run_provenance`

```python
def _absent_breakdown_row(region: Region) -> dict:
    """The method-table row of a donor-absent region: no caller called it.

    A donor-absent region comes from the coverage diagnostic, not from a caller, so every
    caller that ran reads ``no`` for it (and one that could not run, ``not run``).
    """
    return {
        "minor_parent": region.minor_parent,
        "query_start": region.query_start,
        "query_end": region.query_end,
        "methods": (),
        "per_method_support": {},
        "parent_free_support": region.parent_free_support,
        "donor_absent": True,
    }
```

and replace `regions = sorted(regions + absent, key=lambda r: r.msa_start)` (line 485) with

```python
    # One method-table row per region, donor-absent ones included, sorted together so
    # recombination_methods.tsv stays row-for-row with recombination_regions.tsv.
    paired = sorted(
        zip(regions + absent, method_breakdown + [_absent_breakdown_row(r) for r in absent],
            strict=True),
        key=lambda pair: pair[0].msa_start,
    )
    regions = [region for region, _ in paired]
    method_breakdown = [row for _, row in paired]
```

`src/tessera/recomb/report_html.py`, `_method_comparison_html` row loop: before the `phi = (...)` statement add

```python
        donor = (
            "<em>missing</em> (donor absent)" if b.get("donor_absent")
            else html.escape(typed(b["minor_parent"], lineage_map))
        )
```

and replace the first cell `f'<td class="lbl">{html.escape(typed(b["minor_parent"], lineage_map))}</td>'` with `f'<td class="lbl">{donor}</td>'`.

`docs/detection-methods.md`, output-file table, `recombination_methods.tsv` row: "one row per region" becomes "one row per row of `recombination_regions.tsv`, in the same order (a donor-absent region reads `no` for every caller that ran)".

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_methods_table_rows.py tests/unit/test_method_comparison.py tests/unit/test_output_contract.py -q`

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/run.py src/tessera/recomb/report_html.py \
        docs/detection-methods.md tests/unit/test_methods_table_rows.py
git commit -F - <<'EOF'
Keep recombination_methods.tsv row-for-row with recombination_regions.tsv

Donor-absent regions were appended to the region list after the method breakdown was
final, so the methods file had fewer rows than the regions file and a positional join
paired the wrong rows. Reproduced with a called tract plus a tract from a donor absent
from the panel: two region rows, one method row.

The breakdown now gains a row for each donor-absent region and the two lists are
sorted together on the key the regions were already sorted on, so the region list and
its order are unchanged. A donor-absent row reads "no" for every caller that ran, which
is what happened: it comes from the coverage diagnostic, not from a caller.

Audit 4, smaller reporting item (run.py:485).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 9: The pair plot colours a dataset as the rest of the report does  (audit items: B8)

`report_plots.py:153` builds the pair plot's palette from its own two labels, while the HTML, the mosaic and the top-N plot colour from `region_labels(top_datasets, regions)`, and colours are positional. Measured by the audit: one donor is green in the report and orange in the pair plot, where orange denotes another reference. `plot_pairwise` now takes the report's label order.

**Files:**
- Modify: `src/tessera/recomb/report_plots.py:134-153`
- Modify: `src/tessera/recomb/report.py:41`, `:136-137`
- Modify: `tests/unit/test_report_plots.py` (append)

**Interfaces:**
- Produces: `plot_pairwise(result, datasets, regions, output_dir, plot_format, logger, color_labels: list[str] | None = None) -> Path | None`.
- Consumes: `report_plots.region_labels`.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_report_plots.py`:

```python


# --- the pair plot colours a dataset as the rest of the report does ----------

def test_pair_plot_uses_the_report_colours(tmp_path: Path, monkeypatch) -> None:
    """Audit 4, item B8: the pair plot built its palette from its own two labels, so a
    donor ranked third was drawn in the colour the report gives the second reference."""
    import logging

    from matplotlib.axes import Axes
    from matplotlib.colors import to_rgba

    from tessera.recomb.report_plots import plot_pairwise

    drawn: dict[str, tuple] = {}
    original = Axes.plot

    def spy(self, *args, **kwargs):
        if "label" in kwargs:
            drawn[kwargs["label"]] = kwargs.get("color")
        return original(self, *args, **kwargs)

    monkeypatch.setattr(Axes, "plot", spy)
    result = compute_similarity(str(recombinant_msa(tmp_path, recombinant=True)), "query")
    regions = [_region("B")]
    report_order = region_labels(["A", "other", "B"], regions)  # B ranks third
    plot_pairwise(result, ["A", "B"], regions, tmp_path, "png", logging.getLogger("t"),
                  color_labels=report_order)
    html_colours = _color_map(report_order)
    assert drawn["B"] == to_rgba(html_colours["B"])
    assert drawn["A"] == to_rgba(html_colours["A"])


def test_write_reports_gives_the_pair_plot_the_report_colour_order(
    tmp_path: Path, monkeypatch, logger
) -> None:
    from tessera.recomb.run import RecombParams, run_recomb

    seen: dict = {}

    def fake_pairwise(result, datasets, regions, output_dir, plot_format, logger, **kwargs):
        seen.update(kwargs, datasets=datasets)
        return None

    monkeypatch.setattr("tessera.recomb.report.plot_pairwise", fake_pairwise)
    run_recomb(RecombParams(msa=recombinant_msa(tmp_path, recombinant=True),
                            output=tmp_path / "out", query="query", plot_format="png"),
               logger)
    assert "color_labels" in seen
    assert set(seen["datasets"]) <= set(seen["color_labels"])
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_report_plots.py -q`
Expected: 2 failed, 7 passed. Observed on `4b35e39`: `TypeError: plot_pairwise() got an unexpected keyword argument 'color_labels'` and `AssertionError: assert 'color_labels' in {'datasets': ['A', 'B']}`. (The first fails on the interface; the colour mismatch itself was measured by the audit's script and is what the second assertion block pins once the argument exists.)

- [ ] **Step 3: Implement**

`src/tessera/recomb/report_plots.py`, `plot_pairwise`:

```python
def plot_pairwise(
    result: WindowSimilarity,
    datasets: list[str],
    regions: list[Region],
    output_dir: Path,
    plot_format: str,
    logger: logging.Logger,
    color_labels: list[str] | None = None,
) -> Path | None:
    """Static pairwise plot of the two leading datasets, regions shaded.

    ``color_labels`` is the label order the rest of the report colours from (the HTML,
    the mosaic and the top-N plot). Colours are positional, so building the palette from
    the two plotted labels alone gave a donor the colour another reference has elsewhere.
    """
```

and replace `colors = _palette(region_labels([seq1, seq2], regions))` with

```python
    base = list(color_labels) if color_labels is not None else []
    base += [label for label in (seq1, seq2) if label not in base]
    colors = _palette(region_labels(base, regions))
```

`src/tessera/recomb/report.py`: import `region_labels` (`from .report_plots import build_interactive_figure, plot_pairwise, plot_top_n, region_labels`) and replace the pair-plot call with

```python
    # The same label order the HTML colours from, so a dataset is one colour throughout.
    plots.append(plot_pairwise(
        result, pair, regions, output_dir, plot_format, logger,
        color_labels=region_labels(top_datasets, regions),
    ))
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_report_plots.py -q`

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/report_plots.py src/tessera/recomb/report.py \
        tests/unit/test_report_plots.py
git commit -F - <<'EOF'
Colour the pair plot from the report's label order

Colours are assigned by position. The HTML, the mosaic and the top-N plot colour from
the ranked datasets followed by the region parents; the pair plot built its palette from
its own two labels. A donor ranked third was therefore green in the report and orange in
the pair plot, where orange is another reference's colour.

write_reports now passes the report's label order to plot_pairwise, which appends the
two plotted labels if they are missing from it (the fallback pair can lie outside
--top-n). Plot colours only.

Audit 4, item B8.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 10: Reject an output path that is a file, on every command  (audit items: E1)

`main.py:131-133` creates the output directory while opening the run log, and `cmd_detect.py:122`, `cmd_build_panel.py:93`, `cmd_type_lineages.py:52` and `cmd_reassort.py:56` call it before entering `stage_errors`, so an existing file at `--output` gives a raw `FileExistsError` traceback. `fill-references`, `find-references` and `curate-panel` open no run log but never check the path either, so the first writer fails later under "Unexpected error". Only `recomb` checked (`cmd_recomb.py:158`, plan B item B11). Observed while writing this plan: 7 of 7 commands fail the test below on `4b35e39` (four with `FileExistsError(17, 'File exists')`, three reaching the pipeline).

One helper is used everywhere: `_require_output_directory` (`main.py:101-108`). `get_logger(output_dir)` calls it before it creates anything, logging the same message and exiting 1; the three commands that call `get_logger()` without a directory call it first thing inside `stage_errors`. `msa -o` names a file, not a directory, and is not affected.

**Files:**
- Modify: `src/tessera/cli/main.py:131-133`
- Modify: `src/tessera/cli/cmd_fill_references.py:9-19`, `:209`
- Modify: `src/tessera/cli/cmd_find_references.py:9`, `:90`
- Modify: `src/tessera/cli/cmd_curate_panel.py:9`, `:48`
- Create: `tests/unit/test_cli_output_path.py`

**Interfaces:**
- Produces: `get_logger(output_dir: Path | None = None) -> logging.Logger` (same signature; raises `typer.Exit(1)` after logging the error when `output_dir` exists and is not a directory).

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_cli_output_path.py`:

```python
"""Every command that writes into an output directory rejects an existing file there.

``recomb`` checked this (audit B11); five other commands crashed with a raw traceback or
an "Unexpected error" from wherever the directory was first created. None of these tests
may reach the pipeline: the entry points are replaced with a function that fails the test.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

from tessera.cli.main import app

runner = CliRunner()


def _must_not_run(*args, **kwargs):
    raise AssertionError("the command reached the pipeline; it should have been rejected")


@pytest.fixture
def no_pipeline(monkeypatch):
    monkeypatch.setattr("tessera.discover.iterate.fill_references", _must_not_run)
    monkeypatch.setattr("tessera.discover.lineage_assign.assign_lineages", _must_not_run)
    monkeypatch.setattr("tessera.discover.run.find_references", _must_not_run)
    monkeypatch.setattr("tessera.discover.panel.curate_collection", _must_not_run)
    monkeypatch.setattr("tessera.cli.cmd_reassort.assign_segments", _must_not_run)


def _args(command: str, tmp_path: Path, occupied: Path) -> list[str]:
    query = tmp_path / "query.fasta"
    query.write_text(">query\nACGTACGTACGT\n>NA\nACGTACGTACGT\n")
    collection = tmp_path / "collection"
    collection.mkdir()
    (collection / "ref.fasta").write_text(">ref\nACGTACGTACGT\n")
    common = ["--output", str(occupied)]
    if command == "type-lineages":
        return [command, "--collection", str(collection), *common]
    if command == "curate-panel":
        return [command, "--query", str(query), "--collection", str(collection), *common]
    if command == "find-references":
        return [command, "--msa", "example_data/divergent.msa.fasta", "--query", "query",
                *common]
    return [command, "--query", str(query), *common]


@pytest.mark.parametrize("command", [
    "detect", "build-panel", "fill-references", "find-references", "curate-panel",
    "type-lineages", "reassort",
])
def test_command_rejects_an_output_path_that_is_a_file(
    tmp_path: Path, no_pipeline, command: str
) -> None:
    occupied = tmp_path / "out"
    occupied.write_text("not a directory\n")
    result = runner.invoke(app, _args(command, tmp_path, occupied))
    assert result.exit_code == 1, result.output
    assert result.exception is None or isinstance(result.exception, SystemExit), (
        f"raw exception: {result.exception!r}"
    )
    assert "Unexpected error" not in result.output
    assert "Output path exists and is not a directory" in " ".join(result.output.split())
    assert occupied.read_text() == "not a directory\n"  # left untouched


@pytest.mark.parametrize("command", ["detect", "build-panel", "type-lineages", "reassort"])
def test_a_new_nested_output_directory_is_still_created(
    tmp_path: Path, no_pipeline, command: str
) -> None:
    """Review focus 5, guard (passes before and after): the check must not stop a command
    from creating an output directory that does not exist yet (the run log goes there)."""
    out = tmp_path / "new" / "nested" / "out"
    runner.invoke(app, _args(command, tmp_path, out))
    assert out.is_dir()
    assert (out / "tessera.log").exists()
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_cli_output_path.py -q`
Expected: 7 failed, 4 passed (the guard). Observed on `4b35e39`: `detect`, `build-panel`, `type-lineages`, `reassort` fail with `AssertionError: raw exception: FileExistsError(17, 'File exists')`; `fill-references`, `find-references`, `curate-panel` fail with `assert 'Unexpected error' not in '[ERROR] Une...traceback.'` (they reached the pipeline stub).

- [ ] **Step 3: Implement**

`src/tessera/cli/main.py`, replace `get_logger`:

```python
def get_logger(output_dir: Path | None = None) -> logging.Logger:
    """Return the configured logger; with ``output_dir`` also persist a run log there.

    The run log is opened before a command enters :func:`stage_errors`, so an output
    path that is an existing file is checked here. It is rejected with the message every
    other command gives (:func:`_require_output_directory`) instead of a raw
    ``FileExistsError`` traceback.
    """
    level: int = _RUN_STATE["log_level"]  # type: ignore[assignment]
    if output_dir is not None:
        try:
            _require_output_directory(output_dir)
        except UserInputError as exc:
            configure_logging(level=level).error("%s", exc)
            raise typer.Exit(code=1) from exc
    return configure_logging(output_dir=output_dir, level=level)
```

`src/tessera/cli/cmd_fill_references.py`: add `_require_output_directory,` to the `from .main import (...)` list and, inside `stage_errors`, after `_require_lineage_map(lineage_map)`:

```python
        _require_output_directory(output)
```

`src/tessera/cli/cmd_find_references.py`: the import becomes

```python
from .main import (
    _require_file,
    _require_output_directory,
    _require_scan_windows,
    app,
    get_logger,
    stage_errors,
)
```

and after `_require_file(msa, "MSA file")` add `_require_output_directory(output)`.

`src/tessera/cli/cmd_curate_panel.py`: the import becomes

```python
from .main import (
    _require_directory,
    _require_file,
    _require_output_directory,
    app,
    get_logger,
    stage_errors,
)
```

and after `_require_directory(collection, "Collection directory")` add `_require_output_directory(output)`.

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_cli_output_path.py tests/unit/test_cli_input_checks.py tests/unit/test_cli_commands.py -q`

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`

- [ ] **Step 6: Commit**

```bash
git add src/tessera/cli/main.py src/tessera/cli/cmd_fill_references.py \
        src/tessera/cli/cmd_find_references.py src/tessera/cli/cmd_curate_panel.py \
        tests/unit/test_cli_output_path.py
git commit -F - <<'EOF'
Reject an output path that is a file on every command, not only recomb

detect, build-panel, type-lineages and reassort open their run log in the output
directory before entering the error handler, so an existing file at --output gave a raw
FileExistsError traceback. fill-references, find-references and curate-panel never
checked the path and failed later under "Unexpected error". Only recomb rejected it
cleanly (plan B, item B11).

get_logger now applies the same _require_output_directory check before it creates the
directory, and the three commands without a run log call it first inside the error
handler, so every command gives one message and exit code 1. A new or nested output
directory is still created as before.

Audit 4, item E1.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 11: Range-check the discover commands' coverage fractions and search counts  (audit items: E5 -- coverage fractions on `recomb` only; six unbounded count options; `--top-gaps 0`)

`recomb` range-checks `--coverage-floor` and `--coverage-rel-drop` (`cmd_recomb.py:173-175`); `fill-references` and `find-references` take the same two options unchecked, so `--coverage-floor 50` (a percentage) makes every window a gap. Six count options are unbounded: `--max-hits`, `--top-gaps`, `--subtile` (fill and find), `--seed-hits`, `--seed-window` (fill) and `--threads` (fill, detect, build-panel). `--top-gaps 0` searched nothing and exited cleanly, which reads as "no reference is missing". `--fetch-limit` is plan E's (D8) and is left alone.

**Files:**
- Modify: `src/tessera/cli/main.py:117` (two helpers before `_parse_key_values`)
- Modify: `src/tessera/cli/cmd_fill_references.py:9-19`, `:213-214`
- Modify: `src/tessera/cli/cmd_find_references.py` (import list), `:91`
- Modify: `src/tessera/cli/cmd_detect.py:130`, `src/tessera/cli/cmd_build_panel.py:101`
- Modify: `tests/unit/test_cli_input_checks.py` (append)

**Interfaces:**
- Produces: `main._require_coverage_fractions(coverage_floor: float | None, coverage_rel_drop: float) -> None`; `main._require_search_counts(*, max_hits: int, top_gaps: int, subtile: int, threads: int | None = None) -> None`.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_cli_input_checks.py`:

```python


# --- discover commands: coverage fractions and search counts (audit 4, E5) ---

@pytest.mark.parametrize(
    ("option", "value"),
    [("--coverage-floor", "50"), ("--coverage-rel-drop", "5"), ("--top-gaps", "0"),
     ("--max-hits", "0"), ("--subtile", "-1")],
)
@pytest.mark.parametrize("command", ["fill-references", "find-references"])
def test_discover_commands_reject_out_of_range_search_options(
    tmp_path: Path, no_pipeline, command: str, option: str, value: str
) -> None:
    """recomb range-checked the coverage fractions; here --coverage-floor 50 made every
    window a gap, and --top-gaps 0 searched nothing and exited cleanly."""
    head = (["--msa", EXAMPLE, "--query", "query"] if command == "find-references"
            else ["--query", str(_query(tmp_path))])
    result = runner.invoke(app, [command, *head, "--output", str(tmp_path / "out"),
                                 option, value])
    _rejected(result, f"Invalid {option}")


@pytest.mark.parametrize(("option", "value"), [("--seed-hits", "0"), ("--seed-window", "0")])
def test_fill_references_rejects_empty_seed_searches(
    tmp_path: Path, no_pipeline, option: str, value: str
) -> None:
    result = runner.invoke(app, [
        "fill-references", "--query", str(_query(tmp_path)),
        "--output", str(tmp_path / "out"), option, value,
    ])
    _rejected(result, f"Invalid {option}")


@pytest.mark.parametrize("command", ["detect", "fill-references", "build-panel"])
def test_panel_commands_reject_zero_threads(tmp_path: Path, no_pipeline, command: str) -> None:
    result = runner.invoke(app, [
        command, "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        "--threads", "0",
    ])
    _rejected(result, "Invalid --threads")
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_cli_input_checks.py -q -k "search_options or empty_seed or zero_threads"`
Expected: 15 failed. Observed on `4b35e39`: every case fails `assert 'Unexpected error' not in ...` -- the command passed validation and reached the pipeline stub.

- [ ] **Step 3: Implement**

`src/tessera/cli/main.py`, before `_parse_key_values`:

```python
def _require_coverage_fractions(coverage_floor: float | None, coverage_rel_drop: float) -> None:
    """The two coverage options are similarities on a 0-1 scale on every command.

    ``recomb`` checked them; on the discover commands ``--coverage-floor 50`` (a
    percentage) made every window a gap and the run searched the whole genome.
    """
    _require_range(coverage_rel_drop, "--coverage-rel-drop", lo=0.0, hi=1.0)
    if coverage_floor is not None:
        _require_range(coverage_floor, "--coverage-floor", lo=0.0, hi=1.0)


def _require_search_counts(
    *, max_hits: int, top_gaps: int, subtile: int, threads: int | None = None,
) -> None:
    """The per-round search sizes the discover commands share.

    ``--top-gaps 0`` searched nothing and exited cleanly, which reads as "no reference
    is missing"; a zero hit count does the same per gap.
    """
    _require_range(max_hits, "--max-hits", lo=1)
    _require_range(top_gaps, "--top-gaps", lo=1)
    _require_range(subtile, "--subtile", lo=0)
    if threads is not None:
        _require_range(threads, "--threads", lo=1)
```

`src/tessera/cli/cmd_fill_references.py`: add `_require_coverage_fractions,` and `_require_search_counts,` to the `from .main import (...)` list (keep it sorted for ruff `I`), and after `_require_scan_windows(window_size, window_step)`:

```python
        _require_coverage_fractions(coverage_floor, coverage_rel_drop)
        _require_search_counts(
            max_hits=max_hits, top_gaps=top_gaps, subtile=subtile, threads=threads,
        )
        _require_range(seed_hits, "--seed-hits", lo=1)
        _require_range(seed_window, "--seed-window", lo=1)
```

`src/tessera/cli/cmd_find_references.py`: add `_require_coverage_fractions,` and `_require_search_counts,` to the import list (the list from Task 10), and after `_require_scan_windows(window_size, window_step)`:

```python
        _require_coverage_fractions(coverage_floor, coverage_rel_drop)
        _require_search_counts(max_hits=max_hits, top_gaps=top_gaps, subtile=subtile)
```

`src/tessera/cli/cmd_detect.py` and `src/tessera/cli/cmd_build_panel.py`: after `_require_scan_windows(window_size, window_step)` add

```python
        _require_range(threads, "--threads", lo=1)
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_cli_input_checks.py -q`

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`

- [ ] **Step 6: Commit**

```bash
git add src/tessera/cli/main.py src/tessera/cli/cmd_fill_references.py \
        src/tessera/cli/cmd_find_references.py src/tessera/cli/cmd_detect.py \
        src/tessera/cli/cmd_build_panel.py tests/unit/test_cli_input_checks.py
git commit -F - <<'EOF'
Range-check the discover commands' coverage fractions and search counts

recomb rejects a coverage fraction outside 0-1; fill-references and find-references
took the same options unchecked, so --coverage-floor 50 made every window a gap. Six
count options were unbounded, and --top-gaps 0 searched nothing and exited cleanly,
which reads as "no reference is missing" rather than as a search that never ran.

Two small helpers in cli/main.py hold the checks the discover commands share; --threads
is checked on detect, build-panel and fill-references too. Input checks only.

Audit 4, item E5 (coverage fractions, count options, --top-gaps 0).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 12: Option help and interactions on `fill-references`, `find-references`, `msa` and `build-panel`  (audit items: E5 -- `--seed-source` help, `find-references -q`, `msa --reference` with `--query-as-backbone`, `msa` overwrite, `build-panel` options)

Each is a small CLI defect a reviewer would approve or reject together:

- `fill-references --seed-source` accepts `nextclade` (`cmd_fill_references.py:219`) and `docs/reference-panels.md` describes it, but the help (`:38-42`) lists only `blast`, `local` and `ncbi-virus`.
- `find-references -q` is an alignment label declared as `Path` (`cmd_find_references.py:17-19`), so `sample/` reaches the scan as `sample`.
- `msa --reference` silently overrides `--query-as-backbone` (`core/io.py:346-356`, explicit reference first), so the user gets a different coordinate system without a word. Both now are an error.
- `msa` overwrites an existing alignment without notice; it now logs a warning. `msa --threads 0` is also rejected.
- `build-panel` cannot set `--pool-consensus`, `--keep-recombinant-lineages` or `--deep-typing`, three options that shape the panel it exists to build; `detect` and `fill-references` expose them and `FillParams` takes them. They are added with the same help text. (Whether `--pool-consensus` does anything outside a Nextclade pool is plan E's item; the help here says it applies with `--nextclade`.)

**Files:**
- Modify: `src/tessera/cli/cmd_fill_references.py:38-42`
- Modify: `src/tessera/cli/cmd_find_references.py:17-19`, `:92-93`
- Modify: `src/tessera/cli/cmd_msa.py:9`, `:44-46`
- Modify: `src/tessera/cli/cmd_build_panel.py:80-81`, `:103-116`
- Modify: `docs/reference-panels.md:22-29`, `README.md:111-112`
- Modify: `tests/unit/test_cli_input_checks.py` (append)

**Interfaces:** none beyond the CLI options (`build-panel --pool-consensus/--no-pool-consensus`, `--keep-recombinant-lineages/--no-keep-recombinant-lineages`, `--deep-typing/--no-deep-typing`, all default off).

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_cli_input_checks.py`:

```python


def test_find_references_keeps_the_query_label_verbatim(
    tmp_path: Path, monkeypatch
) -> None:
    """-q is an alignment label, but was declared a path: 'sample/' became 'sample'."""
    seen = {}

    def capture(params, logger):
        seen["query"] = params.query

    monkeypatch.setattr("tessera.discover.run.find_references", capture)
    result = runner.invoke(app, [
        "find-references", "--msa", EXAMPLE, "--query", "sample/",
        "--output", str(tmp_path / "out"),
    ])
    assert result.exit_code == 0, result.output
    assert seen["query"] == "sample/"


def test_build_panel_passes_the_panel_shaping_options(tmp_path: Path, monkeypatch) -> None:
    """build-panel could not set three options that shape the panel it builds."""
    seen = {}

    def capture(params, logger):
        seen.update(vars(params))

    monkeypatch.setattr("tessera.discover.iterate.fill_references", capture)
    result = runner.invoke(app, [
        "build-panel", "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        # --nextclade: plan E Task 11 refuses --pool-consensus under any other seed
        # source, so this is the one combination valid before and after plan E.
        "--nextclade",
        "--pool-consensus", "--keep-recombinant-lineages", "--deep-typing",
    ])
    assert result.exit_code == 0, result.output
    assert seen["pool_consensus"] and seen["keep_recombinant"] and seen["deep_typing"]


def test_seed_source_help_names_nextclade() -> None:
    import typer

    cmd = typer.main.get_command(app).commands["fill-references"]
    param = next(p for p in cmd.params if "--seed-source" in getattr(p, "opts", []))
    assert "'nextclade'" in (param.help or "")


# --- msa: backbone choice and overwrite (audit 4, E5) ------------------------

def test_msa_rejects_reference_with_query_as_backbone(tmp_path: Path, monkeypatch) -> None:
    """--reference silently won, so the query was not the backbone the user asked for."""
    monkeypatch.setattr("tessera.msa.build.build_msa", _must_not_run)
    result = runner.invoke(app, [
        "msa", "--query", str(_query(tmp_path)), "--collection", str(_collection(tmp_path)),
        "--output", str(tmp_path / "aln.fasta"), "--aligner", "mafft",
        "--reference", "ref", "--query-as-backbone",
    ])
    _rejected(result, "give one")


def test_msa_says_when_it_overwrites_an_alignment(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr("tessera.msa.build.build_msa", lambda params, logger: None)
    existing = tmp_path / "aln.fasta"
    existing.write_text(">old\nACGT\n")
    result = runner.invoke(app, [
        "msa", "--query", str(_query(tmp_path)), "--collection", str(_collection(tmp_path)),
        "--output", str(existing), "--aligner", "mafft",
    ])
    assert result.exit_code == 0, result.output
    assert "Overwriting the existing alignment" in result.output
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_cli_input_checks.py -q -k "verbatim or shaping or seed_source_help or msa_"`
Expected: 5 failed. Observed on `4b35e39`: `assert 'sample' == 'sample/'`; `build-panel` exits 2 with `Usage: root build-panel [OPTIONS]` (no such option); `assert "'nextclade'" in ...`; `assert 'Unexpected error' not in ...` (the `msa` run reached the `build_msa` stub); `assert 'Overwriting the existing alignment' in ''`.

- [ ] **Step 3: Implement**

`src/tessera/cli/cmd_fill_references.py`, the `--seed-source` help:

```python
        help="Fresh-start source (no --collection): 'blast' (NCBI nt), 'local' (a genome pool "
        "directory, --candidate-pool), 'ncbi-virus' (a taxon-scoped set, --taxon), or "
        "'nextclade' (a Nextclade dataset's reference tree, --nextclade-dataset or "
        "auto-detected from the query).",
```

`src/tessera/cli/cmd_find_references.py`: `query: str = typer.Option(..., "-q", "--query", help="Query label as it appears in the MSA."),` and in `FindRefParams(...)` pass `query=query` (was `query=str(query)`).

`src/tessera/cli/cmd_msa.py`: the import block becomes

```python
from ..core.errors import UserInputError
from .main import (
    DEFAULT_THREADS,
    _parse_key_values,
    _require_choice,
    _require_range,
    app,
    get_logger,
    stage_errors,
)
```

and after `_require_choice(aligner, set(aligner_registry.names()), "--aligner")`:

```python
        _require_range(threads, "--threads", lo=1)
        # Both name the backbone. --reference used to win silently, so a user who asked
        # for the query as backbone got a different coordinate system without a word.
        if reference is not None and query_as_backbone:
            raise UserInputError(
                "--reference and --query-as-backbone both choose the backbone; give one."
            )
        if output.exists():
            logger.warning("Overwriting the existing alignment %s.", output)
```

`src/tessera/cli/cmd_build_panel.py`: after the `threads` option add

```python
    pool_consensus: bool = typer.Option(
        False, "--pool-consensus/--no-pool-consensus",
        help="With a Nextclade pool (--nextclade), use one denoised consensus genome per "
        "clade (a stable per-lineage reference) instead of every tree tip.",
    ),
    keep_recombinant_lineages: bool = typer.Option(
        False, "--keep-recombinant-lineages/--no-keep-recombinant-lineages",
        help="Keep recombinant (CRF/URF/X) lineages in the panel. Off by default: "
        "they carry both parents' segments and mask the true parents. Turn on to "
        "recruit the best-matching genome even if it is a known hybrid. "
        "Applies only when the references carry lineage labels.",
    ),
    deep_typing: bool = typer.Option(
        False, "--deep-typing/--no-deep-typing",
        help="Type the recruited panel with the full lineage ladder (Nextclade "
        "nearest-neighbour + de-novo ANI clustering) instead of header mining alone. "
        "Needs skani; off by default.",
    ),
```

and in `FillParams(...)` after `report=False, lineage_map=lineage_map,`:

```python
            pool_consensus=pool_consensus,
            keep_recombinant=keep_recombinant_lineages,
            deep_typing=deep_typing,
```

`docs/reference-panels.md`, after "Pass `--collection` to grow an existing collection instead of recruiting from scratch." add "`--keep-recombinant-lineages`, `--deep-typing` and `--pool-consensus` shape the panel as they do for `detect`." `README.md:112`, after "`--query-as-backbone` uses a single-contig query as the backbone" add "; it cannot be combined with `--reference`".

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_cli_input_checks.py tests/unit/test_cli_commands.py tests/unit/test_cli_deep_typing.py tests/unit/test_cli_keep_recombinant.py -q`

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`

- [ ] **Step 6: Commit**

```bash
git add src/tessera/cli/cmd_fill_references.py src/tessera/cli/cmd_find_references.py \
        src/tessera/cli/cmd_msa.py src/tessera/cli/cmd_build_panel.py \
        docs/reference-panels.md README.md tests/unit/test_cli_input_checks.py
git commit -F - <<'EOF'
Fix option help and interactions on four commands

--seed-source accepted "nextclade" but its help did not list it. find-references -q is
an alignment label but was declared a path, so "sample/" reached the scan as "sample".
msa --reference silently overrode --query-as-backbone, giving the user a different
backbone and coordinate system than asked for; giving both is now an error. msa
overwrote an existing alignment without notice; it now says so. build-panel could not
set --pool-consensus, --keep-recombinant-lineages or --deep-typing, three options that
shape the panel it exists to build; it now takes them with the same defaults and help
as detect.

Audit 4, item E5 (help and option interactions).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 13: Bound `reassort --margin` and catch an ANI floor given as a fraction  (audit items: E2; E5 -- `--ani-floor 0.8`, the two `--margin` units)

`cmd_reassort.py:62` checks only the lower bound of `--margin`. The near-best set (`constellation.py:31-36`) grows with the margin; once it spans the ANI spread of the candidates every pair shares a strain and every query reads clonal (the audit measured 0.5 giving a reassortant verdict and 40, 100 and 1e6 all giving clonal on one dataset). A margin of zero does the reverse: a 0.01-point tie-break between two tips reads as different parents.

**The bound and why.** `docs/reference-panels.md:341-342` describes the margin as the ANI window within which a strain counts as near-best, default 0.5 points: it exists to absorb near-ties between tips and ANI estimation noise, a fraction of a point. The accepted range is now above 0 and at most 5 points (`MAX_MARGIN = 5.0`), ten times the default; zero is rejected rather than documented, because it turns a tie-break into a call. The cap is a fixed sanity bound, not a calibrated one: it rejects the values the audit used to invert the verdict and leaves the useful range open. The help says what a wide window does.

`--ani-floor` is a percentage (`DEFAULT_ANI_FLOOR = 80.0`); `0.8`, meant as 80 %, is accepted and lets every tip through. A value in (0, 1] is now rejected with the likely intended value. `recomb --margin` is a 0-1 similarity fraction while `reassort --margin` is in ANI percentage points; both help strings now state their unit (recomb's in Task 5).

**Files:**
- Modify: `src/tessera/reassort/constellation.py:15` (add `MAX_MARGIN`)
- Modify: `src/tessera/cli/cmd_reassort.py:20`, `:39-46`, `:59-62`
- Modify: `docs/reference-panels.md:340-342`
- Modify: `tests/unit/test_cli_input_checks.py` (append)

**Interfaces:**
- Produces: `constellation.MAX_MARGIN: float = 5.0`.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_cli_input_checks.py`:

```python


# --- reassort: margin and ANI units (audit 4, E2 and E5) ---------------------

@pytest.mark.parametrize(
    ("option", "value", "needle"),
    [("--margin", "0", "Invalid --margin"),
     ("--margin", "40", "Invalid --margin"),
     ("--ani-floor", "0.8", "Did you mean --ani-floor 80?")],
)
def test_reassort_rejects_a_margin_or_floor_on_the_wrong_scale(
    tmp_path: Path, no_reassort, option: str, value: str, needle: str
) -> None:
    """A margin of 40 made every query clonal and 0 made a 0.01-point tie-break a
    reassortment; an ANI floor of 0.8 (meant as 80 %) accepted every tip."""
    result = runner.invoke(app, [
        "reassort", "--query", str(_segments(tmp_path)), "--output", str(tmp_path / "out"),
        option, value,
    ])
    _rejected(result, needle)


def test_reassort_accepts_the_margin_bound(tmp_path: Path, monkeypatch) -> None:
    """Guard: 5 points is inside the range (the bound is inclusive)."""
    from tessera.reassort.assign import ReassortmentResult

    monkeypatch.setattr("tessera.cli.cmd_reassort.assign_segments",
                        lambda query, **kwargs: ReassortmentResult())
    result = runner.invoke(app, [
        "reassort", "--query", str(_segments(tmp_path)), "--output", str(tmp_path / "out"),
        "--margin", "5",
    ])
    assert result.exit_code == 0, result.output
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_cli_input_checks.py -q -k "wrong_scale or margin_bound"`
Expected: 3 failed, 1 passed (the guard). Observed on `4b35e39`: each of the three fails `assert 'Unexpected error' not in ...` -- the value passed validation and reached the `assign_segments` stub.

- [ ] **Step 3: Implement**

`src/tessera/reassort/constellation.py`, after `DEFAULT_MARGIN`:

```python
# The widest margin accepted. The window absorbs near-ties between tips and ANI estimation
# noise, a fraction of a point; once it spans the ANI spread of the candidates every pair
# shares a strain and every query reads clonal (measured: 0.5 gives reassortant on a case
# that 40 calls clonal).
MAX_MARGIN = 5.0
```

`src/tessera/cli/cmd_reassort.py`: `from ..reassort.constellation import DEFAULT_MARGIN, MAX_MARGIN`; the two options:

```python
    ani_floor: float = typer.Option(
        DEFAULT_ANI_FLOOR, "--ani-floor",
        help="A segment below this ANI (a percentage, 0-100) to every reference tip is "
        "left unassigned.",
    ),
    margin: float = typer.Option(
        DEFAULT_MARGIN, "--margin",
        help="ANI window in percentage points (above 0, at most 5) for a segment's "
        "near-best parents. Not the 0-1 similarity fraction 'recomb --margin' takes. A "
        "wide window lets every pair share a strain, so every query reads clonal.",
    ),
```

and replace the checks at lines 59-62 with

```python
        # ANI is a percentage. Above 100 nothing can be assigned; a value in (0, 1] is
        # almost certainly a fraction meant as a percentage (0.8 for 80), and would
        # accept every tip skani reports.
        _require_range(ani_floor, "--ani-floor", lo=0.0, hi=100.0)
        if 0.0 < ani_floor <= 1.0:
            raise UserInputError(
                f"Invalid --ani-floor {ani_floor:g}: ANI is a percentage (0-100), not a "
                f"fraction. Did you mean --ani-floor {100 * ani_floor:g}?"
            )
        # A margin of zero makes a 0.01-point tie-break between two tips read as
        # different parents; a wide one makes every pair share a strain, so the verdict
        # is clonal whatever the data (0.5 gave reassortant where 40 gave clonal).
        _require_range(margin, "--margin", lo=0.0, hi=MAX_MARGIN, lo_open=True)
```

The existing `test_reassort_rejects_out_of_range_options` (`--margin -3`, `--ani-floor 500/-1`) still passes: the messages keep the `Invalid --margin` / `Invalid --ani-floor` prefix.

`docs/reference-panels.md`, the concordance sentence ("within `--margin` ANI points, default 0.5") becomes "within `--margin` ANI percentage points, default 0.5; the accepted range is above 0 and at most 5, because a window wide enough to span the candidates' ANI spread makes every pair concordant and every query clonal, and a zero window turns a tie-break between two tips into different parents. `--ani-floor` is a percentage too; a value of 1 or below is rejected as a likely fraction."

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_cli_input_checks.py tests/unit/test_cli_reassort.py -q`

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`

- [ ] **Step 6: Commit**

```bash
git add src/tessera/reassort/constellation.py src/tessera/cli/cmd_reassort.py \
        docs/reference-panels.md tests/unit/test_cli_input_checks.py
git commit -F - <<'EOF'
Bound reassort --margin and catch an ANI floor given as a fraction

reassort --margin had no upper bound. The near-best set grows with the margin, and
once it spans the candidates' ANI spread every pair shares a strain and every query
reads clonal: on one dataset 0.5 gave a reassortant verdict and 40, 100 and 1e6 all gave
clonal. A margin of zero did the reverse, making a 0.01-point tie-break between two
tips read as different parents.

The accepted range is now above 0 and at most 5 ANI points, ten times the default. The
window exists to absorb near-ties and ANI noise, a fraction of a point; the cap is a
sanity bound that rejects the inverting values, not a calibrated optimum. --ani-floor
0.8, meant as 80 %, accepted every tip; a value in (0, 1] is now rejected with the
likely intended value. The help states that reassort --margin is in ANI points, unlike
recomb --margin.

Audit 4, items E2 and E5 (ANI units).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 14: A discordant pair inside one parent group, empty candidate lists, and the pairwise evidence  (audit items: E3, E14, E4)

**E3.** `constellation.py:62-68` builds groups from concordant edges only and lets any discordant pair set the verdict. A middle segment concordant with both ends links all three into one group while the ends disagree, so the user sees a reassortant verdict beside a single group with no parent strain, and every row names group 0. The documented explanation of an empty parent group (`docs/reference-panels.md:347-351`) covers only the variant that stays clonal.

**The rule chosen:** a reassortant verdict needs the discordant segments in different parent groups. A discordant pair inside one group gives `undetermined`, with a note naming the pair and saying the parents cannot be separated at this margin. Why not split the group: concordance is not transitive at a fixed margin, so the middle segment is equally concordant with both ends and any split would be arbitrary; and why not keep `reassortant`: the verdict would have no parent groups to report, which is the defect. `undetermined` with the reason is the faithful statement. A discordant pair across two groups still gives `reassortant` (the seven-plus-one test is unchanged).

**E14.** `call_constellation` is public (`reassort/__init__.py`). A segment with an empty candidate list has an empty near-best set; paired with a segment whose near-best strains appear in its universe it was marked discordant, so an empty list produced a reassortant verdict. Unreachable through the CLI (only assigned segments are passed), but unguarded; such a pair is now uninformative.

**E4.** `cmd_reassort.py:110-111` logs the pairwise notes only for an undetermined verdict, so a reassortant call records nowhere which pair disagreed or on which strains. The constellation now returns structured pair records; the command writes them to `segment_pairs.tsv` and logs the notes for every verdict, and a discordant note names the two near-best strain sets.

**Files:**
- Modify: `src/tessera/reassort/constellation.py:1-8` (docstring), after `MAX_MARGIN` (new `PairCall`), `:24-28` (`Constellation`), `:45-69` (`call_constellation` body)
- Modify: `src/tessera/reassort/assign.py:28`, `:47-54`, `:187-190`
- Modify: `src/tessera/reassort/__init__.py`
- Modify: `src/tessera/cli/cmd_reassort.py:108-112`
- Modify: `docs/reference-panels.md:343-351`
- Modify: `tests/unit/test_constellation.py`, `tests/unit/test_cli_input_checks.py` (append)

**Interfaces:**
- Produces: `constellation.PairCall` dataclass (`segment_a: str`, `segment_b: str`, `relation: str` one of `"concordant" | "discordant" | "discordant_within_group" | "uninformative"`, `note: str`); `Constellation.pairs: list[PairCall]`; `ReassortmentResult.pairs: list[PairCall]`; output file `segment_pairs.tsv` (`segment_a`, `segment_b`, `relation`, `note`). `Constellation.pair_notes` is kept (`[p.note for p in pairs]`).

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_constellation.py`:

```python


# --- audit 4: E3, E14 ---------------------------------------------------------

def test_discordant_pair_inside_one_concordant_group_is_undetermined():
    """HA shares Y with NA and X with PB1, so all three form one group, while NA and PB1
    disagree. The verdict was reassortant beside a single group that names no parent."""
    candidates = {
        "HA": [("X", 99.9), ("Y", 99.8)],
        "NA": [("Y", 99.9), ("Z", 99.8)],
        "PB1": [("X", 99.9)],
    }
    universe = {"X", "Y", "Z"}
    c = call_constellation(candidates, {k: set(universe) for k in candidates}, margin=0.5)
    assert len(c.groups) == 1 and c.groups[0].parent_strains == []
    assert c.verdict == "undetermined"
    within = [p for p in c.pairs if p.relation == "discordant_within_group"]
    assert [(p.segment_a, p.segment_b) for p in within] == [("NA", "PB1")]
    assert "cannot be separated" in within[0].note


def test_discordant_pair_across_groups_is_still_reassortant():
    """Guard (passes before and after): the seven-plus-one case keeps its verdict."""
    seven = [f"seg{i}" for i in range(7)]
    candidates = {s: [("A", 99.0)] for s in seven}
    candidates["odd"] = [("B", 99.0)]
    universes = {s: {"A", "B"} for s in [*seven, "odd"]}
    assert call_constellation(candidates, universes, margin=0.5).verdict == "reassortant"


def test_a_segment_with_no_candidates_is_uninformative_not_discordant():
    """The public function made a reassortant call from an empty candidate list."""
    universe = {"A", "B"}
    c = call_constellation({"HA": [], "NA": [("A", 99.0)]},
                           {"HA": set(universe), "NA": set(universe)}, margin=0.5)
    assert c.verdict == "undetermined"
    assert [p.relation for p in c.pairs] == ["uninformative"]
```

Append to `tests/unit/test_cli_input_checks.py`:

```python


def test_reassort_writes_the_pairwise_evidence_for_a_reassortant_call(
    tmp_path: Path, monkeypatch
) -> None:
    """The notes were logged only for an undetermined verdict, so a reassortant call
    recorded nowhere which pair disagreed or on which strains."""
    from tessera.reassort.assign import ReassortmentResult, SegmentAssignment
    from tessera.reassort.constellation import call_constellation

    call = call_constellation({"HA": [("S1", 99.0)], "NA": [("S3", 99.0)]},
                              {"HA": {"S1", "S3"}, "NA": {"S1", "S3"}}, margin=0.5)

    def fake_assign(query, **kwargs):
        return ReassortmentResult(
            segments=[SegmentAssignment("HA", "d", "S1", "c1", 99.0, "assigned"),
                      SegmentAssignment("NA", "d", "S3", "c2", 99.0, "assigned")],
            verdict=call.verdict, groups=call.groups, pair_notes=call.pair_notes,
            pairs=call.pairs,
        )

    monkeypatch.setattr("tessera.cli.cmd_reassort.assign_segments", fake_assign)
    out = tmp_path / "out"
    result = runner.invoke(app, [
        "reassort", "--query", str(_segments(tmp_path)), "--output", str(out),
    ])
    assert result.exit_code == 0, result.output
    assert call.verdict == "reassortant"
    lines = (out / "segment_pairs.tsv").read_text().splitlines()
    assert lines[0] == "segment_a\tsegment_b\trelation\tnote"
    assert lines[1].split("\t")[:3] == ["HA", "NA", "discordant"]
    assert "S1 vs S3" in lines[1]
    assert "Pairwise evidence" in result.output
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_constellation.py tests/unit/test_cli_input_checks.py -q -k "within_group or no_candidates or across_groups or pairwise_evidence"`
Expected: 3 failed, 1 passed (the across-groups guard). Observed on `4b35e39`: `assert 'reassortant' == 'undetermined'` for both constellation tests, and `AssertionError: [ERROR] Unexpected error: 'Constellation' object has no attribute 'pairs'` for the CLI test.

- [ ] **Step 3: Implement**

`src/tessera/reassort/constellation.py`: append to the module docstring

```python
A reassortant verdict needs the discordant segments in different parent groups. Groups are
the connected components of the concordant pairs, and concordance is not transitive: a
middle segment can share a strain with each of two segments that disagree with each other.
Such a discordant pair inside one group gives a verdict of ``undetermined``, with a note,
rather than a reassortant verdict beside a single group that names no parent.
```

after `MAX_MARGIN` add

```python


@dataclass
class PairCall:
    """The relation of one segment pair, with the strains it rests on, in words."""

    segment_a: str
    segment_b: str
    relation: str  # "concordant" | "discordant" | "discordant_within_group" | "uninformative"
    note: str
```

add `pairs: list[PairCall] = field(default_factory=list)` as the last field of `Constellation`, and replace the body of `call_constellation` (from `segments = list(candidates)` to the `return`) with

```python
    segments = list(candidates)
    near = {s: _near_best(candidates[s], margin) for s in segments}

    concordant: list[tuple[str, str]] = []
    discordant: list[tuple[str, str]] = []
    pairs: list[PairCall] = []
    for i, j in combinations(segments, 2):
        shared = near[i] & near[j]
        if not near[i] or not near[j]:
            # A segment with no candidate strain has no parent to compare. Unreachable
            # through the CLI, which passes assigned segments only, but the function is
            # public: an empty list must not read as a different parent.
            pairs.append(PairCall(i, j, "uninformative",
                                  f"{i}?{j} uninformative (no candidate strain for one segment)"))
        elif shared:
            concordant.append((i, j))
            pairs.append(PairCall(i, j, "concordant",
                                  f"{i}~{j} concordant on {', '.join(sorted(shared))}"))
        elif (near[i] & universes[j]) or (near[j] & universes[i]):
            discordant.append((i, j))
            pairs.append(PairCall(
                i, j, "discordant",
                f"{i}!{j} discordant (cross-typed, different parents: "
                f"{', '.join(sorted(near[i]))} vs {', '.join(sorted(near[j]))})",
            ))
        else:
            pairs.append(PairCall(i, j, "uninformative",
                                  f"{i}?{j} uninformative (no shared strain vocabulary)"))

    groups = _components(segments, concordant, near)
    group_of = {s: k for k, g in enumerate(groups) for s in g.segments}
    within = {(i, j) for i, j in discordant if group_of[i] == group_of[j]}
    for pair in pairs:
        if (pair.segment_a, pair.segment_b) in within:
            pair.relation = "discordant_within_group"
            pair.note += (
                "; but both are linked into one group through segments concordant with "
                "each, so the parents cannot be separated at this margin (a smaller "
                "--margin may resolve it)"
            )
    if within:
        verdict = "undetermined"
    elif discordant:
        verdict = "reassortant"
    elif concordant and len(groups) == 1:
        verdict = "clonal"
    else:
        verdict = "undetermined"
    return Constellation(
        groups=groups, verdict=verdict, pair_notes=[p.note for p in pairs], pairs=pairs,
    )
```

`src/tessera/reassort/assign.py`: `from .constellation import DEFAULT_MARGIN, PairCall, ParentGroup, call_constellation`; add `pairs: list[PairCall] = field(default_factory=list)` to `ReassortmentResult` after `pair_notes`; after `result.pair_notes = call.pair_notes` add `result.pairs = call.pairs`.

`src/tessera/reassort/__init__.py`: export `PairCall` (`from .constellation import Constellation, PairCall, ParentGroup, call_constellation` and `"PairCall",` in `__all__`).

`src/tessera/cli/cmd_reassort.py`, replace lines 108-112 with

```python
        # The pairwise evidence behind the verdict, for every verdict: which pair
        # disagreed, and on which strains, is what a reader needs to check a call.
        ptsv = output / "segment_pairs.tsv"
        with open(ptsv, "w") as fo:
            fo.write("segment_a\tsegment_b\trelation\tnote\n")
            for pair in result.pairs:
                fo.write(f"{pair.segment_a}\t{pair.segment_b}\t{pair.relation}\t"
                         f"{' '.join(pair.note.split())}\n")

        logger.info("Reassortment verdict: %s", result.verdict.upper())
        logger.info("Parent constellation: %s", mosaic or "(no segments assigned)")
        if result.pair_notes:
            logger.info("Pairwise evidence: %s", "; ".join(result.pair_notes))
        logger.info("Wrote %s, %s and %s", tsv, ctsv, ptsv)
```

`docs/reference-panels.md`: after "A query is `reassortant` when any pair is discordant" insert "and the discordant segments fall in different parent groups"; after the paragraph on transitively linked groups add

```markdown
The reverse case -- a middle segment concordant with two segments that are discordant
with each other -- links all three into one group although two of them disagree. The
parents then cannot be separated at this margin, and the verdict is `undetermined` with a
note naming the pair, not `reassortant` beside a group with no parent strain; a smaller
`--margin` may resolve it. Every pair's relation and the strains it rests on are written
to `out/segment_pairs.tsv` (`segment_a`, `segment_b`, `relation`, `note`) and to the log,
whatever the verdict.
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_constellation.py tests/unit/test_cli_input_checks.py tests/unit/test_reassort_assign.py tests/unit/test_reassort_benchmark.py -q`

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`

- [ ] **Step 6: Commit**

```bash
git add src/tessera/reassort/constellation.py src/tessera/reassort/assign.py \
        src/tessera/reassort/__init__.py src/tessera/cli/cmd_reassort.py \
        docs/reference-panels.md tests/unit/test_constellation.py \
        tests/unit/test_cli_input_checks.py
git commit -F - <<'EOF'
Call reassortment only across parent groups, and write the evidence

Groups are the connected components of the concordant segment pairs, and concordance
is not transitive at a fixed margin. A middle segment concordant with two segments that
disagree with each other links all three into one group, and the old rule still called
the query reassortant: the user saw a reassortant verdict beside one group with no
parent strain. Such a pair now gives "undetermined" with a note naming it. Splitting the
group instead would be arbitrary, since the middle segment is equally concordant with
both ends. A discordant pair across two groups is still reassortant.

The public call_constellation also turned a segment with an empty candidate list into a
discordant pair, and so into a reassortant verdict; that pair is now uninformative
(unreachable through the CLI, which passes assigned segments only).

The pairwise evidence was logged only for an undetermined verdict, so a reassortant
call recorded nowhere which pair disagreed or on which strains. Every verdict now logs
it and writes segment_pairs.tsv, and a discordant note names both near-best strain sets.

Audit 4, items E3, E14 and E4.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 15: `reassort --threads` for segment scans, and say when several clades tie  (audit items: E5 minor -- no thread option; per-segment clade tie-break)

`reassort` has no thread option, so `--scan-segments` aligns every segment with the `MsaParams` default of one thread (`scan.py:157-158`). It gains `-t/--threads` (default 1, so nothing changes unless asked), passed through `assign_segments` to `scan_segment` and `MsaParams`.

The per-segment clade is the top-ranked tip's (`assign.py:118-131`); tips are ranked by ANI (`assign.py:119`, a stable sort) and an exact tie keeps the pool's order, which `cached_genomes` sorts by path (`core/cache.py:217`), so the first file name wins. That is deterministic for a given cache, so this is an arbitrary-choice and reporting problem, not a reproducibility one: the `clade` column of `reassortment.tsv` (`cmd_reassort.py:95-96`) names one clade while the constellation call uses the whole near-best set, and nothing says several clades sat inside it. A pure helper, `clades_within_margin`, lists them, and the typing step logs a line when there is more than one. The reported clade and the `reassortment.tsv` columns are unchanged.

**Files:**
- Modify: `src/tessera/reassort/scan.py:109-113`, `:157-158`
- Modify: `src/tessera/reassort/assign.py:85` (new `clades_within_margin` before `_type_segment`), `:127-129`, `:135-147`, `:197-199`
- Modify: `src/tessera/cli/cmd_reassort.py:51-53`, `:62`, `:80-85`
- Modify: `docs/reference-panels.md` (the `--scan-segments` paragraph and the typing paragraph)
- Modify: `tests/unit/test_reassort_assign.py` -- **three existing tests change**: the `fake_scan` stubs in `test_scan_segments_scans_assigned_only`, `test_scan_segments_marks_unassigned` and `test_scan_segments_gives_alike_named_segments_separate_directories` declare the old `scan_segment` signature (`dir_name=None):`) and fail with `TypeError: ... got an unexpected keyword argument 'threads'` once `assign_segments` passes it. Each becomes `dir_name=None, threads=1):`; nothing else in them changes.
- Modify: `tests/unit/test_cli_input_checks.py` (append)

**Interfaces:**
- Produces: `scan_segment(segment, seq, dataset, out_dir, *, aligner, cache_dir, logger, dir_name=None, threads: int = 1)`; `assign_segments(..., threads: int = 1, logger)`; `assign.clades_within_margin(ranked: list[Path], ani: dict[Path, tuple[float, float]], margin: float) -> list[str]`.

- [ ] **Step 1: Write the failing test**

In `tests/unit/test_reassort_assign.py`, change the three `fake_scan` signatures named above from `dir_name=None):` to `dir_name=None, threads=1):`, then append:

```python


# --- audit 4, E5: segment-scan threads and clade ties ---------------------------

def test_scan_segments_receive_the_thread_count(tmp_path, monkeypatch):
    """reassort had no thread option, so every segment scan aligned single-threaded."""
    tip = tmp_path / "pool" / "strainA.fasta"
    tip.parent.mkdir(parents=True)
    tip.write_text(">x\nACGT\n")
    _patch(monkeypatch,
           resolve=lambda fasta, override, *, email, logger: _DS("ds"),
           tips_by_path={"ds": [tip]},
           ani_by_path={"pool": {tip: (99.0, 99.0)}})
    monkeypatch.setattr(assign, "require_aligner", lambda aligner: None)
    seen = {}

    def fake_scan(segment, seq, dataset, out_dir, *, aligner, cache_dir, logger,
                  dir_name=None, threads=1):
        from tessera.reassort.scan import SegmentScan
        seen[segment] = threads
        return SegmentScan(segment, True, False, 0, "none")

    monkeypatch.setattr(assign, "scan_segment", fake_scan)
    q = _write_query(tmp_path, [("HA", "AAAA")])
    assign_segments(q, output=tmp_path / "out", scan_segments=True, threads=4, logger=LOG)
    assert seen == {"HA": 4}


def test_clades_within_margin_lists_every_clade_in_the_near_best_window(tmp_path) -> None:
    """The reported clade is the top tip's, ties broken by file name; the near-best
    window can hold several clades, and nothing said so."""
    tips = []
    for name, clade in (("a", "3C.2a"), ("b", "3C.3a"), ("c", "3C.2a"), ("d", "2")):
        tip = tmp_path / f"{name}.fasta"
        tip.write_text(f">{name} {clade}\nACGT\n")
        tips.append(tip)
    ani = {tips[0]: (99.5, 90.0), tips[1]: (99.4, 90.0), tips[2]: (99.2, 90.0),
           tips[3]: (97.0, 90.0)}
    assert assign.clades_within_margin(tips, ani, 0.5) == ["3C.2a", "3C.3a"]
    assert assign.clades_within_margin(tips, ani, 0.05) == ["3C.2a"]
```

Append to `tests/unit/test_cli_input_checks.py`:

```python


def test_reassort_passes_the_thread_count(tmp_path: Path, monkeypatch) -> None:
    from tessera.reassort.assign import ReassortmentResult

    seen = {}

    def fake_assign(query, **kwargs):
        seen.update(kwargs)
        return ReassortmentResult()

    monkeypatch.setattr("tessera.cli.cmd_reassort.assign_segments", fake_assign)
    result = runner.invoke(app, [
        "reassort", "--query", str(_segments(tmp_path)), "--output", str(tmp_path / "out"),
        "--threads", "4",
    ])
    assert result.exit_code == 0, result.output
    assert seen["threads"] == 4
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_reassort_assign.py tests/unit/test_cli_input_checks.py -q -k "thread_count or clades_within"`
Expected: 3 failed. Observed on `4b35e39`: `TypeError: assign_segments() got an unexpected keyword argument 'threads'`; `AttributeError: module 'tessera.reassort.assign' has no attribute 'clades_within_margin'`; `reassort` exits 2 with `Usage: root reassort [OPTIONS]` (no `--threads`). The three edited stubs pass before and after (an unused keyword with a default).

- [ ] **Step 3: Implement**

`src/tessera/reassort/scan.py`: the signature gains `threads: int = 1` after `dir_name: str | None = None,`, and the `build_msa` call becomes

```python
        build_msa(
            MsaParams(query=query, collection=collection, output=msa, aligner=aligner,
                      threads=threads),
            logger,
        )
```

`src/tessera/reassort/assign.py`: add before `_type_segment`

```python
def clades_within_margin(
    ranked: list[Path], ani: dict[Path, tuple[float, float]], margin: float
) -> list[str]:
    """The distinct clades of the tips within ``margin`` ANI of the best, best-first.

    ``ranked`` is best-first. A tip without a clade label is skipped.
    """
    if not ranked:
        return []
    floor = ani[ranked[0]][0] - margin
    clades: list[str] = []
    for tip in ranked:
        if ani[tip][0] < floor:
            break
        clade = _clade_of_tip(tip)
        if clade and clade not in clades:
            clades.append(clade)
    return clades
```

in `_type_segment`, after the `logger.info("[%s] %s -> %s (%s), ANI %.1f%%", ...)` call:

```python
    near_clades = clades_within_margin(ranked, ani, margin)
    if len(near_clades) > 1:
        # The reported clade is the top tip's; a tie in ANI is broken by tip file name.
        # Say when the near-best window holds several clades, so the clade column is not
        # read as a unique assignment.
        logger.info(
            "[%s] %d clades lie within %g ANI points of the best tip (%s); the clade "
            "reported is the highest-ANI tip's, ties broken by tip name.",
            seg, len(near_clades), margin, ", ".join(near_clades),
        )
```

`assign_segments` gains `threads: int = 1,` after `aligner: str = "mafft",`, and its `scan_segment(...)` call passes `dir_name=scan_dirs[s.segment], threads=threads))`.

`src/tessera/cli/cmd_reassort.py`: after the `aligner` option

```python
    threads: int = typer.Option(
        1, "-t", "--threads", help="Aligner worker threads for --scan-segments."),
```

after the margin check `_require_range(threads, "--threads", lo=1)`, and pass `threads=threads` to `assign_segments(...)`.

`docs/reference-panels.md`: in the `--scan-segments` paragraph, after "It needs an aligner (`--aligner`, default `mafft`)" add "and aligns with `--threads` workers (default 1)"; in the typing paragraph, after the sentence ending "gives a `(strain, clade, ANI)`." add "When several clades lie within `--margin` of the best tip the reported clade is the highest-ANI tip's, ties broken by tip name, and the log lists the clades in the window."

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_reassort_assign.py tests/unit/test_reassort_scan.py tests/unit/test_cli_input_checks.py tests/unit/test_cli_reassort.py -q`

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`

- [ ] **Step 6: Commit**

```bash
git add src/tessera/reassort/scan.py src/tessera/reassort/assign.py \
        src/tessera/cli/cmd_reassort.py docs/reference-panels.md \
        tests/unit/test_reassort_assign.py tests/unit/test_cli_input_checks.py
git commit -F - <<'EOF'
Give reassort a thread option, and log when several clades tie

reassort had no thread option, so --scan-segments aligned every segment with one
thread. It now takes -t/--threads (default 1, so a run that does not set it is
unchanged) and passes it to the segment alignments.

Each segment's clade is the top-ranked tip's, and an ANI tie keeps the pool's file-name
order, with nothing to say that several clades sat inside the near-best window. The
typing step now logs the clades in the window when there is more than one; the reported
clade and the output columns are unchanged.

Three existing tests stub scan_segment with its old signature; they gain the threads
parameter and nothing else.

Audit 4, item E5 (reassort threads, clade tie-break).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 16: The self-contained harnesses parse their arguments  (audit items: F1)

`validation/run_specificity.py:221-230` scans `sys.argv` for two flags by hand: `--help` starts a ten-replicate run, an unknown flag is ignored, and a missing or non-numeric value raises `IndexError` / `ValueError`. `run_regimes.py:231-237` and `characterise_false_regions.py:81-82` do the same. Each gets a `parse_args` built on `argparse` with the same options and defaults, `allow_abbrev=False` (so `--rep 3` is an error rather than silently `--reps 3`), and a check that counts are at least 1. `run_benchmark.py`, `run_coalescent_benchmark.py` and `run_caller_benchmark.py` hand-parse too; see Self-review notes.

**Files:**
- Modify: `validation/run_specificity.py:42-46` (import), `:221-230`
- Modify: `validation/run_regimes.py:38-40` (import), `:231-237`
- Modify: `validation/characterise_false_regions.py:19-22` (import), `:81-82`
- Create: `tests/unit/test_harness_arguments.py`

**Interfaces:**
- Produces: `parse_args(argv: list[str]) -> argparse.Namespace` in each of the three scripts (`reps`, `min_methods`; `run_regimes` also `no_cluster_lineages`; `characterise_false_regions` only `reps`). `main(argv)` keeps its signature.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_harness_arguments.py`:

```python
"""The self-contained harnesses parse their arguments (audit 4, item F1).

``run_specificity.py`` scanned ``sys.argv`` for two flags by hand: ``--help`` started a
ten-replicate run, an unknown flag was ignored and a missing value raised an IndexError.
``run_regimes.py`` and ``characterise_false_regions.py`` did the same.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_VALIDATION = Path(__file__).resolve().parents[2] / "validation"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _VALIDATION / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


HARNESSES = ["run_specificity", "run_regimes", "characterise_false_regions"]


def _no_work(monkeypatch, module) -> None:
    """Fail the test if a harness gets past argument parsing into a simulation."""
    def boom(*args, **kwargs):
        raise AssertionError("the harness started work")

    for attr in ("simulate_clonal", "simulate_recombinant", "simulate_near_identical",
                 "_load_specificity", "_scan", "scan"):
        if hasattr(module, attr):
            monkeypatch.setattr(module, attr, boom)


@pytest.mark.parametrize("name", HARNESSES)
def test_help_exits_without_running(name: str, monkeypatch, capsys) -> None:
    module = _load(name)
    _no_work(monkeypatch, module)
    with pytest.raises(SystemExit) as exc:
        module.main(["--help"])
    assert exc.value.code == 0
    assert "--reps" in capsys.readouterr().out


@pytest.mark.parametrize("name", HARNESSES)
@pytest.mark.parametrize("argv", [["--rep", "3"], ["--reps"], ["--reps", "three"],
                                  ["--reps", "0"]])
def test_bad_arguments_are_an_error(name: str, argv: list[str], monkeypatch) -> None:
    module = _load(name)
    _no_work(monkeypatch, module)
    with pytest.raises(SystemExit) as exc:
        module.main(argv)
    assert exc.value.code == 2


def test_specificity_defaults_are_unchanged() -> None:
    """Pins the defaults the hand parser had: 10 replicates, gate 2. It fails before the
    change only because parse_args is new."""
    args = _load("run_specificity").parse_args([])
    assert (args.reps, args.min_methods) == (10, 2)
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_harness_arguments.py -q`
Expected: 16 failed. Observed on `4b35e39`: the `--help` and `--rep 3` cases fail with `AssertionError: the harness started work`; `--reps` with `IndexError: list index out of range`; `--reps three` with `ValueError: invalid literal for int() with base 10: 'three'`; `--reps 0` with `Failed: DID NOT RAISE <class 'SystemExit'>` for `run_specificity` and `run_regimes` (zero replicates run nothing and return 0) and `AssertionError: the harness started work` for `characterise_false_regions`; the defaults test with `AttributeError: module 'run_specificity' has no attribute 'parse_args'`.

- [ ] **Step 3: Implement**

`validation/run_specificity.py`: add `import argparse` to the imports and replace the head of `main` (lines 221-230) with

```python
def parse_args(argv: list[str]) -> argparse.Namespace:
    """The harness options. ``--help`` prints them and exits before any simulation; an
    unknown flag or a non-integer value is an error, not silently ignored."""
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        prog="run_specificity.py",
        description="Count the false regions Tessera reports on simulated clonal data, "
        "with a matching positive control. Needs no aligner, network or data.",
    )
    parser.add_argument("--reps", type=int, default=10,
                        help="replicates per scenario (default 10)")
    parser.add_argument("--min-methods", type=int, default=2,
                        help="agreement gate passed to the scan (default 2, the harness "
                        "default; the CLI default is 1)")
    args = parser.parse_args(argv)
    if args.reps < 1 or args.min_methods < 1:
        parser.error("--reps and --min-methods must be at least 1")
    return args


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    reps, min_methods = args.reps, args.min_methods
```

`validation/run_regimes.py`: add `import argparse` and replace lines 231-237 with

```python
def parse_args(argv: list[str]) -> argparse.Namespace:
    """The harness options; ``--help`` exits before any simulation."""
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        prog="run_regimes.py",
        description="Simulate the near-identical and sibling-tract regimes and score the "
        "scan. Needs no aligner, network or data.",
    )
    parser.add_argument("--reps", type=int, default=10, help="replicates (default 10)")
    parser.add_argument("--min-methods", type=int, default=1,
                        help="agreement gate passed to the scan (default 1)")
    parser.add_argument("--no-cluster-lineages", action="store_true",
                        help="turn lineage clustering off")
    args = parser.parse_args(argv)
    if args.reps < 1 or args.min_methods < 1:
        parser.error("--reps and --min-methods must be at least 1")
    return args


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    reps, min_methods = args.reps, args.min_methods
    cluster = not args.no_cluster_lineages
```

`validation/characterise_false_regions.py`: add `import argparse` and replace lines 81-82 with

```python
def parse_args(argv: list[str]) -> argparse.Namespace:
    """The harness options; ``--help`` exits before any simulation."""
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        prog="characterise_false_regions.py",
        description="List each false region the callers report on the specificity "
        "harness's clonal data. Needs no aligner, network or data.",
    )
    parser.add_argument("--reps", type=int, default=10,
                        help="replicates per scenario (default 10)")
    args = parser.parse_args(argv)
    if args.reps < 1:
        parser.error("--reps must be at least 1")
    return args


def main(argv: list[str]) -> int:
    reps = parse_args(argv).reps
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_harness_arguments.py tests/unit/test_specificity_scoring.py tests/unit/test_regime_scoring.py -q`

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`

- [ ] **Step 6: Commit**

```bash
git add validation/run_specificity.py validation/run_regimes.py \
        validation/characterise_false_regions.py tests/unit/test_harness_arguments.py
git commit -F - <<'EOF'
Parse the self-contained harnesses' arguments with argparse

run_specificity.py, the script CLAUDE.md and CONTRIBUTING.md name for checking the
false-positive rate, scanned sys.argv for two flags by hand: --help started a
ten-replicate run, an unknown flag was ignored silently, and a missing or non-numeric
value raised IndexError or ValueError. run_regimes.py and characterise_false_regions.py
did the same.

Each now has a parse_args with the same options and defaults; --help exits before any
simulation, abbreviations are refused (so --rep 3 is an error, not --reps 3), and a
count below 1 is an error. The documented invocations work unchanged.

Audit 4, item F1.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 17: Development commands and the measured false-region rate in the documentation  (audit items: F2, F3, F4, F5)

- **F2.** `CLAUDE.md:68` and `validation/README.md:115`, `:235`, `:256` prepend the aligner env to `PATH`, which also puts that env's Python first, so the harness runs under an interpreter without Tessera installed; `CONTRIBUTING.md:42` does the same for the `tessera` env. All become the appending form the rest of `validation/README.md` already uses. `CLAUDE.md` is checked-in project instructions; only this one line changes.
- **F3.** `CONTRIBUTING.md:31-32` says mypy is not blocking in CI and has a backlog; CI runs it as a blocking step (`.github/workflows/ci.yml:52-53`) and it is clean over `src/tessera`.
- **F4.** `README.md` (Development) lints `src tests`; CI and `CONTRIBUTING.md` lint `src tests validation`.
- **F5.** The README does not mention `--min-methods` or the false-region rate at the shipped default, although the open decision record (`docs/superpowers/specs/2026-10-01-min-methods-default-decision.md`, option C) says it does. A short section quotes the record's numbers with their replicate counts and intervals. The default is not changed.

**Files:**
- Modify: `CLAUDE.md:68`, `validation/README.md:115`, `:235`, `:256`, `CONTRIBUTING.md:31-32`, `:42`, `README.md:181-198`
- Create: `tests/unit/test_docs_commands.py`

**Interfaces:** none.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_docs_commands.py`:

```python
"""The development commands the documentation gives are the ones that work (audit 4, F2-F4).

Prepending the aligner environment to PATH also puts that environment's Python first,
so the harness ran under an interpreter without Tessera installed. The README lint line
left out ``validation``, which CI lints, and CONTRIBUTING.md said mypy was not blocking.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
DOCS = ["CLAUDE.md", "CONTRIBUTING.md", "README.md", "validation/README.md"]
PREPEND = re.compile(r'export PATH="\$HOME/[^"]*/bin:\$PATH"')


@pytest.mark.parametrize("name", DOCS)
def test_no_document_prepends_an_environment_to_path(name: str) -> None:
    text = (ROOT / name).read_text()
    assert not PREPEND.search(text), PREPEND.search(text).group(0)


def test_readme_lints_what_ci_lints() -> None:
    assert "ruff check src tests validation" in (ROOT / "README.md").read_text()


def test_contributing_states_that_mypy_is_blocking() -> None:
    text = (ROOT / "CONTRIBUTING.md").read_text()
    assert "not yet blocking" not in text
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_docs_commands.py -q`
Expected: 5 failed, 1 passed (`README.md` has no `PATH` line). Observed on `4b35e39`: `AssertionError: export PATH="$HOME/miniforge3/envs/recomfi-aln/bin:$PATH"` (`CLAUDE.md`, `validation/README.md`), `AssertionError: export PATH="$HOME/miniforge3/envs/tessera/bin:$PATH"` (`CONTRIBUTING.md`), the README lint line, and `assert 'not yet blocking' not in ...`.

- [ ] **Step 3: Implement**

`CLAUDE.md:68`:

```
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"   # append; prepending swaps the Python
```

`validation/README.md:115` and `:235`: `export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"`; `:256`: `export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"   # skani on PATH`.

`CONTRIBUTING.md:42`: `export PATH="$PATH:$HOME/miniforge3/envs/tessera/bin"   # append, never prepend`. `CONTRIBUTING.md:31-32`, replace the mypy bullet (which holds two non-ASCII dashes) with

```markdown
- `mypy` is clean. CI runs it as a blocking step over `src/tessera`.
```

`README.md`, Development section:

```
pip install -e ".[dev]"
ruff check src tests validation
mypy
pytest                       # add -m "not requires_binary" to skip aligner-dependent tests
```

`README.md`, Known limitations: add after the existing paragraph

```markdown
**False regions at the default setting.** On simulated clonal data -- four clades about
16 % apart, four tips per clade about 2.4 % apart and a non-recombinant query, so every
reported region is false -- the default `--min-methods 1` reported at least one false
region in 17 of 40 runs (42 %, 95 % CI 29-58 %); most were HMM switches between two tips
of the query's own clade. Requiring two agreeing callers (`--min-methods 2`) lowered this
to 2 of 40 (5 %, CI 1-17 %) and detected 9 of 10 simulated 3 kb tracts, against 10 of 10
at the default. On panels of near-identical references (about 0.2 % apart) the rates were
1 of 20 (CI 1-24 %) and 0 of 20 (CI 0-16 %). The default stays at 1 because a gate of two
cost three true detections on the curated hybrid benchmark (`rsv_a`, `mpox`,
`masksib_rsv`): at very low divergence only the HMM has power. For a panel holding several
near-equidistant relatives of the query, consider `--min-methods 2`, and treat a region
found by one caller as a candidate to confirm. These are simulation results with modest
replicate counts; see `validation/README.md` and
`docs/superpowers/specs/2026-10-01-min-methods-default-decision.md`.
```

(Source of every number: the decision record's "What was measured" tables and option A, measured on 2026-10-02 at 10 replicates on the plan C branch this plan builds on.)

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_docs_commands.py -q`

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`

- [ ] **Step 6: Commit**

```bash
git add CLAUDE.md CONTRIBUTING.md README.md validation/README.md \
        tests/unit/test_docs_commands.py
git commit -F - <<'EOF'
Document development commands that work and the default's false-region rate

CLAUDE.md, CONTRIBUTING.md and three places in validation/README.md prepended the
aligner environment to PATH, which also puts that environment's Python first, so the
harness ran under an interpreter without Tessera installed; they now append, as the
rest of validation/README.md already did. CONTRIBUTING.md said mypy was not blocking in
CI and had a backlog; CI runs it as a blocking step and it is clean. The README lint
line left out validation/, which CI lints.

The README now states the measured false-region rate at the shipped --min-methods 1
(17 of 40 simulated clonal runs, CI 29-58 %) beside the cost of --min-methods 2, with
the replicate counts and intervals from the open decision record. The default is not
changed; the record says documenting the rate is the step it expects to happen either
way.

Audit 4, items F2-F5.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 18: Changelog, the region-call check, pull request  (audit items: none; the Global Constraints check)

**Files:**
- Modify: `CHANGELOG.md` (under `## [Unreleased]`, a new subsection at its top)
- Create (untracked): `validation/data/audit-f/spec_after.txt`

- [ ] **Step 1: Changelog**

Under `## [Unreleased]` in `CHANGELOG.md`, before the existing "Fixes from the post-1.2.0 audit" paragraph, add

```markdown
### Fourth audit, plan F: reporting, CLI and reassortment

From `docs/superpowers/specs/2026-10-12-audit-swot-gap.md` (B1, B3, B5, B7-B9, the smaller
reporting items, E1-E4, E5 in part, E14, F1-F5; help text and documentation for B2, B4,
B6, A4). None of these changes which regions are called or their coordinates: the
specificity harness gives the same output before and after
(`run_specificity.py --reps 3 --min-methods 1`: 7/12 runs, 8 regions, all HMM).

- A test or caller that gave no result is reported as `not run`, with the reason, in the
  report, the TSVs and `run_provenance.json`: the PHI test under `--no-phi`, with too few
  informative sites or when not testable (`parent_free_support` now reads `not run`
  instead of `no`), and an HMM that lineage clustering left with one state.
- A `low_information` coverage gap is no longer announced as a possible missing
  reference.
- `run_provenance.json` gains `coordinates`; under informative-site windowing
  `window / step` says it describes the base-pair scan.
- The report headline names the run's reconciled major parent; the pair plot uses the
  report's colours; `recombination_methods.tsv` has one row per row of
  `recombination_regions.tsv`, donor-absent regions included.
- Every command rejects an output path that is a file. The discover commands range-check
  the coverage fractions and search counts (`--top-gaps 0` is an error). `msa` refuses
  `--reference` with `--query-as-backbone` and warns before overwriting.
  `find-references -q` is kept verbatim. `build-panel` takes `--pool-consensus`,
  `--keep-recombinant-lineages` and `--deep-typing`.
- `reassort --margin` must be above 0 and at most 5 ANI points, and an `--ani-floor`
  of 1 or below is rejected as a likely fraction. A discordant pair inside one parent
  group gives `undetermined` instead of `reassortant`; every pair's relation is written
  to `segment_pairs.tsv`. `reassort` takes `--threads` for `--scan-segments`.
- The help of `--min-region`, `--margin` and `--merge-gap` names the callers that read
  them; the validation harnesses parse their arguments with `argparse`.
```

- [ ] **Step 2: The region-call check (nothing else running)**

```bash
python validation/run_specificity.py --reps 3 --min-methods 1 > validation/data/audit-f/spec_after.txt 2>&1
diff <(sed 1,2d validation/data/audit-f/spec_before.txt) <(sed 1,2d validation/data/audit-f/spec_after.txt) && echo IDENTICAL
```

Expected: `IDENTICAL`. Measured while writing this plan on the full candidate set: identical (`TOTAL 7/12 (58%, CI 32%-81%) 8 hmm=8`; positive control `3/3 | 3/3 | 13 bp`). If it differs, do not open the pull request: bisect the branch with the same command and fix the task that changed a region.

- [ ] **Step 3: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: clean, clean, and 844 passed, 5 deselected (757 at the base plus 87 new tests; measured on the full candidate set while writing this plan -- quote the real number).

- [ ] **Step 4: Commit and open the pull request**

```bash
git add CHANGELOG.md
git commit -F - <<'EOF'
Record the plan F reporting, CLI and reassortment changes in the changelog

The specificity harness gives identical output before and after the branch
(run_specificity.py --reps 3 --min-methods 1: 7/12 runs with a false region, 8 regions,
all HMM; positive control 3/3 detected, 3/3 donor, 13 bp), which is the check that none
of these changes altered which regions are called or their coordinates.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
git push -u origin fix-audit4-reporting
gh pr create --base fix-audit-caller-behaviour --title "Audit 4, plan F: reporting faithfulness, CLI and reassortment" --body-file - <<'EOF'
## Summary
Reporting, CLI and reassortment fixes from the fourth audit (items listed in CHANGELOG.md).
No change to which regions are called or their coordinates.

## Region-call check
`python validation/run_specificity.py --reps 3 --min-methods 1`, before and after: identical
(7/12 runs with a false region, 8 regions, hmm = 8; positive control 3/3, 3/3, 13 bp).

## Not done here
- A one-state HMM does not lower the `--min-methods` gate; `--merge-gap` is not passed to the
  coverage diagnostic. Both would change reported regions (plan G).
- See the plan's Self-review notes for items dropped or deferred.

## Tests
`ruff`, `mypy` clean; `pytest -m "not requires_binary"`: <count> passed.

<the Claude Code attribution line: the robot emoji, then "Generated with [Claude Code](https://claude.com/claude-code)">
EOF
```

The attribution line is written as a description here to keep this file ASCII; put the real line in the body. Base the pull request on `main` instead if `fix-audit-caller-behaviour` has been merged by then.

---

## Self-review notes

**How this plan was verified.** Every task's code was applied, cumulatively, in a detached worktree of `4b35e39`, and every new test was run against both that tree and an untouched `4b35e39` tree. Each "Observed" failure in a Step 2 is the message seen on the untouched tree; guards are marked as such and passed on both. On the full candidate set: `ruff` clean, `mypy` clean over 79 files, `pytest -m "not requires_binary"` 844 passed. The specificity harness (`--reps 3 --min-methods 1`) gave identical output before and after (7/12 runs with a false region, 8 regions, all HMM; positive control 3/3 detected, 3/3 donor, 13 bp). The harnesses that need an aligner or the network were not run. Task order matters only where a later task names an earlier one's interface (Tasks 2 -> 3 -> 8 share `not_run`; Task 10's import list is extended in Task 11; Task 13's `MAX_MARGIN` is used in Task 14's module).

**Deferred to plan G, because each would change which regions are reported:**
- A one-state HMM (Task 3) does not lower the `--min-methods` gate, unlike an unrunnable barcode caller; and an HMM-only run with one state is reported, not refused (a barcode-only run is refused). Whether either should change is a calling decision.
- B6: `run.py:443-447` does not pass `merge_gap` to `CoverageParams`. Wiring it changes the coverage gaps, which become donor-absent regions, whenever `--merge-gap` differs from the window size.
- B4 code half (recomputing or relabelling `n_windows` / `posterior` for a trimmed HMM span), and enforcing `--min-region` across callers (B2, A4). This plan documents both.

**Not planned, with the reason:**
- E5 "a computed return value is discarded" (located by the audit pass after the plan was drafted): `reassort/scan.py:154` unpacks `window, step, _sel = window_params(len(seq))` and drops the selection window that `scan.py:41-47` computes. The segment scan builds its panel from per-clade consensus genomes with no regional selection, so the value has no consumer; the function mirrors `validation/run_hybrids.py:202`, where the selection window is used, and `tests/unit/test_reassort_scan.py:65-66` pins the three-tuple. Dead code with no effect on any output, so no task: removing it would change a tested signature for no user-visible gain. A maintainer who prefers it gone can drop the third element and update that test.
- E5 "`--metric` has one legal value": the value is validated (`similarity.py:239-241` raises a `UserInputError`) and the help says `pdist`. Removing the option breaks existing command lines; that is the maintainer's call (spec G6 says "remove or implement").
- E5 items assigned to plan E: `--candidate-pool` silently overriding `--nextclade`, `--pool-consensus` being inert outside a Nextclade pool, and the candidate-pool reader's extension list. Task 12 adds `--pool-consensus` to `build-panel` with help that says it applies with `--nextclade`; it does not touch the precedence or the reader.
- `validation/run_benchmark.py`, `run_coalescent_benchmark.py` and `run_caller_benchmark.py` hand-parse `sys.argv` the same way (F1 names only `run_specificity.py`; Task 16 also covers `run_regimes.py` and `characterise_false_regions.py`, which the plan's scope asked to check). A follow-up can give them the same `parse_args`.
- Found while planning, not in the spec: a donor-absent region's `parent_free_support` is never computed (`coverage.gaps_as_regions` leaves it `False`), so it reads `no` even when PHI is significant and an Rmin interval overlaps it. Fixing it changes a reported value; left for the maintainer. Task 8's method row copies the region's value, so the two files agree.

**Judgment calls a reviewer may want to revisit:**
- `reassort --margin` upper bound of 5 ANI points (Task 13) is a sanity bound, ten times the default, not a calibrated value; it rejects the values that inverted the verdict (40, 100, 1e6) and leaves a wide useful range.
- E3 rule (Task 14): a discordant pair inside one concordant group gives `undetermined`, not a split into groups, because the split would be arbitrary when the middle segment is concordant with both ends.
- The "six count options" of E5 are not named in the spec; Task 11 bounds `--max-hits`, `--top-gaps`, `--subtile`, `--seed-hits`, `--seed-window` and `--threads`, and leaves `--fetch-limit` to plan E (D8).

**Where the spec is wrong or stale (as checked against `4b35e39`):**
- B9 is reproducible, so CONFIRMED rather than SUSPECTED, but only without a called region: the HMM's major must be a pooled cluster's representative that wins some windows, so `reconcile_major`'s zero-win guard does not replace it (Task 7 has the fixture).
- E1 affects seven commands, not four plus `fill-references`: `find-references` and `curate-panel` never check the output path either (they reach the pipeline). The section 3 test gap "the five uncovered commands" undercounts likewise.
- B1 is wider than `--no-phi`: the same "no" was written when PHI had fewer than ten informative sites or was not testable at its window. Task 2 covers all three.
- F2 also applies to `CONTRIBUTING.md:42` (the `tessera` env), and `validation/README.md` has three prepending lines (`:115`, `:235`, `:256`), not one.
- Line numbers in B2 and B6 are 16 lines off: `min_region` is read at `regions.py:316` and `:442` (spec: 300, 426), `margin` at `:414` and `merge_gap` at `:432` (spec: 398, 416).
