# Audit Plan C: Caller Behaviour Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the caller-behaviour defects C1-C7 of the post-1.2.0 audit, each behind a measured harness gate, and produce the numbers the maintainer needs for the `--min-methods` decision (C8).

**Architecture:** Every task changes what the scan reports, so every task is evaluate-first: a failing test, one candidate implementation, then a harness gate with numeric acceptance criteria and a stated action when the gate fails. Independent fixes come first (3SEQ/MaxChi tract start, p-value underflow, ensemble merge, pool typing, PHI window), then the two that interact (HMM region span, then lineage clustering, which only makes sense once HMM spans are tight), and last a characterisation of the false-positive rate that ends in a decision record rather than a code change. A new self-contained harness (`validation/run_regimes.py`) simulates the two regimes the audit used and no existing harness samples.

**Tech Stack:** Python 3.11+, numpy, pytest, ruff, mypy. No new dependency. Harnesses: `validation/run_specificity.py`, `run_regimes.py` (new), `run_validation.py`, `run_benchmark.py`, `run_hybrids.py`.

**Spec:** `docs/superpowers/specs/2026-10-01-post-1.2.0-audit-design.md` (section "Plan C"). Read it before starting; the "Deviations from the spec" section below lists where this plan departs from it and why.

## Global Constraints

- No new runtime dependency.
- Modest scientific language in code, docs and messages; reported numbers must be faithful (state what passes, what fails, what was skipped).
- A new behaviour needs a test that fails without the change. Where a test in this plan is a guard that passes before and after, the plan says so.
- Anything touching a caller, a default or region calling is validated on `validation/run_specificity.py` at **both** `--min-methods 1` and `--min-methods 2` and on `validation/run_hybrids.py`, before and after, and the numbers go in the PR description.
- "Could not test" must never be reported as "tested, found nothing".
- Work on branch `fix-audit-caller-behaviour`. Commit messages explain the problem, the evidence and the trade-off, and end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`.
- Append the aligner env to `PATH`; never prepend it (prepending swaps the Python interpreter): `export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"`.
- **A failed gate is a result, not an obstacle.** If a task's gate fails: `git revert` that task's commit, write the measured numbers under "Gate results" in the PR description, and move on to the next task. Do not adjust thresholds, constants or the test until the gate passes.
- The `--min-methods` default is not changed in this branch. Task 10 presents it to the maintainer as a decision.
- Do not touch plan A or plan B items. In particular the PHI "not testable" wording and the coverage-gap relabelling belong to plan B.

## What was measured while writing this plan

The candidates below were applied in a scratch worktree of `main` at `457bdfb` and measured. These numbers are the reference the gates compare against; the executor re-measures at 10 replicates in Task 1. Replicate counts are small, so the intervals are wide -- read them as "does the failure mode exist", not as rates.

| measurement | main (`457bdfb`) | with C3, C4, C7, C2 | with all candidates |
|---|---|---|---|
| `run_specificity.py --reps 3 --min-methods 2`: runs with a false region | 0/12 | 0/12 | 0/12 |
| `run_specificity.py --reps 3 --min-methods 1`: runs with a false region | 7/12 (8 regions, hmm = 8) | 7/12 (8, hmm = 8) | 7/12 (8, hmm = 8) |
| positive control (both gates): detected / donor / median breakpoint error | 3/3, 3/3, 55 bp | 3/3, 3/3, 13 bp | 3/3, 3/3, 13 bp |
| `run_validation.py` | 6 PASS, 1 FAIL (`hcv_clonal_1b`), 1 SKIP | same verdicts, same regions, coordinates change | same |
| `run_regimes.py --reps 10`, near-identical: detected / HMM among callers / donor / error | 10/10, 0/10, 10/10, 465 bp | 6/6, 0/6, 6/6, 263 bp (6 replicates; `main` gives the same at 6) | 10/10, 10/10, 10/10, 465 bp |
| near-identical clonal, 20 panels, `--min-methods 1`: runs with a false region | 0/20 | not run at 20 (0/6 in the regime harness) | 1/20 (hmm) |
| near-identical clonal, 20 panels, `--min-methods 2` | 0/20 | not run | 0/20 |
| sibling tract (true difference 600 bp), 10 reps: median reported length | 11700 bp | 1001 bp (6 replicates) | 1001 bp |
| `run_benchmark.py` (PHI): power / specificity | 8/185, 36/36 | -- | 13/185, 36/36 |

**Not measured here:** `run_hybrids.py`. It downloads the Nextclade index, reference and tree for every case, and this session ran without network. It is an executor gate for every task that names it. C3 and C4 were measured together, not separately. C7 was measured only inside the cumulative stack.

## Deviations from the spec

1. **C2 -- the sign test stays on the HMM segment; only the reported span is trimmed.** The spec says to recompute the sign test on the trimmed span. Trimming removes flanking sites that favour the major, so recomputing would lower p-values on a span chosen by looking at the data, in the caller that already produces every false region at `--min-methods 1`. Keeping the test on the segment means the set of called regions cannot change through this step; the measurements above confirm it (same eight false regions, same validation verdicts). The similarities are recomputed on the reported span.
2. **C2 -- a segment touching the first or last window is searched to the alignment end.** Trimming can only shrink a span, so it cannot fix a tract that starts at column 0 and is reported from the first window centre. This small extension does.
3. **C2 -- one alternative was tried and not chosen.** Taking the 3SEQ maximum-descent interval over the full extent of the segment's windows gave a slightly lower breakpoint error on the positive control (9 bp against 13 bp) but changed rows on real data: on `hiv_crf02ag` one HMM region grew by 61 bp, absorbed a coverage gap and removed a donor-absent row. The chosen rule changes coordinates only.
4. **C5 -- the rule is `min(--phi-window, sites // 10)`.** The spec left the rule open. Six rules were compared (Task 7). The benchmark does not separate them strongly: PHI power on it is low under every rule (see Task 7's note on a finding outside this plan's scope).
5. **C6 -- pool labels are passed explicitly, not inferred from the header shape.** The spec says "when a header comes from a pool, read the second token". Cached pools hold multi-word clades (`clade 2 wild-type`, `Clade VI`), so the label is everything after the accession; and a user's `--candidate-pool` holds free-text titles, so whether a header is structured is decided by where the genome came from, carried in a small sidecar.
6. **C1 -- "specificity no worse" holds on `run_specificity.py` but not strictly on the new near-identical clonal panels** (1/20 against 0/20 at `--min-methods 1`). That is the cost of the HMM having a vote there at all; Task 9 states it and Task 10 carries it into the decision.

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `validation/run_regimes.py` | create | simulate and score the near-identical and sibling-tract regimes (opt-in harness) |
| `tests/unit/test_regime_scoring.py` | create | unit tests for that harness's simulation and scoring |
| `tests/integration/test_regimes.py` | create | one-seed pipeline regressions for C1 and C2 |
| `src/tessera/recomb/threeseq.py` | modify | C3 tract start, C4 exact p-value |
| `src/tessera/recomb/maxchi.py` | modify | C3: the permutation null uses the same tract rule |
| `src/tessera/recomb/ensemble.py` | modify | C7 merge edge cases |
| `src/tessera/recomb/typing.py` | modify | C6: read a pool label verbatim |
| `src/tessera/discover/iterate.py` | modify | C6: pass pool labels to selection and panel typing |
| `src/tessera/recomb/diagnostics.py` | modify | C5: cap the PHI window |
| `src/tessera/cli/cmd_recomb.py` | modify | C5: `--phi-window` help text |
| `src/tessera/recomb/regions.py` | modify | C2: reported span of an HMM region |
| `src/tessera/recomb/clusters.py` | modify | C1: pairwise agreement on informative columns |
| `validation/characterise_false_regions.py` | create | C8: one row per false region |
| `docs/superpowers/specs/2026-10-01-min-methods-default-decision.md` | create | C8: decision record |
| `docs/detection-methods.md`, `docs/reference-panels.md`, `validation/README.md`, `CHANGELOG.md` | modify | describe each change where it lands |

## Review Focus

Inputs the spec implies but its own examples do not exercise. Each has a test in the task that owns the code.

1. **A query with gaps upstream of a region.** After the span is trimmed, `query_start` / `query_end` must still be the query coordinates of the reported columns. *Task 8, `test_trimmed_region_reports_query_coordinates_across_a_query_gap`.*
2. **`--min-methods 2` after trimming.** Agreement is counted on overlapping regions, so a trimmed HMM region must still overlap the 3SEQ / MaxChi tract or a real event is dropped. *Task 8, `test_trimmed_hmm_region_still_merges_with_the_site_callers`.*
3. **`--informative-sites` forced on a divergent panel that holds true duplicates.** Clustering then compares on informative columns too; duplicates must still pool. *Task 9, `test_duplicates_still_merge_when_informative_windowing_is_forced`.*
4. **An alignment with exactly ten informative sites.** The capped PHI window is then one rank; the signal must still be computed. *Task 7, `test_signal_at_the_minimum_number_of_informative_sites`.*
5. **A user's `--candidate-pool` whose FASTA headers end in a single word.** That word must not be taken for a clade. *Task 6, `test_select_from_still_mines_a_user_pool`.*

---

### Task 1: Branch and baseline

Nothing is committed in this task. Its output is a directory of harness results that every later gate compares against. `validation/data/` is git-ignored, so the results live there.

**Files:**
- Create (untracked, git-ignored): `validation/data/audit-c/baseline/*.txt`

**Interfaces:**
- Produces: `validation/data/audit-c/baseline/` holding `spec_mm1.txt`, `spec_mm2.txt`, `validation.txt`, `validation_regions.tsv`, `benchmark.txt`, `hybrids.txt`. Later tasks write the same file names under `validation/data/audit-c/<task>/` and diff against these.

- [ ] **Step 1: Branch**

```bash
git switch main && git pull --ff-only
git switch -c fix-audit-caller-behaviour
pip install -e ".[dev]"
```

- [ ] **Step 2: Confirm the starting point is green**

Run: `ruff check src tests validation && pytest -m "not requires_binary" -q && mypy src`
Expected: `All checks passed!`, `586 passed, 1 deselected`, `Success: no issues found in 79 source files`. If any of these differs, stop: the baseline is not the one this plan was written against.

- [ ] **Step 3: Write the gate script**

The same commands run after every task, so put them in one untracked script. Create `validation/data/audit-c/gate.sh`:

```bash
#!/bin/bash
# gate.sh <label>: run the offline-capable harnesses and save their output under
# validation/data/audit-c/<label>/. Run from the repository root.
set -u
out="validation/data/audit-c/$1"
mkdir -p "$out"
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"   # appended, never prepended
python validation/run_specificity.py --reps 10 --min-methods 2 > "$out/spec_mm2.txt" 2>&1
python validation/run_specificity.py --reps 10 --min-methods 1 > "$out/spec_mm1.txt" 2>&1
python validation/run_validation.py > "$out/validation.txt" 2>&1
python - > "$out/validation_regions.tsv" <<'EOF'
import csv
from pathlib import Path
for d in sorted(Path("validation/data").iterdir()):
    f = d / "_run" / "recomb_out" / "recombination_regions.tsv"
    if not f.exists() or d.name == "orthopox_example":
        continue
    for r in csv.DictReader(open(f), delimiter="\t"):
        print(d.name, r["minor_parent"], r["major_parent"], r["query_start"], r["query_end"],
              r.get("methods") or "-", r.get("donor_absent", ""), r.get("donor_undercovered", ""),
              sep="\t")
EOF
python validation/run_benchmark.py > "$out/benchmark.txt" 2>&1
grep -E "TOTAL|detected" "$out/spec_mm2.txt" "$out/spec_mm1.txt"
tail -9 "$out/validation.txt"
tail -2 "$out/benchmark.txt"
```

```bash
mkdir -p validation/data/audit-c && chmod +x validation/data/audit-c/gate.sh
```

- [ ] **Step 4: Capture the baseline**

Run: `validation/data/audit-c/gate.sh baseline` (about 15 minutes).

Expected, at 10 replicates (the 3-replicate values measured for this plan are in the table above):
- `spec_mm2.txt`: `TOTAL 0/40` or close to it; positive control `detected 10/10`.
- `spec_mm1.txt`: roughly half the runs carry a false region, almost all from `hmm`.
- `validation.txt`: `PASS` for `sarscov2_xbb`, `hiv1_crf`, `norovirus_gii`, `enterovirus_e11`, `hiv_crf02ag`, `hcv_2k1b`; `FAIL` for `hcv_clonal_1b` (a known open item: a 12 bp MaxChi-only region); `SKIP` for `orthopox_example`.
- `benchmark.txt`: `power 0.04 (300 recombining) | specificity 1.00 (60 clonal) | 139 untestable`.

What to do when something is absent:
- **No aligner on `PATH`** (`mafft`, `minimap2`): `run_validation.py` prints `SKIP` with `<binary> not on PATH` for every dataset. Create the env (`conda env create -f environment.yml`) and append its `bin` to `PATH`. Do not continue past Task 7 without it.
- **Datasets not fetched**: `run_validation.py` prints `SKIP ... sequences not present`. Run `python validation/fetch.py` (needs `efetch` and network).
- **No benchmark alignments**: `run_benchmark.py` prints `[SKIP] no benchmark alignments`. Download `performance.tar.gz` from Dryad (doi:10.5061/dryad.d7wm37q6f) in a browser and extract the `msa_*.fasta` files into `validation/data/benchmark/`. Without it, Task 7's gate cannot be evaluated; record that in the PR.

- [ ] **Step 5: Capture the hybrids baseline**

```bash
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
python validation/run_hybrids.py > validation/data/audit-c/baseline/hybrids.txt 2>&1
grep -A60 "^case " validation/data/audit-c/baseline/hybrids.txt | tail -60
```

This needs network (Nextclade index, references and trees) and `mafft`, `skani`, `skDER`; allow 30-60 minutes. Expected: a table with one row per case and a final line `N/M passed (S skipped, E error)` with `E = 0`.

If every row reads `ERROR: Could not fetch the Nextclade dataset index`, there is no network. Then the hybrids gate cannot be run for any task: say so in the PR description under "Skipped", and do not merge Tasks 6, 8 or 9 until someone has run it. `CONTRIBUTING.md` requires it for those three.

- [ ] **Step 6: Record the baseline**

Copy the summary lines of the five files into a scratch note; they open the PR description in Task 11. No commit.

---

### Task 2: Regime harness

The audit found two defects with simulations no harness covers. This task commits those simulations so the later gates are reproducible.

**Where it lives, and why both places.** `run_specificity.py` set the pattern: the simulation and scoring sit in an opt-in `validation/` script because rates need replicates and replicates take minutes, and CI unit-tests that script's pure logic by loading it by path. The same split is used here. The harness (`validation/run_regimes.py`) reports rates; the one-seed regressions that must hold on every commit go in the fast suite (`tests/integration/test_regimes.py`, added in Tasks 8 and 9) and import the simulators from the harness, so there is one copy of the simulation.

**Files:**
- Create: `validation/run_regimes.py`
- Create: `tests/unit/test_regime_scoring.py`
- Modify: `validation/README.md` (layout block, and a new section after the specificity section)

**Interfaces:**
- Produces, in `validation/run_regimes.py`:
  - `QUERY = "QUERY"`, `NEAR_REFS = ("R0", "R1", "R2", "R3", "R4")`
  - `simulate_near_identical(seed: int, *, recombinant: bool, length: int = 120_000, distance: float = 0.001, tract: tuple[int, int] = (40_000, 48_000)) -> tuple[dict[str, str], tuple[int, int] | None]`
  - `simulate_sibling_tract(seed: int, *, length: int = 30_000, distance: float = 0.05, tract: tuple[int, int] = (12_000, 18_000), extra: int = 600) -> tuple[dict[str, str], tuple[int, int]]`
  - `score_tract(rows, tract, *, donor) -> dict` with keys `detected`, `hmm`, `donor_ok`, `breakpoint_error`
  - `score_clonal(rows) -> tuple[int, dict[str, int]]`
  - `score_span(rows, span) -> dict` with keys `detected`, `reported_bp`, `excess_bp`
  - `scan(seqs, logger, *, min_methods: int = 1, cluster_lineages: bool = True) -> list[dict]`
  - `_evolve(seq: np.ndarray, distance: float, rng) -> np.ndarray`, `_as_text(seqs) -> dict[str, str]`
- Tests load the module by path as `rg`, exactly as `tests/unit/test_specificity_scoring.py` loads `run_specificity.py`.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_regime_scoring.py`:

```python
"""Unit tests for the regime harness's simulation and scoring (no binaries, no network)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_regimes.py"
_SPEC = importlib.util.spec_from_file_location("run_regimes", _PATH)
rg = importlib.util.module_from_spec(_SPEC)
sys.modules["run_regimes"] = rg
_SPEC.loader.exec_module(rg)

SMALL = {"length": 6000, "distance": 0.004, "tract": (2000, 3000)}


def _identity(a: str, b: str, lo: int, hi: int) -> float:
    return sum(1 for i in range(lo, hi) if a[i] == b[i]) / (hi - lo)


# --- simulation ------------------------------------------------------------

def test_near_identical_simulation_is_deterministic_and_aligned():
    a, tract = rg.simulate_near_identical(3, recombinant=True, **SMALL)
    b, _ = rg.simulate_near_identical(3, recombinant=True, **SMALL)
    assert a == b and tract == SMALL["tract"]
    assert rg.simulate_near_identical(4, recombinant=True, **SMALL)[0] != a
    assert rg.QUERY in a and set(rg.NEAR_REFS) <= set(a)
    assert len({len(s) for s in a.values()}) == 1


def test_near_identical_panel_is_near_identical():
    seqs, _ = rg.simulate_near_identical(1, recombinant=False, **SMALL)
    # two references each 0.4 % from the root are ~0.8 % apart
    assert 0.985 < _identity(seqs["R0"], seqs["R1"], 0, SMALL["length"]) < 0.998


def test_near_identical_recombinant_carries_the_donor_over_the_tract_only():
    seqs, tract = rg.simulate_near_identical(1, recombinant=True, **SMALL)
    lo, hi = tract
    query = seqs[rg.QUERY]
    assert query[lo:hi] == seqs["R1"][lo:hi]
    assert _identity(query, seqs["R0"], 0, lo) > _identity(query, seqs["R1"], 0, lo)


def test_near_identical_clonal_has_no_tract():
    seqs, tract = rg.simulate_near_identical(1, recombinant=False, **SMALL)
    assert tract is None
    assert _identity(seqs[rg.QUERY], seqs["R0"], 0, SMALL["length"]) > 0.998


def test_sibling_tract_differs_from_the_query_only_over_the_span():
    seqs, (lo, hi) = rg.simulate_sibling_tract(
        2, length=6000, distance=0.05, tract=(2000, 3000), extra=300
    )
    query, sibling = seqs[rg.QUERY], seqs["sibling"]
    assert (lo, hi) == (3000, 3300)
    assert query[:lo] == sibling[:lo] and query[hi:] == sibling[hi:]
    assert query[lo:hi] != sibling[lo:hi]
    assert query[lo:hi] == seqs["parent_A"][lo:hi]
    assert sibling[hi:] == seqs["parent_A"][hi:]  # identical downstream of the span


# --- scoring ---------------------------------------------------------------

def _row(**kw):
    row = {"minor_parent": "R1", "query_start": "100", "query_end": "200",
           "donor_absent": "no", "methods": "3seq,maxchi"}
    row.update(kw)
    return row


def test_score_tract_reports_detection_donor_and_error():
    got = rg.score_tract([_row(query_start="90", query_end="230")], (100, 200), donor="R1")
    assert got == {"detected": True, "hmm": False, "donor_ok": True, "breakpoint_error": 20}


def test_score_tract_notes_when_the_hmm_was_among_the_callers():
    rows = [_row(methods="hmm,3seq")]
    assert rg.score_tract(rows, (100, 200), donor="R1")["hmm"] is True


def test_score_tract_miss_and_wrong_donor():
    assert rg.score_tract([_row(query_start="500", query_end="600")], (100, 200),
                          donor="R1")["detected"] is False
    assert rg.score_tract([_row(minor_parent="R3")], (100, 200),
                          donor="R1")["donor_ok"] is False


def test_score_tract_ignores_donor_absent_rows():
    rows = [_row(donor_absent="yes", methods="")]
    assert rg.score_tract(rows, (100, 200), donor="R1")["detected"] is False


def test_score_clonal_counts_present_regions_per_caller():
    rows = [_row(methods="hmm"), _row(methods="hmm,3seq"), _row(donor_absent="yes")]
    assert rg.score_clonal(rows) == (2, {"hmm": 2, "3seq": 1})
    assert rg.score_clonal([]) == (0, {})


def test_score_span_reports_the_excess_over_the_true_difference():
    rows = [_row(query_start="17900", query_end="29600"),
            _row(query_start="18007", query_end="18591")]
    got = rg.score_span(rows, (18000, 18600))
    assert got == {"detected": True, "reported_bp": 11700, "excess_bp": 11100}
    tight = rg.score_span(rows[1:], (18000, 18600))
    assert tight["excess_bp"] == 0
    assert rg.score_span([], (18000, 18600))["detected"] is False
```

- [ ] **Step 2: Run them to see them fail**

Run: `pytest tests/unit/test_regime_scoring.py -q`
Expected: collection error, `FileNotFoundError: ... validation/run_regimes.py`.

- [ ] **Step 3: Write the harness**

Create `validation/run_regimes.py`:

```python
#!/usr/bin/env python
"""Opt-in harness: two regimes the other harnesses do not sample.

``run_specificity.py`` simulates clades ~16 % apart, so it always runs under base-pair
windowing, and its positive control has parents that differ everywhere. Two situations
found by the post-1.2.0 audit fall outside it:

    near_identical   a panel whose references are ~0.2 % apart (the mpox / VZV /
                     within-lineage SARS-CoV-2 regime, analysed under informative-site
                     windowing). Run twice: with a donor tract spliced into the query
                     (is it found, and does the HMM contribute?) and clonal (what is
                     reported when there is nothing to find?).
    sibling_tract    the panel holds a sibling recombinant whose donor tract runs 600 bp
                     further than the query's. Against that sibling the query differs
                     only over those 600 bp; downstream of them the donor and the
                     backbone are identical, so nothing marks where the region ends.

Like ``run_specificity.py`` this needs no aligner, no network and no downloaded data:
the simulated sequences are already aligned. It is opt-in because a full run is minutes
of wall clock; the simulation and scoring logic is unit-tested in CI
(``tests/unit/test_regime_scoring.py``).

    python validation/run_regimes.py                  # 10 replicates, --min-methods 1
    python validation/run_regimes.py --reps 3         # quick look
    python validation/run_regimes.py --min-methods 2  # with the agreement gate
    python validation/run_regimes.py --no-cluster-lineages

``--min-methods`` defaults to 1 here because that is the CLI default; the harness
reports what a user gets.

Caveat: a star tree under JC69 is simpler than real viral evolution, and a handful of
replicates carries real sampling error. These numbers show whether a failure mode
exists and roughly how large it is, not its magnitude on real panels.
"""

from __future__ import annotations

import logging
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np

from tessera.recomb.run import RecombParams, run_recomb

QUERY = "QUERY"
_BASES = np.frombuffer(b"ACGT", dtype=np.uint8)

# near_identical: five references on a star tree, each this far from the root.
NEAR_LENGTH = 120_000
NEAR_DISTANCE = 0.001
NEAR_REFS = ("R0", "R1", "R2", "R3", "R4")
NEAR_TRACT = (40_000, 48_000)  # donor tract (from R1) in the R0-derived query

# sibling_tract: parents 10 % apart; the sibling's tract runs 600 bp past the query's.
SIB_LENGTH = 30_000
SIB_DISTANCE = 0.05
SIB_QUERY_TRACT = (12_000, 18_000)
SIB_EXTRA = 600


# --- simulation ------------------------------------------------------------------

def _evolve(seq: np.ndarray, distance: float, rng) -> np.ndarray:
    """JC69: mutate each site with probability 3/4(1 - exp(-4/3 * d))."""
    p = 0.75 * (1.0 - np.exp(-4.0 / 3.0 * distance))
    hit = rng.random(seq.size) < p
    out = seq.copy()
    n = int(hit.sum())
    if n:  # a mutation always changes the base
        out[hit] = (seq[hit] + rng.integers(1, 4, size=n)) % 4
    return out


def _as_text(seqs: dict[str, np.ndarray]) -> dict[str, str]:
    return {k: _BASES[v].tobytes().decode("ascii") for k, v in seqs.items()}


def simulate_near_identical(
    seed: int, *, recombinant: bool, length: int = NEAR_LENGTH,
    distance: float = NEAR_DISTANCE, tract: tuple[int, int] = NEAR_TRACT,
) -> tuple[dict[str, str], tuple[int, int] | None]:
    """A near-identical panel and a query descended from ``R0``.

    With ``recombinant`` the query carries ``R1`` over ``tract``; otherwise it is a
    plain descendant of ``R0`` and nothing is recombined. Returns ``(sequences, tract)``
    with ``tract`` ``None`` for the clonal case.
    """
    rng = np.random.default_rng(seed)
    root = rng.integers(0, 4, size=length)
    seqs = {label: _evolve(root, distance, rng) for label in NEAR_REFS}
    query = _evolve(seqs["R0"], distance / 10.0, rng)
    if recombinant:
        lo, hi = tract
        query[lo:hi] = seqs["R1"][lo:hi]
    seqs[QUERY] = query
    return _as_text(seqs), (tract if recombinant else None)


def simulate_sibling_tract(
    seed: int, *, length: int = SIB_LENGTH, distance: float = SIB_DISTANCE,
    tract: tuple[int, int] = SIB_QUERY_TRACT, extra: int = SIB_EXTRA,
) -> tuple[dict[str, str], tuple[int, int]]:
    """A query and a sibling recombinant whose donor tract runs ``extra`` bp further.

    Both are a ``parent_A`` backbone carrying ``parent_B``; the query over ``tract``, the
    sibling over ``tract`` extended by ``extra``. Returns ``(sequences, span)`` where
    ``span`` is the only stretch on which the query and the sibling differ -- there the
    query matches ``parent_A``. Everywhere downstream ``parent_A`` and the sibling are
    identical.
    """
    rng = np.random.default_rng(seed)
    root = rng.integers(0, 4, size=length)
    parent_a, parent_b, other = (_evolve(root, distance, rng) for _ in range(3))
    lo, hi = tract
    query = parent_a.copy()
    query[lo:hi] = parent_b[lo:hi]
    sibling = parent_a.copy()
    sibling[lo:hi + extra] = parent_b[lo:hi + extra]
    seqs = {QUERY: query, "sibling": sibling, "parent_A": parent_a,
            "parent_B": parent_b, "other": other}
    return _as_text(seqs), (hi, hi + extra)


# --- scoring ---------------------------------------------------------------------

def _present(rows: list[dict]) -> list[dict]:
    """Rows that claim a donor; a ``donor_absent`` row is a coverage statement."""
    return [r for r in rows if r.get("donor_absent") != "yes"]


def _methods(row: dict) -> list[str]:
    return [m.strip() for m in (row.get("methods") or "").split(",") if m.strip()]


def score_tract(rows: list[dict], tract: tuple[int, int], *, donor: str) -> dict:
    """Detection of a known tract, whether the HMM was among the callers, donor
    attribution and breakpoint error."""
    lo, hi = tract
    overlapping = [
        r for r in _present(rows)
        if int(r["query_start"]) < hi and int(r["query_end"]) > lo
    ]
    if not overlapping:
        return {"detected": False, "hmm": False, "donor_ok": False,
                "breakpoint_error": None}
    best = max(
        overlapping,
        key=lambda r: min(int(r["query_end"]), hi) - max(int(r["query_start"]), lo),
    )
    return {
        "detected": True,
        "hmm": any("hmm" in _methods(r) for r in overlapping),
        "donor_ok": best["minor_parent"] == donor,
        "breakpoint_error": (
            abs(int(best["query_start"]) - lo) + abs(int(best["query_end"]) - hi)
        ) // 2,
    }


def score_clonal(rows: list[dict]) -> tuple[int, dict[str, int]]:
    """``(false regions, per-caller counts)`` for a run that must report nothing."""
    present = _present(rows)
    per_caller: Counter[str] = Counter()
    for row in present:
        per_caller.update(_methods(row))
    return len(present), dict(per_caller)


def score_span(rows: list[dict], span: tuple[int, int]) -> dict:
    """How much longer than the true difference the reported region is.

    ``reported_bp`` is the length of the longest reported region overlapping ``span``
    (``None`` when none does); ``excess_bp`` is that length minus the span's.
    """
    lo, hi = span
    overlapping = [
        r for r in _present(rows)
        if int(r["query_start"]) < hi and int(r["query_end"]) > lo
    ]
    if not overlapping:
        return {"detected": False, "reported_bp": None, "excess_bp": None}
    reported = max(int(r["query_end"]) - int(r["query_start"]) for r in overlapping)
    return {"detected": True, "reported_bp": reported,
            "excess_bp": max(0, reported - (hi - lo))}


# --- running ---------------------------------------------------------------------

def _write_fasta(path: Path, seqs: dict[str, str]) -> None:
    with path.open("w") as fh:
        for label, seq in seqs.items():
            fh.write(f">{label}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i : i + 70] + "\n")


def _read_regions(path: Path) -> list[dict]:
    if not path.exists():
        return []
    lines = [x for x in path.read_text().splitlines() if x.strip()]
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, x.split("\t"), strict=False)) for x in lines[1:]]


def scan(
    seqs: dict[str, str], logger: logging.Logger, *, min_methods: int = 1,
    cluster_lineages: bool = True,
) -> list[dict]:
    """Run the shipped pipeline on one simulated alignment; return its region rows."""
    with tempfile.TemporaryDirectory() as td:
        msa = Path(td) / "aln.fasta"
        _write_fasta(msa, seqs)
        out = Path(td) / "out"
        run_recomb(
            RecombParams(msa=msa, output=out, query=QUERY, plot_format="png",
                         min_methods=min_methods, cluster_lineages=cluster_lineages),
            logger,
        )
        return _read_regions(out / "recombination_regions.tsv")


def _median(values: list[int]) -> str:
    return f"{int(np.median(values))} bp" if values else "n/a"


def main(argv: list[str]) -> int:
    reps, min_methods = 10, 1
    if "--reps" in argv:
        reps = int(argv[argv.index("--reps") + 1])
    if "--min-methods" in argv:
        min_methods = int(argv[argv.index("--min-methods") + 1])
    cluster = "--no-cluster-lineages" not in argv

    logger = logging.getLogger("tessera.regimes")
    logger.addHandler(logging.NullHandler())
    logger.propagate = False  # the scan is chatty; the table below is the output

    def run(seqs: dict[str, str]) -> list[dict]:
        return scan(seqs, logger, min_methods=min_methods, cluster_lineages=cluster)

    print(f"Tessera regime harness -- {reps} replicate(s), --min-methods {min_methods}, "
          f"lineage clustering {'on' if cluster else 'off'}\n")

    detected = hmm = donor_ok = 0
    errors: list[int] = []
    for rep in range(reps):
        seqs, tract = simulate_near_identical(3000 + rep, recombinant=True)
        assert tract is not None
        got = score_tract(run(seqs), tract, donor="R1")
        detected += got["detected"]
        hmm += got["hmm"]
        donor_ok += got["donor_ok"]
        if got["breakpoint_error"] is not None:
            errors.append(got["breakpoint_error"])
    print(f"near_identical, tract {NEAR_TRACT[0]}-{NEAR_TRACT[1]} from R1:")
    print(f"  detected {detected}/{reps} | called by the HMM {hmm}/{reps} | "
          f"correct donor {donor_ok}/{reps} | median breakpoint error {_median(errors)}")

    bad = regions = 0
    callers: Counter[str] = Counter()
    for rep in range(reps):
        seqs, _ = simulate_near_identical(4000 + rep, recombinant=False)
        n_false, per_caller = score_clonal(run(seqs))
        bad += bool(n_false)
        regions += n_false
        callers.update(per_caller)
    print("near_identical, clonal (every region is a false positive):")
    print(f"  runs with a false region {bad}/{reps} | false regions {regions} | "
          + (", ".join(f"{k}={v}" for k, v in sorted(callers.items())) or "--"))

    found = 0
    reported: list[int] = []
    excess: list[int] = []
    for rep in range(reps):
        seqs, span = simulate_sibling_tract(5000 + rep)
        got = score_span(run(seqs), span)
        found += got["detected"]
        if got["reported_bp"] is not None:
            reported.append(got["reported_bp"])
            excess.append(got["excess_bp"])
    print(f"sibling_tract, true difference {SIB_EXTRA} bp:")
    print(f"  detected {found}/{reps} | median reported length {_median(reported)} | "
          f"median excess {_median(excess)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
```

- [ ] **Step 4: Run the tests**

Run: `pytest tests/unit/test_regime_scoring.py -q && ruff check validation tests`
Expected: `11 passed`, `All checks passed!`

- [ ] **Step 5: Run the harness on the unchanged callers**

```bash
python validation/run_regimes.py --reps 10 | tee validation/data/audit-c/baseline/regimes_mm1.txt
python validation/run_regimes.py --reps 10 --min-methods 2 | tee validation/data/audit-c/baseline/regimes_mm2.txt
```

Expected (measured for this plan, both gates):

```
near_identical, tract 40000-48000 from R1:
  detected 10/10 | called by the HMM 0/10 | correct donor 10/10 | median breakpoint error 465 bp
near_identical, clonal (every region is a false positive):
  runs with a false region 0/10 | false regions 0 | --
sibling_tract, true difference 600 bp:
  detected 10/10 | median reported length 11700 bp | median excess 11100 bp
```

`called by the HMM 0/10` is defect C1 and `median reported length 11700 bp` is defect C2. If these two lines do not reproduce, stop and report: the defects this plan fixes are not present on your checkout.

- [ ] **Step 6: Document the harness**

In `validation/README.md`, in the layout block, add after the `run_specificity.py` line:

```
  run_regimes.py       near-identical panels and sibling tracts, simulated (no aligner needed)
```

Then add this section immediately before the `## Prerequisites` heading:

````markdown
## Regimes the other harnesses do not sample

`run_regimes.py` simulates two situations that fall outside `run_specificity.py`, whose
clades are about 16 % apart and so always run under base-pair windowing:

| regime | what it asks |
|---|---|
| `near_identical` with a tract | references about 0.2 % apart (the mpox / VZV regime, analysed under informative-site windowing): is the tract found, and is the HMM among the callers? |
| `near_identical` clonal | the same panel with nothing recombined: every region is a false positive |
| `sibling_tract` | the panel holds a sibling recombinant whose tract runs 600 bp past the query's; downstream of those 600 bp the donor and the backbone are identical, so nothing marks where the region ends |

It defaults to `--min-methods 1`, the CLI default, so it reports what a user gets; pass
`--min-methods 2` for the agreement gate and `--no-cluster-lineages` to switch lineage
clustering off. Like the specificity harness it needs no aligner, network or data.

```
python validation/run_regimes.py --reps 3      # quick look
```

The same caveat applies: a star tree under JC69 and a handful of replicates show whether a
failure mode exists, not how large it is on real panels.
````

- [ ] **Step 7: Commit**

```bash
git add validation/run_regimes.py tests/unit/test_regime_scoring.py validation/README.md
git commit -m "$(cat <<'EOF'
Add a harness for near-identical panels and sibling tracts

The post-1.2.0 audit found two defects with simulations that no harness
covers: on a panel ~0.2 % divergent the HMM never contributes a call, and
when a sibling recombinant is in the panel an HMM region runs on through
columns where donor and backbone are identical. run_specificity.py cannot
see either -- its clades are ~16 % apart, so it always runs under base-pair
windowing, and its positive control's parents differ everywhere.

run_regimes.py simulates both, plus a clonal near-identical panel so that
any fix which gives the HMM a vote there is measured against its
false-positive cost. It follows run_specificity.py: self-contained, opt-in,
with the simulation and scoring unit-tested in CI.

On main: HMM among the callers 0/10; sibling tract of 600 bp reported as
11700 bp (median, 10 replicates).

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: 3SEQ and MaxChi tract start on ties (C3)

`threeseq.max_descent` finds the trough of the walk and then takes the **first** maximum before it. If the walk reaches its maximum, dips, returns to the same height and only then descends, the reported tract starts at the first maximum and includes the balanced excursion. MaxChi scores the same tract, and its permutation null (`maxchi._exceedances`) applies the same rule to every permuted walk, so the two must change together or observed and null statistics stop being comparable.

**Files:**
- Modify: `src/tessera/recomb/threeseq.py` (`max_descent`, the `peak = ...` line)
- Modify: `src/tessera/recomb/maxchi.py` (`_exceedances`, the `peak = ...` line)
- Test: `tests/unit/test_threeseq.py`, `tests/unit/test_maxchi.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `max_descent(steps).start_site` is the last index at which the walk is at its maximum before the trough. Signature unchanged. Task 8 calls `triplet_steps` from this module; it does not depend on this change.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_threeseq.py`:

```python


def test_max_descent_starts_at_the_last_maximum_before_the_trough() -> None:
    # The walk reaches its maximum (2) at site 2, dips, returns to 2 at site 4 and only
    # then descends. The donor run is sites [4, 8); sites 2-3 are a balanced excursion.
    steps = np.array([1, 1, -1, 1, -1, -1, -1, -1])
    d = max_descent(steps)
    assert d.depth == 4
    assert (d.start_site, d.end_site) == (4, 8)
    assert set(steps[d.start_site:d.end_site]) == {-1}
```

- [ ] **Step 2: Add the consistency guard for MaxChi**

In `tests/unit/test_maxchi.py`, change the import block to include `_exceedances`:

```python
from tessera.recomb.maxchi import (
    _chi2_2x2,
    _exceedances,
    maxchi_pvalue,
    maxchi_statistic,
)
```

and append:

```python


def test_permutation_null_uses_the_same_tract_rule_as_the_observed_walk() -> None:
    """The vectorised null must score each permuted walk on the tract ``max_descent``
    would report for it; otherwise observed and null statistics are not comparable.

    A consistency guard: it holds before and after the tie rule changes, and fails if
    only one of ``threeseq.max_descent`` / ``maxchi._exceedances`` is changed.
    """
    steps = np.array([1] * 30 + [-1] * 14)
    np.random.default_rng(5).shuffle(steps)
    k, seed = 400, 11
    order = np.argsort(np.random.default_rng(seed).random((k, steps.size)), axis=1)
    stats = []
    for perm in steps[order]:
        d = max_descent(perm)
        stats.append(maxchi_statistic(perm, d.start_site, d.end_site))
    observed = float(np.median(stats))
    expected = sum(s >= observed for s in stats)
    assert 0 < expected < k  # the threshold splits the permutations
    assert _exceedances(steps, observed, k, np.random.default_rng(seed)) == expected
```

- [ ] **Step 3: Run the tests**

Run: `pytest tests/unit/test_threeseq.py tests/unit/test_maxchi.py -q`
Expected: 1 failed -- `test_max_descent_starts_at_the_last_maximum_before_the_trough` with `assert (2, 8) == (4, 8)`. The MaxChi guard passes (both sides still use the first maximum).

- [ ] **Step 4: Change the rule in `max_descent`**

In `src/tessera/recomb/threeseq.py`, replace

```python
    peak = int(np.argmax(cumulative[: trough + 1])) if trough > 0 else 0
```

with

```python
    # The tract starts at the LAST maximum before the trough. The first one (plain
    # argmax) would pull in any balanced excursion that returns to the same height,
    # reporting the tract as starting earlier than the donor run does.
    peak = trough - int(np.argmax(cumulative[trough::-1])) if trough > 0 else 0
```

- [ ] **Step 5: Run the tests again**

Run: `pytest tests/unit/test_threeseq.py tests/unit/test_maxchi.py -q`
Expected: 1 failed -- now the MaxChi guard, because the null still uses the first maximum. This is the guard doing its job.

- [ ] **Step 6: Change the rule in the MaxChi null**

In `src/tessera/recomb/maxchi.py`, in `_exceedances`, replace

```python
    # peak = where the running max was reached, on or before the trough
    masked = np.where(cols[None, :] <= trough[:, None], cum, np.iinfo(np.int64).min)
    peak = np.argmax(masked, axis=1)
```

with

```python
    # peak = the LAST place the running max was reached, on or before the trough --
    # the same rule threeseq.max_descent applies to the observed walk, so the null
    # and the observed statistic are computed on the same kind of tract.
    masked = np.where(cols[None, :] <= trough[:, None], cum, np.iinfo(np.int64).min)
    peak = length - np.argmax(masked[:, ::-1], axis=1)
```

- [ ] **Step 7: Run the tests and the fast suite**

Run: `pytest tests/unit/test_threeseq.py tests/unit/test_maxchi.py -q && pytest -m "not requires_binary" -q && ruff check src tests`
Expected: all pass.

- [ ] **Step 8: Harness gate**

Run: `validation/data/audit-c/gate.sh task3` and, if Task 1 Step 5 produced a hybrids baseline, `python validation/run_hybrids.py > validation/data/audit-c/task3/hybrids.txt 2>&1`.

Acceptance (all must hold):
- `spec_mm2.txt` and `spec_mm1.txt`: the `TOTAL` line (runs with a false region, false regions, per-caller counts) is identical to the baseline.
- Positive control: `detected` and `correct donor` identical to the baseline; median breakpoint error not higher than the baseline. (Measured at 3 replicates: 55 bp -> 50 bp.)
- `validation.txt`: the same PASS / FAIL / SKIP per dataset as the baseline.
- `diff validation/data/audit-c/baseline/validation_regions.tsv validation/data/audit-c/task3/validation_regions.tsv`: the same rows; only 3SEQ / MaxChi start coordinates may move, and only later. (Measured: one row, `hiv1_crf` donor `C`, 8433 -> 8439.)
- Hybrids: no case changes from PASS to FAIL.

The 3SEQ p-value depends only on the depth of the descent, which this does not change. The MaxChi p-value can change, because the tract's chi-square is computed on a different interval; that is why the false-region counts are part of the gate.

If the gate fails: `git checkout -- src/tessera/recomb/threeseq.py src/tessera/recomb/maxchi.py`, keep the two tests out of the commit, record the numbers, go to Task 4.

- [ ] **Step 9: Commit**

```bash
git add src/tessera/recomb/threeseq.py src/tessera/recomb/maxchi.py tests/unit/test_threeseq.py tests/unit/test_maxchi.py
git commit -m "$(cat <<'EOF'
Start a 3SEQ/MaxChi tract at the last maximum before the trough

max_descent took the first maximum of the walk before the trough. When the
walk reaches its maximum, dips, and returns to the same height before
descending, that reports the tract from the earlier point and includes a
balanced stretch that is not donor tract: for steps +1 +1 -1 +1 -1 -1 -1 -1
it gave sites 2-8 (support 5/6) where the donor run is 4-8 (4/4).

MaxChi scores the same tract and its permutation null applied the same rule
to each permuted walk, so both change together; a guard test now fails if
only one of them does.

The 3SEQ p-value is a function of the descent depth and is unchanged.
Harness: <paste the spec_mm1 / spec_mm2 TOTAL lines, positive-control line
and the validation_regions diff from validation/data/audit-c/task3/>.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

Replace the `<paste ...>` line with the measured lines before committing.

---

### Task 4: Exact 3SEQ p-values keep their magnitude (C4)

`descent_pvalue_exact` propagates the probability that the drawdown never reaches the depth and returns one minus it. Below about 1e-16 that subtraction is exactly `0.0`, which is what `recombination_regions.tsv` then holds as `pvalue` and `qvalue` for the strongest regions. Propagating the complementary probability directly keeps every term a sum of non-negative products.

**Files:**
- Modify: `src/tessera/recomb/threeseq.py` (`descent_pvalue_exact`, whole function)
- Test: `tests/unit/test_threeseq.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `descent_pvalue_exact(m: int, n: int, depth: int) -> float`, same signature and same values to within floating-point error, except that values below ~1e-16 are no longer `0.0`. Still returns exactly `0.0` when `n < depth` (the event is impossible) and `1.0` when `depth <= 0`.

- [ ] **Step 1: Write the failing test**

In `tests/unit/test_threeseq.py`, add `import math` after `import itertools`, and append:

```python


def test_descent_pvalue_exact_keeps_small_values() -> None:
    # depth == n: every down-step must be consecutive, so the count of favourable
    # arrangements is the m + 1 positions of that block among the up-steps.
    for m, n in ((200, 100), (50, 30), (12, 9)):
        expected = (m + 1) / math.comb(m + n, n)
        got = descent_pvalue_exact(m, n, n)
        assert got > 0.0
        assert math.isclose(got, expected, rel_tol=1e-9)
```

- [ ] **Step 2: Run it to see it fail**

Run: `pytest tests/unit/test_threeseq.py::test_descent_pvalue_exact_keeps_small_values -q`
Expected: FAIL with `assert 0.0 > 0.0` (the true value for `(200, 100)` is about 4.8e-80).

- [ ] **Step 3: Rewrite the recursion**

In `src/tessera/recomb/threeseq.py`, replace the whole of `descent_pvalue_exact` with:

```python
def descent_pvalue_exact(m: int, n: int, depth: int) -> float:
    """Exact ``P(max drawdown >= depth)`` over the uniform arrangements of ``m`` +1 and
    ``n`` -1 steps.

    Dynamic program over the walk's current drawdown ``delta`` in ``[0, depth)``: an
    up-step moves ``delta -> max(0, delta - 1)``, a down-step ``delta -> delta + 1``.
    ``g[delta]`` is the probability that the drawdown *reaches* ``depth`` from that
    state; a down-step from ``delta == depth - 1`` reaches it with probability 1.

    The reaching probability is propagated directly rather than as one minus the
    probability of never reaching it: every term is a sum of non-negative products, so
    a very small p-value keeps its magnitude instead of cancelling to exactly ``0.0``
    below about 1e-16. Probabilities (not path counts) are propagated, so there is no
    big-integer overflow.
    """
    if depth <= 0:
        return 1.0
    if n < depth:  # need at least ``depth`` consecutive down-steps to reach it
        return 0.0

    def after_down(row_below: np.ndarray) -> np.ndarray:
        down = np.empty(depth)
        down[: depth - 1] = row_below[1:depth]  # delta -> delta + 1
        down[depth - 1] = 1.0  # ... and from depth-1 the walk reaches ``depth``
        return down

    # prev[j] is the length-``depth`` vector g(i-1, j, .). Build row i from row i-1.
    # i == 0 row: only down-steps remain, a single forced order; nothing left -> 0.
    prev = [np.zeros(depth)]
    for j in range(1, n + 1):
        prev.append(after_down(prev[j - 1]))
    for i in range(1, m + 1):
        cur: list[np.ndarray] = [np.zeros(depth)]  # j == 0: only up-steps remain
        for j in range(1, n + 1):
            total = i + j
            up = np.empty(depth)
            up[0] = prev[j][0]
            up[1:] = prev[j][: depth - 1]  # delta -> max(0, delta - 1)
            cur.append((i / total) * up + (j / total) * after_down(cur[j - 1]))
        prev = cur
    return float(prev[n][0])
```

- [ ] **Step 4: Run the 3SEQ tests**

Run: `pytest tests/unit/test_threeseq.py -q`
Expected: all pass, including `test_descent_pvalue_exact_matches_brute_force` (every `m, n < 7` against enumeration) and `test_descent_pvalue_edges`.

- [ ] **Step 5: Confirm the region table no longer holds a zero**

```bash
tessera recomb --msa example_data/divergent.msa.fasta --query query --output /tmp/tessera_c4 \
    --window-size 300 --window-step 30 >/dev/null 2>&1
cut -f14,15,16 /tmp/tessera_c4/recombination_regions.tsv
```

Expected: a `pvalue` and `qvalue` that are small and non-zero (before this task both read `0.0`), with `test` = `3SEQ max-descent (exact)`.

- [ ] **Step 6: Harness gate**

Run: `validation/data/audit-c/gate.sh task4`

Acceptance: every line of the gate identical to `task3` (or to `baseline`, if Task 3 was reverted). The value of an exact p-value changes only below 1e-16; no comparison against `--alpha` can flip. `diff validation/data/audit-c/task3/validation_regions.tsv validation/data/audit-c/task4/validation_regions.tsv` must be empty.

If it is not empty, something other than underflow changed: revert and record.

- [ ] **Step 7: Commit**

```bash
git add src/tessera/recomb/threeseq.py tests/unit/test_threeseq.py
git commit -m "$(cat <<'EOF'
Keep the magnitude of very small exact 3SEQ p-values

descent_pvalue_exact propagated the probability of never reaching the
descent depth and returned one minus it. Below about 1e-16 that subtraction
is exactly 0.0, so the strongest regions were written to
recombination_regions.tsv with pvalue = 0.0 and qvalue = 0.0 -- a value no
permutation or exact test can produce.

The recursion now propagates the probability of reaching the depth. Every
term is a sum of non-negative products, so nothing cancels: m=200, n=100,
depth=100 gives 4.8e-80, equal to the closed form 201 / C(300, 100). Values
above 1e-16 are unchanged (the brute-force test over all m, n < 7 still
passes), so no call can change.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: Ensemble merge edge cases (C7)

Two defects in `recomb/ensemble.py`, both shown at unit level only; neither was reached end to end in the audit, so this task has unit tests and a regression gate, not a detection gate.

1. `_merge_group` relabels every region's major parent to the canonical backbone. When the callers disagree on the backbone, a site caller's region can name the canonical backbone as its *donor*, and after relabelling the row reads "H donated into H".
2. `_group` sorts by `(minor_parent, query_start)` and puts each region into the first group it overlaps. On a typed panel two labels of one lineage arrive out of coordinate order, so a region can overlap two existing groups; it joins one, and the same donor lineage is reported as two overlapping regions.

**Files:**
- Modify: `src/tessera/recomb/ensemble.py` (`_merge_group`, `_group`)
- Test: `tests/unit/test_ensemble.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `consensus_regions(...)` unchanged in signature. A merged `Region` never has `minor_parent == major_parent`; regions of one donor (genome, or lineage when typed) that overlap transitively are one region.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_ensemble.py`:

```python


def test_region_is_not_relabelled_as_a_donor_into_itself() -> None:
    # The callers disagree on the backbone: the HMM reports H (it set the sibling S
    # aside), 3SEQ kept S and found an H tract in it. H is the canonical backbone, so
    # relabelling the 3SEQ region's major to H would read "H donated into H".
    hmm = mk("X", 100, 200, "hmm", major="H", qvalue=1e-6)
    seq = mk("H", 500, 600, "3seq", major="S", qvalue=1e-9)
    merged, _ = consensus_regions({"hmm": [hmm], "3seq": [seq]}, major="H")
    assert all(r.minor_parent != r.major_parent for r in merged)
    by_minor = {r.minor_parent: r.major_parent for r in merged}
    assert by_minor == {"X": "H", "H": "S"}  # the tested pair is kept for the odd one out


def test_same_lineage_regions_merge_transitively_across_labels() -> None:
    # A1 and A2 are one lineage. A1 has two separate regions; the A2 region overlaps
    # both. All three are one event; sorting by label must not leave two groups.
    lmap = {"A1": "L", "A2": "L"}
    first = mk("A1", 0, 10, "hmm", qvalue=1e-6)
    second = mk("A1", 20, 30, "hmm", qvalue=1e-6)
    bridge = mk("A2", 5, 25, "3seq", qvalue=1e-9)
    merged, breakdown = consensus_regions(
        {"hmm": [first, second], "3seq": [bridge]}, major="M", lineage_map=lmap
    )
    assert len(merged) == 1
    assert (merged[0].query_start, merged[0].query_end) == (0, 30)
    assert merged[0].methods == ("hmm", "3seq")
    assert len(breakdown) == 1
```

- [ ] **Step 2: Run them to see them fail**

Run: `pytest tests/unit/test_ensemble.py -q`
Expected: 2 failed. The first fails on `assert all(r.minor_parent != r.major_parent ...)` (the 3SEQ region is relabelled `H` -> `H`); the second on `assert 2 == 1`.

- [ ] **Step 3: Keep the tested pair when the canonical backbone is the donor**

In `src/tessera/recomb/ensemble.py`, in `_merge_group`, replace

```python
    return Region(
        minor_parent=best.minor_parent,
        major_parent=major or best.major_parent,
```

with

```python
    # The canonical backbone replaces each caller's own, so the table reads against one
    # backbone -- except when it IS this region's donor. That happens when the callers
    # disagree on the backbone (the HMM set a sibling aside, a site caller did not): the
    # region was called against the other backbone, and relabelling it would report a
    # genome as a donor into itself. Keep the pair the caller actually tested.
    relabel = bool(major) and major != best.minor_parent
    return Region(
        minor_parent=best.minor_parent,
        major_parent=major if relabel and major else best.major_parent,
```

- [ ] **Step 4: Fold every group a region touches into one**

In the same file, in `_group`, replace

```python
    for region in sorted(regions, key=lambda r: (r.minor_parent, r.query_start)):
        placed = False
        for group in groups:
            if _same_donor(group[0], region, lineage_map) and any(
                _overlap(region, member) for member in group
            ):
                group.append(region)
                placed = True
                break
        if not placed:
            groups.append([region])
    return groups
```

with

```python
    for region in sorted(regions, key=lambda r: (r.minor_parent, r.query_start)):
        hits = [
            group for group in groups
            if _same_donor(group[0], region, lineage_map)
            and any(_overlap(region, member) for member in group)
        ]
        if not hits:
            groups.append([region])
            continue
        # A region can bridge groups that did not overlap each other -- on a typed panel
        # the sort is by genome label, so two labels of one lineage arrive out of
        # coordinate order. Fold every group it touches into the first, or the same
        # donor would be reported as two overlapping regions.
        target = hits[0]
        target.append(region)
        for other in hits[1:]:
            target.extend(other)
        groups = [g for g in groups if all(g is not other for other in hits[1:])]
    return groups
```

The groups are compared by identity (`is not`), not equality: two groups can hold equal regions.

- [ ] **Step 5: Run the tests**

Run: `pytest tests/unit/test_ensemble.py -q && pytest -m "not requires_binary" -q && ruff check src tests && mypy src`
Expected: all pass, no new mypy error.

- [ ] **Step 6: Regression gate**

Run: `validation/data/audit-c/gate.sh task5`

Acceptance: identical to the previous task's output in every file, including an empty `diff` of `validation_regions.tsv`. Neither edge case occurs in these harnesses (untyped panels; callers agree on the backbone), so any difference is a regression in the ordinary path: revert and record.

- [ ] **Step 7: Commit**

```bash
git add src/tessera/recomb/ensemble.py tests/unit/test_ensemble.py
git commit -m "$(cat <<'EOF'
Close two edge cases in the ensemble merge

The merge relabels every region's major parent to the canonical backbone.
When the callers disagree on the backbone -- the HMM set a sibling aside, a
site caller kept it -- a site caller's region can name the canonical
backbone as its donor, and the relabelled row then read as a genome donating
into itself. Such a region now keeps the pair its caller tested.

Grouping put each region into the first group it overlapped. On a typed
panel the sort is by genome label, so a region of one label can overlap two
groups of another label of the same lineage; it joined one and the lineage
was reported as two overlapping regions. A region now folds every group it
touches into one.

Both were found by reading and reproduced with constructed regions only;
neither was reached end to end. The harness outputs are unchanged.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Read pool clade labels as written (C6)

Tessera writes pool genomes itself with the defline `>{accession} {label}` (`discover/nextclade.py`, `_write_genome`; `discover/pool.py`, the NCBI Virus split). `_select_from` and `_type_panel` then pass that defline to `genotype_from_title`, a free-text miner that wants a digit-bearing token of at least four characters. Short and multi-word clades are lost, and with them lineage-aware selection and the recombinant-lineage exclusion for those genomes.

Measured on the 33 pools cached on the audit machine (genomes carrying a label, typed correctly by mining, typed when the label is read as written):

| pool | labelled | mined | read |
|---|---|---|---|
| measles | 846 | 0 | 846 |
| mpox (all clades) | 1556 | 2 | 1556 |
| VZV | 256 | 0 | 256 |
| HIV-1 (hxb2) | 1052 | 639 | 1052 |
| RSV-A | 1677 | 1170 | 1677 |
| SARS-CoV-2 (XBB dataset) | 2803 | 2672 (421 recombinant-named) | 2803 (549 recombinant-named) |

And the effect on selection from three of those pools (one tip held out as the query, local skani only): measles 1 genome selected -> 4 (4 clades); HIV-1 9 genomes in 5 clades -> 9 genomes in 9 clades; mpox 5 -> 8.

**This changes which references a panel contains for most pathogens.** `run_hybrids.py` does not exercise `_select_from`: it types its pools from the Nextclade tree and calls `select_regional` directly, which is the behaviour this task gives the product. So the hybrids harness has been validating lineage-typed selection that the CLI did not perform, and it cannot serve as this task's gate. The gate below is an end-to-end panel build.

**Files:**
- Modify: `src/tessera/recomb/typing.py` (new `label_from_pool_header`, `pool_labels`; `build_lineage_map` gains `pool_labels`)
- Modify: `src/tessera/discover/iterate.py` (`_seed_from_pool`, `_select_from`, `_type_panel`; new `_record_pool_labels`, `_read_pool_labels`, `POOL_LABELS_TSV`)
- Modify: `docs/reference-panels.md`
- Test: `tests/unit/test_typing.py`, `tests/unit/test_select_from_lineage.py`, `tests/unit/test_iterate.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces:
  - `typing.label_from_pool_header(title: str) -> str | None`
  - `typing.pool_labels(files: Iterable[Path]) -> dict[str, str]` (file label -> clade)
  - `typing.build_lineage_map(*, user_tsv=None, datasets_rows=None, title_by_label=None, organism=None, pool_labels: Mapping[str, str] | None = None)`; rows from pool labels carry the source string `"pool-header"`
  - `iterate._select_from(params, genomes, logger, *, pool_headers: bool = False)`
  - `iterate.POOL_LABELS_TSV = "pool_labels.tsv"`; `<output>/pool_labels.tsv` is `label<TAB>clade`, one row per pool-seeded genome

- [ ] **Step 1: Write the failing tests for the typing helpers**

In `tests/unit/test_typing.py`, replace the import block with:

```python
from tessera.recomb.typing import (
    build_lineage_map,
    dominant_lineage_token,
    genotype_from_title,
    label_from_pool_header,
    lineage_map_from_rows,
    lineage_of,
    load_lineage_map,
    pool_labels,
    read_lineage_rows,
    titles_from_collection,
    typed,
    write_lineage_map,
)
```

and append:

```python


# --- structured pool labels -------------------------------------------------

def test_label_from_pool_header_takes_the_label_verbatim():
    # the short and multi-word clades that title mining drops
    assert label_from_pool_header("MN908947.3 B.1") == "B.1"
    assert label_from_pool_header("OX1 XBB") == "XBB"
    assert label_from_pool_header("K03455 B") == "B"
    assert label_from_pool_header("MV1 D8") == "D8"
    assert label_from_pool_header("NC_063383 IIb") == "IIb"
    assert label_from_pool_header("VZV1 clade 2 wild-type") == "clade 2 wild-type"
    assert genotype_from_title("MN908947.3 B.1") is None  # why mining is not enough


def test_label_from_pool_header_none_for_unlabelled_genomes():
    for header in ("ACC1", "ACC1 ", "ACC1 NA", "ACC1 example", ""):
        assert label_from_pool_header(header) is None


def test_pool_labels_reads_each_file_header(tmp_path: Path):
    (tmp_path / "ACC1.fasta").write_text(">ACC1 B\nACGT\n")
    (tmp_path / "ACC2.fasta").write_text(">ACC2 NA\nACGT\n")
    (tmp_path / "ACC3.2.fasta").write_text(">ACC3.2 clade 5\nACGT\n")
    files = sorted(tmp_path.iterdir())
    assert pool_labels(files) == {"ACC1": "B", "ACC3.2": "clade 5"}


def test_build_lineage_map_pool_label_beats_title_and_yields_to_datasets_and_user(
    tmp_path: Path,
):
    user = tmp_path / "user.tsv"
    user.write_text("U1\tUSER\n")
    rows = build_lineage_map(
        user_tsv=user,
        datasets_rows=[("D1", "DATASETS")],
        title_by_label={"P1": "P1 CRF01_AE", "T1": "T1 Norovirus GII.4 isolate x"},
        pool_labels={"P1": "A1", "D1": "pool", "U1": "pool"},
    )
    assert dict((label, (g, src)) for label, g, src in rows) == {
        "P1": ("A1", "pool-header"),      # not the token mined from its title
        "T1": ("GII.4", "title"),         # no pool label: mined as before
        "D1": ("DATASETS", "ncbi-datasets"),
        "U1": ("USER", "user"),
    }
```

- [ ] **Step 2: Write the failing tests for selection**

Append to `tests/unit/test_select_from_lineage.py`:

```python


def _captured_lineage_map(monkeypatch, tmp_path, *, pool_headers: bool):
    pool = tmp_path / "pool"
    pool.mkdir()
    genomes = [_write(pool, "ACC1", "B"), _write(pool, "ACC2", "XBB"),
               _write(pool, "ACC3", "NA")]
    captured = {}

    def fake_select_regional(query, genomes, **kwargs):
        captured.update(kwargs)
        return PoolSelection(selected=list(genomes), table=[])

    monkeypatch.setattr("tessera.discover.pool.select_regional", fake_select_regional)
    params = FillParams(query=tmp_path / "q.fasta", collection=None, output=tmp_path / "out")
    it._select_from(params, genomes, _LOG, pool_headers=pool_headers)
    return captured["lineage_of"]


def test_select_from_reads_structured_pool_labels(monkeypatch, tmp_path):
    # Short clade labels as a Nextclade pool writes them. Mined as free text they are all
    # dropped, so lineage selection and the recombinant exclusion never see them.
    lineage_of = _captured_lineage_map(monkeypatch, tmp_path, pool_headers=True)
    assert lineage_of == {"ACC1": "B", "ACC2": "XBB"}  # "NA" is not a clade


def test_select_from_still_mines_a_user_pool(monkeypatch, tmp_path):
    # A --candidate-pool holds the user's own files; their headers are free text and a
    # trailing word must not be taken for a clade.
    assert _captured_lineage_map(monkeypatch, tmp_path, pool_headers=False) is None
```

- [ ] **Step 3: Extend the existing seeding test**

In `tests/unit/test_iterate.py`, in `test_seed_source_nextclade_routes_through_pool_selection`, replace

```python
    def fake_select(params, genomes, logger):
        from tessera.discover.pool import PoolSelection
        return PoolSelection(selected=list(genomes))
```

with

```python
    def fake_select(params, genomes, logger, *, pool_headers=False):
        from tessera.discover.pool import PoolSelection
        captured["pool_headers"] = pool_headers
        return PoolSelection(selected=list(genomes))
```

and after the line `assert captured["dataset"] == "nextstrain/sars-cov-2/XBB"` add:

```python
    # a Nextclade pool's headers hold structured clade labels; they are read as such
    assert captured["pool_headers"] is True
    # ... and carried to the panel's lineage table, although "A1" is too short to mine
    assert (out / "pool_labels.tsv").read_text() == "REF1\tA1\n"
    assert "REF1\tA1\tpool-header" in (out / "lineages.tsv").read_text()
```

(The test's pool genome is written as `>REF1 A1`.)

- [ ] **Step 4: Run them to see them fail**

Run: `pytest tests/unit/test_typing.py tests/unit/test_select_from_lineage.py tests/unit/test_iterate.py -q`
Expected: `test_typing.py` fails at import (`cannot import name 'label_from_pool_header'`); the two new selection tests fail with `TypeError: _select_from() got an unexpected keyword argument 'pool_headers'`; the seeding test fails on `assert captured["pool_headers"] is True` (the flag is never passed).

- [ ] **Step 5: Add the typing helpers**

In `src/tessera/recomb/typing.py`, insert immediately before `def first_header(`:

```python
# What the pool builders write in the label slot when a genome carries no clade: an
# empty note, Nextclade's "NA", and the marker given to a dataset's example sequences.
_NO_POOL_LABEL = frozenset({"", "na", "example"})


def label_from_pool_header(title: str) -> str | None:
    """The clade/lineage label of a pool genome's defline, verbatim.

    Tessera writes pool genomes itself -- Nextclade tree tips and NCBI Virus genomes --
    with the defline ``{accession} {label}``, so the label is everything after the
    accession and needs no mining. :func:`genotype_from_title` is for free-text GenBank
    titles: it wants a digit-bearing token of four characters or more, which drops
    ``B.1``, ``XBB``, ``D8``, ``IIb``, the HIV pure subtypes and any multi-word clade
    (``clade 2 wild-type``). Returns ``None`` for a genome that carries no label.
    """
    parts = title.split(None, 1)
    label = parts[1].strip() if len(parts) == 2 else ""
    return None if label.lower() in _NO_POOL_LABEL else label


def pool_labels(files: Iterable[Path]) -> dict[str, str]:
    """``label -> clade`` for pool genome files, read from their deflines."""
    out: dict[str, str] = {}
    for path in files:
        label = label_from_pool_header(first_header(path))
        if label:
            out[strip_sequence_extension(path.name)] = label
    return out


```

- [ ] **Step 6: Give `build_lineage_map` a pool-label tier**

In the same file, replace the signature, docstring and first loop of `build_lineage_map`:

```python
    title_by_label: Mapping[str, str] | None = None,
    organism: str | None = None,
) -> list[tuple[str, str, str]]:
    """Merge typed names from all sources by priority into sorted sidecar rows.

    Priority (highest wins): user map > NCBI datasets lineage > mined title token.
    Keys are normalized with ``strip_sequence_extension`` so they match MSA leaves.
    """
    merged: dict[str, tuple[str, str]] = {}  # label -> (genotype, source)
    # Source 3 (lowest): a token mined from each reference's title/header note.
    for label, title in (title_by_label or {}).items():
        genotype = genotype_from_title(title, organism)
        if genotype:
            merged[strip_sequence_extension(label)] = (genotype, "title")
```

with

```python
    title_by_label: Mapping[str, str] | None = None,
    organism: str | None = None,
    pool_labels: Mapping[str, str] | None = None,
) -> list[tuple[str, str, str]]:
    """Merge typed names from all sources by priority into sorted sidecar rows.

    Priority (highest wins): user map > NCBI datasets lineage > pool label > mined
    title token. ``pool_labels`` are the structured labels of genomes that came from a
    pool Tessera built (see :func:`label_from_pool_header`); a genome that has one is
    never typed from a mined token. Keys are normalized with
    ``strip_sequence_extension`` so they match MSA leaves.
    """
    merged: dict[str, tuple[str, str]] = {}  # label -> (genotype, source)
    # Source 4 (lowest): a token mined from each reference's title/header note.
    for label, title in (title_by_label or {}).items():
        genotype = genotype_from_title(title, organism)
        if genotype:
            merged[strip_sequence_extension(label)] = (genotype, "title")
    # Source 3: the label a pool builder wrote into the header, taken verbatim.
    for label, clade in (pool_labels or {}).items():
        if clade:
            merged[strip_sequence_extension(label)] = (clade, "pool-header")
```

The parameter shadows the module-level function `pool_labels` inside `build_lineage_map` only; the function body does not call it.

- [ ] **Step 7: Pass pool labels through recruitment**

In `src/tessera/discover/iterate.py`:

(a) add `pool_labels,` to the `from ..recomb.typing import (...)` block, between `organism_from_title,` and `titles_from_collection,`;

(b) after the `from .run import (...)` block and before `@dataclass class FillParams`, add:

```python
POOL_LABELS_TSV = "pool_labels.tsv"  # sidecar: label<TAB>clade of pool-seeded genomes
```

(c) insert immediately before `def _type_panel(`:

```python
def _record_pool_labels(output: Path, genomes: list[Path]) -> None:
    """Add the pool labels of ``genomes`` to ``<output>/pool_labels.tsv``.

    Seeded genomes are copied into the collection next to downloaded ones whose headers
    are GenBank titles; once there, nothing says which header holds a structured label.
    The sidecar carries that across to :func:`_type_panel`.
    """
    merged = _read_pool_labels(output)
    merged.update(pool_labels(genomes))
    if not merged:
        return
    output.mkdir(parents=True, exist_ok=True)
    with open(output / POOL_LABELS_TSV, "w") as fo:
        for label, clade in sorted(merged.items()):
            fo.write(f"{label}\t{clade}\n")


def _read_pool_labels(output: Path) -> dict[str, str]:
    """``label -> clade`` from ``<output>/pool_labels.tsv`` (empty when absent)."""
    path = output / POOL_LABELS_TSV
    if not path.exists():
        return {}
    labels: dict[str, str] = {}
    for line in path.read_text().splitlines():
        label, _, clade = line.partition("\t")
        if label and clade:
            labels[label] = clade
    return labels


```

(d) in `_type_panel`, replace

```python
        lineage_rows = build_lineage_map(
            user_tsv=params.lineage_map,
            title_by_label=titles_from_collection(coll_files),
            organism=params.taxon,
        )
```

with

```python
        panel_labels = {strip_sequence_extension(p.name) for p in coll_files}
        lineage_rows = build_lineage_map(
            user_tsv=params.lineage_map,
            title_by_label=titles_from_collection(coll_files),
            organism=params.taxon,
            pool_labels={
                label: clade
                for label, clade in _read_pool_labels(params.output).items()
                if label in panel_labels
            },
        )
```

(e) in `_seed_from_pool`, replace the last line

```python
    _copy_into(_select_from(params, genomes, logger).selected, collection, logger)
```

with

```python
    # A pool Tessera built carries each genome's clade in its header; a user's own
    # --candidate-pool carries whatever titles its files happen to have.
    structured = force_ncbi or params.seed_source != "local"
    selected = _select_from(params, genomes, logger, pool_headers=structured).selected
    if structured:
        _record_pool_labels(params.output, selected)
    _copy_into(selected, collection, logger)
```

(f) replace the head of `_select_from`

```python
def _select_from(params: FillParams, genomes: list[Path], logger: logging.Logger):
    from .pool import select_regional

    # Type the pool from its headers (user map > NCBI datasets lineage > mined title
    # token) so the panel is reduced by lineage rather than clade-blind ANI. An empty
    # map (untyped pool) yields the pre-lineage global behaviour.
    rows = build_lineage_map(
        user_tsv=params.lineage_map,
        title_by_label=titles_from_collection(genomes),
        organism=params.taxon,
    )
```

with

```python
def _select_from(
    params: FillParams, genomes: list[Path], logger: logging.Logger,
    *, pool_headers: bool = False,
):
    from .pool import select_regional

    # Type the pool from its headers (user map > pool label > mined title token) so the
    # panel is reduced by lineage rather than clade-blind ANI. With ``pool_headers`` the
    # header's label is read as written; mining it as free text loses short clades
    # (B.1, XBB, D8, IIb, the HIV pure subtypes), which then skip lineage selection and
    # the recombinant-lineage exclusion. An empty map (untyped pool) yields the
    # pre-lineage global behaviour.
    rows = build_lineage_map(
        user_tsv=params.lineage_map,
        title_by_label=titles_from_collection(genomes),
        organism=params.taxon,
        pool_labels=pool_labels(genomes) if pool_headers else None,
    )
```

The `--deep-typing` branch of `_type_panel` is left as it is: there, genomes the miner leaves untyped are typed against the Nextclade tips by nearest reference.

- [ ] **Step 8: Run the tests**

Run: `pytest tests/unit/test_typing.py tests/unit/test_select_from_lineage.py tests/unit/test_iterate.py tests/unit/test_pool.py tests/unit/test_lineage_select.py -q && ruff check src tests && mypy src`
Expected: all pass.

- [ ] **Step 9: Document it**

In `docs/reference-panels.md`, in the `nextclade` bullet, replace

```
dataset reference plus its mutations and labelled by clade, so the report names
parents by clade. Fetched pools are cached per dataset version.
```

with

```
dataset reference plus its mutations and labelled by clade, so the report names
parents by clade. The label is read from the pool as written (short clades such as
`B.1`, `D8` or `IIb` included), drives lineage-aware selection, and is recorded for the
seeded genomes in `<output>/pool_labels.tsv`. Fetched pools are cached per dataset
version.
```

(If the line breaks in the file differ, keep the wording and re-wrap to 100 columns.)

- [ ] **Step 10: Offline measurement of the selection change**

This needs a populated Nextclade cache (`~/.cache/tessera/nextclade`) and `skani`; skip it and say so if either is absent.

```bash
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
python - <<'EOF'
import logging, tempfile
from collections import Counter
from pathlib import Path
from tessera.core.cache import cached_genomes
from tessera.core.io import strip_sequence_extension
from tessera.discover import iterate as it
from tessera.discover.iterate import FillParams
from tessera.recomb.typing import pool_labels

log = logging.getLogger("sel"); log.addHandler(logging.NullHandler()); log.propagate = False
root = Path.home() / ".cache/tessera/nextclade"
for pool_dir in sorted(p for p in root.iterdir() if p.is_dir() and not p.name.endswith("_consensus")):
    genomes = cached_genomes(pool_dir)
    if len(genomes) < 20:
        continue
    truth = pool_labels(genomes)
    query = genomes[len(genomes) // 2]          # hold one tip out as the query
    pool = [g for g in genomes if g != query]
    for mode in (False, True):
        with tempfile.TemporaryDirectory() as td:
            params = FillParams(query=query, collection=None, output=Path(td))
            try:
                selected = it._select_from(params, pool, log, pool_headers=mode).selected
            except Exception as exc:  # skani rejects very short gene datasets
                print(pool_dir.name[:44], mode, f"skipped: {type(exc).__name__}", sep="\t")
                continue
        clades = Counter(truth.get(strip_sequence_extension(g.name), "-") for g in selected)
        print(pool_dir.name[:44], f"pool_headers={mode}", f"selected={len(selected)}",
              f"clades={len(clades)}", sep="\t")
EOF
```

Save the output as `validation/data/audit-c/task6/selection.txt`. Expected direction: with `pool_headers=True` the number of clades represented in the selection is equal or higher for every pool. Measured for this plan: measles 1 -> 4, HIV-1 5 -> 9, mpox 5 -> 8 clades.

- [ ] **Step 11: End-to-end gate (needs network, mafft, skani, skDER, nextclade data)**

For measles, mpox and HIV-1, take one genome from the cached pool as a stand-in query and run the shipped command on `main` and on this branch:

```bash
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
for ds in nextstrain/measles/genome/WHO-2012 nextstrain/mpox/all-clades community/neherlab/hiv-1/hxb2; do
  name=$(echo "$ds" | tr '/-' '__')
  query=$(ls ~/.cache/tessera/nextclade/${name}_*/*.fasta | head -200 | tail -1)
  mkdir -p validation/data/audit-c/task6
  tessera detect --query "$query" --nextclade --nextclade-dataset "$ds" \
      --output "validation/data/audit-c/task6/$name" 2>&1 | tee "validation/data/audit-c/task6/$name.log" \
      | grep -E "Lineage selection|Seeding|Caller\(s\)"
done
```

Run the same loop against `main` first, writing to `validation/data/audit-c/baseline-task6/` instead: `git stash`, run the loop with that output directory, `git stash pop`. (The package is installed editable, so the checked-out source is what runs.)

Acceptance:
- Each query is an unrecombined tree tip. This branch must not report a present-donor region (a row with `donor_absent = no` in `recombination_regions.tsv`) that `main` does not report. A new region here is a false positive introduced by the panel change.
- The `Lineage selection:` log line shows at least as many lineages represented as on `main`.
- `lineages.tsv` names the seeded references by clade with source `pool-header`.

If a query gains a region: revert the commit, record the three logs, and report it. Do not adjust `select_regional`.

If there is no network, this step cannot run. Record "Task 6 end-to-end gate: not run" in the PR and leave the commit in place only if the maintainer agrees; the unit tests and Step 10 alone do not validate a panel-composition change.

- [ ] **Step 12: Commit**

```bash
git add src/tessera/recomb/typing.py src/tessera/discover/iterate.py docs/reference-panels.md \
    tests/unit/test_typing.py tests/unit/test_select_from_lineage.py tests/unit/test_iterate.py
git commit -m "$(cat <<'EOF'
Read pool clade labels as written instead of mining them

Tessera writes pool genomes with the defline "accession label", then passed
that defline to genotype_from_title, which is built for free-text GenBank
titles and wants a digit-bearing token of four or more characters. Short and
multi-word clades did not survive: on cached pools every measles tip (846),
all but two mpox tips (1556) and 413 of 1052 HIV-1 tips -- all the pure
subtypes -- came back untyped, and 128 recombinant-named SARS-CoV-2 tips
were never excluded because they were never typed.

For those genomes lineage-aware selection fell back to ANI dereplication.
run_hybrids.py did not see it: it types its pools from the tree and calls
select_regional itself, so the harness was validating a selection the CLI
did not perform.

The label is now read verbatim for pools Tessera built (Nextclade, NCBI
Virus) and carried to the panel's lineage table through a pool_labels.tsv
sidecar. A user's --candidate-pool is still mined: its headers are free
text.

This changes panel composition. Selection from cached pools:
<paste validation/data/audit-c/task6/selection.txt summary>.
End-to-end: <paste the three before/after lines, or "not run: no network">.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

---

### Task 7: Cap the PHI window (C5)

The PHI window is counted in informative-site ranks and defaults to 100. PHI compares the incompatibility of nearby site pairs with a random reordering; once the window spans most of the informative sites nearly every pair is "nearby" and the test cannot reject. On `example_data/divergent.msa.fasta` (159 informative sites; four callers agree on the tract) the default gives p = 1.0, and a window of 20 gives p = 0.001.

Six rules were compared while writing this plan (significant / testable; "testable" excludes alignments where the window covers every pair):

| rule | benchmark power (rec > 0) | benchmark specificity (rec = 0) | clonal simulations, false | `divergent` example | `cryptic_insert` example | near-identical recombinant |
|---|---|---|---|---|---|---|
| 100 (current) | 8/73 | 16/16 | 0/20 | 1.0 | not testable | 0/5 |
| `sites // 2` | 8/185 | 36/36 | 0/20 | 1.0 | 1.0 | 0/5 |
| `sites // 4` | 11/185 | 36/36 | 0/20 | 0.001 | 0.81 | 0/5 |
| **`sites // 10`** | **13/185** | **36/36** | **0/20** | **0.001** | **0.026** | **4/5** |
| fixed 20 | 21/150 | 29/30 | 3/20 | 0.001 | 1.0 | 2/5 |
| fixed 10 | 16/164 | 30/30 | 2/20 | 0.001 | 0.48 | 4/5 |

The candidate is `min(--phi-window, sites // 10)`: it has no false positive on any clonal set, it is the only rule that recovers both shipped examples, and it leaves the window at 100 for any alignment with 1000 or more informative sites, so divergent panels are unaffected.

**A finding outside this plan, to be reported and not fixed here.** PHI power on the Jaya 2023 benchmark is low under every rule (at most 14 %). The benchmark's most divergent class (`mut = 0.1`) yields zero biallelic informative columns, and `mut = 0.01` about 65: `biallelic_columns` discards every column with a third allele, which at high divergence is most of them. The window is not what limits PHI there. Record this in the PR under "Found, not fixed".

**Interaction with plan B.** With this cap the effective window is always below `sites - 1`, so the "PHI not testable (window covers every pair)" case that plan B item B7 words can no longer occur through the default path. Plan B's wording stays correct but becomes unreachable; tell whoever executes plan B.

**Files:**
- Modify: `src/tessera/recomb/diagnostics.py` (new `effective_phi_window`, `_WINDOW_FRACTION`; `recombination_signal`)
- Modify: `src/tessera/cli/cmd_recomb.py` (`--phi-window` help)
- Modify: `docs/detection-methods.md`
- Test: `tests/unit/test_diagnostics.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `diagnostics.effective_phi_window(window: int, n_informative: int) -> int`. `RecombinationSignal.phi_window` holds the window actually used.

- [ ] **Step 1: Write the failing tests**

In `tests/unit/test_diagnostics.py`, add `effective_phi_window,` to the `from tessera.recomb.diagnostics import (...)` block (after `corroborating_intervals,`), add the line `from tessera.recomb.similarity import _read_alignment` after that block, and insert before the `# --- corroborating intervals` comment:

```python
def test_effective_phi_window_is_capped_at_a_tenth_of_the_sites() -> None:
    assert effective_phi_window(100, 159) == 15
    assert effective_phi_window(100, 54) == 5
    assert effective_phi_window(100, 10) == 1  # never below one rank
    # unchanged once the alignment holds ten windows' worth of sites
    assert effective_phi_window(100, 1000) == 100
    assert effective_phi_window(100, 5000) == 100
    assert effective_phi_window(8, 5000) == 8  # a smaller request is honoured


def test_signal_at_the_minimum_number_of_informative_sites() -> None:
    # Ten informative sites is the floor for reporting a signal at all; the window is
    # then a single rank and the test must still return a valid p-value.
    rows = _block_alignment(recombinant=True, per_block=5)
    signal = recombination_signal(rows, "s0", lambda c: c)
    assert signal is not None and signal.n_informative == 10
    assert signal.phi_window == 1
    assert 0.0 < signal.phi_p <= 1.0


def test_default_window_has_power_on_the_shipped_example(example_data) -> None:
    # 159 informative sites: a 100-rank window spans most pairs, and PHI was 1.0 on an
    # alignment whose recombinant tract four callers agree on.
    rows = _read_alignment(str(example_data / "divergent.msa.fasta"))
    signal = recombination_signal(rows, "query", lambda c: c)
    assert signal is not None
    assert signal.n_informative == 159
    assert signal.phi_window == 15  # the window used, not the one asked for
    assert signal.phi_p < 0.05


```

- [ ] **Step 2: Run them to see them fail**

Run: `pytest tests/unit/test_diagnostics.py -q`
Expected: collection error, `ImportError: cannot import name 'effective_phi_window'`.

- [ ] **Step 3: Implement the cap**

In `src/tessera/recomb/diagnostics.py`:

(a) after `_PERMUTATIONS = 1000` add:

```python
# The PHI window is capped at 1/_WINDOW_FRACTION of the informative sites; see
# effective_phi_window.
_WINDOW_FRACTION = 10
```

(b) insert immediately before `def recombination_signal(`:

```python
def effective_phi_window(window: int, n_informative: int) -> int:
    """The PHI window actually used: ``window``, capped at a tenth of the sites.

    PHI compares the incompatibility of *nearby* site pairs with that of a random
    reordering. "Nearby" only means something while the window is small against the
    number of informative sites: once it spans most of them nearly every pair is
    "nearby", the statistic barely changes under permutation, and the test cannot
    reject whatever the data hold (with ``window >= n - 1`` it is exactly invariant and
    p is always 1). Capping at ``n // 10`` keeps at least ten windows' worth of sites,
    and leaves the window untouched on any alignment with ``10 * window`` sites or more.
    """
    return max(1, min(int(window), n_informative // _WINDOW_FRACTION))


```

(c) in `recombination_signal`, replace

```python
    incompatible = incompatibility_matrix(allele1, allele0)
    p, observed = phi_pvalue(incompatible, window, seed=seed)
```

with

```python
    incompatible = incompatibility_matrix(allele1, allele0)
    # `window` is an upper bound. The reported `phi_window` is the one actually used,
    # so the report never states a window the test did not run with.
    window = effective_phi_window(window, int(z))
    p, observed = phi_pvalue(incompatible, window, seed=seed)
```

The rest of the function already passes `window` to `phi_profile` and stores it as `phi_window`.

- [ ] **Step 4: Run the tests**

Run: `pytest tests/unit/test_diagnostics.py -q`
Expected: 12 passed (the nine existing tests use a 24-site alignment with `window=8`, which the cap turns into 2; they still pass).

- [ ] **Step 5: Say so in the CLI and the docs**

In `src/tessera/cli/cmd_recomb.py`, replace

```python
        help="PHI test window width, in informative-site ranks.",
```

with

```python
        help="PHI test window width, in informative-site ranks. An upper bound: the "
        "window is capped at a tenth of the informative sites, and the report states "
        "the one used.",
```

In `docs/detection-methods.md`, replace the sentence

```
parent-free outputs. The diagnostic runs for every `--method`; disable with
`--no-phi`, or widen its window with `--phi-window`.
```

with

```
parent-free outputs. The diagnostic runs for every `--method`; disable with
`--no-phi`.

The PHI window (`--phi-window`, default 100) is counted in informative-site ranks and
is an **upper bound**: the window used is capped at a tenth of the informative sites.
PHI asks whether *nearby* sites are more compatible than a random reordering, and
"nearby" stops meaning anything once the window spans most of the sites -- at 159
informative sites a 100-rank window gave p = 1 on an alignment whose tract every caller
found, where a 15-rank window gives p = 0.001. Alignments with a thousand informative
sites or more are unaffected. The report and `recombination_profile.tsv` state the
window that was used.
```

- [ ] **Step 6: Harness gate**

Run: `validation/data/audit-c/gate.sh task7`

Acceptance:
- `benchmark.txt`: specificity `1.00`; the power count not lower than the baseline's. (Measured: 8 -> 13 significant of 185 testable recombining alignments; specificity 36/36 before and after.)
- `spec_mm1.txt`, `spec_mm2.txt`, `validation.txt`, `validation_regions.tsv`: identical to the previous task. PHI does not gate any region; it only sets the `parent_free_support` flag, which `validation_regions.tsv` does not list.
- Additionally count how many clonal runs now carry a significant PHI, which the gate script does not print:

```bash
python - <<'EOF'
import importlib.util, sys
import numpy as np
from pathlib import Path
from tessera.recomb.diagnostics import recombination_signal
spec = importlib.util.spec_from_file_location("run_specificity", Path("validation/run_specificity.py"))
rs = importlib.util.module_from_spec(spec); sys.modules["run_specificity"] = rs; spec.loader.exec_module(rs)
sig = n = 0
for scenario in rs.SCENARIOS:
    for rep in range(10):
        seqs = rs.simulate_clonal(scenario, seed=1000 + rep)
        rows = {k: np.frombuffer(v.encode(), dtype=np.uint8) for k, v in seqs.items()}
        s = recombination_signal(rows, rs.QUERY, lambda c: c)
        n += 1
        sig += s is not None and s.phi_p < 0.05
print(f"clonal alignments with PHI p < 0.05: {sig}/{n}")
EOF
```

  Acceptance: at most 4 of 40 (10 %; the documented rate on clonal data is 2.5-7.5 % against a nominal 5 %). Measured for this plan: 2/40. These panels have about 2600 informative sites, so the window is 100 before and after and the count should equal the one on `main`.

If the benchmark specificity drops below 1.00 or the clonal count exceeds 4/40: revert and record.

- [ ] **Step 7: Commit**

```bash
git add src/tessera/recomb/diagnostics.py src/tessera/cli/cmd_recomb.py docs/detection-methods.md tests/unit/test_diagnostics.py
git commit -m "$(cat <<'EOF'
Cap the PHI window at a tenth of the informative sites

The PHI window is counted in informative-site ranks and defaulted to 100
whatever the alignment held. PHI asks whether nearby sites are more
compatible than a random reordering; once the window spans most of the
sites nearly every pair is nearby and the test cannot reject. With 101 sites
or fewer it is exactly invariant under permutation (p = 1 always). On the
shipped divergent example, 159 sites, it gave p = 1.0 for a tract all four
callers find.

The window is now min(--phi-window, sites // 10). Of six rules compared,
this is the only one with no false positive on clonal simulations that also
recovers both shipped examples (p = 0.001 and 0.026). Alignments with 1000+
informative sites keep a window of 100.

Jaya 2023 benchmark: power <before> -> <after> of 185, specificity 36/36
before and after. Power there stays low for a reason the window does not
touch -- the biallelic-column filter leaves few sites at high divergence --
which is recorded in the PR, not addressed here.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

---

### Task 8: Report the span an HMM region's donor actually explains (C2)

An HMM segment is where the path sat in the donor state. Two things make that wider than the donor tract. Once in the donor state the path has no emission reason to leave while donor and major are identical, and the jump penalty keeps it there: a 600 bp difference was reported as 11.7 kb. And a segment's ends sit on the window grid (`positions[lo] - step` to `positions[hi] + step`), so a tract starting at column 0 is reported from the first window centre.

**Design.** The reported span is from the first to the last site inside the segment at which the query matches the donor and not the major. Sites matching both parents or neither cannot hold a region open. A segment that includes the first (last) window is searched from the alignment start (to its end). **The sign test is still computed on the segment**, so the set of called regions does not depend on this step; `mean_sim_minor` / `mean_sim_major` are identity over the reported columns (the quantity the site callers already report), under either windowing.

**Tried and not chosen:** the 3SEQ maximum-descent interval over the full extent of the segment's windows. Positive-control breakpoint error 9 bp against 13 bp, but on `hiv_crf02ag` an HMM region grew from 3744-4584 to 3834-4645, overlapped the coverage gap at 4624-7309 and removed that donor-absent row. A change that only tightens coordinates is easier to defend.

**Files:**
- Modify: `src/tessera/recomb/regions.py` (imports; `_call_regions_hmm` pass 2; new `_tract_span`)
- Create: `tests/integration/test_regimes.py`
- Modify: `docs/detection-methods.md`

**Interfaces:**
- Consumes: `validation/run_regimes.py` from Task 2 (`simulate_sibling_tract`, `scan`, `_evolve`, `_as_text`, `QUERY`); `threeseq.triplet_steps(rows, query, major, minor) -> (steps, columns)` where `steps` is `+1` for a site matching the major and `-1` for a site matching the minor.
- Produces: `regions._tract_span(work: WindowSimilarity, seg: Segment, major: str) -> tuple[int, int]`. HMM `Region.msa_start/msa_end/query_start/query_end` are the trimmed span; `pvalue`, `qvalue`, `support`, `n_windows`, `breakpoint_lo/hi` are unchanged. Task 9 appends tests to `tests/integration/test_regimes.py`.

- [ ] **Step 1: Write the failing tests**

Create `tests/integration/test_regimes.py`:

```python
"""Pipeline-level regressions for the regimes simulated by validation/run_regimes.py.

One seed each at reduced size, so they run in the fast suite; the harness itself
reports rates over replicates.
"""

from __future__ import annotations

import importlib.util
import logging
import sys
from pathlib import Path

import numpy as np

from tessera.recomb.analyze import analyze
from tessera.recomb.regions import RegionParams, call_regions
from tessera.recomb.similarity import compute_similarity

from ..conftest import write_fasta

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_regimes.py"
_SPEC = importlib.util.spec_from_file_location("run_regimes", _PATH)
rg = importlib.util.module_from_spec(_SPEC)
sys.modules["run_regimes"] = rg
_SPEC.loader.exec_module(rg)

WINDOW, STEP = 1000, 100
_LOG = logging.getLogger("tessera.test")


def _hmm_regions(result, window: int = WINDOW):
    params = RegionParams.with_defaults(window, method="hmm")
    regions, major, _ = call_regions(result, analyze(result), window, params)
    return regions, major


# --- HMM region span -------------------------------------------------------

def test_hmm_region_stops_where_donor_and_major_become_identical(tmp_path) -> None:
    # The query differs from the sibling only over `span`; downstream of it parent_A
    # (the donor relative to the sibling) and the sibling are identical.
    seqs, (lo, hi) = rg.simulate_sibling_tract(5000)
    msa = write_fasta(tmp_path / "sib.fasta", seqs)
    result = compute_similarity(str(msa), rg.QUERY, WINDOW, STEP)
    regions, major = _hmm_regions(result)
    assert major == "sibling"
    assert [r.minor_parent for r in regions] == ["parent_A"]
    region = regions[0]
    assert abs(region.msa_start - lo) <= STEP and abs(region.msa_end - hi) <= STEP
    assert region.length_bp <= (hi - lo) + STEP  # was 11.7 kb for a 600 bp difference


def test_hmm_region_reaches_the_alignment_start_for_a_terminal_tract(tmp_path) -> None:
    rng = np.random.default_rng(1)
    root = rng.integers(0, 4, size=3000)
    parent_a, parent_b, other = (rg._evolve(root, 0.05, rng) for _ in range(3))
    query = parent_a.copy()
    query[:800] = parent_b[:800]  # the donor tract starts at column 0
    seqs = rg._as_text({rg.QUERY: query, "parent_A": parent_a, "parent_B": parent_b,
                        "other": other})
    msa = write_fasta(tmp_path / "term.fasta", seqs)
    result = compute_similarity(str(msa), rg.QUERY, 300, 30)
    regions, _ = _hmm_regions(result, window=300)
    assert [r.minor_parent for r in regions] == ["parent_B"]
    # the first window centre is column 150; the tract must not be reported from there
    assert regions[0].msa_start <= 30
    assert abs(regions[0].msa_end - 800) <= 30


def test_trimmed_region_reports_query_coordinates_across_a_query_gap(tmp_path) -> None:
    # 300 alignment columns the query does not have (a backbone insertion), upstream of
    # the region: MSA and query coordinates then differ by exactly those 300 columns.
    seqs, (lo, hi) = rg.simulate_sibling_tract(5000)
    gap_at, gap = 3000, 300
    padded = {}
    for label, seq in seqs.items():
        insert = "-" * gap if label == rg.QUERY else seqs["parent_A"][gap_at:gap_at + gap]
        padded[label] = seq[:gap_at] + insert + seq[gap_at:]
    msa = write_fasta(tmp_path / "gapped.fasta", padded)
    result = compute_similarity(str(msa), rg.QUERY, WINDOW, STEP)
    regions, _ = _hmm_regions(result)
    assert [r.minor_parent for r in regions] == ["parent_A"]
    region = regions[0]
    assert abs(region.msa_start - (lo + gap)) <= STEP
    assert region.query_start == region.msa_start - gap
    assert region.query_end == region.msa_end - gap
    assert region.length_bp == region.length_msa


def test_trimmed_hmm_region_still_merges_with_the_site_callers(tmp_path) -> None:
    # Agreement is counted on overlapping regions. A trimmed HMM region must still
    # overlap the 3SEQ / MaxChi tract, or --min-methods 2 would drop a real event.
    seqs, (lo, hi) = rg.simulate_sibling_tract(5000)
    rows = rg.scan(seqs, _LOG, min_methods=2)
    hits = [r for r in rows if int(r["query_start"]) < hi and int(r["query_end"]) > lo]
    assert len(hits) == 1
    assert {"hmm", "3seq", "maxchi"} <= set(hits[0]["methods"].split(","))
```

The last two tests are guards for the review-focus items: they pass before and after this task and pin behaviour the trimming must not break.

- [ ] **Step 2: Run them to see them fail**

Run: `pytest tests/integration/test_regimes.py -q`
Expected: 2 failed, 2 passed.
- `test_hmm_region_stops_where_donor_and_major_become_identical`: `assert (100 <= 100 and 11000 <= 100)` -- the region is 17900-29600.
- `test_hmm_region_reaches_the_alignment_start_for_a_terminal_tract`: `assert 120 <= 30`.

- [ ] **Step 3: Implement `_tract_span`**

In `src/tessera/recomb/regions.py`:

(a) change two import lines:

```python
from .hmm import DEFAULT_JUMP_RATE, Segment, segment_query
```

and add after `from .stats import benjamini_hochberg, sign_test_pvalue`:

```python
from .threeseq import triplet_steps
```

(`threeseq` imports `regions` only inside its caller function, so there is no import cycle.)

(b) insert immediately before `def _call_regions_heuristic(`:

```python
def _tract_span(work: WindowSimilarity, seg: Segment, major: str) -> tuple[int, int]:
    """MSA columns ``[start, end)`` of an HMM segment that its donor actually explains.

    The span from the first to the last site where the query matches the donor and not
    the major. Columns where the query matches both parents or neither carry no
    information about which one it was copied from, so they cannot hold a region open:
    a run of windows in which donor and major are identical is left out, and so is the
    padding the window grid adds at either end.

    A segment that includes the first (last) window is searched from the alignment
    start (to its end): window centres stop half a window short of the ends, so a tract
    reaching an end would otherwise be reported from the first centre.

    Falls back to the segment's own span when no site in it favours the donor.
    """
    lo = 0 if seg.start_window == 0 else seg.msa_start
    hi = work.width if seg.end_window == len(work.positions) - 1 else seg.msa_end
    steps, cols = triplet_steps(work.rows, work.query, major, seg.state)
    donor_cols = cols[(steps == -1) & (cols >= lo) & (cols < hi)]
    if donor_cols.size == 0:
        return seg.msa_start, seg.msa_end
    return int(donor_cols[0]), int(donor_cols[-1]) + 1


```

(c) in `_call_regions_hmm`, replace everything from `        support = favor_minor / (favor_minor + favor_major)` down to and including the line `                query_start=seg.query_start, query_end=seg.query_end,` -- that is, this block:

```python
        support = favor_minor / (favor_minor + favor_major)
        if work.window_spans:
            # Informative-site windowing: the per-window values are identity at
            # polymorphic columns only, far below the identity over all columns (a
            # 99 %-identical donor can sit near 0.8). Reporting their mean as the
            # region's similarity would misstate it, and the coverage check compares
            # this value against a base-pair threshold. Use identity over the region's
            # columns, the same quantity the site-based callers report.
            mean_minor = region_identity(
                work.rows, work.query, seg.state, seg.msa_start, seg.msa_end
            )
            mean_major = region_identity(
                work.rows, work.query, major, seg.msa_start, seg.msa_end
            )
        else:
            idx = range(seg.start_window, seg.end_window + 1)
            minor_sims = [work.similarities[seg.state][i] for i in idx
                          if not isnan(work.similarities[seg.state][i])]
            major_sims = [work.similarities[major][i] for i in idx
                          if not isnan(work.similarities[major][i])]
            mean_minor = mean(minor_sims) if minor_sims else float("nan")
            mean_major = mean(major_sims) if major_sims else float("nan")
        regions.append(
            Region(
                minor_parent=seg.state, major_parent=major,
                msa_start=seg.msa_start, msa_end=seg.msa_end,
                query_start=seg.query_start, query_end=seg.query_end,
```

with

```python
        support = favor_minor / (favor_minor + favor_major)
        # The segment is where the HMM path sat in the donor state, which is not the
        # same as where the donor is the better match: once in the donor state the path
        # has no reason to leave while donor and major are identical, and its ends are
        # only as fine as the window grid. The reported span is the part of the segment
        # that the distinguishing sites support (see _tract_span). The sign test above
        # stays on the segment, so which regions are called does not depend on this
        # step -- only their coordinates do.
        lo, hi = _tract_span(work, seg, major)
        # Identity over the reported columns, the same quantity the site-based callers
        # report, on the same scale under either windowing (the per-window values of
        # informative-site windowing are identity at polymorphic columns only).
        mean_minor = region_identity(work.rows, work.query, seg.state, lo, hi)
        mean_major = region_identity(work.rows, work.query, major, lo, hi)
        regions.append(
            Region(
                minor_parent=seg.state, major_parent=major,
                msa_start=lo, msa_end=hi,
                query_start=work.column_to_query(lo), query_end=work.column_to_query(hi),
```

`isnan` and `mean` are still used by `_call_regions_heuristic`; leave their imports.

- [ ] **Step 4: Run the tests**

Run: `pytest tests/integration/test_regimes.py tests/unit/test_hmm.py tests/unit/test_clusters.py tests/integration -q && pytest -m "not requires_binary" -q && ruff check src tests && mypy src`
Expected: all pass. No existing test pins an HMM region's coordinates to the window grid.

- [ ] **Step 5: Document it**

In `docs/detection-methods.md`, in the "HMM caller" section, add after the paragraph that ends "...is flagged as marginal rather than dropped. The legacy `--method heuristic` (margin / merge-gap / min-region) is kept for comparison.":

```markdown
**The reported span is not the segment.** A segment is where the HMM path sat in the
donor state. Once there, the path has no reason to leave while donor and major are
identical, and its ends are only as fine as the window grid. The coordinates reported
for an HMM region therefore run from the first to the last site inside the segment at
which the query matches the donor and not the major; columns matching both parents or
neither cannot hold a region open. A segment that includes the first or last window is
searched to the end of the alignment, so a tract reaching a genome end is reported from
that end. The sign test, its p-value and `support` are computed on the whole segment, so
this step changes where a region is drawn and never whether it is called.
`mean_sim_minor` / `mean_sim_major` are identity over the reported columns.
```

- [ ] **Step 6: Harness gate**

```bash
validation/data/audit-c/gate.sh task8
python validation/run_regimes.py --reps 10 | tee validation/data/audit-c/task8/regimes_mm1.txt
python validation/run_regimes.py --reps 10 --min-methods 2 | tee validation/data/audit-c/task8/regimes_mm2.txt
python validation/run_hybrids.py > validation/data/audit-c/task8/hybrids.txt 2>&1
```

Acceptance (the spec's criteria for C2, made numeric):
- `sibling_tract`: median reported length **at most 1100 bp** at both gates (baseline 11700 bp; measured 1001 bp). The remaining excess over 600 bp is the Bootscan member of the merged region, whose span runs between window centres; it is not addressed here.
- `spec_mm1.txt`, `spec_mm2.txt`: the `TOTAL` lines identical to the previous task (measured: identical -- the same eight false regions, with tighter coordinates).
- Positive control: detection and donor unchanged; **median breakpoint error not higher** than the previous task's (measured 50 bp -> 13 bp).
- `validation.txt`: same verdict per dataset. `validation_regions.tsv`: **the same rows in the same order**; only coordinates of rows whose `methods` include `hmm` may differ. Measured differences, for comparison:

  | dataset | donor | before | after |
  |---|---|---|---|
  | `hcv_2k1b` | GT2a_JFH1 | 33-3127 | 33-3113 |
  | `hiv_crf02ag` | G | 1587-2624 | 1606-2614 |
  | `hiv_crf02ag` | G | 3744-4584 | 3744-4539 |
  | `hiv_crf02ag` | G | 7644-8989 | 7644-9058 (segment includes the last window) |
  | `norovirus_gii` | GII.1_Hawaii | 5040-7447 | 5070-7484 (last window) |
  | `sarscov2_xbb` | BA.2.75 | 21883-26888 | 22007-26112 |

  `hcv_2k1b`'s published breakpoint is near nt 3187: the reported end moves 14 bp further from it (3127 -> 3113). State that in the PR.
- Hybrids: no case changes from PASS to FAIL; where the table prints a breakpoint or donor-agreement column, no value worse than the baseline's.

If a row appears or disappears in `validation_regions.tsv`, or a hybrids case flips to FAIL: revert, record the diff, go to Task 9 -- and note that Task 9's gate was written assuming this task is in place (see its Step 6).

- [ ] **Step 7: Commit**

```bash
git add src/tessera/recomb/regions.py tests/integration/test_regimes.py docs/detection-methods.md
git commit -m "$(cat <<'EOF'
Report the span an HMM region's donor explains, not the whole segment

An HMM segment is where the path sat in the donor state. After a tract the
path has no emission reason to return to the major while the two are
identical downstream, and the jump penalty keeps it where it is: with a
sibling recombinant in the panel a 600 bp difference was reported as 11.7 kb
(39 % of the query), and the ensemble's union span carried it into the
merged region. Segment ends also sit on the window grid, so a tract starting
at column 0 was reported from the first window centre.

The reported span now runs from the first to the last site in the segment
where the query matches the donor and not the major, and a segment that
includes the first or last window is searched to the alignment end. Sites
matching both parents or neither cannot hold a region open.

The sign test stays on the segment. Recomputing it on the trimmed span would
lower p-values on a span chosen from the data, in the caller that already
produces every false region at --min-methods 1; keeping it means the set of
called regions cannot change here. mean_sim_minor / mean_sim_major are now
identity over the reported columns under either windowing.

Not chosen: the 3SEQ max-descent interval over the segment's full window
extent. 9 bp against 13 bp on the positive control, but it changed rows on
hiv_crf02ag (a region grew into a coverage gap and removed a donor-absent
row).

Harness: <paste sibling_tract lines, spec TOTAL lines, positive-control line,
validation_regions diff, hybrids summary or "hybrids not run: no network">.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

---

### Task 9: Lineage clustering on informative columns under informative-site windowing (C1)

`cluster_references` merges two references when their per-window identity over **all** columns stays at or above 0.985 with no region-sized run below it. On a panel less than about 1.5 % divergent every pair qualifies, the whole panel pools into one lineage, and the HMM has one state. That is the regime informative-site windowing exists for, so in that regime the HMM never contributes a call when there are 4-200 references.

**Design.** When the scan is under informative-site windowing (`result.window_spans` is non-empty), compare each pair on the informative columns only. There the 0.985 floor again separates true duplicates, which agree, from distinct lineages, which differ at a large share of the columns where the panel varies. Under base-pair windowing nothing changes.

**The cost, measured.** With clustering no longer pooling the panel, the HMM competes the references individually, as it does with `--no-cluster-lineages`. On 20 clonal near-identical panels that produced one HMM-only false region at `--min-methods 1` (1/20, against 0/20 on `main`) and none at `--min-methods 2`. The spec's criterion "specificity no worse" holds on `run_specificity.py` (which never enters informative-site windowing) and does not strictly hold on these panels. This is a trade the maintainer should see: it is carried into Task 10.

**Depends on Task 8.** Without the span trimming, an HMM region on such a panel is wide (median breakpoint error of the merged region 1087 bp with clustering off on `main`, against 263 bp without the HMM, at 6 replicates). With Task 8 in place the merged error is the same with and without the HMM (465 bp at 10 replicates). If Task 8 was reverted, do not start this task; record "C1 not attempted: depends on C2".

**Files:**
- Modify: `src/tessera/recomb/clusters.py` (imports; `cluster_references`)
- Modify: `tests/integration/test_regimes.py` (imports; three appended tests)
- Modify: `docs/detection-methods.md`

**Interfaces:**
- Consumes: `tests/integration/test_regimes.py` and its `rg`, `_hmm_regions`, `WINDOW` from Task 8; `rg.simulate_near_identical` from Task 2; `similarity._informative_column_mask(rows, ref_labels) -> np.ndarray` (boolean mask of columns polymorphic among the references).
- Produces: `cluster_references(result, window_size, params)` unchanged in signature.

- [ ] **Step 1: Write the failing tests**

In `tests/integration/test_regimes.py`, add to the imports

```python
from tessera.recomb.clusters import cluster_references
```

(after the `analyze` import) and change the similarity import to

```python
from tessera.recomb.similarity import compute_similarity, compute_similarity_informative
```

then append to the file:

```python
# --- lineage clustering on a near-identical panel ---------------------------

NEAR = {"length": 30_000, "distance": 0.004, "tract": (12_000, 18_000)}


def test_near_identical_references_are_not_pooled_into_one_lineage(tmp_path) -> None:
    seqs, _ = rg.simulate_near_identical(3, recombinant=True, **NEAR)
    seqs["R0_copy"] = seqs["R0"]  # a true duplicate: this is what clustering is for
    msa = write_fasta(tmp_path / "near.fasta", seqs)
    result = compute_similarity_informative(str(msa), rg.QUERY)
    params = RegionParams.with_defaults(WINDOW, method="hmm")
    clusters = {frozenset(c) for c in cluster_references(result, WINDOW, params)}
    assert frozenset({"R0", "R0_copy"}) in clusters
    assert len(clusters) == len(rg.NEAR_REFS)  # every other reference on its own


def test_duplicates_still_merge_when_informative_windowing_is_forced(tmp_path) -> None:
    # --informative-sites on a divergent panel: the comparison is on informative columns
    # there too, and true duplicates must still pool while distinct parents stay apart.
    rng = np.random.default_rng(4)
    root = rng.integers(0, 4, size=6000)
    parent_a, parent_b = rg._evolve(root, 0.05, rng), rg._evolve(root, 0.05, rng)
    seqs = rg._as_text({rg.QUERY: parent_a, "pa0": parent_a, "pa1": parent_a,
                        "pb0": parent_b, "pb1": parent_b})
    msa = write_fasta(tmp_path / "dups.fasta", seqs)
    result = compute_similarity_informative(str(msa), rg.QUERY)
    params = RegionParams.with_defaults(WINDOW, method="hmm")
    clusters = {frozenset(c) for c in cluster_references(result, WINDOW, params)}
    assert clusters == {frozenset({"pa0", "pa1"}), frozenset({"pb0", "pb1"})}


def test_hmm_calls_the_tract_on_a_near_identical_panel(tmp_path) -> None:
    seqs, (lo, hi) = rg.simulate_near_identical(3, recombinant=True, **NEAR)
    msa = write_fasta(tmp_path / "near.fasta", seqs)
    result = compute_similarity_informative(str(msa), rg.QUERY)
    regions, major = _hmm_regions(result)
    assert major == "R0"
    donors = [r for r in regions if r.minor_parent == "R1"]
    assert len(donors) == 1
    assert donors[0].msa_start < hi and donors[0].msa_end > lo
```

`test_duplicates_still_merge_when_informative_windowing_is_forced` is a guard: it passes before and after.

- [ ] **Step 2: Run them to see them fail**

Run: `pytest tests/integration/test_regimes.py -q`
Expected: 2 failed, 5 passed.
- `test_near_identical_references_are_not_pooled_into_one_lineage`: the only cluster is `{R0, R0_copy, R1, R2, R3, R4}`.
- `test_hmm_calls_the_tract_on_a_near_identical_panel`: `assert 0 == 1` -- no region.

- [ ] **Step 3: Compare on informative columns**

In `src/tessera/recomb/clusters.py`, replace the import

```python
from .similarity import WindowSimilarity, _best_per_window, _canonical_mask
```

with

```python
from .similarity import (
    WindowSimilarity,
    _best_per_window,
    _canonical_mask,
    _informative_column_mask,
)
```

and in `cluster_references` replace

```python
    canon = {label: _canonical_mask(result.rows[label]) for label in labels}
    uf = _UnionFind(labels)
    for i, a in enumerate(labels):
        a_row, a_canon = result.rows[a], canon[a]
        for b in labels[i + 1:]:
            comp = a_canon & canon[b]
```

with

```python
    canon = {label: _canonical_mask(result.rows[label]) for label in labels}
    # Informative-site windowing is what a near-identical panel is scanned under, and
    # on such a panel identity over all columns is above the floor for every pair: the
    # whole panel would pool into one lineage and leave the HMM a single state. Compare
    # the pair where the panel varies instead. There the floor again separates true
    # duplicates, which agree, from distinct lineages, which differ at a large share of
    # the informative columns.
    informative = (
        _informative_column_mask(result.rows, labels) if result.window_spans else None
    )
    uf = _UnionFind(labels)
    for i, a in enumerate(labels):
        a_row, a_canon = result.rows[a], canon[a]
        for b in labels[i + 1:]:
            comp = a_canon & canon[b]
            if informative is not None:
                comp = comp & informative
```

The next line (`match = comp & (a_row == result.rows[b])`) and the rest of the loop are unchanged.

- [ ] **Step 4: Run the tests**

Run: `pytest tests/integration/test_regimes.py tests/unit/test_clusters.py -q && pytest -m "not requires_binary" -q && ruff check src tests && mypy src`
Expected: all pass. The existing clustering tests run under base-pair windowing and are unaffected.

- [ ] **Step 5: Document it**

In `docs/detection-methods.md`, in "Low-divergence panels", add after the paragraph ending "...you cannot localise a switch more finely than the spacing of the discriminating sites.":

```markdown
Before the HMM runs, near-duplicate references are pooled into lineages so that
duplicates do not tie every window (`--cluster-lineages`, on by default for panels of 4
to 200 references). Two references pool when their agreement stays at or above 98.5 % in
every window. Under informative-site windowing that agreement is measured **on the
informative columns**: on a near-identical panel every pair is above 98.5 % over all
columns, and measuring it there would pool the whole panel and leave the HMM nothing to
choose between. On the informative columns only true duplicates agree.
```

- [ ] **Step 6: Harness gate**

```bash
validation/data/audit-c/gate.sh task9
python validation/run_regimes.py --reps 10 | tee validation/data/audit-c/task9/regimes_mm1.txt
python validation/run_regimes.py --reps 10 --min-methods 2 | tee validation/data/audit-c/task9/regimes_mm2.txt
python - <<'EOF' | tee validation/data/audit-c/task9/near_clonal.txt
import importlib.util, logging, sys
from collections import Counter
from pathlib import Path
spec = importlib.util.spec_from_file_location("run_regimes", Path("validation/run_regimes.py"))
rg = importlib.util.module_from_spec(spec); sys.modules["run_regimes"] = rg; spec.loader.exec_module(rg)
log = logging.getLogger("x"); log.addHandler(logging.NullHandler()); log.propagate = False
for min_methods in (1, 2):
    bad = 0
    callers = Counter()
    for rep in range(20):
        seqs, _ = rg.simulate_near_identical(4000 + rep, recombinant=False)
        n, per = rg.score_clonal(rg.scan(seqs, log, min_methods=min_methods))
        bad += bool(n)
        callers.update(per)
    print(f"near-identical clonal, --min-methods {min_methods}: "
          f"runs with a false region {bad}/20 {dict(callers)}")
EOF
python validation/run_hybrids.py > validation/data/audit-c/task9/hybrids.txt 2>&1
```

Acceptance (the spec's criteria for C1):
- **The HMM calls the tract.** `near_identical`: `called by the HMM` at least 9/10 at both gates (baseline 0/10; measured 10/10), `detected` 10/10, `correct donor` 10/10.
- Median breakpoint error on `near_identical` not higher than Task 8's (measured 465 bp before and after).
- **Specificity no worse on `run_specificity.py`:** `spec_mm1.txt` and `spec_mm2.txt` `TOTAL` lines identical to Task 8's. (Those panels run under base-pair windowing, so this is expected by construction; a difference means the change leaked into base-pair mode.)
- `validation.txt` and `validation_regions.tsv` identical to Task 8's (measured: identical).
- **Near-identical clonal:** 0/20 at `--min-methods 2`. At `--min-methods 1` record the count. Measured: 1/20 (one HMM-only region, q = 0.021) against 0/20 on `main`. A count of 0 or 1 is consistent with what was measured; **2 or more of 20 is beyond it -- stop, revert, and take the number to the maintainer.**
- **Hybrids unchanged or better:** no case from PASS to FAIL. Look in particular at the cases analysed under informative-site windowing (`mode` column: `mpox`, `vzv`, `ebola`, `lowdiv_rsv`, `neg_sarscov2`) and at every `neg_*` case: a negative control that gains an HMM region fails the gate.

If the gate fails: revert, record all five outputs. The recorded 1/20 goes into Task 10's decision record whether or not this task is kept.

- [ ] **Step 7: Commit**

```bash
git add src/tessera/recomb/clusters.py tests/integration/test_regimes.py docs/detection-methods.md
git commit -m "$(cat <<'EOF'
Cluster on informative columns under informative-site windowing

Lineage clustering merges two references whose per-window identity over all
columns stays at or above 98.5 %. On a panel less than ~1.5 % divergent every
pair qualifies, the whole panel pools into one lineage and the HMM has a
single state: with 4-200 references it could not call anything, in exactly
the regime informative-site windowing was added for. In 10 of 10 simulated
panels ~0.2 % divergent the HMM was absent from a call that 3SEQ and MaxChi
made; with --no-cluster-lineages it was present in all 10.

Under informative-site windowing the pairwise agreement is now measured on
the informative columns, where the floor again separates true duplicates
from distinct lineages. Base-pair windowing is unchanged.

The cost: the HMM now competes those references individually and can be
wrong on its own. On 20 clonal near-identical panels: <N>/20 runs with an
HMM-only false region at --min-methods 1 (0/20 before), 0/20 at
--min-methods 2. run_specificity.py is unchanged at both gates.

Harness: <paste near_identical lines, near_clonal lines, spec TOTAL lines,
hybrids summary or "hybrids not run: no network">.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

---

### Task 10: Characterise the false regions and record the decision (C8)

At the CLI default (`--min-methods 1`) the specificity harness reports a false region in about half of its clonal runs, all from the HMM. `run_specificity.py` defaults to `--min-methods 2`, so its clean result does not describe the shipped default. This task does not change the default. It commits a script that shows what the false regions are, and a decision record that puts the options and their measured costs in front of the maintainer.

**Files:**
- Create: `validation/characterise_false_regions.py`
- Create: `docs/superpowers/specs/2026-10-01-min-methods-default-decision.md`
- Modify: `validation/README.md`

**Interfaces:**
- Consumes: `validation/run_specificity.py` (`SCENARIOS`, `QUERY`, `simulate_clonal`).
- Produces: `characterise_false_regions.false_region_rows(scenario, seed, seqs, query) -> list[dict]` and a TSV on stdout with the columns in `COLUMNS`.

- [ ] **Step 1: Write the script**

Create `validation/characterise_false_regions.py`:

```python
#!/usr/bin/env python
"""Describe each false region the callers report on the specificity harness's clonal data.

``run_specificity.py`` counts false regions; this lists them, one row per region and
per caller, with the quantities needed to see *why* each was called: its length, how
many sites distinguish the two parents and which way they lean, its p- and q-value,
whether donor and major belong to the same simulated clade, and whether lineage
clustering pooled anything. It exists to inform a decision about the ``--min-methods``
default and the HMM's significance gate, not to gate a change.

Needs no aligner, network or downloaded data.

    python validation/characterise_false_regions.py             # 10 replicates
    python validation/characterise_false_regions.py --reps 3
"""

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

from tessera.recomb.analyze import analyze
from tessera.recomb.clusters import all_singletons, cluster_references
from tessera.recomb.regions import DEFAULT_METHODS, RegionParams, call_regions
from tessera.recomb.similarity import compute_similarity, discordant_counts

WINDOW, STEP = 1000, 100  # the CLI defaults, as run_specificity.py uses
COLUMNS = (
    "scenario", "seed", "caller", "minor", "major", "msa_start", "msa_end", "length",
    "favor_minor", "favor_major", "pvalue", "qvalue", "sim_minor", "sim_major",
    "same_clade", "panel_clustered",
)


def _load_specificity():
    path = Path(__file__).resolve().parent / "run_specificity.py"
    spec = importlib.util.spec_from_file_location("run_specificity", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["run_specificity"] = module
    spec.loader.exec_module(module)
    return module


def false_region_rows(scenario: str, seed: int, seqs: dict[str, str], query: str) -> list[dict]:
    """One row per region any default caller reports on a clonal alignment."""
    rows: list[dict] = []
    with tempfile.TemporaryDirectory() as td:
        msa = Path(td) / "aln.fasta"
        msa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
        result = compute_similarity(str(msa), query, WINDOW, STEP)
    analysis = analyze(result)
    clustered = not all_singletons(
        cluster_references(result, WINDOW, RegionParams.with_defaults(WINDOW))
    )
    for method in DEFAULT_METHODS:
        params = RegionParams.with_defaults(WINDOW, method=method)
        regions, _major, _siblings = call_regions(result, analysis, WINDOW, params)
        for r in regions:
            favor_minor, favor_major = discordant_counts(
                result.rows, query, r.major_parent, r.minor_parent, r.msa_start, r.msa_end
            )
            rows.append({
                "scenario": scenario, "seed": seed, "caller": method,
                "minor": r.minor_parent, "major": r.major_parent,
                "msa_start": r.msa_start, "msa_end": r.msa_end,
                "length": r.msa_end - r.msa_start,
                "favor_minor": favor_minor, "favor_major": favor_major,
                "pvalue": "" if r.pvalue is None else f"{r.pvalue:.3g}",
                "qvalue": "" if r.qvalue is None else f"{r.qvalue:.3g}",
                "sim_minor": r.mean_sim_minor, "sim_major": r.mean_sim_major,
                # tips are named <clade letter><index>, e.g. A2
                "same_clade": r.minor_parent[0] == r.major_parent[0],
                "panel_clustered": clustered,
            })
    return rows


def main(argv: list[str]) -> int:
    reps = int(argv[argv.index("--reps") + 1]) if "--reps" in argv else 10
    rs = _load_specificity()
    print("\t".join(COLUMNS))
    total = 0
    for scenario in rs.SCENARIOS:
        for rep in range(reps):
            seed = 1000 + rep
            for row in false_region_rows(scenario, seed,
                                         rs.simulate_clonal(scenario, seed=seed), rs.QUERY):
                print("\t".join(str(row[c]) for c in COLUMNS))
                total += 1
    print(f"# {total} false region(s) over {reps * len(rs.SCENARIOS)} clonal alignment(s); "
          "single-caller rows, before the ensemble merge and any --min-methods gate",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
```

- [ ] **Step 2: Run it**

```bash
mkdir -p validation/data/audit-c/task10
python validation/characterise_false_regions.py --reps 10 > validation/data/audit-c/task10/false_regions.tsv
column -t -s$'\t' validation/data/audit-c/task10/false_regions.tsv | head -40
ruff check validation
```

The table it must produce has one row per region and these columns:

`scenario  seed  caller  minor  major  msa_start  msa_end  length  favor_minor  favor_major  pvalue  qvalue  sim_minor  sim_major  same_clade  panel_clustered`

For reference, on `main` at 3 replicates it produced eight rows, all `caller = hmm`:

| scenario | seed | minor | major | length | favor_minor | favor_major | pvalue | qvalue | same_clade | panel_clustered |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | 1000 | A2 | A1 | 1400 | 20 | 8 | 0.0178 | 0.0357 | True | False |
| clean | 1001 | A1 | A3 | 4400 | 59 | 40 | 0.035 | 0.035 | True | False |
| asrv | 1000 | A2 | A3 | 2100 | 29 | 14 | 0.0158 | 0.0315 | True | False |
| asrv | 1002 | A1 | A0 | 1800 | 29 | 14 | 0.0158 | 0.0394 | True | False |
| asrv | 1002 | A2 | A0 | 1000 | 15 | 3 | 0.00377 | 0.0188 | True | False |
| lineage_rate | 1000 | A2 | A1 | 1600 | 20 | 8 | 0.0178 | 0.0357 | True | False |
| lineage_rate | 1001 | A1 | A0 | 1400 | 24 | 5 | 0.000273 | 0.00164 | True | False |
| rate_shift | 1002 | A3 | A0 | 1300 | 43 | 23 | 0.00933 | 0.0373 | True | False |

(With Task 8 in place the same eight regions appear with lengths 894-4049.)

- [ ] **Step 3: Summarise what the table shows**

Compute, from your 10-replicate table: the number of rows per caller; the share with `same_clade = True`; the share with `panel_clustered = True`; the median `length`; the median `favor_minor + favor_major`; and how many rows have `qvalue > 0.01`. These six numbers go into the decision record in Step 5.

What the 3-replicate table showed, to compare against:
- every false region is an HMM call;
- every one is a switch between two tips of the query's **own clade** (the references about 2.4 % apart), never across clades;
- none of the panels was pooled by lineage clustering (the tips are below the 98.5 % floor), so the tips competed individually;
- the regions are 1-4.4 kb and rest on 18-99 distinguishing sites;
- seven of eight have a q-value between 0.01 and 0.04.

The mechanism this points to: the HMM chooses a segment because the donor wins there, and the sign test is then run on the same sites. The test does not know a segment was searched for, so it is anti-conservative; 3SEQ and MaxChi, whose nulls account for the scan, called none of the eight.

- [ ] **Step 4: Measure the two options that do not need a design**

(a) The agreement gate: already in `validation/data/audit-c/task9/spec_mm2.txt` and `regimes_mm2.txt`.

(b) `--alpha 0.01` for the whole scan, which the table suggests would remove most rows. `run_specificity.py` has no alpha flag, so:

```bash
python - <<'EOF' | tee validation/data/audit-c/task10/alpha_001.txt
import importlib.util, logging, sys, tempfile
from pathlib import Path
from tessera.recomb.run import RecombParams, run_recomb
spec = importlib.util.spec_from_file_location("run_specificity", Path("validation/run_specificity.py"))
rs = importlib.util.module_from_spec(spec); sys.modules["run_specificity"] = rs; spec.loader.exec_module(rs)
log = logging.getLogger("x"); log.addHandler(logging.NullHandler()); log.propagate = False

def scan(seqs):
    with tempfile.TemporaryDirectory() as td:
        msa = Path(td) / "aln.fasta"
        rs._write_fasta(msa, seqs)
        run_recomb(RecombParams(msa=msa, output=Path(td) / "out", query=rs.QUERY,
                                plot_format="png", min_methods=1, alpha=0.01), log)
        return rs._read_regions(Path(td) / "out" / "recombination_regions.tsv")

bad = runs = 0
for scenario in rs.SCENARIOS:
    for rep in range(10):
        n, _ = rs.score_negative(scan(rs.simulate_clonal(scenario, seed=1000 + rep)))
        bad += bool(n)
        runs += 1
hits = 0
for rep in range(10):
    seqs, tract = rs.simulate_recombinant(seed=2000 + rep)
    hits += rs.score_positive(scan(seqs), tract, donor_prefix="B")["detected"]
print(f"--alpha 0.01, --min-methods 1: runs with a false region {bad}/{runs}; "
      f"positive control detected {hits}/10")
EOF
```

Neither option is implemented as a default here; these are numbers for the record.

- [ ] **Step 5: Write the decision record**

Create `docs/superpowers/specs/2026-10-01-min-methods-default-decision.md` with the text below. The numbers in it were measured while this plan was written, at the replicate counts stated (option B's on a scratch copy that is not part of this branch -- keep those as they are, labelled as measured on 2026-10-01); **replace every other number with your own 10-replicate measurement** (the file names are given beside each) and keep the replicate count beside every number.

```markdown
# Decision record: false regions at the default agreement gate

**Status:** open -- for the maintainer. **Date opened:** 2026-10-01.
**Question:** at the CLI default (`--min-methods 1`) the specificity harness reports a
false region in roughly half of its clonal runs. Should the default change, should the
HMM's significance gate change, or should the default stay and be documented?

## What was measured

All on simulated clonal data (`validation/run_specificity.py`: four clades about 16 %
apart, four tips each about 2.4 % apart, a clonal query in clade A). Every region is a
false positive by construction. Replicate counts are small; the intervals are wide.

| setting | runs with a false region | source | positive control (3 kb tract) |
|---|---|---|---|
| `--min-methods 1` (CLI default) | 7/12 (58 %, CI 32-81 %) | hmm = 8 of 8 regions | 3/3 detected, 3/3 donor |
| `--min-methods 2` (harness default) | 0/12 (CI 0-24 %) | -- | 3/3 detected, 3/3 donor |

Source files: `validation/data/audit-c/task9/spec_mm1.txt`, `spec_mm2.txt`.

On near-identical panels (`validation/run_regimes.py`, references about 0.2 % apart),
after lineage clustering stopped pooling the whole panel (audit item C1):

| setting | clonal runs with a false region | tract detected | HMM among the callers |
|---|---|---|---|
| `--min-methods 1` | 1/20 (hmm) | 10/10 | 10/10 |
| `--min-methods 2` | 0/20 | 10/10 | 10/10 |

Source files: `validation/data/audit-c/task9/near_clonal.txt`, `regimes_mm1.txt`,
`regimes_mm2.txt`.

## What the false regions are

From `validation/characterise_false_regions.py` (3 replicates, 8 regions):

- all eight are HMM calls; 3SEQ, MaxChi and Bootscan called none of them;
- all eight are switches between two tips of the query's own clade, about 2.4 % apart;
  none crosses a clade;
- none of the panels was pooled by lineage clustering (the tips sit below its 98.5 %
  floor), so the tips competed as individual genomes;
- lengths 1.0-4.4 kb, resting on 18-99 distinguishing sites;
- seven of eight have a q-value between 0.01 and 0.04.

The likely mechanism: the HMM picks a segment because the donor wins there, and the sign
test is then run on the same sites. It is not aware that a segment was searched for.
3SEQ and MaxChi test against nulls that account for the scan.

## Options

### A. Change the default to `--min-methods 2`

- For: 0/12 and 0/20 false regions here, with the positive controls intact.
- Against: this was tried and reverted (PR #47). On the hybrid harness it lost three
  true detections (`rsv_a`, `mpox`, `masksib_rsv`) and gained nothing, because every
  negative control already passed. At very low divergence only the HMM has power, so a
  gate of two is not reachable there.
- Evidence still needed: `run_hybrids.py` at `--min-methods 2` on the current branch.

### B. Keep `--min-methods 1`; make the HMM's gate aware of the scan

One concrete form was measured, not implemented: an HMM segment must pass both the sign
test and the 3SEQ max-descent test for its (major, donor) pair, taking the larger
p-value.

| measurement (3 or 10 replicates as stated) | current gate | scan-aware gate |
|---|---|---|
| `run_specificity.py --reps 3 --min-methods 1`: runs with a false region | 7/12 | 2/12 |
| positive control | 3/3, 3/3, 13 bp | 3/3, 3/3, 13 bp |
| near-identical, 10 reps: HMM among the callers | 10/10 | 10/10 |
| near-identical clonal, 10 reps: runs with a false region | 1/10 | 1/10 |
| `run_validation.py` verdicts | 6 PASS, 1 FAIL, 1 SKIP | the same |
| `hiv1_crf`, AE_env region 5569-8177 | `hmm,bootscan` | `bootscan` (the HMM vote is lost) |
| p-value of the four strongest regions | 1e-32 to 1e-67 | 5e-05 (permutation floor) |

- For: removes about two thirds of the false regions without an agreement gate.
- Against: it removed a real HMM call on HIV-1 CRF01_AE; it replaces exact small
  p-values by the permutation floor wherever the triplet is too large for the exact
  test; it was not run on the hybrid harness, where the HMM-only detections live.
- Evidence still needed: `run_hybrids.py`; a design that keeps the exact p-value.

### C. Keep everything; document the rate and recommend `--min-methods 2` for redundant panels

- For: no behaviour change; the README and `--min-methods` help already describe the
  trade.
- Against: a default that reports a false region in about half of clonal runs on a
  redundant panel is a weak default for a first-time user, and `hcv_clonal_1b` (a real
  clonal control) currently fails on a single-caller region.
- Needed: state the measured rate in `README.md` and `validation/README.md`, and make
  the harness and the CLI use the same default or say plainly that they do not.

### D. Lower `--alpha` for the HMM only

Seven of the eight false regions have q between 0.01 and 0.04.

Measured with `--alpha 0.01` for the whole scan, 3 replicates
(`validation/data/audit-c/task10/alpha_001.txt`): 1/12 clonal runs with a false region,
positive control 3/3 detected. It was not measured on real data or on the hybrid harness,
where weak true regions (several real HMM calls have q between 0.01 and 0.05) would be
lost, and an alpha for one caller that differs from the others needs a reason beyond "it
removes these rows".

## Recommendation from the audit

Not A on its own -- the hybrid harness already answered that. B is the option that
addresses the mechanism, but the form measured here is not ready: it cost a real HMM
call and coarsened p-values. The audit's suggestion is to take C now (it is
documentation only) and treat B as a design task of its own, gated on the hybrid
harness.

## Decision

_To be filled in by the maintainer: option chosen, date, and the harness numbers it was
taken on._
```

The final section is left for the maintainer on purpose: this is the one line in the plan that the executor must not fill in.

- [ ] **Step 6: List the script in the validation README**

In `validation/README.md`, in the layout block, add after the `run_regimes.py` line:

```
  characterise_false_regions.py  one row per false region on the clonal simulations
```

- [ ] **Step 7: Commit**

```bash
git add validation/characterise_false_regions.py validation/README.md \
    docs/superpowers/specs/2026-10-01-min-methods-default-decision.md
git commit -m "$(cat <<'EOF'
Characterise the false regions at the default gate; open a decision record

run_specificity.py defaults to --min-methods 2 and is clean there. The CLI
defaults to 1, where the same harness reports a false region in about half
of its clonal runs. The clean result therefore does not describe what a
user gets, and nothing in the repository said so.

characterise_false_regions.py lists each false region with the quantities
needed to see why it was called. On the current callers every one is an HMM
call between two tips of the query's own clade, resting on a few dozen
distinguishing sites, with q mostly between 0.01 and 0.04 -- consistent with
a sign test run on a segment the HMM chose because the donor wins there.

The default is not changed here. The decision record sets out four options
with what was measured for each and what is still missing (above all the
hybrid harness), and leaves the choice to the maintainer.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 8: Ask the maintainer**

Post the decision record's "What was measured", "Options" and "Recommendation" sections in the PR description and ask for a decision. Do not implement any option in this branch.

---

### Task 11: Changelog, final gate, pull request

**Files:**
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes: the gate outputs under `validation/data/audit-c/` from every earlier task.

- [ ] **Step 1: Write the changelog**

In `CHANGELOG.md`, under `## [Unreleased]`, add the block below. **Delete the bullet of any task whose gate failed and was reverted**, and replace each number with the one you measured.

```markdown
### Changed

- **HMM region coordinates are the span the donor explains, not the HMM segment.** A
  segment ran on through columns where donor and major are identical (a 600 bp difference
  was reported as 11.7 kb when a sibling recombinant was in the panel) and its ends sat on
  the window grid. The reported span now runs from the first to the last site in the
  segment where the query matches the donor and not the major, and reaches the alignment
  end for a terminal tract. **`msa_start`, `msa_end`, `query_start`, `query_end`,
  `length_bp`, `mean_sim_minor` and `mean_sim_major` change for regions the HMM called**;
  which regions are called, and their p-values, do not.
- **On near-identical panels the HMM can contribute again.** Lineage clustering pooled an
  entire panel less than about 1.5 % divergent into one lineage, leaving the HMM a single
  state. Under informative-site windowing the pairwise agreement is now measured on the
  informative columns. On simulated clonal panels of this kind the HMM alone then reported
  a false region in 1 of 20 runs at `--min-methods 1` (0 of 20 before, and 0 of 20 at
  `--min-methods 2`).
- **The PHI window is capped at a tenth of the informative sites.** `--phi-window` is an
  upper bound; the report states the window used. With the fixed 100-rank window the test
  could not reject on alignments with few informative sites (p = 1 on the shipped
  `divergent` example, now 0.001). Alignments with 1000 or more informative sites are
  unaffected.
- **Pool clade labels are read as written.** Short and multi-word clades (`B.1`, `XBB`,
  `D8`, `IIb`, the HIV pure subtypes, `clade 2 wild-type`) were dropped when a Nextclade
  or NCBI Virus pool was typed, so lineage-aware selection and the recombinant-lineage
  exclusion skipped those genomes. **Panels seeded from such pools change**; the labels
  of seeded genomes are recorded in `<output>/pool_labels.tsv`.

### Fixed

- A 3SEQ or MaxChi tract started at the first maximum of the discriminating-site walk
  before its trough, not the last, so it could begin a few sites early and include a
  stretch that is not donor tract. MaxChi's permutation null follows the same rule.
- Exact 3SEQ p-values below about 1e-16 were written as `0.0`. They now keep their
  magnitude.
- The ensemble merge could label a region with the same genome as donor and major when
  the callers disagreed on the backbone, and could report one donor lineage as two
  overlapping regions on a typed panel.

### Added

- `validation/run_regimes.py`: a self-contained harness for near-identical panels and
  sibling tracts, the two regimes the specificity harness does not sample.
- `validation/characterise_false_regions.py`, and a decision record on the false-positive
  rate at the default agreement gate
  (`docs/superpowers/specs/2026-10-01-min-methods-default-decision.md`).
```

- [ ] **Step 2: Final checks**

Run: `ruff check src tests validation && pytest -m "not requires_binary" -q && mypy src`
Expected with every task kept: `All checks passed!`, `618 passed, 1 deselected`, `Success: no issues found in 79 source files`. (586 at the start; this plan adds 32 tests.) With tasks reverted the count is lower by that task's tests; it must not be below 586.

- [ ] **Step 3: Final gate**

```bash
validation/data/audit-c/gate.sh final
python validation/run_regimes.py --reps 10 | tee validation/data/audit-c/final/regimes_mm1.txt
python validation/run_regimes.py --reps 10 --min-methods 2 | tee validation/data/audit-c/final/regimes_mm2.txt
python validation/run_hybrids.py > validation/data/audit-c/final/hybrids.txt 2>&1
diff validation/data/audit-c/baseline/validation_regions.tsv validation/data/audit-c/final/validation_regions.tsv
```

Acceptance: `final` equals the last kept task's outputs. Build the before/after table for the PR from `baseline/` and `final/`.

- [ ] **Step 4: Commit and open the pull request**

```bash
git add CHANGELOG.md
git commit -m "$(cat <<'EOF'
Record the caller-behaviour changes in the changelog

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
git push -u origin fix-audit-caller-behaviour
gh pr create --title "Caller behaviour: audit plan C" --body "$(cat <<'EOF'
Implements plan C of the post-1.2.0 audit
(docs/superpowers/specs/2026-10-01-post-1.2.0-audit-design.md, items C1-C8).

## Harness results, before and after

<table built from validation/data/audit-c/baseline and /final: specificity at both
gates, positive control, run_validation verdicts and the validation_regions diff,
run_regimes at both gates, near-identical clonal, PHI benchmark, hybrids>

## Gate results per task

<one line per task: kept or reverted, with the numbers that decided it>

## Skipped

<every gate that could not be run and why -- e.g. "run_hybrids.py: no network";
"Task 6 end-to-end: not run">

## Found, not fixed

- PHI has low power on the Jaya 2023 benchmark under every window rule: the
  biallelic-column filter leaves no informative sites at mut = 0.1 and about 65 at
  mut = 0.01.
- The merged region on the sibling-tract regime is still about 400 bp longer than the
  true difference: that is the Bootscan member, whose span runs between window centres.
- `hcv_clonal_1b` still fails (a 12 bp MaxChi-only region); `validation/README.md` still
  lists it as 0 regions. Plan B owns that text.

## Decision needed

<the "What was measured", "Options" and "Recommendation" sections of
docs/superpowers/specs/2026-10-01-min-methods-default-decision.md>

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

Fill every `<...>` in the body from the saved outputs before running `gh pr create`. Do not merge; hand the PR to the maintainer.
