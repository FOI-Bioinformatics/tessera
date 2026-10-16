# Audit Plan G: Caller Statistics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the caller-statistics defects A1-A14 of the fourth audit, make `--min-region` bound every reported region (B2), restate HMM and merged-row window counts in one unit (B4 code half, A8), close the trailing-window and incompatible-pair items, and add the caller-side tests the audit names -- each behind a measured harness gate.

**Architecture:** Plan C's evaluate-first discipline throughout: a failing test, one candidate implementation, then a harness gate with numeric acceptance criteria and a stated action on failure. A small extension of `validation/run_regimes.py` comes first, because no harness samples a redundant panel and A1/A2 cannot be gated without one. Then independent numerical fixes whose harness output must be identical to the baseline (sign test, HMM segment spacing, merged-row statistics, Bootscan null, MaxChi overflow, memory caps, terminal breakpoint interval), then the behaviour-changing ones (cluster-size emission, shared lineage clusters, the descent floor, `--min-region`, the trailing window), and last a documentation correction and decision-record addendum for the 3SEQ/MaxChi dependence (A12) and the barcode vote (A14).

**Tech Stack:** Python 3.11+, numpy, pytest, ruff, mypy. No new dependency.

**Spec:** `docs/superpowers/specs/2026-10-12-audit-swot-gap.md` (items A1-A14 in section 1.1; B2; the code half of B4; the two "further items from this pass" -- trailing unscanned columns and the incompatible-pair list; the "Caller side" test gaps in section 3; gap G2). Read sections 1.1 and 3 before starting. Read also `docs/superpowers/plans/2026-10-01-audit-c-caller-behaviour.md` (the discipline this plan follows) and `docs/superpowers/specs/2026-10-01-min-methods-default-decision.md` (Task 15 appends to it).

## Global Constraints

- No new runtime dependency.
- Modest scientific language in code, docs and messages; reported numbers must be faithful (state what passes, what fails, what was skipped).
- "Could not test" must never be reported as "tested, found nothing".
- A new behaviour needs a test that fails without the change; if a test is a guard that passes before and after, the plan says so.
- Branch: create `fix-audit4-callers` from `fix-audit-caller-behaviour` (or from `main` once that branch is merged). Commit messages explain the problem, the evidence and the trade-off, and end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Append the aligner env to PATH, never prepend (prepending swaps the Python interpreter): `export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"`.
- The dev Mac has 18 GB RAM: run one harness at a time.
- Each task ends with `ruff check src tests validation`, `mypy`, and `pytest -m "not requires_binary" -q` green before its commit.
- No Unicode in code or docs beyond what files already contain; ASCII punctuation.
- Anything touching a caller, a default or region calling is validated on `validation/run_specificity.py` at **both** `--min-methods 1` and `--min-methods 2`, on `validation/run_regimes.py` at both gates, and on `validation/run_hybrids.py`, before and after; the numbers go in the PR description.
- **A failed gate is a result, not an obstacle.** If a task's gate fails: `git revert` that task's commit (or do not commit it), write the measured numbers under "Gate results" in the PR description, and move on to the next task. Do not adjust thresholds, constants or the test until the gate passes. A later task that names the reverted one under "Depends on" is then not attempted; record "not attempted: depends on Task N".
- The `--min-methods` default is not changed in this branch. Task 15 records the evidence for the open decision; it does not decide it.
- Do not touch items owned by other audit-4 plans (B1, B3, B5-B9, C*, D*, E*, F*). B4 is split: this plan changes the code (`n_windows`, `posterior` after trimming); wording in reports that are not region columns belongs to the reporting plan.

## What was measured while writing this plan

The candidates below were applied cumulatively in a scratch worktree of `4b35e39` and measured with `validation/run_specificity.py --reps 3` at both gates, `validation/run_regimes.py --reps 3` (at `--min-methods 2` for the first and last columns only), and the redundant regimes at 10 replicates with a scratch driver over the same simulator and seeds as Task 2. These are the reference values the gates compare against; the executor re-measures at 10 replicates in Tasks 1 and 2. Replicate counts are small: read them as "does the failure mode exist", not as rates.

| measurement | `4b35e39` | + Tasks 3-9 (numerical) | + Tasks 10-11 (A1, A2) | + Tasks 12-14 (all candidates) |
|---|---|---|---|---|
| `run_specificity` `--min-methods 1`: runs with a false region | 7/12 (8, hmm = 8) | 7/12 (8, hmm = 8) | 7/12 (8, hmm = 8) | 7/12 (8, hmm = 8) |
| `run_specificity` `--min-methods 2` | 0/12 | 0/12 | 0/12 | 0/12 |
| positive control (both gates): detected / donor / median breakpoint error | 3/3, 3/3, 13 bp | same | same | same |
| `near_identical`: detected / HMM among callers / donor / error | 3/3, 3/3, 3/3, 378 bp | same | same | same (also at `--min-methods 2`) |
| `near_identical` clonal: runs with a false region | 0/3 | 0/3 | 0/3 | 0/3 |
| `sibling_tract`: median reported length | 1001 bp | 1001 bp | 1001 bp | 1001 bp |
| `redundant_backbone` (10 reps, `--min-methods 1`): HMM among callers | 3/10 | -- | 6/10 (Task 10 alone: 6/10) | 6/10 |
| `redundant_donor` (10 reps, `--min-methods 1`): one row / HMM / Bootscan | 7/10, 0/10, 1/10 | -- | 10/10, 10/10, 10/10 (Task 10 alone: 7/10, 10/10, 1/10) | 10/10, 10/10, 10/10 |
| `redundant_donor` (10 reps, `--min-methods 2`): detected | 7/10 | -- | -- | 10/10 |
| control, one genome per lineage (10 reps): HMM among callers | 8/10 | -- | -- | 8/10 |
| one 16-genome specificity scan, wall clock | 4.3 s | -- | 4.3 s | 4.5 s (45 s with the first, rejected form of Task 12) |
| the three 3-replicate harnesses together | 6:46 | -- | -- | 6:53 |
| fast suite | 757 passed | -- | -- | 796 passed |

**Evaluated and not chosen -- padding HMM segments by half the window width (spec A10, literal reading).** Applied on top of all candidates (base-pair segments padded by `window_size // 2`, informative-site segments taken to their window spans):

| measurement | chosen (median step) | half the window width |
|---|---|---|
| `run_specificity` `--min-methods 1`: runs with a false region | 7/12 (8 hmm) | 1/12 (2 hmm) |
| positive control: median breakpoint error | 13 bp | 208 bp |
| `near_identical`: HMM among callers | 3/3 | 0/3 |
| `redundant_backbone` (3 reps): HMM among callers | 2/3 | 1/3 |

Wider segments add flanking sites that favour the major to the HMM's sign test, which removes most of the HMM's false regions on the specificity panels -- and the HMM's true calls on near-identical panels, while the breakpoints move 200 bp. That is a change to the HMM's significance gate, the subject of option B in the open decision record, not a coordinate fix; it is reported to the maintainer (Task 17, "Decision needed") and not implemented here.

**Not measured here:** `run_validation.py` and `run_hybrids.py` (no aligner binaries and no network in this session). Both are executor gates for every task that names them. A1 and A2 were measured with the harness's redundant regimes and a scratch driver only; no real panel with near-duplicates was scanned.

## Review Focus

Inputs or failure modes the spec implies but none of its own reproductions exercise, most likely to bite a user first. Each has a test in the task that owns the code.

1. **A near-duplicate of the backbone reported as a donor.** With the backbone lineage redundant, a backbone copy could win Bootscan windows (or any caller's ranking) as if it were a donor. *Task 11, `test_no_backbone_duplicate_is_reported_as_a_donor` (guard: passes before and after; it pins what the cluster grouping must not break).*
2. **`--min-region` together with `--min-methods 2`.** The length filter must act on the merged span, so an event called by several callers survives a floor that one caller's own shorter tract would not. *Task 13, `test_min_region_is_applied_to_the_merged_span_before_the_agreement_gate` (guard).*
3. **A duplicated donor on a near-identical panel.** The shared clusters are computed on the informative-site scan the HMM uses; Bootscan and the merge work on base-pair windows. The event must still be one row. *Task 11, `test_duplicate_donor_on_a_near_identical_panel_is_one_row` (guard; see Self-review note 3 for the case it does not cover).*
4. **The sign test at its edges far above 1000 sites** (every site against the donor, every site for it, an even split). *Task 3, `test_sign_test_pvalue_one_sided_edges_at_large_n`.*
5. **The trailing gap under informative-site windowing**, where the windows are counted in informative sites, not columns. *Task 14, `test_informative_windows_reach_the_last_informative_site`.*

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `validation/run_regimes.py` | modify | two redundant-lineage regimes (backbone, donor) |
| `tests/unit/test_regime_scoring.py` | modify | unit tests for the new simulator and scorer |
| `src/tessera/recomb/stats.py` | modify | A9: exact sign-test tail in log space |
| `src/tessera/recomb/hmm.py` | modify | A10: segment padding from the median window step; B4: per-window posteriors on `Segment` |
| `src/tessera/recomb/regions.py` | modify | B4: HMM `n_windows` / `posterior_support` over the trimmed span; A2: shared clusters |
| `src/tessera/recomb/ensemble.py` | modify | A11/A8: merged identity and window count over the merged span; A2: cluster-aware grouping; B2: `filter_by_length` |
| `src/tessera/recomb/run.py` | modify | pass the scan and the shared clusters to the merge and callers; B2 filter; a one-state HMM leaves the gate count (B7); `--merge-gap` to `CoverageParams` (B6) |
| `src/tessera/recomb/bootscan.py` | modify | A5: block permutation that keeps every win; A2: support counted per lineage cluster |
| `src/tessera/recomb/maxchi.py` | modify | A7: float chi-square in the null; A3/A6: shallow donors counted in the BH family |
| `src/tessera/recomb/threeseq.py` | modify | A3/A6: no depth floor; anti-diagonal exact DP; A14: terminal breakpoint interval |
| `src/tessera/recomb/geneconv.py` | modify | A13: batched permutations; A14: terminal breakpoint interval |
| `src/tessera/recomb/diagnostics.py` | modify | A13: one gamete matrix at a time; Hudson-Kaplan without a pair list |
| `src/tessera/recomb/clusters.py` | modify | A1: mean, not sum, within a cluster; A2: `lineage_clusters` |
| `src/tessera/recomb/similarity.py` | modify | trailing columns: a final window flush with the end |
| `src/tessera/cli/cmd_recomb.py` | modify | B2/A4: `--min-region` help text; B6: `--merge-gap` help text |
| `tests/unit/test_untestable_callers.py` | create | B7 gate half, B6 coverage half (Task 16) |
| `tests/unit/test_stats.py`, `test_hmm.py`, `test_ensemble.py`, `test_bootscan.py`, `test_maxchi.py`, `test_scaling_guards.py`, `test_diagnostics.py`, `test_threeseq.py`, `test_clusters.py`, `test_similarity.py` | modify | the tests below |
| `tests/integration/test_regimes.py` | modify | pipeline-level regressions on the regime simulators |
| `docs/detection-methods.md`, `validation/README.md`, `CHANGELOG.md` | modify | describe each change where it lands |
| `docs/superpowers/specs/2026-10-01-min-methods-default-decision.md` | modify | A12/A14 addendum |

## Tasks

### Task 1: Branch and baseline  (no audit item; infrastructure)

Nothing is committed in this task. Its output is a directory of harness results that every later gate compares against. `validation/data/` is git-ignored, so the results live there.

**Files:**
- Create (untracked, git-ignored): `validation/data/audit-g/gate.sh`, `validation/data/audit-g/baseline/*.txt`

**Interfaces:**
- Produces: `validation/data/audit-g/<label>/` holding `spec_mm1.txt`, `spec_mm2.txt`, `regimes_mm1.txt`, `regimes_mm2.txt`, `validation.txt`, `validation_regions.tsv`, `hybrids.txt`. Later tasks write the same names under their own label and diff against the previous task's.

- [ ] **Step 1: Branch**

```bash
git switch fix-audit-caller-behaviour && git pull --ff-only
git switch -c fix-audit4-callers
pip install -e ".[dev]"
```

- [ ] **Step 2: Confirm the starting point is green**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: `All checks passed!`, `Success: no issues found in 79 source files`, `757 passed, 5 deselected`. If any of these differs, stop: the baseline is not the one this plan was written against (`4b35e39`).

- [ ] **Step 3: Write the gate script**

Create `validation/data/audit-g/gate.sh`:

```bash
#!/bin/bash
# gate.sh <label> [reps]: run the offline harnesses and save their output under
# validation/data/audit-g/<label>/. Run from the repository root, one at a time.
set -u
out="validation/data/audit-g/$1"
reps="${2:-10}"
mkdir -p "$out"
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"   # appended, never prepended
python validation/run_specificity.py --reps "$reps" --min-methods 1 > "$out/spec_mm1.txt" 2>&1
python validation/run_specificity.py --reps "$reps" --min-methods 2 > "$out/spec_mm2.txt" 2>&1
python validation/run_regimes.py --reps "$reps" > "$out/regimes_mm1.txt" 2>&1
python validation/run_regimes.py --reps "$reps" --min-methods 2 > "$out/regimes_mm2.txt" 2>&1
python validation/run_validation.py > "$out/validation.txt" 2>&1
python - > "$out/validation_regions.tsv" <<'PY'
import csv
from pathlib import Path
for d in sorted(Path("validation/data").iterdir()):
    f = d / "_run" / "recomb_out" / "recombination_regions.tsv"
    if not f.exists() or d.name == "orthopox_example":
        continue
    for r in csv.DictReader(open(f), delimiter="\t"):
        print(d.name, r["minor_parent"], r["major_parent"], r["query_start"], r["query_end"],
              r.get("methods") or "-", r.get("n_windows", ""), r.get("mean_sim_minor", ""),
              r.get("pvalue", ""), r.get("donor_undercovered", ""), r.get("donor_absent", ""),
              sep="\t")
PY
grep -E "TOTAL|detected" "$out/spec_mm1.txt" "$out/spec_mm2.txt"
grep -E "detected|runs with" "$out/regimes_mm1.txt" "$out/regimes_mm2.txt"
tail -9 "$out/validation.txt"
```

```bash
chmod +x validation/data/audit-g/gate.sh
```

The script prints `n_windows`, `mean_sim_minor` and `pvalue` per region because Tasks 3-5 change exactly those columns and nothing else; a diff that touches any other column is a regression.

- [ ] **Step 4: Capture the baseline**

Run: `validation/data/audit-g/gate.sh baseline` (about 25 minutes at 10 replicates).

Expected, at 3 replicates (measured while writing this plan; 10-replicate values will differ within sampling error):
- `spec_mm1.txt`: `TOTAL 7/12 ... 8 hmm=8`; positive control `detected 3/3 | correct donor 3/3 | median breakpoint error 13 bp`.
- `spec_mm2.txt`: `TOTAL 0/12`; positive control `3/3, 3/3, 13 bp`.
- `regimes_mm1.txt`: near_identical `3/3 detected, HMM 3/3, donor 3/3, 378 bp`; near_identical clonal `0/3`; sibling_tract `1001 bp`.
- `validation.txt`: the verdicts recorded on this branch by plan C (6 PASS, 1 FAIL `hcv_clonal_1b`, 1 SKIP), if the datasets and aligners are present.

What to do when something is absent:
- **No aligner on `PATH`**: `run_validation.py` prints `SKIP ... not on PATH` for every dataset. Create the env (`conda env create -f environment.yml`) and append its `bin`. Without it, record "run_validation: not run" for every gate; do not report it as unchanged.
- **Datasets not fetched**: `SKIP ... sequences not present`; run `python validation/fetch.py` (needs `efetch` and network).

- [ ] **Step 5: Capture the hybrids baseline**

```bash
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
python validation/run_hybrids.py > validation/data/audit-g/baseline/hybrids.txt 2>&1
tail -5 validation/data/audit-g/baseline/hybrids.txt
```

Needs network and `mafft`, `skani`, `skDER`; 30-60 minutes. If every row reads `ERROR: Could not fetch the Nextclade dataset index`, there is no network: the hybrids gate cannot be run for any task. Say so in the PR description under "Skipped", and do not merge Tasks 10, 11, 12, 14 or 16 until someone has run it (they change which regions are called).

- [ ] **Step 6: Record the baseline**

Copy the summary lines of each file into a scratch note; they open the PR description in Task 17. No commit.

---

### Task 2: Regime harness: redundant lineages  (infrastructure for A1, A2; gap G1's "cheapest next step")

No harness samples a panel in which one lineage is represented by several near-duplicates, the situation of an untyped recruited panel and of A1 and A2. This task adds two regimes to `validation/run_regimes.py` so that Tasks 10 and 11 have a gate, and so that the measurement is reproducible. It changes no source under `src/`.

**Files:**
- Modify: `validation/run_regimes.py` (module docstring lines 1-34; constants after `SIB_EXTRA = 600` at line 61; new `simulate_redundant` before `# --- scoring` at line 127; new `score_lineage_tract` before `score_clonal` at line 163; `main`, before `return 0` at line 294)
- Modify: `tests/unit/test_regime_scoring.py` (append)
- Modify: `validation/README.md` ("Regimes the other harnesses do not sample", the table at lines 91-95)

**Interfaces:**
- Consumes: `_evolve`, `_as_text`, `_present`, `_methods`, `QUERY`, `scan` (all existing in `run_regimes.py`).
- Produces, in `validation/run_regimes.py`:
  - `RED_LENGTH = 12_000`, `RED_DISTANCE = 0.03`, `RED_COPY_DISTANCE = 0.002`, `RED_COPIES = 3`, `RED_BACKBONE_TRACT = (5_500, 6_500)`, `RED_DONOR_TRACT = (4_000, 7_000)`
  - `simulate_redundant(seed: int, *, backbone_copies: int = 1, donor_copies: int = 1, length: int = RED_LENGTH, distance: float = RED_DISTANCE, tract: tuple[int, int] = RED_BACKBONE_TRACT) -> tuple[dict[str, str], tuple[int, int]]` -- labels `QUERY`, `E`, `F`, `A0..`, `D0..`
  - `score_lineage_tract(rows: list[dict], tract: tuple[int, int], *, lineage: str) -> dict` with keys `detected`, `rows`, `donor_ok`, `hmm`, `bootscan`
  - Tasks 10, 11 and 13 call `simulate_redundant` and `score_lineage_tract` from `tests/integration/test_regimes.py`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_regime_scoring.py`:

```python
# --- redundant lineages ------------------------------------------------------

RED_SMALL = {"length": 3000, "tract": (1000, 1500)}


def test_redundant_simulation_labels_and_copies():
    seqs, tract = rg.simulate_redundant(1, backbone_copies=3, donor_copies=2, **RED_SMALL)
    assert tract == RED_SMALL["tract"]
    assert set(seqs) == {rg.QUERY, "E", "F", "A0", "A1", "A2", "D0", "D1"}
    assert len({len(s) for s in seqs.values()}) == 1
    assert rg.simulate_redundant(1, backbone_copies=3, donor_copies=2, **RED_SMALL)[0] == seqs


def test_redundant_copies_are_near_duplicates_and_lineages_are_distinct():
    seqs, _ = rg.simulate_redundant(2, backbone_copies=2, donor_copies=2, length=20_000)
    n = 20_000
    assert _identity(seqs["A0"], seqs["A1"], 0, n) > 0.99  # ~0.4 % apart
    assert _identity(seqs["D0"], seqs["D1"], 0, n) > 0.99
    assert _identity(seqs["A0"], seqs["D0"], 0, n) < 0.96  # lineages ~6 % apart


def test_redundant_query_carries_the_donor_lineage_over_the_tract():
    seqs, (lo, hi) = rg.simulate_redundant(3, donor_copies=1, length=20_000,
                                           tract=(8000, 12_000))
    query = seqs[rg.QUERY]
    assert _identity(query, seqs["D0"], lo, hi) > _identity(query, seqs["A0"], lo, hi)
    assert _identity(query, seqs["A0"], 0, lo) > _identity(query, seqs["D0"], 0, lo)


def test_score_lineage_tract_counts_rows_and_callers():
    rows = [_row(minor_parent="D1", methods="hmm,bootscan"),
            _row(minor_parent="D0", methods="3seq", query_start="150", query_end="190")]
    got = rg.score_lineage_tract(rows, (100, 200), lineage="D")
    assert got == {"detected": True, "rows": 2, "donor_ok": True, "hmm": True,
                   "bootscan": True}


def test_score_lineage_tract_wrong_lineage_and_miss():
    assert rg.score_lineage_tract([_row(minor_parent="A1")], (100, 200),
                                  lineage="D")["donor_ok"] is False
    miss = rg.score_lineage_tract([_row(query_start="500", query_end="600")], (100, 200),
                                  lineage="D")
    assert miss == {"detected": False, "rows": 0, "donor_ok": False, "hmm": False,
                    "bootscan": False}
```

- [ ] **Step 2: Run them, expected failure**

Run: `pytest tests/unit/test_regime_scoring.py -q`
Expected: `5 failed, 11 passed`, each with `AttributeError: module 'run_regimes' has no attribute 'simulate_redundant'` or `... 'score_lineage_tract'` (observed).

- [ ] **Step 3: Implement**

(a) In `validation/run_regimes.py`, after `SIB_EXTRA = 600` insert:

```python
# redundant: a divergent panel (lineages ~6 % apart) in which one lineage is represented
# by several near-duplicates (~0.4 % apart), as on an untyped recruited panel. Either the
# backbone lineage (a short tract from a single donor genome) or the donor lineage (a
# longer tract, so Bootscan has windows to call) is the redundant one.
RED_LENGTH = 12_000
RED_DISTANCE = 0.03
RED_COPY_DISTANCE = 0.002
RED_COPIES = 3
RED_BACKBONE_TRACT = (5_500, 6_500)
RED_DONOR_TRACT = (4_000, 7_000)
```

(b) Immediately before the line `# --- scoring ---...` insert:

```python
def simulate_redundant(
    seed: int, *, backbone_copies: int = 1, donor_copies: int = 1,
    length: int = RED_LENGTH, distance: float = RED_DISTANCE,
    tract: tuple[int, int] = RED_BACKBONE_TRACT,
) -> tuple[dict[str, str], tuple[int, int]]:
    """A query from backbone lineage ``A`` carrying lineage ``D`` over ``tract``.

    Four lineages (``A``, ``D``, ``E``, ``F``) on a star tree, each ``distance`` from the
    root. ``A`` is represented by ``backbone_copies`` genomes (``A0``, ``A1``, ...) and
    ``D`` by ``donor_copies`` (``D0``, ...), each ``RED_COPY_DISTANCE`` from its lineage.
    The query is ``RED_COPY_DISTANCE`` from ``A`` with ``D``'s sequence over ``tract``.
    Returns ``(sequences, tract)``.
    """
    rng = np.random.default_rng(seed)
    root = rng.integers(0, 4, size=length)
    lineage_a, lineage_d, lineage_e, lineage_f = (_evolve(root, distance, rng) for _ in range(4))
    query = _evolve(lineage_a, RED_COPY_DISTANCE, rng)
    lo, hi = tract
    query[lo:hi] = lineage_d[lo:hi]
    seqs = {QUERY: query, "E": lineage_e, "F": lineage_f}
    for i in range(backbone_copies):
        seqs[f"A{i}"] = _evolve(lineage_a, RED_COPY_DISTANCE, rng)
    for i in range(donor_copies):
        seqs[f"D{i}"] = _evolve(lineage_d, RED_COPY_DISTANCE, rng)
    return _as_text(seqs), tract
```

(c) Immediately before `def score_clonal(` insert:

```python
def score_lineage_tract(rows: list[dict], tract: tuple[int, int], *, lineage: str) -> dict:
    """Detection of a tract whose donor is any genome of ``lineage`` (a label prefix).

    ``rows`` counts the reported regions overlapping the tract (one event should be one
    row); ``donor_ok`` is true when every such row names a genome of the lineage.
    """
    lo, hi = tract
    overlapping = [
        r for r in _present(rows)
        if int(r["query_start"]) < hi and int(r["query_end"]) > lo
    ]
    methods = {m for r in overlapping for m in _methods(r)}
    return {
        "detected": bool(overlapping),
        "rows": len(overlapping),
        "donor_ok": bool(overlapping)
        and all(r["minor_parent"].startswith(lineage) for r in overlapping),
        "hmm": "hmm" in methods,
        "bootscan": "bootscan" in methods,
    }
```

(d) In `main`, replace the final `return 0` (after the `sibling_tract` print) with:

```python
    for name, kwargs, tract in (
        ("redundant_backbone", {"backbone_copies": RED_COPIES}, RED_BACKBONE_TRACT),
        ("redundant_donor", {"donor_copies": RED_COPIES}, RED_DONOR_TRACT),
    ):
        tally: Counter[str] = Counter()
        for rep in range(reps):
            seqs, _ = simulate_redundant(6000 + rep, tract=tract, **kwargs)
            got = score_lineage_tract(run(seqs), tract, lineage="D")
            tally.update(k for k in ("detected", "donor_ok", "hmm", "bootscan") if got[k])
            tally["one_row"] += got["rows"] == 1
        print(f"{name}, {RED_COPIES} near-duplicates, tract {tract[0]}-{tract[1]} from D:")
        print(f"  detected {tally['detected']}/{reps} | one row {tally['one_row']}/{reps} | "
              f"donor lineage {tally['donor_ok']}/{reps} | called by the HMM "
              f"{tally['hmm']}/{reps} | Bootscan among callers {tally['bootscan']}/{reps}")
    return 0
```

(e) In the module docstring, after the `sibling_tract` paragraph (ending "...so nothing marks where the region ends."), add:

```
    redundant_*      a divergent panel (lineages ~6 % apart) in which the backbone
                     lineage, or the donor lineage, is represented by three
                     near-duplicates ~0.4 % apart, as on an untyped recruited panel.
                     Is the tract found, is it one row, and are the HMM and Bootscan
                     among its callers?
```

and change the first line of the docstring from `Opt-in harness: two regimes the other harnesses do not sample.` to `Opt-in harness: regimes the other harnesses do not sample.`

(f) In `validation/README.md`, change "simulates two situations" to "simulates situations" and append two rows to the regime table:

```markdown
| `redundant_backbone` | the backbone lineage is three near-duplicates (~0.4 % apart) and a 1 kb tract comes from a single donor genome: is the HMM still among the callers? (audit 4, item A1) |
| `redundant_donor` | the donor lineage is three near-duplicates and the tract is 3 kb: is the event one row, and are the HMM and Bootscan among its callers? (audit 4, item A2) |
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_regime_scoring.py -q`
Expected: `16 passed`.

- [ ] **Step 5: Full checks, and the new baseline for the redundant regimes**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: all green, `762 passed`.

Then, because the harness now prints two more regimes, re-capture the regimes baseline on the unchanged source:

```bash
python validation/run_regimes.py --reps 10 > validation/data/audit-g/baseline/regimes_mm1.txt 2>&1
python validation/run_regimes.py --reps 10 --min-methods 2 > validation/data/audit-g/baseline/regimes_mm2.txt 2>&1
```

Expected for the new lines (measured at 3 replicates on `4b35e39`, and at 10 replicates of the same simulator with a scratch driver):

| regime (10 reps, `--min-methods 1`) | detected | one row | donor lineage | HMM among callers | Bootscan among callers |
|---|---|---|---|---|---|
| `redundant_backbone` | 10/10 | 10/10 | 10/10 | 3/10 | 1/10 |
| `redundant_donor` | 10/10 | 7/10 | 10/10 | 0/10 | 1/10 |
| control: one genome per lineage, 1 kb tract | 10/10 | 10/10 | 10/10 | 8/10 | 1/10 |

The existing `near_identical` and `sibling_tract` lines must be identical to `baseline` from Task 1 (the simulators and seeds are unchanged).

- [ ] **Step 6: Commit**

```bash
git add validation/run_regimes.py tests/unit/test_regime_scoring.py validation/README.md
git commit -m "$(cat <<'MSG'
Add redundant-lineage regimes to the regime harness

The fourth audit found two defects that appear only when one lineage of
the panel is represented by several near-duplicates (A1, A2), the shape of
an untyped recruited panel, and no harness samples that. run_regimes.py
now simulates it twice: a redundant backbone with a 1 kb tract from a
single donor genome, and a redundant donor lineage with a 3 kb tract.

On the unchanged callers (10 replicates): with three copies of the
backbone the HMM is among the callers in 3/10 runs, against 8/10 when each
lineage is one genome; with three copies of the donor the event is split
into two rows in 3/10 runs, and the HMM and Bootscan are among its callers
in 0/10 and 1/10. These lines are the gate for the fixes that follow.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

---

### Task 3: Exact sign-test tail in log space  (audit item: A9; test gap "the sign test above a thousand trials")

`sign_test_pvalue` (`src/tessera/recomb/stats.py:67-79`) uses the exact binomial tail up to n = 1000 and a continuity-corrected normal approximation above it, whose `1 - cdf` is exactly `0.0` beyond z of about 8.3. Reproduced on `4b35e39`: 900:100 gives 6.7e-162; 909:101, 1500:300 and 1827:200 give `0.0`. So the HMM's strongest regions report p = q = 0, and zeros tie in `ensemble._best`.

**Deviation from the spec's remedy.** The spec suggests lgamma for the exact branch and a Mills-ratio asymptotic for the normal branch. There is no need for a normal branch: the exact tail in log space costs one vectorised pass over `n - a` terms (`log C(n, a)` from `lgamma`, the rest by the ratio `(n - i) / (i + 1)`, summed with log-sum-exp), which is exact at every n and fast at genome scale. A value below the smallest positive normal double (2.2e-308; 1827:200 is about 1e-322) is returned as that value, `P_FLOOR`, never as 0. `sign_test_log_pvalue` returns the unfloored log for callers that need it.

**Files:**
- Modify: `src/tessera/recomb/stats.py` (imports at lines 18-19; `sign_test_pvalue`, lines 67-79)
- Test: `tests/unit/test_stats.py` (append)

**Interfaces:**
- Produces: `stats.P_FLOOR: float` (`sys.float_info.min`); `stats.sign_test_log_pvalue(a: int, b: int) -> float` (natural log, exact, may be below `log(P_FLOOR)`); `stats.sign_test_pvalue(a: int, b: int) -> float` unchanged in signature, now exact for all n and never `0.0` (`1.0` when `a + b == 0`).

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_stats.py`:

```python
def _exact_log_tail(a: int, b: int) -> float:
    """log P(X >= a), X ~ Bin(a + b, 0.5), from exact integer binomial coefficients."""
    n = a + b
    return math.log(sum(math.comb(n, i) for i in range(a, n + 1))) - n * math.log(2.0)


def test_sign_test_pvalue_is_exact_above_a_thousand_sites():
    # Above n = 1000 the test used a normal approximation whose 1 - cdf is exactly 0.0
    # beyond z ~ 8.3, so 909:101 reported p = 0 while 900:100 reported 6.7e-162.
    for a, b in ((909, 101), (1500, 300), (600, 500), (5200, 4800)):
        got = sign_test_pvalue(a, b)
        assert got > 0.0
        assert math.isclose(math.log(got), _exact_log_tail(a, b), rel_tol=1e-9)


def test_sign_test_pvalue_is_monotone_across_the_old_cliff():
    # More discordant sites in the donor's favour can only lower the p-value.
    assert 0.0 < sign_test_pvalue(909, 101) < sign_test_pvalue(900, 100)
    assert 0.0 < sign_test_pvalue(1001, 0) < sign_test_pvalue(1000, 0)


def test_sign_test_pvalue_below_the_float_range_is_floored_not_zero():
    from tessera.recomb import stats

    # The true value, about 1e-322, is below the smallest normal double.
    assert stats.sign_test_pvalue(1827, 200) == stats.P_FLOOR > 0.0
    assert stats.sign_test_log_pvalue(1827, 200) < math.log(stats.P_FLOOR)
    assert math.isclose(stats.sign_test_log_pvalue(1827, 200), _exact_log_tail(1827, 200),
                        rel_tol=1e-9)


def test_sign_test_pvalue_one_sided_edges_at_large_n():
    # Review focus: every site against the donor, or every site for it, far above 1000.
    from tessera.recomb import stats

    assert sign_test_pvalue(0, 5000) == 1.0
    assert stats.sign_test_pvalue(5000, 0) == stats.P_FLOOR
    assert math.isclose(sign_test_pvalue(1, 0), 0.5)
    assert math.isclose(math.log(sign_test_pvalue(2500, 2500)), _exact_log_tail(2500, 2500),
                        rel_tol=1e-9)
```

- [ ] **Step 2: Run them, expected failure**

Run: `pytest tests/unit/test_stats.py -q`
Expected (observed): `4 failed, 6 passed`. `test_sign_test_pvalue_is_exact_above_a_thousand_sites` fails with `assert 0.0 > 0.0` (909:101); `test_sign_test_pvalue_is_monotone_across_the_old_cliff` with `assert 0.0 < 0.0`; the two that use `stats.P_FLOOR` with `AttributeError: module 'tessera.recomb.stats' has no attribute 'P_FLOOR'`.

- [ ] **Step 3: Implement**

In `src/tessera/recomb/stats.py` change the imports

```python
import math
from statistics import NormalDist
```

to

```python
import math
import sys
from statistics import NormalDist

import numpy as np
```

and replace the whole of `sign_test_pvalue` with:

```python
# The smallest positive normal double. A sign-test p-value below it is reported as this
# value rather than as 0.0: zero is a value no test can produce, and zeros tie.
P_FLOOR = sys.float_info.min


def sign_test_log_pvalue(a: int, b: int) -> float:
    """Natural log of the one-sided sign-test p-value ``P(X >= a)``, ``X ~ Bin(a+b, 0.5)``.

    Exact for every ``n``, computed in log space: ``log C(n, a)`` from ``lgamma``, the
    following terms by the ratio ``C(n, i+1) / C(n, i) = (n - i) / (i + 1)``, summed with
    log-sum-exp. Nothing is subtracted from one, so the value keeps its magnitude however
    small it is.
    """
    n = a + b
    if n == 0 or a <= 0:
        return 0.0
    if a > n:
        return float("-inf")
    i = np.arange(a, n, dtype=np.float64)  # ratios for i = a .. n-1 give terms a+1 .. n
    first = math.lgamma(n + 1) - math.lgamma(a + 1) - math.lgamma(n - a + 1)
    log_terms = np.concatenate(([first], first + np.cumsum(np.log((n - i) / (i + 1)))))
    top = float(log_terms.max())
    return top + math.log(float(np.exp(log_terms - top).sum())) - n * math.log(2.0)


def sign_test_pvalue(a: int, b: int) -> float:
    """One-sided sign-test p-value: ``P(X >= a)`` for ``X ~ Binomial(a+b, 0.5)``.

    The probability of seeing at least ``a`` discordant sites favouring the donor
    if the two parents were equally good (the null). Small means a real preference.

    Exact at every ``n`` (see :func:`sign_test_log_pvalue`). An earlier version switched
    to a normal approximation above 1000 sites, whose ``1 - cdf`` is exactly ``0.0`` beyond
    z of about 8.3, so the strongest regions were reported with p = 0. A value below
    :data:`P_FLOOR` is returned as ``P_FLOOR``, never as zero.
    """
    if a + b == 0:
        return 1.0
    log_p = sign_test_log_pvalue(a, b)
    return min(1.0, max(P_FLOOR, math.exp(log_p))) if log_p > float("-inf") else P_FLOOR
```

`_NORMAL` is still used by `_z_two_sided` / `_z_one_sided`; leave it.

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_stats.py tests/unit/test_regions.py tests/unit/test_hmm.py -q`
Expected: all pass. The existing `test_sign_test_pvalue` (`0.0 < sign_test_pvalue(900, 700) < 0.001`) still holds; its value moves from the normal approximation to the exact tail.

- [ ] **Step 5: Full checks and regression gate**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q` (expected `766 passed`), then `validation/data/audit-g/gate.sh task3`.

Acceptance: every `TOTAL`, positive-control and regime line identical to `baseline` (measured: identical at 3 replicates). In `diff validation/data/audit-g/baseline/validation_regions.tsv validation/data/audit-g/task3/validation_regions.tsv` only the `pvalue` column of rows whose methods include `hmm` may change, and only for regions with more than 1000 discordant sites. A p-value that moved from the normal approximation to the exact value can in principle cross `--alpha` after BH (600:500 moves from 0.001418 to 0.001410); if a row appears or disappears, record it in the PR -- it is the exact value deciding, not a regression -- and keep the commit only if the hybrids gate shows no PASS-to-FAIL.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/stats.py tests/unit/test_stats.py
git commit -m "$(cat <<'MSG'
Compute the HMM sign-test tail exactly, in log space

Above 1000 discordant sites sign_test_pvalue switched to a normal
approximation whose 1 - cdf is exactly 0.0 beyond z ~ 8.3. A region with
909:101 discordant sites therefore reported p = q = 0 while one with
900:100 reported 6.7e-162: the strongest evidence carried the least
informative number, and zero q-values tie in the ensemble's ordering.

The tail is now exact at every n: log C(n, a) from lgamma, the following
terms by the ratio (n - i) / (i + 1), summed with log-sum-exp. Nothing is
subtracted from one. Below the smallest normal double (2.2e-308) the value
is reported as that floor, never as zero; sign_test_log_pvalue keeps the
unfloored log. Values that were already exact (n <= 1000) agree to 1e-12.

Harness: <paste the spec / regime lines and the validation_regions diff
from validation/data/audit-g/task3/>.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

Replace the `<paste ...>` line with the measured lines before committing.

---

### Task 4: HMM segment spacing and the statistics of a trimmed span  (audit items: A10; B4, code half; test gap "the HMM segment span under a window step smaller than the window width")

Two defects in how an HMM segment is turned into a region.

1. `hmm._path_to_segments` (`src/tessera/recomb/hmm.py:182`) pads a segment's ends by `positions[1] - positions[0]`. Under informative-site windowing that is one unrepresentative gap applied genome-wide (33 columns on the audit's panel where the median is about 237). `clusters.median_window_step` exists for exactly this and `tests/unit/test_clusters.py:156` names this expression as what the old code used; the HMM was never migrated. Under base-pair windowing every gap equals the step, so nothing changes there.
2. After plan C's trimming (`regions.py:327`), `n_windows` and `posterior_support` (`regions.py:338,341`) still describe the untrimmed segment: on the sibling-tract regime a 585-column region reports `n_windows = 116`. Both are now restated over the segment's windows that overlap the reported span. The sign test and its p-value stay on the segment, as plan C decided (its deviation 1).

**Evaluated and not chosen: padding by half the window width.** The spec states the padding "is the window step, not half the window width". A candidate that padded base-pair segments by `window_size // 2` and informative-site segments by their window spans was measured (see "What was measured"); it is not part of this plan for the reason given there. This task only replaces the single sampled gap by the median one.

**Files:**
- Modify: `src/tessera/recomb/hmm.py` (imports at lines 19-24; `Segment`, lines 31-45; `_path_to_segments`, lines 179-211)
- Modify: `src/tessera/recomb/regions.py` (imports at line 30; `_call_regions_hmm` lines 331-341; new `_span_windows` before `_tract_span`)
- Modify: `docs/detection-methods.md` (the "reported span is not the segment" paragraph, lines 73-82)
- Test: `tests/unit/test_hmm.py` (append), `tests/integration/test_regimes.py` (append)

**Interfaces:**
- Consumes: `clusters.median_window_step(result) -> int`, `clusters._window_bounds(result, window_size) -> list[tuple[int, int]]`.
- Produces: `Segment.window_posteriors: list[float]` (posterior of the segment's state at each of its windows, default `[]`); `regions._span_windows(work, window_size, seg, lo, hi) -> tuple[int, float]`. HMM `Region.n_windows` = segment windows overlapping the reported span; `Region.posterior_support` = their mean posterior. Task 5 restates `n_windows` again for merged rows.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_hmm.py`:

```python
def _uneven_two_state_result():
    """As _two_state_result, but the first two windows sit 3 columns apart, as the dense
    5' end of an informative-site scan does; every later gap is 50 columns."""
    result = _two_state_result()
    positions = [22, 25] + [25 + 50 * i for i in range(1, len(result.positions) - 1)]
    result.positions = positions
    result.query_positions = positions
    return result


def test_segment_ends_use_the_median_window_spacing():
    from tessera.recomb.clusters import median_window_step

    result = _uneven_two_state_result()
    segments, _ = segment_query(result, jump_rate=1e-3)
    assert [s.state for s in segments] == ["refA", "refB"]
    second = segments[1]
    step = median_window_step(result)
    assert step == 50
    # The first gap (3 columns) is not the spacing of the scan.
    assert second.msa_start == result.positions[second.start_window] - step
    assert segments[0].msa_end == result.positions[segments[0].end_window] + step


def test_segment_carries_its_per_window_posteriors():
    segments, _ = segment_query(_two_state_result(), jump_rate=1e-3)
    for seg in segments:
        assert len(seg.window_posteriors) == seg.n_windows
        assert abs(sum(seg.window_posteriors) / seg.n_windows - seg.mean_posterior) < 1e-4
```

Append to `tests/integration/test_regimes.py`:

```python
def test_hmm_window_count_and_posterior_describe_the_trimmed_span(tmp_path) -> None:
    # The sibling tract's segment runs ~11 kb; the reported span is ~600 bp. n_windows
    # counted the segment's windows (over 100) beside a 600 bp region.
    seqs, (lo, hi) = rg.simulate_sibling_tract(5000)
    msa = write_fasta(tmp_path / "sib.fasta", seqs)
    result = compute_similarity(str(msa), rg.QUERY, WINDOW, STEP)
    regions, _ = _hmm_regions(result)
    region = regions[0]
    # windows overlapping a span of L columns: at most (L + window) / step + 1
    assert region.n_windows <= (region.length_msa + WINDOW) // STEP + 1
    assert 0.0 <= region.posterior_support <= 1.0
```

- [ ] **Step 2: Run them, expected failure**

Run: `pytest tests/unit/test_hmm.py tests/integration/test_regimes.py -q`
Expected (observed): 3 failed.
- `test_segment_ends_use_the_median_window_spacing`: `assert 472 == (475 - 50)` (padded by the first gap, 3).
- `test_segment_carries_its_per_window_posteriors`: `AttributeError: 'Segment' object has no attribute 'window_posteriors'`.
- `test_hmm_window_count_and_posterior_describe_the_trimmed_span`: `assert 116 <= (((585 + 1000) // 100) + 1)`.

- [ ] **Step 3: Implement**

(a) `src/tessera/recomb/hmm.py`: change `from dataclasses import dataclass` to `from dataclasses import dataclass, field`, and add `from .clusters import median_window_step` above `from .similarity import WindowSimilarity` (`clusters` imports only `similarity`, so there is no cycle). In `Segment`, after `breakpoint_hi: int | None` add:

```python
    # Posterior of this segment's state at each of its windows (start_window ..
    # end_window), so a statistic can be restated over part of the segment.
    window_posteriors: list[float] = field(default_factory=list)
```

In `_path_to_segments` replace

```python
    half = (result.positions[1] - result.positions[0]) if len(result.positions) > 1 else 0
```

with

```python
    # A segment's ends are padded by the typical window spacing. Under informative-site
    # windowing the gap between the first two windows is not typical (it tracks local
    # site density), so the median gap is used, as everywhere else that needs the step.
    half = median_window_step(result) if len(result.positions) > 1 else 0
```

and in the `Segment(...)` call replace `            breakpoint_lo=bp_lo, breakpoint_hi=bp_hi,` with

```python
            breakpoint_lo=bp_lo, breakpoint_hi=bp_hi,
            window_posteriors=[float(x) for x in post[lo_i : hi_i + 1, state_idx]],
```

(b) `src/tessera/recomb/regions.py`: change the import to `from .clusters import _window_bounds, all_singletons, cluster_references, clustered_view`. In `_call_regions_hmm`, after the two `region_identity` lines add

```python
        # The window count and the posterior describe the reported span too, not the
        # segment it was trimmed from: the segment's windows that overlap the span.
        n_windows, posterior = _span_windows(work, window_size, seg, lo, hi)
```

and in the `Region(...)` call replace `n_windows=seg.n_windows,` by `n_windows=n_windows,` and `posterior_support=seg.mean_posterior,` by `posterior_support=posterior,`. Insert before `def _tract_span(`:

```python
def _span_windows(
    work: WindowSimilarity, window_size: int, seg: Segment, lo: int, hi: int
) -> tuple[int, float]:
    """``(n_windows, mean posterior)`` over the segment's windows overlapping ``[lo, hi)``.

    Falls back to the whole segment when no window of it overlaps the span (a span
    can only lie outside its segment's windows at an alignment end).
    """
    bounds = _window_bounds(work, window_size)
    idx = [
        i for i in range(seg.start_window, seg.end_window + 1)
        if bounds[i][0] < hi and bounds[i][1] > lo
    ]
    if not idx or not seg.window_posteriors:
        return seg.n_windows, seg.mean_posterior
    post = [seg.window_posteriors[i - seg.start_window] for i in idx]
    return len(idx), round(float(sum(post) / len(post)), 4)
```

(c) `docs/detection-methods.md`: at the end of the paragraph that begins "**The reported span is not the segment.**" (it ends "`mean_sim_minor` / `mean_sim_major` are identity over the reported columns."), add:

```markdown
`n_windows` and `posterior_support` are restated over the reported span as well: the
segment's windows that overlap it, and their mean posterior for the donor state. A
segment's ends are padded by the median spacing between windows, which under
informative-site windowing is not the spacing of the first two windows.
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_hmm.py tests/unit/test_clusters.py tests/integration/test_regimes.py -q`
Expected: all pass.

- [ ] **Step 5: Full checks and regression gate**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q` (expected `769 passed`), then `validation/data/audit-g/gate.sh task4`.

Acceptance (measured identical at 3 replicates, cumulatively with Tasks 3, 5-9):
- `spec_mm1.txt`, `spec_mm2.txt`, the positive control and every `regimes_*` line identical to `task3`. The specificity panels run under base-pair windowing, where the median gap is the step: a difference there means the change leaked into base-pair mode. `near_identical` runs under informative-site windowing; its breakpoint error must not rise (measured 378 bp before and after).
- `validation_regions.tsv`: the same rows; only `n_windows` of `hmm` rows may change, and coordinates of `hmm` rows only on datasets analysed under informative-site windowing.
- Hybrids: no case from PASS to FAIL.

If a row appears or disappears: the median padding changed which segments clear the `--min-region` segment filter. Revert, record the rows, go to Task 5.

**If plan F Task 5 has landed:** it adds to the same paragraph "`n_windows` and `posterior` describe the segment too, not the reported span ...". Replace that sentence with the one above rather than keeping both. If plan F lands later, its Task 5 executor drops that sentence.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/hmm.py src/tessera/recomb/regions.py docs/detection-methods.md tests/unit/test_hmm.py tests/integration/test_regimes.py
git commit -m "$(cat <<'MSG'
Pad HMM segments by the median window step; restate n_windows and posterior

_path_to_segments padded a segment's ends by the gap between the first two
windows. Under informative-site windowing that gap tracks local site
density and is not the scan's spacing (33 columns against a median of 237
on one audit panel); median_window_step was written for exactly this and
the HMM had not been migrated. Base-pair windowing is unchanged.

After an HMM region's span is trimmed, n_windows and posterior_support
still described the untrimmed segment: a 585-column region reported 116
windows. Both are now taken over the segment's windows that overlap the
reported span. The sign test stays on the segment, as before.

Not chosen: padding by half the window width (see the plan's
measurements).

Harness: <paste the lines from validation/data/audit-g/task4/>.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

---

### Task 5: A merged row's identity and window count over its own span  (audit items: A11, A8; test gap "a merged region's identity against its own span")

`ensemble._merge_group` (`src/tessera/recomb/ensemble.py:88-131`) takes the union of its members' spans (lines 107-110) but keeps the most significant member's `mean_sim_minor`, `mean_sim_major` and `margin` (lines 112-114): an HMM member over 1000-2000 at 0.99 merged with a 3SEQ member over 1900-6000 reports 0.99 over 1000-6000. And `n_windows = max(...)` (line 111) is a maximum over incommensurable units: the HMM and Bootscan count windows, 3SEQ and MaxChi count discriminating sites, GENECONV a run length. Both are restated over the merged span from the base-pair scan. `consensus_regions` is called for every run, single-caller ones included (`run.py:384`), so every row of `recombination_regions.tsv` then carries the same two definitions.

`flag_undercovered_regions` (`run.py`, after the merge) reads `mean_sim_minor`, so `donor_undercovered` can change on a merged row whose union span is less similar to the donor than its best member was. That is the intended consequence and is why the gate diffs that column.

**Files:**
- Modify: `src/tessera/recomb/ensemble.py` (imports at lines 19-20; `_merge_group`, lines 88-131; `consensus_regions`, lines 181-203)
- Modify: `src/tessera/recomb/run.py` (the `consensus_regions(...)` call, lines 384-388)
- Modify: `docs/detection-methods.md` (output table, the `recombination_regions.tsv` row at line 338)
- Test: `tests/unit/test_ensemble.py` (append)

**Interfaces:**
- Consumes: `clusters._window_bounds`, `similarity.region_identity(rows, query, ref, start, end) -> float`.
- Produces: `consensus_regions(per_method, *, major, rmin_intervals=None, lineage_map=None, scan: tuple[WindowSimilarity, int] | None = None)`; `_merge_group(group, major, scan=None)`. With `scan` the merged row's `mean_sim_minor` / `mean_sim_major` are identity over `[msa_start, msa_end)`, `margin` their difference, and `n_windows` the number of base-pair windows overlapping that span. Without `scan`, the old behaviour (guarded by a test). Task 11 adds a `clusters` keyword to the same function.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_ensemble.py`:

```python
def _scan_result():
    """A 100-column alignment: B matches the query on [10, 40), A on [40, 100). Five
    20-column base-pair windows centred at 10, 30, 50, 70, 90."""
    import numpy as np

    from tessera.recomb.similarity import WindowSimilarity

    width = 100
    query = np.frombuffer(b"A" * width, dtype=np.uint8)
    minor = np.frombuffer(b"C" * 10 + b"A" * 30 + b"C" * 60, dtype=np.uint8)
    major = np.frombuffer(b"A" * 10 + b"C" * 30 + b"A" * 60, dtype=np.uint8)
    positions = [10, 30, 50, 70, 90]
    return WindowSimilarity(
        positions=positions, query_positions=positions,
        similarities={"A": [1.0] * 5, "B": [1.0] * 5}, query="q", width=width,
        query_cumulative=np.arange(width + 1),
        rows={"q": query, "A": major, "B": minor},
    )


def test_merged_identity_is_measured_over_the_merged_span() -> None:
    # The HMM member covers [10, 40) where B is identical to the query; the 3SEQ member
    # covers [30, 80). The merged row spans [10, 80) and must not report the HMM
    # member's 1.0 over a span of which it describes less than half.
    hmm = mk("B", 10, 40, "hmm", qvalue=1e-9)
    hmm.mean_sim_minor, hmm.mean_sim_major = 1.0, 0.0
    seq = mk("B", 30, 80, "3seq", qvalue=1e-3)
    merged, _ = consensus_regions({"hmm": [hmm], "3seq": [seq]}, major="A",
                                  scan=(_scan_result(), 20))
    assert len(merged) == 1
    region = merged[0]
    assert (region.msa_start, region.msa_end) == (10, 80)
    assert region.mean_sim_minor == round(30 / 70, 4)
    assert region.mean_sim_major == round(40 / 70, 4)
    assert region.margin == round(30 / 70 - 40 / 70, 4)


def test_merged_window_count_is_in_one_unit() -> None:
    # 3SEQ counts discriminating sites in n_windows, the HMM counts windows; the merged
    # row restates it as base-pair windows overlapping the merged span.
    hmm = mk("B", 10, 40, "hmm", qvalue=1e-9)
    seq = mk("B", 30, 80, "3seq", qvalue=1e-3)
    seq.n_windows = 37  # discriminating sites, not windows
    merged, _ = consensus_regions({"hmm": [hmm], "3seq": [seq]}, major="A",
                                  scan=(_scan_result(), 20))
    # windows [0,20) [20,40) [40,60) [60,80) [80,100): four overlap [10, 80)
    assert merged[0].n_windows == 4


def test_merge_without_a_scan_keeps_the_best_members_values() -> None:
    # Guard: callers that do not pass the scan (tests, the API) get the old behaviour.
    hmm = mk("B", 10, 40, "hmm", qvalue=1e-9)
    seq = mk("B", 30, 80, "3seq", qvalue=1e-3)
    merged, _ = consensus_regions({"hmm": [hmm], "3seq": [seq]}, major="A")
    assert merged[0].mean_sim_minor == hmm.mean_sim_minor
```

`test_merge_without_a_scan_keeps_the_best_members_values` is a guard: it passes before and after.

- [ ] **Step 2: Run them, expected failure**

Run: `pytest tests/unit/test_ensemble.py -q`
Expected (observed): 2 failed, both `TypeError: consensus_regions() got an unexpected keyword argument 'scan'`.

- [ ] **Step 3: Implement**

(a) In `src/tessera/recomb/ensemble.py`, change the imports to

```python
from .clusters import _window_bounds
from .regions import CALLERS, Region
from .similarity import WindowSimilarity, region_identity
from .typing import LineageMap, lineage_of
```

(`clusters` and `similarity` do not import `ensemble`; no cycle.)

(b) Replace the signature and docstring of `_merge_group`

```python
def _merge_group(group: list[Region], major: str | None) -> Region:
    """One consensus region from a set of overlapping same-minor-parent regions."""
```

with

```python
def _merge_group(
    group: list[Region], major: str | None,
    scan: tuple[WindowSimilarity, int] | None = None,
) -> Region:
    """One consensus region from a set of overlapping same-minor-parent regions.

    With ``scan`` (the base-pair similarity result and its window size) the merged row's
    identities and window count are restated over the merged span; without it they are
    the most significant member's, as before.
    """
```

(c) In `_merge_group`, replace from `    relabel = bool(major) and major != best.minor_parent` down to and including `        margin=best.margin,` with:

```python
    relabel = bool(major) and major != best.minor_parent
    major_label = major if relabel and major else best.major_parent
    msa_start = min(r.msa_start for r in group)
    msa_end = max(r.msa_end for r in group)
    n_windows = max(r.n_windows for r in group)
    sim_minor, sim_major, margin = best.mean_sim_minor, best.mean_sim_major, best.margin
    if scan is not None:
        # The span is the union of the members, so the identities must be measured over
        # that union: the best member's own values describe its own (often much shorter)
        # span. The window count is restated in one unit -- base-pair scan windows
        # overlapping the span -- because the callers count different things (windows,
        # discriminating sites, run length) and a maximum across them has no unit.
        result, window_size = scan
        q = result.query
        raw_minor = region_identity(result.rows, q, best.minor_parent, msa_start, msa_end)
        raw_major = region_identity(result.rows, q, major_label, msa_start, msa_end)
        sim_minor, sim_major = round(raw_minor, 4), round(raw_major, 4)
        margin = round(raw_minor - raw_major, 4)
        n_windows = sum(
            1 for lo, hi in _window_bounds(result, window_size) if lo < msa_end and hi > msa_start
        )
    return Region(
        minor_parent=best.minor_parent,
        major_parent=major_label,
        msa_start=msa_start,
        msa_end=msa_end,
        query_start=q_start,
        query_end=q_end,
        n_windows=n_windows,
        mean_sim_minor=sim_minor,
        mean_sim_major=sim_major,
        margin=margin,
```

(the rest of the `Region(...)` call, from `posterior_support=best.posterior_support,`, is unchanged).

(d) In `consensus_regions`, add the keyword after `lineage_map: LineageMap | None = None,`:

```python
    scan: tuple[WindowSimilarity, int] | None = None,
```

append to its docstring

```
    ``scan`` is ``(base-pair result, window size)``. When given, each merged row's
    ``mean_sim_minor`` / ``mean_sim_major`` / ``margin`` are identity over the merged span
    and ``n_windows`` counts the base-pair windows overlapping it.
```

and change `        region = _merge_group(group, major)` to `        region = _merge_group(group, major, scan)`.

(e) In `src/tessera/recomb/run.py`, in the `consensus_regions(` call add after `lineage_map=lineage_map,`:

```python
        scan=(bp_result, params.window_size),
```

(f) In `docs/detection-methods.md`, in the output table row for `recombination_regions.tsv`, after "the `test` and `statistic` that produced them," insert: "`n_windows` (base-pair scan windows overlapping the region, the same unit for every caller), mean similarities (identity over the region's own span, for a merged row the whole merged span),". Remove the now-duplicated "mean similarities," that follows.

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_ensemble.py tests/integration -q`
Expected: all pass.

- [ ] **Step 5: Full checks and regression gate**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q` (expected `772 passed`), then `validation/data/audit-g/gate.sh task5`.

Acceptance (measured identical at 3 replicates, cumulatively with Tasks 3, 4, 6-9):
- every `TOTAL`, positive-control and regime line identical to `task4` (`score_*` read coordinates, donors and methods only).
- `validation_regions.tsv`: the same rows in the same order with the same coordinates; `n_windows` and `mean_sim_minor` may change on every row; `donor_undercovered` may change only on rows with more than one method. Write each `donor_undercovered` change into the PR with the old and new `mean_sim_minor`.
- Hybrids: no case from PASS to FAIL.

If coordinates, donors or methods change: revert and record; this task must not change which regions are reported.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/ensemble.py src/tessera/recomb/run.py docs/detection-methods.md tests/unit/test_ensemble.py
git commit -m "$(cat <<'MSG'
Restate a merged region's identity and window count over its own span

The ensemble merge took the union of its members' spans but kept the most
significant member's identities: an HMM member at 0.99 over 1000-2000
merged with a 3SEQ member over 1900-6000 reported 0.99 over 1000-6000,
an identity measured over a fifth of the span. It also took the maximum
n_windows across callers that count different things (windows,
discriminating sites, run length).

With the base-pair scan passed in, a merged row now reports identity over
the merged span and the number of base-pair windows overlapping it. Every
run goes through the merge, so every row uses the same definitions.
donor_undercovered reads mean_sim_minor and can change on merged rows.

Harness: <paste the lines and the validation_regions diff, with each
donor_undercovered change, from validation/data/audit-g/task5/>.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

---

### Task 6: Bootscan's block permutation keeps every win  (audit item: A5; test gap "the permuted win count when the window count is not a multiple of the block")

`bootscan.run_pvalue` (`src/tessera/recomb/bootscan.py:49-80`) pads the win vector to whole blocks, shuffles the blocks and then cuts each row back to `n_win` (line 73). Whenever `n_win` is not a multiple of `block`, wins that a shuffle moves past the end are dropped, so the null has fewer wins than the observed vector and shorter runs: an anti-conservative p-value, although the docstring says the win count is held fixed. The clearest case: 11 windows, block 10, every window a win -- every arrangement has a run of 11, so p must be 1, and the current code returns 0.49. The fix removes the padding cells wherever they land instead of cutting the row.

**Files:**
- Modify: `src/tessera/recomb/bootscan.py` (`run_pvalue` lines 62-73; new `_block_permutations` after `run_pvalue`)
- Test: `tests/unit/test_bootscan.py` (append)

**Interfaces:**
- Produces: `bootscan._block_permutations(wins: np.ndarray, block: int, perms: int, rng) -> np.ndarray` of shape `(perms, wins.size)`, every row holding `wins.sum()` wins. `run_pvalue` unchanged in signature.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_bootscan.py`:

```python
# --- the run-length null and lineage clusters --------------------------------

def test_block_permutation_keeps_every_win_when_the_last_block_is_partial() -> None:
    # 11 windows, block 10: the second block is one window long. Every window is a win,
    # so every arrangement has a run of 11 and p must be 1. Cutting the shuffled vector
    # back to 11 cells dropped ten wins whenever the short block came first.
    from tessera.recomb.bootscan import run_pvalue

    assert run_pvalue(np.ones(11, dtype=bool), 11, block=10) == 1.0


def test_block_permutations_hold_the_win_count_and_length() -> None:
    from tessera.recomb.bootscan import _block_permutations

    wins = np.zeros(103, dtype=bool)
    wins[[0, 1, 2, 50, 51, 100, 101, 102]] = True
    perms = _block_permutations(wins, 10, 400, np.random.default_rng(0))
    assert perms.shape == (400, 103)
    assert set(perms.sum(axis=1)) == {8}
```

- [ ] **Step 2: Run them, expected failure**

Run: `pytest tests/unit/test_bootscan.py -q`
Expected (observed): 2 failed -- `assert 0.49325337331334335 == 1.0`, and `ImportError: cannot import name '_block_permutations'`.

- [ ] **Step 3: Implement**

In `run_pvalue`, delete the line `    n_blocks = int(np.ceil(n_win / block))` and replace

```python
    rng = np.random.default_rng(seed)
    padded = np.zeros(n_blocks * block, dtype=bool)
    padded[:n_win] = wins
    blocks = padded.reshape(n_blocks, block)
    order = np.argsort(rng.random((perms, n_blocks)), axis=1)
    shuffled = blocks[order].reshape(perms, n_blocks * block)[:, :n_win]

```

with

```python
    shuffled = _block_permutations(wins, block, perms, np.random.default_rng(seed))
```

and insert after `run_pvalue`:

```python
def _block_permutations(wins: np.ndarray, block: int, perms: int, rng) -> np.ndarray:
    """``(perms, n_win)`` block permutations of ``wins``, every row holding its wins.

    The vector is padded to whole blocks, the blocks are shuffled, and the padding cells
    are then removed wherever they landed. Cutting the shuffled row back to ``n_win``
    instead (the earlier approach) dropped any wins that a shuffle moved past the end
    whenever ``n_win`` was not a multiple of ``block``, so the null had fewer wins than
    the observed vector and its runs were shorter than they should be.
    """
    n_win = int(wins.size)
    n_blocks = int(np.ceil(n_win / block))
    padded = np.zeros(n_blocks * block, dtype=bool)
    padded[:n_win] = wins
    real = np.zeros(n_blocks * block, dtype=bool)
    real[:n_win] = True
    order = np.argsort(rng.random((perms, n_blocks)), axis=1)
    values = padded.reshape(n_blocks, block)[order].reshape(perms, -1)
    keep = real.reshape(n_blocks, block)[order].reshape(perms, -1)
    # Every row has exactly n_win real cells, so the boolean selection reshapes cleanly.
    return values[keep].reshape(perms, n_win)
```

When `n_win` is a multiple of `block` there is no padding and the rows are exactly what the old code produced from the same random stream, so those p-values are unchanged.

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_bootscan.py -q`
Expected: all pass.

- [ ] **Step 5: Full checks and regression gate**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q` (expected `774 passed`), then `validation/data/audit-g/gate.sh task6`.

Acceptance: every line identical to `task5` (measured identical at 3 replicates, cumulatively). Bootscan p-values rise slightly where `n_win % block != 0`; a Bootscan-only row could therefore disappear from `validation_regions.tsv`. If one does, record it with its old p and q -- it was called on an anti-conservative null -- and keep the commit unless the hybrids gate shows a PASS-to-FAIL.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/bootscan.py tests/unit/test_bootscan.py
git commit -m "$(cat <<'MSG'
Hold the win count fixed in Bootscan's block permutation

run_pvalue padded the window vector to whole blocks, shuffled the blocks
and cut each row back to n_win. Whenever n_win is not a multiple of the
block, wins moved past the end were dropped, so the null had fewer wins
and shorter runs than the data: an anti-conservative p-value, although the
docstring states the win count is fixed. With 11 windows, block 10 and
every window a win, p must be 1; it was 0.49.

The padding cells are now removed wherever they land. When n_win is a
multiple of the block the rows are unchanged.

Harness: <paste from validation/data/audit-g/task6/>.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

---

### Task 7: MaxChi's permutation null without int64 overflow  (audit item: A7, SUSPECTED)

`maxchi._exceedances` (`src/tessera/recomb/maxchi.py:80-87`) computes the chi-square denominator `r1 * r2 * tot_major * tot_minor` in int64, which wraps above roughly 110 000 discriminating sites with balanced counts; the guard `denom > 0` then sets that permutation's statistic to zero, while the observed statistic (Python integers, `_chi2_2x2`) does not wrap. The overflow is confirmed in isolation; it was not reached through `_exceedances`, so this is a guard plus a widening, as the spec allows for a SUSPECTED item. The computation moves into a small function so it can be tested at the sizes where it matters.

**Files:**
- Modify: `src/tessera/recomb/maxchi.py` (`_exceedances` lines 80-88; new `_chi2_batch` after it)
- Test: `tests/unit/test_maxchi.py` (append)

**Interfaces:**
- Produces: `maxchi._chi2_batch(in_major, in_minor, tot_major, tot_minor) -> np.ndarray` (float64; equals `_chi2_2x2` element-wise). `_exceedances` unchanged in signature and in value below the overflow.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_maxchi.py`:

```python
def test_vectorised_chi_square_does_not_wrap_at_bacterial_scale() -> None:
    # 120 000 discriminating sites, balanced: the product of the four margins is about
    # 1.3e19, above the int64 maximum. Computed in int64 it wrapped negative and the
    # permutation's statistic was forced to zero.
    from tessera.recomb.maxchi import _chi2_batch

    in_major = np.array([40_000, 30_000], dtype=np.int64)
    in_minor = np.array([20_000, 30_000], dtype=np.int64)
    tot = np.array([60_000, 60_000], dtype=np.int64)
    got = _chi2_batch(in_major, in_minor, tot, tot)
    want = [_chi2_2x2(40_000, 20_000, 20_000, 40_000), _chi2_2x2(30_000, 30_000, 30_000, 30_000)]
    assert np.allclose(got, want, rtol=1e-12)
    assert got[0] > 0.0
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_maxchi.py -q`
Expected (observed): 1 failed, `ImportError: cannot import name '_chi2_batch'`. (The arithmetic it pins: `60_000 * 60_000 * 60_000 * 60_000` as int64 is negative.)

- [ ] **Step 3: Implement**

In `_exceedances`, replace from `    tot_major, tot_minor = cum_major[:, -1], cum_minor[:, -1]` to the end of the function with:

```python
    chi2 = _chi2_batch(in_major, in_minor, cum_major[:, -1], cum_minor[:, -1])
    return int(np.count_nonzero(chi2 >= observed))
```

and insert after `_exceedances`:

```python
def _chi2_batch(
    in_major: np.ndarray, in_minor: np.ndarray, tot_major: np.ndarray, tot_minor: np.ndarray
) -> np.ndarray:
    """Vectorised :func:`_chi2_2x2` for tract-versus-rest tables.

    Every factor is converted to float64 before it is multiplied. The denominator is a
    product of four counts; in int64 it wraps above about 110 000 discriminating sites
    with balanced counts, turns negative, and the positivity guard below then set that
    permutation's statistic to zero -- while the observed statistic, computed with Python
    integers, did not wrap. float64 holds the product without wrapping at any size a
    genome alignment reaches.
    """
    a = in_major.astype(np.float64)
    b = in_minor.astype(np.float64)
    c = tot_major.astype(np.float64) - a
    d = tot_minor.astype(np.float64) - b
    denom = (a + b) * (c + d) * (a + c) * (b + d)
    det = a * d - b * c
    n = a + b + c + d
    return np.where(denom > 0, n * det**2 / np.where(denom > 0, denom, 1.0), 0.0)
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_maxchi.py tests/unit/test_scaling_guards.py -q`
Expected: all pass, including the existing `test_permutation_null_uses_the_same_tract_rule_as_the_observed_walk` (the null's values are unchanged below the overflow).

- [ ] **Step 5: Full checks and regression gate**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q` (expected `775 passed`), then `validation/data/audit-g/gate.sh task7`.

Acceptance: every file identical to `task6`, including `validation_regions.tsv` (no harness alignment approaches 110 000 discriminating sites; measured identical). Any difference: revert and record.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/maxchi.py tests/unit/test_maxchi.py
git commit -m "$(cat <<'MSG'
Compute MaxChi's permutation chi-square in float64

The vectorised null multiplied four int64 margins. Above about 110 000
discriminating sites with balanced counts the product wraps negative, the
positivity guard then set that permutation's statistic to zero, and the
observed statistic, computed with Python integers, did not wrap: the two
sides of the test would disagree and p would fall to the permutation
floor. Confirmed in isolation, not reached on any alignment tested, so
this is a widening and a guard test at 120 000 sites. Values below the
overflow are unchanged.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

---

### Task 8: Memory: gamete matrices, the Hudson-Kaplan pair list, GENECONV permutations  (audit item: A13; "the incompatible-pair list"; test gaps "the Hudson-Kaplan and GENECONV memory paths", "a dense incompatibility matrix", "the PHI window a test names")

Three memory sites outside the audited caps, measured on `4b35e39` with `tracemalloc` (numpy reports its allocations to it):

| site | input | peak before | peak after |
|---|---|---|---|
| `diagnostics.incompatibility_matrix` (lines 131-143): four float64 z x z gamete matrices live together | z = 2000, 40 sequences | 136 MB | 25 MB |
| `diagnostics.hudson_kaplan_rmin` (line 213): every incompatible pair as a Python tuple | dense 2000 x 2000 | 381 MB | 0.25 MB |
| `geneconv.run_pvalue_perm` (lines 52-69): the full permutation matrix, unbatched | 20 000 sites, 2000 permutations | 680 MB | 250 MB |

At the PHI column cap (z = 5000) the first two scale to about 900 MB and 2 GB on a dense matrix. The fixes: build the four gamete matrices one at a time in float32 (every count is at most the number of sequences, far below 2**24, so float32 is exact) and reduce each to a boolean before the next; run the Hudson-Kaplan greedy column by column (the same choice as sorting all pairs by `(j, i)`, proved by a test against that algorithm); batch the GENECONV permutations as MaxChi does, drawing the same random stream so the p-value is identical. The GENECONV batch is 10 million cells: a 5-million-cell batch kept peak memory at 125 MB but took 4.2 s against 2.7 s unbatched on the 20 000-site input; 10 million takes 3.4 s.

The test gap "the PHI window a test names": `test_phi_clonal_not_significant_recombinant_significant` (`tests/unit/test_diagnostics.py:102-112`) asks for `window=8` on 24 informative sites, and `effective_phi_window` runs it with 2. The test is changed to ask for 2 and to assert the window it ran with. That is a test-only change; it passes before and after.

**Files:**
- Modify: `src/tessera/recomb/diagnostics.py` (`incompatibility_matrix` lines 131-143; `hudson_kaplan_rmin` lines 203-225)
- Modify: `src/tessera/recomb/geneconv.py` (constants after line 33; `run_pvalue_perm` lines 52-69)
- Test: `tests/unit/test_scaling_guards.py` (append), `tests/unit/test_diagnostics.py` (lines 102-108)

**Interfaces:**
- Produces: `incompatibility_matrix(allele1, allele0) -> np.ndarray[bool]`, `hudson_kaplan_rmin(incompatible, positions) -> tuple[int, list[tuple[int, int]]]` and `run_pvalue_perm(steps, run_len, *, seed=0) -> float`, all unchanged in signature and value. New constant `geneconv._BATCH_CELLS = 10_000_000`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_scaling_guards.py`:

```python
# --- PHI, Hudson-Kaplan and GENECONV: peak memory -------------------------

def _random_gametes(z: int, n: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    allele1 = (rng.random((z, n)) < 0.4).astype(np.float64)
    return allele1, 1.0 - allele1


def _peak_mb(fn, *args):
    import tracemalloc

    tracemalloc.start()
    try:
        out = fn(*args)
        return out, tracemalloc.get_traced_memory()[1] / 1e6
    finally:
        tracemalloc.stop()


def test_incompatibility_matrix_does_not_hold_four_float64_matrices() -> None:
    # z = 2000: one float64 z*z matrix is 32 MB. Four live together peaked at 136 MB;
    # one float32 matrix at a time stays near 25 MB.
    from tessera.recomb.diagnostics import incompatibility_matrix

    a1, a0 = _random_gametes(2000, 40)
    matrix, peak = _peak_mb(incompatibility_matrix, a1, a0)
    assert matrix.dtype == bool and matrix.shape == (2000, 2000)
    assert peak < 60


def _rmin_by_sorted_pairs(incompatible: np.ndarray, positions: np.ndarray):
    """The textbook greedy: every incompatible pair, sorted by right endpoint."""
    pairs = [(int(i), int(j)) for i, j in zip(*np.where(np.triu(incompatible, 1)), strict=True)]
    pairs.sort(key=lambda p: (p[1], p[0]))
    intervals, last = [], -1
    for i, j in pairs:
        if i >= last:
            intervals.append((int(positions[i]), int(positions[j])))
            last = j
    return len(intervals), intervals


def test_hudson_kaplan_matches_the_sorted_pair_greedy() -> None:
    from tessera.recomb.diagnostics import hudson_kaplan_rmin

    rng = np.random.default_rng(1)
    for density in (0.0, 0.02, 0.3, 1.0):
        upper = np.triu(rng.random((60, 60)) < density, 1)
        matrix = upper | upper.T
        positions = np.arange(60) * 7
        assert hudson_kaplan_rmin(matrix, positions) == _rmin_by_sorted_pairs(matrix, positions)


def test_hudson_kaplan_does_not_materialise_the_pairs() -> None:
    # A dense 2000 x 2000 matrix holds ~2 million incompatible pairs; as Python tuples
    # they cost ~380 MB. The column-wise greedy needs almost nothing beyond the matrix.
    from tessera.recomb.diagnostics import hudson_kaplan_rmin

    matrix = np.ones((2000, 2000), dtype=bool)
    np.fill_diagonal(matrix, False)
    (rmin, intervals), peak = _peak_mb(hudson_kaplan_rmin, matrix, np.arange(2000))
    assert rmin == len(intervals) == 1999  # every adjacent pair is its own interval
    assert peak < 20


def test_geneconv_permutations_are_batched() -> None:
    # 10 000 discriminating sites: 4000 permutations, 40 million cells. Unbatched, the
    # random draw, the argsort and the mask peaked at about 680 MB per candidate donor.
    from tessera.recomb import geneconv

    steps = np.random.default_rng(2).choice([1, -1], size=10_000).astype(np.int64)
    p, peak = _peak_mb(geneconv.run_pvalue_perm, steps, 12)
    assert 0.0 < p <= 1.0
    assert peak < 300


def test_geneconv_batching_does_not_change_the_p_value(monkeypatch) -> None:
    # The batches draw the same random stream in the same order as one full draw, so
    # the p-value is identical. Shrink the batch so a small input takes several.
    from tessera.recomb import geneconv

    steps = np.random.default_rng(3).choice([1, -1], size=200).astype(np.int64)
    perms = geneconv._PERMUTATIONS
    order = np.argsort(np.random.default_rng(0).random((perms, steps.size)), axis=1)
    is_minor = steps[order] == -1
    longest = np.zeros(perms, dtype=np.int32)
    cur = np.zeros(perms, dtype=np.int32)
    for j in range(steps.size):
        cur = np.where(is_minor[:, j], cur + 1, 0)
        longest = np.maximum(longest, cur)
    unbatched = float((np.count_nonzero(longest >= 6) + 1) / (perms + 1))

    monkeypatch.setattr(geneconv, "_BATCH_CELLS", 50_000, raising=False)
    assert geneconv.run_pvalue_perm(steps, 6) == unbatched
```

`test_hudson_kaplan_matches_the_sorted_pair_greedy` and `test_geneconv_batching_does_not_change_the_p_value` are guards (they pass before: the reference algorithm is the current one, and `raising=False` lets the monkeypatch set an attribute the old module ignores); they pin that the rewrite changes no value.

In `tests/unit/test_diagnostics.py`, replace

```python
def test_phi_clonal_not_significant_recombinant_significant() -> None:
    identity = lambda c: c  # noqa: E731
    clonal = recombination_signal(_block_alignment(recombinant=False), "s0", identity,
                                  window=8, seed=1)
    recomb = recombination_signal(_block_alignment(recombinant=True), "s0", identity,
                                  window=8, seed=1)
    assert clonal is not None and recomb is not None
```

with

```python
def test_phi_clonal_not_significant_recombinant_significant() -> None:
    identity = lambda c: c  # noqa: E731
    # 24 informative sites: effective_phi_window caps any request at 24 // 10 = 2, so
    # ask for the window the test actually runs with (it used to ask for 8 and get 2).
    clonal = recombination_signal(_block_alignment(recombinant=False), "s0", identity,
                                  window=2, seed=1)
    recomb = recombination_signal(_block_alignment(recombinant=True), "s0", identity,
                                  window=2, seed=1)
    assert clonal is not None and recomb is not None
    assert clonal.n_informative == recomb.n_informative == 24
    assert clonal.phi_window == recomb.phi_window == 2
```

- [ ] **Step 2: Run them, expected failure**

Run: `pytest tests/unit/test_scaling_guards.py tests/unit/test_diagnostics.py -q`
Expected (observed): 3 failed -- the incompatibility-matrix peak (`assert 136.0... < 60`), the Hudson-Kaplan peak (`assert 381.4... < 20`), and the GENECONV peak (unbatched, about 680 MB at 10 000 sites, `assert ... < 300`). The two guards and the PHI test pass.

- [ ] **Step 3: Implement**

(a) Replace `incompatibility_matrix` in `src/tessera/recomb/diagnostics.py` with:

```python
def incompatibility_matrix(allele1: np.ndarray, allele0: np.ndarray) -> np.ndarray:
    """The ``z x z`` four-gamete incompatibility matrix (True = all four gametes present).

    Gamete counts come from indicator dot products, so sequences missing at either
    column never contribute a gamete. The diagonal is False.

    The four gamete-count matrices are formed one at a time, in float32, and each is
    reduced to a boolean before the next is built. Holding all four float64 matrices at
    once peaked at about 900 MB at the column cap (z = 5000), four times what the cap was
    meant to allow. float32 represents every count exactly: a count is at most the number
    of sequences, far below 2**24.
    """
    a1 = allele1.astype(np.float32, copy=False)
    a0 = allele0.astype(np.float32, copy=False)
    incompatible = (a1 @ a1.T) > 0
    incompatible &= (a1 @ a0.T) > 0
    incompatible &= (a0 @ a1.T) > 0
    incompatible &= (a0 @ a0.T) > 0
    np.fill_diagonal(incompatible, False)
    return incompatible
```

(b) Replace `hudson_kaplan_rmin` with:

```python
def hudson_kaplan_rmin(
    incompatible: np.ndarray, positions: np.ndarray
) -> tuple[int, list[tuple[int, int]]]:
    """Hudson-Kaplan ``(Rmin, intervals)`` -- the maximum set of disjoint
    four-gamete-incompatible intervals (a lower bound on recombination events).

    Each incompatible column pair ``(i, j)`` requires a recombination in the rank
    interval ``[i, j]``; the greedy maximum set of non-overlapping intervals (by right
    endpoint) is Rmin. ``intervals`` are the chosen ones in MSA column coordinates.

    The greedy walks the right endpoint ``j`` in order and, at each, selects the pair
    with the smallest left endpoint ``i >= last_right`` if there is one -- the same
    choice as sorting every incompatible pair by ``(j, i)``, without materialising the
    pairs. A dense matrix at the column cap holds millions of them, and a list of Python
    tuples cost hundreds of megabytes on top of the matrix itself.
    """
    z = incompatible.shape[0]
    rmin = 0
    intervals: list[tuple[int, int]] = []
    last_right = 0  # left endpoints are >= 0, so 0 admits every first interval
    for j in range(1, z):
        hits = np.flatnonzero(incompatible[last_right:j, j])
        if hits.size:
            i = last_right + int(hits[0])
            rmin += 1
            intervals.append((int(positions[i]), int(positions[j])))
            last_right = j
    return rmin, intervals
```

Why it is the same greedy: with pairs sorted by `(j, i)`, at a given right endpoint `j` the first pair with `i >= last_right` is selected and `last_right` becomes `j`; every later pair with the same `j` has `i < j = last_right` and is skipped. So at each `j` the choice is "the smallest incompatible `i` in `[last_right, j)`, if any" -- one slice of column `j`.

(c) In `src/tessera/recomb/geneconv.py` add after `_MIN_PERMUTATIONS = 2_000`:

```python
# Cells per permutation batch; bounds peak memory without changing the permutation count.
_BATCH_CELLS = 10_000_000
```

and replace `run_pvalue_perm` with:

```python
def run_pvalue_perm(steps: np.ndarray, run_len: int, *, seed: int = 0) -> float:
    """Permutation ``P(longest minor-run >= run_len)`` over random arrangements of the
    same major/minor counts (exact-in-the-limit; mirrors the 3SEQ permutation fallback)."""
    if run_len <= 0:
        return 1.0
    n_steps = steps.size
    perms = _PERMUTATIONS
    if n_steps * perms > _MAX_PERM_CELLS:
        perms = max(_MIN_PERMUTATIONS, _MAX_PERM_CELLS // max(1, n_steps))
    rng = np.random.default_rng(seed)
    # Batched, as in MaxChi: the full (perms, n_steps) draw is up to 40 million cells,
    # about 730 MB across its random, order and mask arrays, per candidate donor. The
    # batches draw the same random stream in the same order, so the p-value is identical
    # to the unbatched one; only the peak memory changes.
    batch = max(1, min(perms, _BATCH_CELLS // max(1, n_steps)))
    exceed = 0
    remaining = perms
    while remaining > 0:
        k = min(batch, remaining)
        order = np.argsort(rng.random((k, n_steps)), axis=1)
        is_minor = steps[order] == -1
        longest = np.zeros(k, dtype=np.int32)
        cur = np.zeros(k, dtype=np.int32)
        for j in range(n_steps):  # vectorised over permutations; reset-cumulative run
            cur = np.where(is_minor[:, j], cur + 1, 0)
            longest = np.maximum(longest, cur)
        exceed += int(np.count_nonzero(longest >= run_len))
        remaining -= k
    return float((exceed + 1) / (perms + 1))
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_scaling_guards.py tests/unit/test_diagnostics.py tests/unit/test_geneconv.py -q`
Expected: all pass, including the existing `test_hudson_kaplan_rmin_textbook_cases`.

- [ ] **Step 5: Full checks and regression gate**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q` (expected `780 passed`), then `validation/data/audit-g/gate.sh task8`.

Acceptance: every file identical to `task7` (measured identical). PHI and Rmin only set `parent_free_support`, which the gate does not print; additionally check that `recombination_profile.tsv` of every `run_validation.py` dataset is byte-identical to `task7`'s:

```bash
for d in validation/data/*/_run/recomb_out; do shasum "$d/recombination_profile.tsv"; done
```

Any difference: revert and record.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/diagnostics.py src/tessera/recomb/geneconv.py tests/unit/test_scaling_guards.py tests/unit/test_diagnostics.py
git commit -m "$(cat <<'MSG'
Bound three memory sites outside the scaling caps

incompatibility_matrix held four float64 z x z gamete matrices at once
(136 MB at z = 2000, about 900 MB at the 5000-column cap, against the
~225 MB the cap's comment budgets); it now builds them one at a time in
float32, which represents every count exactly, and reduces each to a
boolean before the next (25 MB at z = 2000).

hudson_kaplan_rmin listed every incompatible pair as a Python tuple (381
MB on a dense 2000 x 2000 matrix). The same greedy now runs column by
column; a test checks it against the sorted-pair algorithm.

GENECONV's permutation null allocated its full matrix unbatched (680 MB at
20 000 sites); it is now batched at 10 million cells, drawing the same
random stream, so the p-value is identical (250 MB).

The PHI test that asked for window 8 and ran with 2 now asks for 2 and
asserts it. No harness output changes.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

---

### Task 9: A terminal tract's breakpoint interval; the exact/permutation switch  (audit item: A14, second half; test gap "the exact-versus-permutation switch")

A tract that starts at the first discriminating site gets a zero-width breakpoint interval: `bp_lo = column_to_query(cols[max(0, peak - 1)])` is the same column as `bp_hi = column_to_query(cols[peak])` when `peak == 0` (`threeseq.py:218`, `geneconv.py:121`, and the same line in `maxchi.py:170`, which the spec does not list). There is no discriminating site before such a tract, so the breakpoint may lie anywhere from the alignment start to the first site: the lower bound is query position 0.

The guard test for the switch between the exact descent p-value and its permutation fallback (`threeseq.descent_pvalue`, `_MAX_EXACT_CELLS`) is added here because it belongs to the same module and no task changes that switch.

**Files:**
- Modify: `src/tessera/recomb/threeseq.py` (line 218), `src/tessera/recomb/maxchi.py` (line 170), `src/tessera/recomb/geneconv.py` (line 121)
- Test: `tests/unit/test_threeseq.py` (imports; append)

**Interfaces:**
- Produces: for 3SEQ, MaxChi and GENECONV regions whose tract starts at the first discriminating site, `breakpoint_lo == 0`. Otherwise unchanged.

- [ ] **Step 1: Write the failing test and the guard**

In `tests/unit/test_threeseq.py` add `import pytest` after `import numpy as np`, then append:

```python
def _call(path: Path, method: str):
    result = compute_similarity(str(path), "query", window_size=500, window_step=50)
    return call_regions(result, analyze(result), 500,
                        RegionParams.with_defaults(500, method=method))


def test_descent_pvalue_switches_to_permutation_above_the_cell_cap(monkeypatch) -> None:
    # Guard: which branch answered is reported in the `test` column, and the two agree.
    from tessera.recomb import threeseq

    steps = np.array([1] * 60 + [-1] * 12 + [1] * 60, dtype=np.int64)
    d = max_descent(steps)
    exact, kind = descent_pvalue(steps, d)
    assert kind == "exact"
    monkeypatch.setattr(threeseq, "_MAX_EXACT_CELLS", d.m * d.n * d.depth - 1)
    perm, kind = descent_pvalue(steps, d)
    assert kind == "permutation"
    assert perm == pytest.approx(exact, abs=3 / threeseq._PERMUTATIONS + exact * 0.2)


def test_a_tract_at_the_first_site_has_a_breakpoint_interval_from_the_start(tmp_path) -> None:
    # The tract starts at the first discriminating site, so there is no site before it;
    # the interval ran from that site to itself (zero width) instead of from column 0.
    import random

    from ..factories import _mutate

    rng = random.Random(5)
    a = "".join(rng.choice("ACGT") for _ in range(6000))
    b, other = _mutate(a, 0.05, rng), _mutate(a, 0.15, rng)
    query = b[:1500] + a[1500:]
    path = write_fasta(tmp_path / "terminal.fasta",
                       {"query": query, "A": a, "B": b, "other": other})
    for method in ("3seq", "maxchi", "geneconv"):
        regions, _, _ = _call(path, method)
        terminal = [r for r in regions if r.minor_parent == "B"]
        assert len(terminal) == 1, method
        assert terminal[0].breakpoint_lo == 0, method
        assert terminal[0].breakpoint_hi > 0, method
```

`_call` is shared with Task 12's tests. `test_descent_pvalue_switches_to_permutation_above_the_cell_cap` is a guard (passes before and after).

- [ ] **Step 2: Run, expected failure**

Run: `pytest tests/unit/test_threeseq.py -q`
Expected (observed): 1 failed, `test_a_tract_at_the_first_site_has_a_breakpoint_interval_from_the_start`, with `AssertionError: 3seq` / `assert 25 == 0`: the region starts at column 25, the first discriminating site, and its interval was `[25, 25]`.

- [ ] **Step 3: Implement**

In each of the three files replace the line

```python
        bp_lo = result.column_to_query(int(cols[max(0, peak - 1)]))
```

(`start` instead of `peak` in `geneconv.py`) with

```python
        # The breakpoint lies between the last site before the tract and its first
        # site. A tract starting at the first discriminating site has no site before
        # it, so the interval runs from the alignment start, not zero-width.
        bp_lo = result.column_to_query(int(cols[peak - 1])) if peak > 0 else 0
```

(`start` for `peak` in `geneconv.py`).

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_threeseq.py tests/unit/test_maxchi.py tests/unit/test_geneconv.py -q`
Expected: all pass.

- [ ] **Step 5: Full checks and regression gate**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q` (expected `782 passed`), then `validation/data/audit-g/gate.sh task9`.

Acceptance: every file identical to `task8` (the gate does not print breakpoint intervals; coordinates, p-values and rows cannot change). Any difference: revert and record.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/threeseq.py src/tessera/recomb/maxchi.py src/tessera/recomb/geneconv.py tests/unit/test_threeseq.py
git commit -m "$(cat <<'MSG'
Give a tract at the first discriminating site a breakpoint interval

3SEQ, MaxChi and GENECONV bound a tract's left breakpoint by the site
before the tract and its first site. For a tract starting at the first
discriminating site there is no site before it, and max(0, peak - 1)
made both bounds the same column: a zero-width interval. It now runs from
the alignment start. Adds a guard test for the exact-to-permutation
switch of the 3SEQ p-value.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

---

### Task 10: A pooled lineage carries one genome's evidence  (audit item: A1; test gap "`clustered_view` with clusters of unequal size feeding the emission matrix")

**Evaluate-first.** `clusters.clustered_view` (`src/tessera/recomb/clusters.py:185-188`) sums each cluster's per-window numerators and denominators, and `hmm._emission_matrix` (`hmm.py:126-134`) pads every state to `n_max = den.max(axis=0)` with `log(0.5)` per missing site. A cluster of k genomes therefore carries k times the evidence of a singleton, and each singleton state is charged for the shortfall. Reproduced on `4b35e39`: three near-duplicates of the backbone stop the HMM calling a 1 kb tract from a single donor genome; `--no-cluster-lineages` finds it again.

**Two candidates were compared** (scratch driver, 20 seeds of a 12 kb panel A, D, E, F, ~6 % apart, query = A with a 1 kb tract from D; "dups" = three extra copies of A ~0.4 % apart; HMM caller only):

| emission | HMM calls the D tract, no duplicates | ... three backbone duplicates |
|---|---|---|
| current (sum within a cluster) | 14/20 | 1/20 |
| mean within a cluster | 14/20 | 14/20 |
| sum, padded to the largest single-genome count | 14/20 | 1/20 |

Padding to the single-genome maximum leaves the pooled state's own k-fold likelihood in place and turns its padding into a bonus (`(n_single - k * n) * log 0.5 > 0`), so it cannot help; the mean removes the weighting. The mean is rounded to whole sites so the counts stay integers for every consumer (`WindowSimilarity.numerators` is `dict[str, list[int]]`; widening it to float ripples through `siblings.py` and `stats.proportion_diff_significant`). Rounding moves a count by at most half a site per window. The similarity is computed from the unrounded means, so it is unchanged (the pooled ratio).

The sibling test (`siblings._profiles`) reads the same counts as binomial sample sizes; with the mean, a pooled lineage is no longer k times "more significant" than a singleton in the per-window proportion test. That is the intended direction and is covered by the gate.

**Files:**
- Modify: `src/tessera/recomb/clusters.py` (module docstring lines 5-9; `clustered_view`, lines 166-201)
- Modify: `docs/detection-methods.md` ("Near-duplicate references (lineage clustering)", lines 92-100)
- Test: `tests/unit/test_clusters.py` (imports; `test_clustered_view_pools_counts_and_is_identity_for_singletons`, lines 92-115, **changed**; one appended test), `tests/integration/test_regimes.py` (append)

**Interfaces:**
- Consumes: `rg.simulate_redundant` from Task 2.
- Produces: `clustered_view(result, clusters) -> (view, cluster_size)` unchanged in signature; the view's per-window counts for a cluster are the rounded member means.

**Existing test that changes, and why.** `test_clustered_view_pools_counts_and_is_identity_for_singletons` asserts that pooled counts are the sum of the members' counts -- the behaviour this task removes. It is changed to assert the rounded mean and the pooled ratio as the similarity; its singleton half is unchanged.

- [ ] **Step 1: Write the failing tests**

In `tests/unit/test_clusters.py` add `import pytest` after `from pathlib import Path` (with a blank line), replace

```python
    # Pooling two members sums their per-window counts.
    view, sizes = clustered_view(result, [["sib0", "sib1"], ["pa0"]])
    rep = next(r for r, n in sizes.items() if n == 2)
    for i in range(len(result.positions)):
        assert view.numerators[rep][i] == (
            result.numerators["sib0"][i] + result.numerators["sib1"][i]
        )
        assert view.denominators[rep][i] == (
            result.denominators["sib0"][i] + result.denominators["sib1"][i]
        )
```

with

```python
    # Pooling two members averages their per-window counts (rounded to whole sites), so
    # the pooled lineage carries one genome's worth of evidence; the similarity is the
    # pooled ratio.
    view, sizes = clustered_view(result, [["sib0", "sib1"], ["pa0"]])
    rep = next(r for r, n in sizes.items() if n == 2)
    for i in range(len(result.positions)):
        num = result.numerators["sib0"][i] + result.numerators["sib1"][i]
        den = result.denominators["sib0"][i] + result.denominators["sib1"][i]
        assert view.numerators[rep][i] == round(num / 2)
        assert view.denominators[rep][i] == round(den / 2)
        assert view.similarities[rep][i] == pytest.approx(num / den)
```

and append:

```python
def test_a_pooled_cluster_carries_one_genomes_evidence(tmp_path) -> None:
    # Four exact copies of one genome pooled into a lineage must emit as that genome
    # does. Summed, the copies carried four times its evidence and every singleton state
    # was charged log(0.5) for each site of the difference.
    from tessera.recomb.hmm import _emission_matrix

    lane_a, lane_b = _lanes()
    query = lane_a[:R1] + lane_b[R1:]
    rng = random.Random(12)
    a = _mutate(lane_a, 0.02, rng)
    records = {"query": query, "A": a, "A1": a, "A2": a, "A3": a,
               "D": _mutate(lane_b, 0.02, rng)}
    result = _result(records, tmp_path)
    view, sizes = clustered_view(result, [["A", "A1", "A2", "A3"], ["D"]])
    assert sizes == {"A": 4, "D": 1}
    assert view.numerators["A"] == result.numerators["A"]
    assert view.denominators["A"] == result.denominators["A"]
    pooled = _emission_matrix(view, ["A", "D"], 0.95)
    single = _emission_matrix(result, ["A", "D"], 0.95)
    assert pooled == pytest.approx(single)
```

Append to `tests/integration/test_regimes.py`:

```python
# --- redundant lineages (near-duplicates on a divergent panel) ---------------

def test_redundant_backbone_does_not_silence_the_hmm(tmp_path) -> None:
    # Three near-duplicates of the backbone pooled into one lineage carried three times
    # the evidence of the single donor genome D0, and the HMM lost a 1 kb tract it finds
    # when the backbone is one genome.
    seqs, (lo, hi) = rg.simulate_redundant(6001, backbone_copies=3)
    result = compute_similarity(str(write_fasta(tmp_path / "rb.fasta", seqs)), rg.QUERY,
                                WINDOW, STEP)
    regions, major = _hmm_regions(result)
    assert major.startswith("A")
    donors = [r for r in regions if r.minor_parent == "D0"]
    assert len(donors) == 1
    assert donors[0].msa_start < hi and donors[0].msa_end > lo
```

- [ ] **Step 2: Run them, expected failure**

Run: `pytest tests/unit/test_clusters.py tests/integration/test_regimes.py -q`
Expected (observed): 3 failed -- `assert 599 == 300` (the sum where the mean is now expected), `assert [1172, 1180, ...] == [293, 295, ...]` (four copies carry four times one genome's counts), and `test_redundant_backbone_does_not_silence_the_hmm` with `assert 0 == 1` (no D0 region).

- [ ] **Step 3: Implement**

Replace `clustered_view` in `src/tessera/recomb/clusters.py` with:

```python
def clustered_view(
    result: WindowSimilarity, clusters: list[list[str]]
) -> tuple[WindowSimilarity, dict[str, int]]:
    """A representative-labelled view of ``result`` with per-cluster pooled evidence.

    Each cluster becomes one label (its representative) whose per-window
    ``numerators``/``denominators`` are the element-wise **mean** over the cluster's
    members, rounded to whole sites, so a pooled lineage carries one genome's worth of
    evidence. The similarity is the ratio of the unrounded means, which equals the
    pooled ratio of the sums. ``rows`` is left untouched (all member
    genomes remain), so the downstream discordant-site test stays at genome level.
    Returns ``(view, cluster_size)`` mapping each representative to its member count.

    The mean, not the sum: the HMM emission is a log-likelihood over these counts and
    pads every state to the largest count in the window. Summed counts gave a cluster of
    k genomes k times the evidence of a single genome, and charged every singleton state
    a padding penalty for the difference; three near-duplicates of the backbone were
    enough to stop the HMM calling a 1 kb tract from a singleton donor.
    """
    n_win = len(result.positions)
    pooled_num: dict[str, list[int]] = {}
    pooled_den: dict[str, list[int]] = {}
    pooled_sim: dict[str, list[float]] = {}
    cluster_size: dict[str, int] = {}
    for members in clusters:
        rep = _representative(result, members)
        cluster_size[rep] = len(members)
        num = np.mean([result.numerators[m] for m in members], axis=0)
        den = np.mean([result.denominators[m] for m in members], axis=0)
        # Rounding moves a count by at most half a site per window, well inside the
        # binomial noise of the window, and keeps the counts integers for every consumer.
        pooled_num[rep] = [int(round(float(x))) for x in num]
        pooled_den[rep] = [int(round(float(x))) for x in den]
        pooled_sim[rep] = [
            float(num[i] / den[i]) if den[i] > 0 else float("nan") for i in range(n_win)
        ]
    best_sim, best_label = _best_per_window(pooled_sim, n_win)
    view = replace(
        result,
        similarities=pooled_sim,
        numerators=pooled_num,
        denominators=pooled_den,
        best_sim=best_sim,
        best_label=best_label,
    )
    return view, cluster_size
```

In the module docstring replace "with each cluster's per-window evidence summed." by "with each cluster's per-window evidence averaged over its members (so a pooled lineage weighs as one genome)." and leave the next sentence ("with one genome per cluster the pooled view equals the input") as it is -- it is now true for any cluster of identical members too.

In `docs/detection-methods.md`, in "Near-duplicate references (lineage clustering)", after "The pooled lineage competes as one state under the label of its best-covering member, and a region names that member." add: "Its evidence in each window is the mean over its members, not the sum: a lineage of five genomes weighs as one genome in the emission, so adding copies of a reference cannot outweigh a single-genome donor."

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_clusters.py tests/unit/test_siblings.py tests/integration/test_regimes.py -q`
Expected: all pass, including `test_cluster_caller_excludes_sibling_lineage_and_calls_one_region` (the sibling lineage is still excluded; the pooled minor is still reported with cluster size 3).

- [ ] **Step 5: Full checks and harness gate**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q` (expected `784 passed`), then:

```bash
validation/data/audit-g/gate.sh task10
python validation/run_hybrids.py > validation/data/audit-g/task10/hybrids.txt 2>&1
```

Acceptance:
- `redundant_backbone`, `--min-methods 1`: HMM among the callers **at least 5/10** (baseline 3/10; measured 6/10 with this task applied alone to `4b35e39`, against 8/10 for the one-genome control). Detected and donor lineage 10/10 at both gates.
- `redundant_donor`: HMM among the callers **at least 8/10** (baseline 0/10; measured 10/10 with this task applied alone to `4b35e39`; Tasks 3-9 leave every regime line unchanged, so the cumulative value should match).
- `spec_mm1.txt`, `spec_mm2.txt`, positive control, `near_identical`, `near_identical` clonal and `sibling_tract`: identical to `task9` (measured identical: the specificity tips are about 2.4 % apart and are never pooled, and the near-identical panels hold no duplicates).
- `validation_regions.tsv`: report every changed row. A changed row is expected only on a dataset whose panel holds near-duplicates; a row on any other dataset that changes is a regression.
- Hybrids: no case from PASS to FAIL; no `neg_*` case gains a region.

If the gate fails: revert, record all outputs, and record "Task 11 not attempted: depends on Task 10" -- Task 11's pipeline test expects the HMM among the callers of the redundant donor, which this task provides. Go to Task 12.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/clusters.py docs/detection-methods.md tests/unit/test_clusters.py tests/integration/test_regimes.py
git commit -m "$(cat <<'MSG'
Average, not sum, a lineage cluster's evidence in the HMM

clustered_view summed each cluster's per-window counts and the emission
pads every state to the window's largest count, so a cluster of k
near-duplicates carried k times a singleton's evidence and each singleton
was charged log(0.5) per site of the shortfall. Three copies of the
backbone stopped the HMM calling a 1 kb tract from a single donor genome
that --no-cluster-lineages found.

The pooled counts are now the member means, rounded to whole sites; the
similarity is the pooled ratio as before. Padding to the largest
single-genome count was also measured and does not help (1/20 against
14/20 with the mean, on 20 scratch panels).

Harness: redundant_backbone HMM among callers <n>/10 (baseline 3/10);
redundant_donor <n>/10 (baseline 0/10); specificity and the other regimes
<identical or the diff>; hybrids <summary or "not run: no network">.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

Replace the `<...>` fields with the measured values before committing.

---

### Task 11: One set of lineage clusters for every caller and the merge  (audit item: A2; test gap "two near-duplicates of the donor in an untyped panel yielding one region")

**Depends on Task 10** (the pipeline test below expects the HMM among the callers of a redundant donor).

**Evaluate-first.** Lineage clustering runs only inside `_call_regions_hmm` (`regions.py:280-283`). The HMM then names the cluster representative while 3SEQ and MaxChi name whichever genome their own ranking picks, and `ensemble._same_donor` (`ensemble.py:32-41`) groups by genome label unless a lineage map exists: one event becomes two rows. Separately `bootscan.py:143` requires a single genome to reach 0.70 bootstrap support, which two near-duplicates split. Reproduced on `4b35e39` with a 3 kb tract from a donor lineage of one, two and three genomes ~0.4 % apart (8 seeds): one row naming a D genome 8/8, 6/8, 8/8; Bootscan among the callers 8/8, 4/8, 4/8.

The clusters are now computed once per run on the HMM's scan (`clusters.lineage_clusters`, the same gate the HMM applied) and passed to every caller through `RegionParams.reference_clusters` and to the merge through a `clusters` keyword. The HMM pools them as before; Bootscan sums bootstrap support over a cluster's candidates and names the cluster's best-supported member; the merge treats two members of one cluster as the same donor. With no multi-member cluster every path is unchanged, which the gate checks.

**Interaction with plan C's item C1 (commit `a53ea53`).** Under informative-site windowing, clustering compares a pair on informative columns only, in windows of 40 informative sites, at the 98.5 % floor -- so a pair that differs at even one informative site in a window is not pooled. On a near-identical panel only exact (or nearly exact) duplicates are therefore clustered, and this task groups only those. A donor copy 0.03 % apart on a 0.4 %-divergent panel is not pooled, and the event is still reported as two rows (one naming each copy). That is recorded in Self-review note 3 and not changed here: the remedy is a choice about C1's window and floor, not about A2.

**Files:**
- Modify: `src/tessera/recomb/clusters.py` (new `lineage_clusters` after `all_singletons`, line 154)
- Modify: `src/tessera/recomb/regions.py` (import, line 30; `RegionParams` and `with_defaults`, lines 64-130; `_call_regions_hmm` lines 277-283; new `_with_singletons` before `_span_windows`)
- Modify: `src/tessera/recomb/ensemble.py` (`_same_donor` lines 32-41; `_group` lines 134-156; `consensus_regions`)
- Modify: `src/tessera/recomb/bootscan.py` (new `_candidate_groups`; `call_regions_bootscan` lines 132-165)
- Modify: `src/tessera/recomb/run.py` (imports; `_region_params` lines 329-341; the caller loop at line 350; the `consensus_regions(...)` call)
- Modify: `docs/detection-methods.md` (ensemble paragraph lines 20-26; lineage-clustering paragraph lines 92-100; Bootscan section lines 178-184)
- Test: `tests/unit/test_ensemble.py`, `tests/unit/test_bootscan.py`, `tests/integration/test_regimes.py` (append)

**Interfaces:**
- Consumes: `cluster_references`, `clustered_view` (Task 10's averaging), `rg.simulate_redundant`, `rg.score_lineage_tract`, `NEAR` (already in `tests/integration/test_regimes.py`).
- Produces:
  - `clusters.lineage_clusters(result: WindowSimilarity, window_size: int, params) -> list[list[str]]` -- multi-member clusters only, `[]` when none or when the gate does not apply.
  - `RegionParams.reference_clusters: list[list[str]] | None = None` (and the same keyword on `with_defaults`); `None` = "not computed", so a direct `call_regions` caller gets the old behaviour.
  - `regions._with_singletons(result, pooled) -> list[list[str]]`.
  - `bootscan._candidate_groups(candidates: list[str], clusters: list[list[str]] | None) -> np.ndarray`.
  - `ensemble._same_donor(a, b, lineage_map, cluster_of=None)`, `ensemble._group(regions, lineage_map, cluster_of=None)`, `consensus_regions(..., clusters: list[list[str]] | None = None)`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_ensemble.py`:

```python
def test_regions_naming_two_members_of_one_cluster_are_one_event() -> None:
    # The HMM reports a pooled lineage under its representative (D1); 3SEQ ranks single
    # genomes and names D0. Without the cluster membership they were two rows.
    hmm = mk("D1", 100, 200, "hmm", qvalue=1e-6)
    seq = mk("D0", 110, 190, "3seq", qvalue=1e-9)
    split, _ = consensus_regions({"hmm": [hmm], "3seq": [seq]}, major="A")
    assert len(split) == 2
    merged, breakdown = consensus_regions(
        {"hmm": [hmm], "3seq": [seq]}, major="A", clusters=[["D0", "D1"]]
    )
    assert len(merged) == 1
    assert merged[0].methods == ("hmm", "3seq")
    assert merged[0].minor_parent == "D0"  # the most significant member names the donor
    assert len(breakdown) == 1
```

Append to `tests/unit/test_bootscan.py`:

```python
def _redundant_donor_msa(tmp_path: Path) -> Path:
    """recombinant_msa with the donor lineage B represented by two genomes ~0.4 % apart,
    each 0.2 % from the B sequence the query carries (neither is closer by design)."""
    from tessera.core.io import read_fasta

    from ..factories import _mutate

    rng = random.Random(3)
    records = dict(read_fasta(recombinant_msa(tmp_path, recombinant=True)))
    original_b = records["B"]
    records["B"] = _mutate(original_b, 0.002, rng)
    records["B2"] = _mutate(original_b, 0.002, rng)
    return write_fasta(tmp_path / "redundant.fasta", records)


def test_bootscan_counts_support_per_lineage_cluster(tmp_path: Path) -> None:
    # Two near-identical genomes of the donor lineage split the bootstrap support between
    # them, so neither reaches 70 % and the tract is lost. Counted per cluster, it is not.
    result = compute_similarity(str(_redundant_donor_msa(tmp_path)), "query",
                                window_size=500, window_step=50)
    split, _, _ = call_regions(
        result, analyze(result), 500, RegionParams.with_defaults(500, method="bootscan")
    )
    assert split == []
    pooled, major, _ = call_regions(
        result, analyze(result), 500,
        RegionParams.with_defaults(500, method="bootscan", reference_clusters=[["B", "B2"]]),
    )
    assert major == "A"
    assert len(pooled) == 1
    assert pooled[0].minor_parent in {"B", "B2"}
    assert pooled[0].query_start < 2300 and pooled[0].query_end > 3700
```

Append to `tests/integration/test_regimes.py`:

```python
def test_redundant_donor_is_one_region_with_every_caller(tmp_path) -> None:
    seqs, tract = rg.simulate_redundant(6001, donor_copies=3, tract=rg.RED_DONOR_TRACT)
    got = rg.score_lineage_tract(rg.scan(seqs, _LOG), tract, lineage="D")
    assert got == {"detected": True, "rows": 1, "donor_ok": True, "hmm": True,
                   "bootscan": True}


def test_no_backbone_duplicate_is_reported_as_a_donor(tmp_path) -> None:
    # Review focus: with the backbone lineage redundant, a near-duplicate of the
    # backbone must not win Bootscan windows, or any caller's, as a "donor".
    seqs, _ = rg.simulate_redundant(6002, backbone_copies=3)
    rows = [r for r in rg.scan(seqs, _LOG) if r.get("donor_absent") != "yes"]
    assert rows
    assert not [r for r in rows if r["minor_parent"].startswith("A")]


def test_duplicate_donor_on_a_near_identical_panel_is_one_row(tmp_path) -> None:
    # Review focus: the shared clusters are computed on the informative-site scan the
    # HMM uses, while Bootscan and the merge work on base-pair windows. An exact copy of
    # the donor is pooled there, and the event must still be one row. (A copy differing
    # at a few informative sites is not pooled under informative-site windowing; see the
    # plan's self-review notes.)
    seqs, (lo, hi) = rg.simulate_near_identical(3, recombinant=True, **NEAR)
    seqs["R1_copy"] = seqs["R1"]
    rows = [r for r in rg.scan(seqs, _LOG) if r.get("donor_absent") != "yes"]
    hits = [r for r in rows if int(r["query_start"]) < hi and int(r["query_end"]) > lo]
    assert len(hits) == 1
    assert hits[0]["minor_parent"] in {"R1", "R1_copy"}
    assert "hmm" in hits[0]["methods"].split(",")
```

The last two are review-focus guards: they pass before and after this task.

- [ ] **Step 2: Run them, expected failure**

Run: `pytest tests/unit/test_ensemble.py tests/unit/test_bootscan.py tests/integration/test_regimes.py -q`
Expected (observed, with Tasks 1-10 in place): 3 failed --
- `test_regions_naming_two_members_of_one_cluster_are_one_event`: `TypeError: consensus_regions() got an unexpected keyword argument 'clusters'`;
- `test_bootscan_counts_support_per_lineage_cluster`: `TypeError: RegionParams.with_defaults() got an unexpected keyword argument 'reference_clusters'`;
- `test_redundant_donor_is_one_region_with_every_caller`: `{'bootscan': False} != {'bootscan': True}`, `{'rows': 3} != {'rows': 1}`.

- [ ] **Step 3: Implement**

(a) `src/tessera/recomb/clusters.py`, after `all_singletons`:

```python
def lineage_clusters(result: WindowSimilarity, window_size: int, params) -> list[list[str]]:
    """The panel's lineage clusters with more than one member; ``[]`` when none apply.

    Applies the same gate as the HMM caller (``cluster_lineages`` on, and between
    ``cluster_min_refs`` and ``cluster_max_refs`` references). Computed once per run and
    shared, so that every caller and the ensemble merge see the same lineages: the HMM
    pools them, Bootscan sums their bootstrap support, and the merge treats two members
    of one cluster as the same donor.
    """
    n_refs = len(result.similarities)
    in_range = params.cluster_min_refs <= n_refs <= params.cluster_max_refs
    if not (params.cluster_lineages and in_range):
        return []
    return [c for c in cluster_references(result, window_size, params) if len(c) > 1]
```

(b) `src/tessera/recomb/regions.py`: change the import to `from .clusters import _window_bounds, clustered_view, lineage_clusters`. In `RegionParams`, after `lineage_map: dict[str, str] | None = None` add

```python
    # Lineage clusters computed once for the run (clusters.lineage_clusters): multi-member
    # clusters only, ``[]`` for none. ``None`` means "not computed": the HMM caller then
    # computes them itself, as it did before they were shared.
    reference_clusters: list[list[str]] | None = None
```

add `reference_clusters: list[list[str]] | None = None,` as the last keyword of `with_defaults` and `reference_clusters=reference_clusters,` as the last argument of its `cls(...)`. In `_call_regions_hmm` replace

```python
    n_refs = len(result.similarities)
    if params.cluster_lineages and params.cluster_min_refs <= n_refs <= params.cluster_max_refs:
        clusters = cluster_references(result, window_size, params)
        if not all_singletons(clusters):
            work, cluster_size = clustered_view(result, clusters)
```

with

```python
    pooled = params.reference_clusters
    if pooled is None:
        pooled = lineage_clusters(result, window_size, params)
    if pooled:
        work, cluster_size = clustered_view(result, _with_singletons(result, pooled))
```

and insert before `def _span_windows(`:

```python
def _with_singletons(result: WindowSimilarity, pooled: list[list[str]]) -> list[list[str]]:
    """Every reference in one cluster: the multi-member ``pooled`` clusters plus each
    remaining label alone, in the order the labels first appear in ``result`` (the
    order ``cluster_references`` returns, so the HMM's state order is unchanged)."""
    index = {label: i for i, members in enumerate(pooled) for label in members}
    out: list[list[str]] = []
    seen: set[int] = set()
    for label in result.similarities:
        i = index.get(label)
        if i is None:
            out.append([label])
        elif i not in seen:
            seen.add(i)
            out.append(pooled[i])
    return out
```

(c) `src/tessera/recomb/ensemble.py`: replace `_same_donor` with

```python
def _same_donor(
    a: Region, b: Region, lineage_map: LineageMap | None,
    cluster_of: dict[str, int] | None = None,
) -> bool:
    """Do two regions name the same donor? Equal genome labels; or members of one
    lineage cluster (near-duplicates the HMM pools; see ``clusters.lineage_clusters``);
    or -- when references are typed -- the same lineage. Two callers that pick different
    genomes of one lineage then still count as the same event (and as agreement)."""
    if a.minor_parent == b.minor_parent:
        return True
    if cluster_of:
        ca = cluster_of.get(a.minor_parent)
        if ca is not None and ca == cluster_of.get(b.minor_parent):
            return True
    if lineage_map:
        la = lineage_of(a.minor_parent, lineage_map)
        return la is not None and la == lineage_of(b.minor_parent, lineage_map)
    return False
```

replace the signature and docstring of `_group` with

```python
def _group(
    regions: list[Region], lineage_map: LineageMap | None,
    cluster_of: dict[str, int] | None = None,
) -> list[list[Region]]:
    """Group regions that overlap (transitively) and name the same donor (genome,
    lineage cluster or, when typed, lineage)."""
```

and in its body change `_same_donor(group[0], region, lineage_map)` to `_same_donor(group[0], region, lineage_map, cluster_of)`. In `consensus_regions` add the keyword `clusters: list[list[str]] | None = None,` after `scan`, append to the docstring

```
    ``clusters`` are the run's lineage clusters (``clusters.lineage_clusters``). Regions
    naming two members of one cluster are the same event: the HMM reports a cluster under
    its representative, while 3SEQ and MaxChi rank individual genomes and may name
    another member.
```

and replace `    for group in _group(all_regions, lineage_map):` with

```python
    cluster_of = {label: i for i, members in enumerate(clusters or []) for label in members}
    for group in _group(all_regions, lineage_map, cluster_of):
```

(d) `src/tessera/recomb/bootscan.py`: insert before `def call_regions_bootscan(`:

```python
def _candidate_groups(candidates: list[str], clusters: list[list[str]] | None) -> np.ndarray:
    """Group index per candidate: candidates in one lineage cluster share an index."""
    of = {label: i for i, members in enumerate(clusters or []) for label in members}
    keys: dict[object, int] = {}
    out = np.empty(len(candidates), dtype=np.int64)
    for k, label in enumerate(candidates):
        key: object = ("cluster", of[label]) if label in of else ("genome", label)
        out[k] = keys.setdefault(key, len(keys))
    return out
```

and in `call_regions_bootscan` replace everything from `    bounds = _window_bounds(result, window_size)  # window column spans (per position)` down to (not including) `    if not runs:` with:

```python
    bounds = _window_bounds(result, window_size)  # window column spans (per position)
    width = match_stack.shape[1]
    # Near-duplicates of one lineage split the bootstrap support between them: two
    # genomes 0.4 % apart each win about half the resamples, and neither reaches the
    # 70 % threshold the lineage would. Support is therefore counted per lineage cluster
    # (the run's shared clusters; a genome outside every cluster is its own), and the
    # cluster's best-supported member names the donor.
    group = _candidate_groups(candidates, getattr(params, "reference_clusters", None))
    n_groups = int(group.max()) + 1
    major_group = int(group[major_idx])
    # winner[w] = the candidate index with top bootstrap support at window w (or major)
    winner = np.full(len(bounds), major_idx)
    support = np.zeros(len(bounds))
    for w, (start, end) in enumerate(bounds):
        s, e = max(0, start), min(width, end)
        if e - s < 1:
            continue
        sup = bootstrap_support(match_stack[:, s:e], comp_stack[:, s:e], seed=w)
        by_group = np.bincount(group, weights=sup, minlength=n_groups)
        best_group = int(np.argmax(by_group))
        if (best_group != major_group and by_group[best_group] >= _SUPPORT
                and by_group[best_group] > by_group[major_group]):
            members = np.flatnonzero(group == best_group)
            winner[w] = int(members[np.argmax(sup[members])])
            support[w] = by_group[best_group]

    # Contiguous runs of the same non-major winner -> candidate regions, each scored
    # against the block-permutation null so a run of chance flips cannot pass as an
    # event. Bootstrap support says which reference is closest here; it does not say
    # that the query recombined, which is what the run length is tested for.
    block = _overlap_block(bounds)
    # Runs are of one lineage, not one genome: the member with the top support can change
    # from window to window within a cluster.
    won = np.where(winner == major_idx, -1, group[winner])
    runs: list[tuple[int, int, int, float]] = []  # (donor, start_w, end_w, p)
    w = 0
    while w < len(won):
        if won[w] < 0:
            w += 1
            continue
        lineage = won[w]
        j = w
        while j < len(won) and won[j] == lineage:
            j += 1
        if j - w >= _MIN_WINDOWS:
            # The run's donor is the member that won most of its windows.
            donor = int(np.bincount(winner[w:j]).argmax())
            p = run_pvalue(won == lineage, j - w, block=block, seed=donor)
            runs.append((donor, w, j, p))
        w = j
```

With no clusters every candidate is its own group, numbered in candidate order, so `by_group == sup`, the winner and the run are exactly as before.

(e) `src/tessera/recomb/run.py`: add `from .clusters import lineage_clusters` (after `from .analyze import ...`). Change `def _region_params(method: str) -> RegionParams:` to `def _region_params(method: str, clusters: list[list[str]] | None = None) -> RegionParams:` and add `reference_clusters=clusters,` after `lineage_map=lineage_map,` in its body. After the closing of `_region_params` insert

```python
    # Lineage clusters (near-duplicate references), computed once on the HMM's scan and
    # shared: the HMM pools them, Bootscan counts their support together, and the merge
    # treats two members of one cluster as the same donor.
    reference_clusters = lineage_clusters(result, params.window_size, _region_params("hmm"))
```

change the caller line to

```python
        regs, major_m, sibs = call_regions(
            res, ana, params.window_size, _region_params(method, reference_clusters)
        )
```

and add `clusters=reference_clusters,` to the `consensus_regions(...)` call after `scan=...`.

Clustering now runs even when the HMM is not among `--method`; with `--no-cluster-lineages` it does not run at all (`lineage_clusters` returns `[]`).

(f) `docs/detection-methods.md`:
- in the ensemble paragraph, change "the same minor parent genome, or, when the references are typed (a lineage map), the same **lineage**," to "the same minor parent genome, two members of one lineage cluster (near-duplicates pooled for the HMM, below), or, when the references are typed (a lineage map), the same **lineage**,";
- in "Near-duplicate references", change "so the HMM caller first pools them (`--cluster-lineages`, on by default; hmm only)" to "so they are first grouped into lineage clusters (`--cluster-lineages`, on by default). The clusters are computed once and shared: the HMM pools them, Bootscan counts their bootstrap support together, and the merge treats two members of one cluster as one donor";
- in the Bootscan section, after "carrying that support as a confidence the other callers express only as a p-value." add: "Support is counted per lineage cluster: two near-duplicates of the donor each win about half the resamples, and neither alone would reach 70 %."

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_ensemble.py tests/unit/test_bootscan.py tests/unit/test_clusters.py tests/unit/test_regions.py tests/integration -q`
Expected: all pass.

- [ ] **Step 5: Full checks and harness gate**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q` (expected `789 passed`), then:

```bash
validation/data/audit-g/gate.sh task11
python validation/run_hybrids.py > validation/data/audit-g/task11/hybrids.txt 2>&1
```

Acceptance:
- `redundant_donor`, both gates: one row **at least 9/10**, Bootscan among the callers **at least 8/10**, HMM among the callers not below `task10` (baseline 7/10, 1/10 and 0/10; with Task 10 alone 7/10, 1/10 and 10/10; measured with this task 10/10, 10/10 and 10/10 at `--min-methods 1`). At `--min-methods 2`, **detected at least 9/10**: on the baseline the split rows can each fall below the gate and the event disappears (detected 7/10 at 10 replicates, 0/3 for the harness's first three seeds), with this task 10/10.
- `redundant_backbone`: not worse than `task10` on any column, and no row names an `A` genome as the donor.
- every other line identical to `task10` (measured identical).
- `validation_regions.tsv`: report every changed row; a change is expected only where a panel holds near-duplicates of a donor.
- Hybrids: no case from PASS to FAIL; no `neg_*` case gains a region.

If the gate fails: revert and record.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/clusters.py src/tessera/recomb/regions.py src/tessera/recomb/ensemble.py src/tessera/recomb/bootscan.py src/tessera/recomb/run.py docs/detection-methods.md tests/unit/test_ensemble.py tests/unit/test_bootscan.py tests/integration/test_regimes.py
git commit -m "$(cat <<'MSG'
Share one set of lineage clusters between the callers and the merge

Lineage clustering ran inside the HMM caller only. The HMM then named a
cluster's representative while 3SEQ and MaxChi named whichever genome
ranked first, and the merge grouped by genome label: with the donor
lineage represented by near-duplicates one event became two rows, and
under --min-methods 2 the HMM's row could be suppressed. Bootscan required
a single genome to reach 70 % bootstrap support, which two near-duplicates
split, so it fell silent.

The clusters are now computed once per run and passed to every caller and
to the merge. The HMM pools them as before; Bootscan counts support per
cluster and names its best-supported member; the merge treats two members
of one cluster as one donor. Without multi-member clusters nothing changes.

Harness: redundant_donor one row <n>/10, Bootscan <n>/10 (baseline 7/10,
1/10); other lines <identical or the diff>; hybrids <summary>.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

---

### Task 12: The descent floor and the Benjamini-Hochberg family  (audit items: A3, A6; test gap "what the pre-Benjamini-Hochberg floors exclude, pinned by an exact p-value")

**Evaluate-first.** `threeseq.py:157` sets `_MIN_DESCENT = 4`, and `threeseq.py:189-190` and `maxchi.py:142-143` drop a donor whose descent is shallower before any p-value is computed, although the comment at `threeseq.py:154-155` says the p-value gate does the real work. Two consequences:
- **A3:** a short tract can be significant and is never tested. Three consecutive donor matches among 189 discriminating sites have an exact 3SEQ p of 187 / C(189, 3) = 1.7e-4 and are discarded.
- **A6:** the BH family is built from the survivors only (`threeseq.py:191-198`, `maxchi.py:143-150`). Depth is monotone in p at fixed counts, so the least significant donors are removed before the correction and the rest are inflated (q = p with one survivor where the honest family has two members and q = 2p).

**What was tried first, and its cost.** Removing the floor in both callers left every harness line unchanged but made a 16-genome specificity scan **ten times slower** (4.3 s to 45 s). On the 12 clonal specificity panels at 3 replicates there were 144 donors with a descent shallower than 4; every one had an exact 3SEQ p of 0.21 or more, and their MaxChi nulls (20 000 permutations of walks of 1000-2000 sites each) took 34 s of the 45. The 3SEQ exact DP took the other 7 s: it recursed cell by cell over `m x n` in Python, 0.2-0.4 s per long shallow walk.

**The chosen candidate:**
1. **3SEQ:** no depth floor; every donor with at least 20 discriminating sites is tested and is one member of the BH family (a donor with depth 0 enters with p = 1). The exact DP is rewritten to fill its table one anti-diagonal at a time -- each cell depends only on the previous diagonal, so each diagonal is one vectorised step. The values are identical to the cell-by-cell recursion (checked for every `m, n <= 8`, `depth <= 6`, and on the large cases below to the last printed digit) and it is 10-40 times faster:

   | m, n, depth | p | before | after |
   |---|---|---|---|
   | 186, 3, 3 | 1.69e-4 | <0.01 s | <0.01 s |
   | 1217, 74, 3 | 0.210 | 0.16 s | 0.010 s |
   | 1795, 116, 3 | 0.343 | 0.40 s | 0.016 s |
   | 3000, 300, 40 | 1.26e-38 | 1.70 s | 0.088 s |

2. **MaxChi:** a donor with a descent shallower than 4 is **counted in the BH family with p = 1 and not permuted** (`_MIN_PERMUTED_DEPTH = 4`). This fixes A6 for MaxChi: p = 1 is a valid upper bound, and replacing p-values above alpha by 1 cannot change which donors pass BH, because every such term exceeds alpha (`p * m / rank > alpha` when `p > alpha`). It does not fix A3 for MaxChi, deliberately: MaxChi scores the same tract as 3SEQ (Task 15, A12), 3SEQ now tests the shallow tracts exactly, and the permutation cost is the ten-fold slowdown above. With both, the specificity scan takes 4.5 s (4.3 s before).

GENECONV has the same kind of floor (`geneconv._MIN_RUN = 4`) and the same family pruning; it is opt-in and not named by the spec, so it is left alone (Self-review note 5).

**Files:**
- Modify: `src/tessera/recomb/threeseq.py` (`descent_pvalue_exact` lines 89-129; the floor constants and comment lines 154-157; `call_regions_3seq` lines 175-190)
- Modify: `src/tessera/recomb/maxchi.py` (import of `_MIN_DESCENT` line 26; constants after line 38; `call_regions_maxchi` lines 130-146)
- Modify: `docs/detection-methods.md` (3SEQ section line 163; MaxChi section)
- Modify: `validation/README.md` (line 565-566, the sentence on the pre-BH floors)
- Test: `tests/unit/test_threeseq.py` (append)

**Interfaces:**
- Consumes: `_call` from Task 9 (in `tests/unit/test_threeseq.py`).
- Produces: `threeseq._MIN_DESCENT` is removed (only `maxchi` imported it); `maxchi._MIN_PERMUTED_DEPTH = 4`; `descent_pvalue_exact(m, n, depth) -> float` same signature and values.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_threeseq.py` (before `_call`, or anywhere after the imports):

```python
# --- pre-test floors, the BH family, the exact/permutation switch ------------

def _three_site_alignment(tmp_path: Path, *, decoy: bool = False, sites: int = 3) -> Path:
    """Major P and minor Q differ at 189 columns; the query is P except at ``sites``
    consecutive ones, where it carries Q. With ``decoy`` a third reference R differs from
    P at 30 other columns, where the query matches P (a donor with no descent at all)."""
    rng = np.random.default_rng(0)
    width = 4000
    p = rng.integers(0, 4, size=width)
    q = p.copy()
    q_cols = np.sort(rng.choice(np.arange(0, width, 2), size=189, replace=False))
    q[q_cols] = (p[q_cols] + 1) % 4
    query = p.copy()
    query[q_cols[100:100 + sites]] = q[q_cols[100:100 + sites]]
    seqs = {"query": query, "P": p, "Q": q}
    if decoy:
        r = p.copy()
        r_cols = rng.choice(np.arange(1, width, 2), size=30, replace=False)  # not Q's
        r[r_cols] = (p[r_cols] + 2) % 4
        seqs["R"] = r
    bases = np.frombuffer(b"ACGT", dtype=np.uint8)
    text = {k: bases[v].tobytes().decode() for k, v in seqs.items()}
    return write_fasta(tmp_path / "three.fasta", text)


def test_a_three_site_tract_is_tested_not_discarded_before_the_p_value(tmp_path) -> None:
    # Three consecutive donor matches among 189 sites: exact p = 187 / C(189, 3) = 1.7e-4.
    # The depth floor of 4 discarded the triplet before any p-value was computed.
    regions, major, _ = _call(_three_site_alignment(tmp_path), "3seq")
    assert major == "P"
    assert [r.minor_parent for r in regions] == ["Q"]
    assert math.isclose(regions[0].pvalue, 187 / math.comb(189, 3), rel_tol=1e-9)


def test_the_bh_family_counts_every_tested_donor(tmp_path) -> None:
    # R has 30 discriminating sites and no descent (p = 1). It was pruned by the depth
    # floor before Benjamini-Hochberg, so the family had one member and q equalled p.
    path = _three_site_alignment(tmp_path, decoy=True, sites=8)
    for method in ("3seq", "maxchi"):
        regions, _, _ = _call(path, method)
        donors = [r for r in regions if r.minor_parent == "Q"]
        assert len(donors) == 1, method
        assert math.isclose(donors[0].qvalue, 2 * donors[0].pvalue, rel_tol=1e-9), method
```

- [ ] **Step 2: Run them, expected failure**

Run: `pytest tests/unit/test_threeseq.py -q`
Expected (observed): 2 failed -- `test_a_three_site_tract_is_tested_not_discarded_before_the_p_value` with `assert [] == ['Q']`; `test_the_bh_family_counts_every_tested_donor` with `assert math.isclose(<q>, 2 * <p>)` false for `3seq` (q equals p: the decoy R was pruned before BH).

- [ ] **Step 3: Implement**

(a) Replace `descent_pvalue_exact` in `src/tessera/recomb/threeseq.py` with:

```python
def descent_pvalue_exact(m: int, n: int, depth: int) -> float:
    """Exact ``P(max drawdown >= depth)`` over the uniform arrangements of ``m`` +1 and
    ``n`` -1 steps.

    Dynamic program over the walk's current drawdown ``delta`` in ``[0, depth)``: with
    ``i`` up-steps and ``j`` down-steps still to come, an up-step moves
    ``delta -> max(0, delta - 1)`` with probability ``i / (i + j)`` and a down-step
    ``delta -> delta + 1`` with probability ``j / (i + j)``; a down-step from
    ``delta == depth - 1`` reaches the depth. ``g(i, j, delta)`` is the probability of
    reaching it, and the answer is ``g(m, n, 0)``.

    The reaching probability is propagated directly rather than as one minus the
    probability of never reaching it: every term is a sum of non-negative products, so
    a very small p-value keeps its magnitude instead of cancelling to exactly ``0.0``.
    Probabilities (not path counts) are propagated, so there is no big-integer overflow.

    The table is filled one anti-diagonal ``i + j = t`` at a time: every cell on it
    depends only on the previous diagonal, so each diagonal is one vectorised step. The
    values are those of the cell-by-cell recursion; only the order of evaluation
    changed, which made a long walk with a shallow descent (``m`` in the thousands) ten
    to forty times faster.
    """
    if depth <= 0:
        return 1.0
    if n < depth:  # need at least ``depth`` consecutive down-steps to reach it
        return 0.0
    prev = np.zeros((m + 1, depth))  # diagonal t-1, indexed by i (ups remaining)
    cur = np.zeros((m + 1, depth))
    for t in range(1, m + n + 1):
        lo, hi = max(0, t - n), min(m, t)
        i = np.arange(lo, hi + 1, dtype=np.float64)
        j = t - i
        block = np.zeros((hi - lo + 1, depth))
        # an up-step: (i, j, d) -> (i - 1, j, max(0, d - 1)), from row i - 1 of diagonal t-1
        a = 1 if lo == 0 else 0  # row i = 0 has no up-step
        if hi >= lo + a:
            src = prev[lo + a - 1 : hi]
            up = np.empty_like(src)
            up[:, 0] = src[:, 0]
            up[:, 1:] = src[:, :-1]
            block[a:] += (i[a:] / t)[:, None] * up
        # a down-step: (i, j, d) -> (i, j - 1, d + 1), or the depth is reached
        b = 1 if hi == t else 0  # row j = 0 has no down-step
        if hi - b >= lo:
            src = prev[lo : hi - b + 1]
            down = np.empty_like(src)
            down[:, :-1] = src[:, 1:]
            down[:, -1] = 1.0
            block[: hi - lo + 1 - b] += (j[: hi - lo + 1 - b] / t)[:, None] * down
        cur[lo : hi + 1] = block
        prev, cur = cur, prev
    return float(prev[m][0])
```

(b) Replace

```python
# The triplet test needs enough discriminating sites to be meaningful, and a drawdown
# deeper than trivial noise; the p-value gate does the real work above these floors.
_MIN_DISCRIMINATING = 20
_MIN_DESCENT = 4
```

with

```python
# The triplet test needs enough discriminating sites to be meaningful. There is no floor
# on the descent depth: a shallow descent can still be significant when the donor
# matches are few (three consecutive ones among 189 sites give p = 1.7e-4), so the
# p-value decides, and every donor tested counts towards the Benjamini-Hochberg family.
_MIN_DISCRIMINATING = 20
```

(c) In `call_regions_3seq` delete

```python
        if descent.depth < _MIN_DESCENT:
            continue
```

and replace the comment lines "win no window yet still carry a strong discriminating-site signal. The pre-BH floors / below drop candidates without a real drawdown, so the extra donors cost little power." with "win no window yet still carry a strong discriminating-site signal. Each donor tested / is one member of the Benjamini-Hochberg family below, whatever its p-value."

(d) In `src/tessera/recomb/maxchi.py` remove `    _MIN_DESCENT,` from the `from .threeseq import (...)` block, add after `_MAX_PERM_CELLS = 20_000_000`:

```python
# Below this descent depth a donor enters the BH family with p = 1 instead of being
# permuted; see call_regions_maxchi.
_MIN_PERMUTED_DEPTH = 4
```

and in `call_regions_maxchi` replace

```python
        descent = max_descent(steps)
        if descent.depth < _MIN_DESCENT:
            continue
        chi2 = maxchi_statistic(steps, descent.start_site, descent.end_site)
        p = maxchi_pvalue(steps, chi2)
```

with

```python
        descent = max_descent(steps)
        chi2 = maxchi_statistic(steps, descent.start_site, descent.end_site)
        # A donor with a shallow descent is still a member of the Benjamini-Hochberg
        # family: it is counted, with p = 1, rather than dropped before the correction,
        # which shrank the family to the survivors and inflated their significance. Its
        # permutation null is not run: 20 000 permutations of a long walk per such donor
        # made a 16-genome scan ten times slower, and 3SEQ, which scores the same tract
        # exactly and cheaply, tests these shallow tracts (see threeseq).
        p = maxchi_pvalue(steps, chi2) if descent.depth >= _MIN_PERMUTED_DEPTH else 1.0
```

and change the comment "window yet carry a strong discriminating-site signal (the pre-BH floors drop the rest)." to "window yet carry a strong discriminating-site signal.".

(e) `docs/detection-methods.md`, 3SEQ section, after "Benjamini-Hochberg-corrected across the donors tested." add: "Every donor with at least 20 discriminating sites is tested and counted in that family, however shallow its descent: three consecutive donor matches among 189 sites are significant (p = 1.7e-4)." In the MaxChi section add at the end: "A donor whose descent is shallower than four sites is counted in MaxChi's Benjamini-Hochberg family with p = 1 rather than permuted; 3SEQ tests those short tracts exactly."

(f) `validation/README.md`, replace "The pre-Benjamini-Hochberg floors (`_MIN_DISCRIMINATING` / `_MIN_DESCENT` / min run) drop candidates without a real signal, so the extra donors cost little power." with "Candidates with fewer than 20 discriminating sites (`_MIN_DISCRIMINATING`) are not tested; since audit 4, 3SEQ tests every other donor and MaxChi counts shallow ones in its family with p = 1, so the extra donors enlarge the Benjamini-Hochberg family rather than escape it."

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_threeseq.py tests/unit/test_maxchi.py tests/unit/test_geneconv.py -q`
Expected: all pass, including `test_descent_pvalue_exact_matches_brute_force` and `test_descent_pvalue_exact_keeps_small_values` (the rewritten DP against enumeration and the closed form).

- [ ] **Step 5: Full checks and harness gate**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q` (expected `791 passed`), then:

```bash
validation/data/audit-g/gate.sh task12
python validation/run_hybrids.py > validation/data/audit-g/task12/hybrids.txt 2>&1
```

and time one harness against `task11`:

```bash
/usr/bin/time -p python validation/run_specificity.py --reps 3 --min-methods 2 > /dev/null
```

Acceptance:
- `spec_mm1.txt`, `spec_mm2.txt`: `TOTAL` lines no worse than `task11` (measured identical: 7/12 with 8 hmm regions, 0/12). **A new 3SEQ- or MaxChi-only false region in either file fails the gate.**
- Positive control, every regime line: identical to `task11` (measured identical).
- Runtime of the timed command within 1.5x of the same command on `task11` (measured 6:53 against 6:46 for the three harnesses together, 3 replicates).
- `validation_regions.tsv`: report every changed row with its p and q. Rows can appear (a shallow tract now tested) or have larger q (larger family); a row that disappears because its q rose above alpha is A6 working as intended, and is recorded, not reverted.
- Hybrids: no case from PASS to FAIL; no `neg_*` case gains a region.

If the gate fails: revert and record.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/threeseq.py src/tessera/recomb/maxchi.py docs/detection-methods.md validation/README.md tests/unit/test_threeseq.py
git commit -m "$(cat <<'MSG'
Test every 3SEQ donor and keep shallow donors in the BH family

3SEQ and MaxChi dropped a donor whose descent was shallower than four sites
before computing any p-value. A short tract could be significant and never
tested (three consecutive donor matches among 189 sites: exact p 1.7e-4),
and the Benjamini-Hochberg family held only the survivors, so the least
significant donors were removed before the correction and the rest were
inflated.

3SEQ now tests every donor with at least 20 discriminating sites. Its exact
DP fills the table one anti-diagonal at a time (identical values, 10-40x
faster on long shallow walks), which keeps this affordable: removing the
floor with the old DP and MaxChi's permutations made a 16-genome scan ten
times slower. MaxChi counts a shallow donor in its family with p = 1
instead of permuting it; that cannot change which donors pass BH, and 3SEQ
tests those tracts exactly.

Harness: <paste from validation/data/audit-g/task12/, and the timing>.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

---

### Task 13: `--min-region` bounds every reported region  (audit items: B2, A4; test gaps "`--min-region` against a site caller", "any relation between `--min-region` and a reported HMM region length")

`--min-region` is read only by the HMM (`regions.py:300`, on the untrimmed segment), the heuristic caller (`regions.py:426`), clustering (`clusters.py:122`) and the sibling test (`siblings.py:131,154`). 3SEQ, MaxChi, Bootscan, GENECONV and barcode never read it: with `--min-region 5000` a 175 bp 3SEQ region and a 49 bp MaxChi region are reported. And since plan C's trimming (`regions.py:327`), the HMM's own regions are not bounded either: an 11 200 bp segment passes the filter and is reported as 132 bp. The help text ("Minimum recombinant region length in bp (default: the window size)", `cli/cmd_recomb.py:70-73`) is false for both.

**Decision: apply it, when given, to every merged region; leave the default alone.** Applying the default (the window size, 1000 bp) to every caller would remove the short tracts 3SEQ and MaxChi exist to find -- the 1 kb `redundant_backbone` tract and every sub-window hybrid case recorded in `validation/README.md` -- and change the shipped behaviour for every user. So:
- `--min-region N` given: every region in `recombination_regions.tsv` is at least N query bases long, whichever caller found it. The filter acts on the merged regions, before the agreement gate, so an event called by several callers is judged on its merged span (Review Focus 2).
- `--min-region` not given: unchanged. The HMM keeps its segment floor at the window size (plan C's decision to test the segment, not the trimmed span, stands); the other callers are not bounded. The help text and the documentation now say so.

This changes no harness output (none passes `--min-region`); the gate checks that.

**Files:**
- Modify: `src/tessera/recomb/ensemble.py` (new `filter_by_length` before `consensus_regions`)
- Modify: `src/tessera/recomb/run.py` (import from `.ensemble`; after the `consensus_regions(...)` call, before the agreement gate at line 389)
- Modify: `src/tessera/cli/cmd_recomb.py` (help text, lines 70-73)
- Modify: `docs/detection-methods.md` (HMM section, after the "reported span" paragraph)
- Test: `tests/integration/test_regimes.py` (append)

**Interfaces:**
- Consumes: `rg.simulate_redundant`, `rg._read_regions` (Task 2 / existing).
- Produces: `ensemble.filter_by_length(regions: list[Region], breakdown: list[dict], min_length: int) -> tuple[list[Region], list[dict], int]`.

- [ ] **Step 1: Write the failing test and the guard**

Append to `tests/integration/test_regimes.py`:

```python
def test_min_region_bounds_every_callers_reported_length(tmp_path) -> None:
    from tessera.recomb.run import RecombParams, run_recomb

    seqs, (lo, hi) = rg.simulate_redundant(6001)  # a 1 kb tract, called by 3SEQ/MaxChi
    msa = write_fasta(tmp_path / "mr.fasta", seqs)
    lengths = {}
    for label, min_region in (("default", None), ("bounded", 5000)):
        out = tmp_path / label
        run_recomb(RecombParams(msa=msa, output=out, query=rg.QUERY, plot_format="png",
                                min_region=min_region), _LOG)
        rows = rg._read_regions(out / "recombination_regions.tsv")
        lengths[label] = [int(r["query_end"]) - int(r["query_start"]) for r in rows
                          if r.get("donor_absent") != "yes"]
    assert lengths["default"] and min(lengths["default"]) < 5000  # short tracts reported
    assert all(n >= 5000 for n in lengths["bounded"])


def test_min_region_is_applied_to_the_merged_span_before_the_agreement_gate(tmp_path) -> None:
    # Review focus: --min-region with --min-methods 2. The filter is on merged regions,
    # so a 3 kb event called by several callers survives a 2 kb floor even when one
    # caller's own tract is shorter than the floor.
    from tessera.recomb.run import RecombParams, run_recomb

    seqs, (lo, hi) = rg.simulate_redundant(6001, tract=rg.RED_DONOR_TRACT)
    msa = write_fasta(tmp_path / "mm.fasta", seqs)
    out = tmp_path / "out"
    run_recomb(RecombParams(msa=msa, output=out, query=rg.QUERY, plot_format="png",
                            min_region=2000, min_methods=2), _LOG)
    rows = [r for r in rg._read_regions(out / "recombination_regions.tsv")
            if r.get("donor_absent") != "yes"]
    assert len(rows) == 1
    assert int(rows[0]["query_start"]) < hi and int(rows[0]["query_end"]) > lo
    assert len(rows[0]["methods"].split(",")) >= 2
```

The second test is a review-focus guard (passes before and after).

- [ ] **Step 2: Run, expected failure**

Run: `pytest tests/integration/test_regimes.py -q -k min_region`
Expected (observed): 1 failed, `test_min_region_bounds_every_callers_reported_length` with `assert False` on `all(n >= 5000 for n in lengths["bounded"])` (3SEQ and MaxChi report the 1 kb tract under `--min-region 5000`).

- [ ] **Step 3: Implement**

(a) In `src/tessera/recomb/ensemble.py`, insert before `def consensus_regions(`:

```python
def filter_by_length(
    regions: list[Region], breakdown: list[dict], min_length: int
) -> tuple[list[Region], list[dict], int]:
    """Keep only regions at least ``min_length`` query bases long (``--min-region``).

    Applied to the merged regions, so it bounds what is reported whichever caller found
    a region; ``regions`` and ``breakdown`` are parallel and filtered together. Returns
    the kept regions, their breakdown rows and how many were removed.
    """
    kept = [
        (region, row) for region, row in zip(regions, breakdown, strict=True)
        if region.length_bp >= min_length
    ]
    return [r for r, _ in kept], [b for _, b in kept], len(regions) - len(kept)
```

(b) In `src/tessera/recomb/run.py` change the import to

```python
from .ensemble import (
    consensus_regions,
    filter_by_agreement,
    filter_by_length,
    reconcile_major,
)
```

and insert immediately before the comment `    # Agreement gate: a region resting on a single caller ...`:

```python
    # --min-region, when given, bounds the length of every reported region. Left unset,
    # only the HMM applies a floor (to its segments, before trimming); the site callers
    # exist to find tracts shorter than a window and are not bounded by default.
    if params.min_region is not None:
        regions, method_breakdown, too_short = filter_by_length(
            regions, method_breakdown, params.min_region
        )
        if too_short:
            logger.info(
                "Removed %d region(s) shorter than --min-region %d bp.",
                too_short, params.min_region,
            )
```

(c) In `src/tessera/cli/cmd_recomb.py` replace the help string with

```python
        help="Minimum length, in query bp, of every reported region, whichever caller "
        "found it. When unset, only the HMM applies a floor: it drops segments shorter "
        "than the window size before testing them, and the other callers report "
        "regions of any length.",
```

(d) In `docs/detection-methods.md`, after the paragraph that begins "**The reported span is not the segment.**", add:

```markdown
**`--min-region`.** When given, it bounds every reported region: a region shorter than
it, in query bases, is removed after the callers' regions are merged and before the
agreement gate, whichever caller found it. When it is not given, only the HMM applies a
length floor, and only to its segments before they are trimmed (the window size): a
trimmed HMM region can be shorter than a window, and 3SEQ, MaxChi, Bootscan and GENECONV
report regions of any length -- finding tracts shorter than a window is what the
site-based callers are for.
```

**If plan F Task 5 has landed** (it rewrites the same `--min-region` help, the `RecombParams.min_region` comment and adds a "Which options apply to which caller" paragraph saying 3SEQ, MaxChi, Bootscan, GENECONV and barcode do not apply `--min-region`): use the help text above, change the `RecombParams` comment to "Region-length floor. When given, every merged region shorter than it is removed; when None, the HMM drops segments shorter than the window size and the site callers are not bounded.", replace F's `--min-region` sentences in that paragraph with the `--min-region` paragraph above, and in F's `tests/unit/test_cli_option_scope.py` replace the `--min-region` assertions with `assert "whichever caller" in _help("recomb", "--min-region")`. If plan F lands later, its Task 5 executor makes the same substitutions.

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/integration/test_regimes.py tests/unit/test_ensemble.py tests/unit/test_cli_commands.py -q`
Expected: all pass.

- [ ] **Step 5: Full checks and regression gate**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q` (expected `793 passed`), then `validation/data/audit-g/gate.sh task13`.

Acceptance: every file identical to `task12` (no harness passes `--min-region`; measured identical). Any difference: revert and record.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/ensemble.py src/tessera/recomb/run.py src/tessera/cli/cmd_recomb.py docs/detection-methods.md tests/integration/test_regimes.py
git commit -m "$(cat <<'MSG'
Make --min-region bound every reported region when it is given

--min-region was read by the HMM, on its untrimmed segment, and by no
other caller: with --min-region 5000 a 175 bp 3SEQ region and a 49 bp
MaxChi region were reported, and since HMM spans are trimmed an 11.2 kb
segment could be reported as 132 bp. The help text said otherwise.

When given, it now removes every merged region shorter than it, before the
agreement gate. The default is unchanged: the HMM keeps its segment floor
at the window size and the site callers are not bounded, because finding
sub-window tracts is their purpose; applying the window size to them by
default would remove those tracts for every user. The help text and the
documentation now state both cases. No harness output changes.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

---

### Task 14: The last columns fall in a window  (audit "further item": trailing unscanned columns; test gap "the trailing unscanned columns")

`similarity._window_starts` (`src/tessera/recomb/similarity.py:192-193`) is `range(0, width - window + 1, step)`, so the last `(width - window) % step` columns are in no window: 499 of 5999 columns at window 1000, step 500. The HMM mitigates it through `_tract_span` (a segment touching the last window is searched to the alignment end), but Bootscan, the barcode caller and the coverage diagnostic see nothing there, and a tract or a coverage gap in those columns is invisible to them. Informative-site windowing (`similarity.py:354` and `411`) has the same gap, counted in informative sites. The fix adds one window flush with the end when the grid does not reach it. 3SEQ, MaxChi and GENECONV work on all discriminating sites and are not affected.

**Files:**
- Modify: `src/tessera/recomb/similarity.py` (`_window_starts`, lines 192-193; `compute_similarity`, line 265; `_informative_windowed_similarity`, line 354; `compute_similarity_informative`, line 411)
- Test: `tests/unit/test_similarity.py` (append)

**Interfaces:**
- Produces: `_window_starts(width: int, window: int, step: int) -> list[int]` (was `range`); the last start is always `width - window`. Every consumer of `positions`, `window_spans`, the per-window counts and `informative` sees one more window when the grid does not fit exactly; the last gap between positions is then shorter than the step. `median_window_step` and `bootscan._overlap_block` use medians and are unaffected.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_similarity.py`:

```python
def test_the_last_columns_fall_in_a_window(tmp_path: Path) -> None:
    # width 5999, window 1000, step 500: starts 0..4500 end at column 5500, leaving 499
    # columns in no window. A final window flush with the end covers them.
    import random

    rng = random.Random(1)
    seq = "".join(rng.choice("ACGT") for _ in range(5999))
    msa = write_fasta(tmp_path / "w.fasta", {"q": seq, "r1": seq, "r2": seq})
    result = compute_similarity(str(msa), "q", window_size=1000, window_step=500)
    assert result.positions[-1] + 500 == 5999  # the last window ends at the last column
    assert result.positions[:-1] == [500 + 500 * i for i in range(10)]
    assert len(result.similarities["r1"]) == len(result.positions) == 11
    assert len(result.informative) == 11


def test_an_exact_fit_adds_no_extra_window(tmp_path: Path) -> None:
    # Guard: when (width - window) is a multiple of step nothing changes.
    import random

    rng = random.Random(2)
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    msa = write_fasta(tmp_path / "x.fasta", {"q": seq, "r1": seq, "r2": seq})
    result = compute_similarity(str(msa), "q", window_size=1000, window_step=500)
    assert result.positions == [500 + 500 * i for i in range(11)]


def test_informative_windows_reach_the_last_informative_site(tmp_path: Path) -> None:
    # Review focus: the same trailing gap under informative-site windowing. 103
    # informative sites, 40 per window, step 5: starts 0..60 end at site 100.
    import random

    import numpy as np

    from tessera.recomb.similarity import (
        _informative_column_mask,
        compute_similarity_informative,
    )

    rng = random.Random(3)
    base = "".join(rng.choice("ACGT") for _ in range(3000))
    other = list(base)
    cols = sorted(rng.sample(range(3000), 103))
    for c in cols:
        other[c] = "A" if base[c] != "A" else "C"
    msa = write_fasta(tmp_path / "i.fasta", {"q": base, "r1": base, "r2": "".join(other)})
    result = compute_similarity_informative(str(msa), "q", info_window=40, info_step=5)
    info = np.flatnonzero(_informative_column_mask(result.rows, ["r1", "r2"])).tolist()
    assert len(info) == 103
    assert result.window_spans[-1] == (info[63], info[102] + 1)
```

`test_an_exact_fit_adds_no_extra_window` is a guard (passes before and after).

- [ ] **Step 2: Run them, expected failure**

Run: `pytest tests/unit/test_similarity.py -q`
Expected (observed): 2 failed -- `assert (5000 + 500) == 5999`, and for the informative case `assert (1936, 2905) == (2031, 2930)` (the last window spans informative sites 60-99, not 63-102).

- [ ] **Step 3: Implement**

(a) Replace `_window_starts` with:

```python
def _window_starts(width: int, window: int, step: int) -> list[int]:
    """Window start columns: every ``step``, plus one window flush with the end.

    Without the final window the last ``(width - window) % step`` columns fell in no
    window at all (499 of 5999 columns at window 1000, step 500), so a tract there was
    invisible to every window-based caller and to the coverage diagnostic.
    """
    starts = list(range(0, width - window + 1, step))
    if starts and starts[-1] != width - window:
        starts.append(width - window)
    return starts
```

(b) In `compute_similarity` change `    starts = list(_window_starts(width, window_size, window_step))` to `    starts = _window_starts(width, window_size, window_step)`.

(c) In `_informative_windowed_similarity` change `    for start in range(0, n_info - info_window + 1, info_step):` to `    for start in _window_starts(n_info, info_window, info_step):`.

(d) In `compute_similarity_informative` change `    starts = list(range(0, n_info - info_window + 1, info_step))` to `    starts = _window_starts(n_info, info_window, info_step)`.

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_similarity.py tests/unit/test_coverage.py tests/unit/test_hmm.py tests/integration -q`
Expected: all pass. No existing test pins the window count of an alignment whose width does not fit the grid.

- [ ] **Step 5: Full checks and harness gate**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q` (expected `796 passed`), then:

```bash
validation/data/audit-g/gate.sh task14
python validation/run_hybrids.py > validation/data/audit-g/task14/hybrids.txt 2>&1
```

Acceptance:
- every `TOTAL`, positive-control and regime line identical to `task13` (measured identical; every simulated panel -- 12 000 columns in `run_specificity.py`, 120 000, 30 000 and 12 000 in `run_regimes.py` -- fits the 1000/100 grid exactly, so this task cannot change them; the informative-site scans of `near_identical` can gain a window).
- `validation_regions.tsv`: report every changed row. A change is expected only where a region, a coverage gap or a donor-absent row touches the last `(width - window) % step` columns of that dataset; anything else is a regression.
- Hybrids: no case from PASS to FAIL; no `neg_*` case gains a region (a new window at a genome end is new evidence for every window-based caller and for the coverage diagnostic).

If the gate fails: revert and record.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/similarity.py tests/unit/test_similarity.py
git commit -m "$(cat <<'MSG'
Add a final window flush with the end of the alignment

Windows started every step from column 0 while a full window still fit,
so the last (width - window) % step columns fell in no window: 499 of
5999 at window 1000, step 500. The HMM searches a segment touching the
last window to the alignment end, but Bootscan, the barcode caller and the
coverage diagnostic saw nothing there. Informative-site windowing had the
same gap in informative sites.

One window flush with the end is added when the grid does not reach it.
An exact fit is unchanged.

Harness: <paste from validation/data/audit-g/task14/>.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

---

### Task 15: 3SEQ and MaxChi are not independent votes; the barcode vote  (audit items: A12, A14 first half; gap G2, third step)

A documentation correction and an addendum to the open decision record. No code changes behaviour. 

**A12.** Both triplet callers score the same tract -- MaxChi calls `threeseq.max_descent` on the same walk and computes its chi-square over that tract -- yet `docs/detection-methods.md` (lines 12-14 and 172-174) and the `maxchi.py` module docstring (lines 3-10) call them "genuinely independent votes" and say a region significant under both is "corroborated by two tests". Re-measured while writing this plan on 300 simulated clonal walks (40-200 discriminating sites, 5 to half of them donor matches, random order): Spearman correlation of the two p-values 0.95 (the audit measured 0.97); at alpha 0.05, 3SEQ rejects 0.033, MaxChi 0.037, both 0.030 -- independence would give 0.0012 -- and MaxChi rejects in 90 % of the walks where 3SEQ does. Each is individually calibrated. This bears on the open `--min-methods` decision: the two false regions that remained at the agreement gate (`--min-methods 2`, 10 replicates) were exactly a 3SEQ + MaxChi pair.

**A14, first half: the barcode caller votes without a significance test.** `barcode.py:155` sets no p-value, q-value or test, yet a barcode region counts toward `--min-methods` like any other. Whether it should is a maintainer decision (it changes what the agreement gate means and is not exercised by any harness: the barcode caller is opt-in and needs a typed panel). This task records the options; it does not change the count.

**Files:**
- Modify: `docs/detection-methods.md` (lines 12-14; lines 166-174; the `--min-methods` paragraph, lines 28-41; the barcode section, line 213-226)
- Modify: `src/tessera/recomb/maxchi.py` (module docstring, lines 1-17 -- text only)
- Modify: `docs/superpowers/specs/2026-10-01-min-methods-default-decision.md` (append, before "## Decision")
- Test: none. This task changes text; Step 2 re-measures the numbers it records.

**Interfaces:** none.

- [ ] **Step 1: Write the failing test**

Not applicable: no behaviour changes. The documentation claim being corrected is checked by Step 2's measurement instead.

- [ ] **Step 2: Re-measure the dependence**

```bash
python - <<'PY'
"""Dependence of the 3SEQ and MaxChi p-values under the clonal null (300 walks)."""
import numpy as np
from tessera.recomb.maxchi import maxchi_pvalue, maxchi_statistic
from tessera.recomb.threeseq import descent_pvalue, max_descent

rng = np.random.default_rng(0)
p3, pm = [], []
for k in range(300):
    n = int(rng.integers(40, 200))
    minor = int(rng.integers(5, n // 2))
    steps = np.array([1] * (n - minor) + [-1] * minor, dtype=np.int64)
    rng.shuffle(steps)
    d = max_descent(steps)
    p3.append(descent_pvalue(steps, d)[0])
    pm.append(maxchi_pvalue(steps, maxchi_statistic(steps, d.start_site, d.end_site), seed=k))
p3, pm = np.array(p3), np.array(pm)
rank = lambda x: np.argsort(np.argsort(x))  # noqa: E731
print(f"Spearman {np.corrcoef(rank(p3), rank(pm))[0, 1]:.2f}")
for alpha in (0.05, 0.01):
    a, b = p3 <= alpha, pm <= alpha
    print(f"alpha {alpha}: 3SEQ {a.mean():.3f} MaxChi {b.mean():.3f} both {(a & b).mean():.3f} "
          f"independent {a.mean() * b.mean():.4f} P(MaxChi | 3SEQ) {(a & b).sum() / max(1, a.sum()):.2f}")
PY
```

Expected (measured on the plan's candidate code; Tasks 7 and 12 do not change these values for walks this short): `Spearman 0.95`; `alpha 0.05: 3SEQ 0.033 MaxChi 0.037 both 0.030 independent 0.0012 P(MaxChi | 3SEQ) 0.90`; `alpha 0.01: 3SEQ 0.003 MaxChi 0.007 both 0.003 ...`. Use the printed numbers in Step 3.

- [ ] **Step 3: Correct the documentation**

(a) `docs/detection-methods.md`, lines 12-14: replace "MaxChi and Bootscan add two more independent votes." with "MaxChi and Bootscan add two more votes, though not independent ones (see the MaxChi section)."

(b) Replace the MaxChi paragraph's last sentence, from "Drawdown and proportion are sensitive to different mosaics" to "exactly as RDP treats them as separate methods.", with:

```markdown
MaxChi scores the tract 3SEQ finds -- the same maximum descent of the same walk -- so the
two are **not independent votes**. On 300 simulated clonal walks their p-values have a
Spearman correlation of 0.95, and at alpha 0.05 both reject in 3.0 % of walks where
independence would give 0.12 %: MaxChi rejects in 90 % of the walks where 3SEQ does. Each
is calibrated on its own (3.3 % and 3.7 % at a nominal 5 %). A region found by 3SEQ and
MaxChi together is therefore one triplet signal tested twice, not corroborated twice, and
`--min-methods 2` satisfied by that pair alone adds little to either caller. The callers
that look at different evidence are the HMM (per-window copying likelihood over the whole
panel, confirmed by a sign test on the segment), Bootscan (per-window bootstrap of the
closest reference, with a run-length null) and the 3SEQ/MaxChi pair (one tract of the
triplet walk); how far those three are independent of each other has not been measured.
```

(c) In the `--min-methods` paragraph (lines 28-41), after "The threshold is clamped to the number of callers actually run, and any suppressed count is logged rather than silently dropped." add: "Two of the default callers, 3SEQ and MaxChi, score the same tract (see the MaxChi section), so a region found by those two alone meets a gate of two without independent corroboration."

(d) In the barcode section, after "...and a run that selected only `barcode` is refused." add: "A barcode region carries a support but no p-value or q-value: the caller has no significance test, and its regions count toward `--min-methods` like any other caller's. Whether they should is an open decision (`docs/superpowers/specs/2026-10-01-min-methods-default-decision.md`)."

(e) `src/tessera/recomb/maxchi.py`, module docstring: replace from "but with an independent statistic." to "exactly as RDP treats MaxChi and 3SEQ as separate methods." with:

```
but with a different statistic. 3SEQ asks whether the minor-matching sites form a deep
contiguous *run* (maximum drawdown); MaxChi asks whether the **proportion** of minor
matches inside that tract differs from the rest -- a 2x2 chi-square. Both score the same
tract (``max_descent``), so they are not independent tests: their p-values are strongly
correlated under the clonal null (see docs/detection-methods.md), and a region found by
both rests on one triplet signal.
```

- [ ] **Step 4: Append to the decision record**

In `docs/superpowers/specs/2026-10-01-min-methods-default-decision.md`, insert before `## Decision`:

```markdown
## Addendum (2026-10-12, audit 4): what a vote is worth

Two findings of the fourth audit bear on this decision. Neither changes a default.

**3SEQ and MaxChi are one signal, not two (audit item A12).** MaxChi computes its
chi-square over the tract 3SEQ's maximum descent finds. On 300 simulated clonal walks
the two p-values have a Spearman correlation of <Step 2 value, measured 0.95>; at alpha
0.05 both reject in <0.030> of walks against <0.0012> expected under independence, and
MaxChi rejects in <0.90> of the walks where 3SEQ does. Both are individually calibrated
(<0.033> and <0.037> at a nominal 0.05). The two false regions that remained at
`--min-methods 2` in the table above were exactly such a pair. For option A this means
the agreement gate of two is weaker than it reads whenever the two agreeing callers are
3SEQ and MaxChi; a gate that counted them as one vote is a further option (A'), not
measured here.

**The barcode caller votes without a significance test (audit item A14).** A barcode
region has a support and no p-value, q-value or test, yet counts toward `--min-methods`.
The options, for the maintainer:

- keep it counting (status quo): a typed panel's barcode call corroborates like any
  other caller;
- record barcode in `methods` but do not count it toward the gate until it has a null
  (for example a permutation of marker positions), stated in the `--min-methods` help;
- give it a null first, then count it.

Nothing was measured for this: no harness selects the barcode caller, which is opt-in and
needs a typed panel. Measuring it needs `run_hybrids.py` with `--method all` on typed
panels.
```

Replace each `<...>` with the value printed in Step 2.

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: all green, `796 passed`. No harness gate: no code path changes.

- [ ] **Step 6: Commit**

```bash
git add docs/detection-methods.md src/tessera/recomb/maxchi.py docs/superpowers/specs/2026-10-01-min-methods-default-decision.md
git commit -m "$(cat <<'MSG'
Correct the claim that 3SEQ and MaxChi are independent votes

MaxChi computes its chi-square over the tract 3SEQ's maximum descent
finds, yet the documentation called the two genuinely independent votes
and a region significant under both corroborated twice. On 300 simulated
clonal walks their p-values correlate at 0.95 and both reject in 3.0 % of
walks against 0.12 % under independence; each is calibrated on its own.
The documentation and MaxChi's docstring now say so, and name the callers
that look at different evidence.

The open --min-methods decision record gains an addendum with these
numbers and with the options for the barcode caller, which counts toward
the agreement gate without a significance test (audit A14). No default
changes.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

---

### Task 16: Callers that could not test do not count toward the agreement gate; `--merge-gap` reaches the coverage diagnostic  (audit items: B7, gate half; B6, coverage half -- both handed over by plan F)

Plan F (`docs/superpowers/plans/2026-10-12-audit-f-reporting-cli.md`) reports these two cases faithfully but leaves alone the part that changes which regions are reported (its Global Constraints, the two "does not do" items). This task decides both.

**1. A one-state HMM stops counting toward `--min-methods`.** When lineage clustering pools every reference into one lineage, the HMM decodes one state and cannot call a region whatever the data; plan F Task 3 reports it as "not run" but keeps it in the count. Then `--method hmm,3seq --min-methods 2` can only be met by a caller that tested nothing: every 3SEQ region is suppressed and the run reads as a negative. Reproduced on plan G's candidate code with plan F's own panel (five references 0.25 % apart, a real 6 kb tract, `--no-informative-sites`): no region at all; with this task, the 3SEQ region is reported and the provenance reads `1 (requested 2; 1 caller ran)`. **Decision: it stops counting**, as an unrunnable barcode caller already does (plan B, item B6), through the same clamp and the same `gate_lowered` warning and provenance line. The effect is confined to that clamp: it changes a result only when the callers that can test number fewer than `--min-methods` (with the default four callers and `--min-methods 2` it changes nothing). A one-reference panel needs no handling: `run_recomb` refuses fewer than two references before any caller runs (`run.py:280`).

**2. `--merge-gap` is passed to the coverage diagnostic.** `run.py:443-447` builds `CoverageParams` without it, so coverage gaps are always merged within one window size whatever the user asks. **Decision: wire it through.** The default is unchanged (`None` means the window size in both places), so no run that does not set the option changes; a run that sets it now merges coverage gaps -- and so donor-absent regions, which are built from them -- with the gap the user gave, which is what the option's name says. Keeping it heuristic-only would leave one option meaning two things in one run.

**Depends on:** plan G Task 11 (`reference_clusters` in `run_recomb`). It is written to apply **with or without plan F**:
- **Plan F not yet landed** (the code as at Task 15 of this plan): the barcode tuple is `not_run` and `gate_lowered = bool(not_run) and ...` (`run.py:398`). Apply Step 3 as written.
- **Plan F Tasks 2 and 3 landed first**: F renames the barcode tuple to `callers_not_run` (F Task 2) and makes `not_run` a `dict[str, str]`; F Task 3 records `not_run["hmm"] = regions.HMM_ONE_STATE` from inside `_call_regions_hmm` under the same condition this task tests (the pooled view has fewer than two states), logs its own warning, and adds a comment "Only the barcode caller lowers the agreement gate ..." before `n_ran`. Then in Step 3 (a) use `callers_not_run` where the code below says `not_run`, compute `hmm_one_state = "hmm" in not_run` instead of from the clusters, drop this task's `logger.warning` (F's warning already names the caller), and delete F's comment. Also change, in F's files: the sentence "It still counts toward `--min-methods`." in `docs/detection-methods.md` to "It does not count toward `--min-methods`, as an unrunnable barcode caller does not."; in `tests/unit/test_not_run_reporting.py::test_barcode_and_one_state_hmm_are_both_reported_not_run` the docstring's "Only the barcode lowers the agreement gate (a one-state HMM still counts; see plan G)." to "Neither counts toward the agreement gate." and the assertion `assert run["min methods (agreement gate)"] == "2"` to `assert run["min methods (agreement gate)"] == "1 (requested 2; 1 caller ran)"`.
- **Plan F lands after this task**: the executor of plan F Task 3 makes the same three changes to F's text, test and comment, and keeps this task's code (reading `"hmm" in not_run` is then equivalent and may replace the cluster count).
- **Plan F Task 5** (help text) and this task touch the same `--merge-gap` help string. Whichever lands second uses this wording: `help="Merge same-parent regions (heuristic caller) and coverage gaps separated by at most this many bp. Default: the window size."`, and in `docs/detection-methods.md` F's sentence "the coverage diagnostic merges gaps within one window size whatever `--merge-gap` says" becomes "`--merge-gap` also sets how far apart two coverage gaps may be and still be merged (default: one window size)"; F's `RecombParams` comment "merge_gap is not passed to the coverage diagnostic, which merges gaps within one window size" becomes "merge_gap also merges coverage gaps"; and F's test `test_margin_and_merge_gap_are_marked_heuristic_only` drops its `--merge-gap` assertion.

**Files:**
- Modify: `src/tessera/recomb/run.py` (after `n_ran = ...`, line 378 before plan G; `gate_lowered`, line 398; `CoverageParams.with_defaults(...)`, lines 443-447)
- Modify: `src/tessera/cli/cmd_recomb.py` (`--merge-gap` help, lines 79-83)
- Modify: `docs/detection-methods.md` (the `--min-methods` paragraph, lines 28-41; the lineage-clustering limitation paragraph, lines 102-111; the `coverage_gaps.tsv` row of the output table)
- Create: `tests/unit/test_untestable_callers.py`

**Interfaces:**
- Consumes: `reference_clusters` (local in `run_recomb`, Task 11); `RecombParams.merge_gap`; `CoverageParams.with_defaults(window_size, floor=None, rel_drop=..., merge_gap: int | None = None)` (exists). With plan F: `ReportContext.not_run`'s source dict `not_run` and `callers_not_run` in `run_recomb` (F Task 2), `regions.HMM_ONE_STATE` (F Task 3).
- Produces: in `run_recomb`, a local `hmm_one_state: bool`; `n_ran` excludes it; `gate_lowered` covers it. No new public name.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_untestable_callers.py`:

```python
"""Callers that cannot test do not count toward the agreement gate; --merge-gap reaches the
coverage diagnostic.

Audit 4, items B7 (the gate half, handed over by plan F) and B6 (the coverage half).
"""

from __future__ import annotations

import json
import random
from pathlib import Path

from tessera.recomb import run as run_module
from tessera.recomb.run import RecombParams, run_recomb


def _rows(path: Path) -> list[dict[str, str]]:
    lines = [x for x in path.read_text().splitlines() if x.strip()]
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, x.split("\t"), strict=False)) for x in lines[1:]]


def _one_lineage_panel(tmp_path: Path) -> Path:
    """Five references 0.25 % apart and a query carrying a real 6 kb tract from one.

    Under base-pair windowing every pair stays above the 98.5 % pooling floor in every
    window, so lineage clustering pools them into a single HMM state. (The same panel
    as plan F's ``test_not_run_reporting._one_lineage_panel``.)
    """
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


def test_a_one_state_hmm_does_not_count_toward_the_agreement_gate(tmp_path, logger) -> None:
    # hmm + 3seq at --min-methods 2. The HMM is left one state and can call nothing, so
    # a gate of two could only be met by a caller that tested nothing: every 3SEQ region
    # was suppressed and the run read as a negative.
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=_one_lineage_panel(tmp_path), output=out, query="query",
                     plot_format="png", methods=("hmm", "3seq"), min_methods=2,
                     informative_sites=False),
        logger,
    )
    regions = _rows(out / "recombination_regions.tsv")
    called = [r for r in regions if r.get("donor_absent") != "yes"]
    assert called, "3SEQ's region is reported once the gate is lowered to the callers that ran"
    assert all(r["methods"] == "3seq" for r in called)
    gate = json.loads((out / "run_provenance.json").read_text())["run"][
        "min methods (agreement gate)"]
    assert gate == "1 (requested 2; 1 caller ran)"


def test_a_testable_hmm_still_counts(tmp_path, logger) -> None:
    # Guard: with clustering off the HMM has five states and counts as before.
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=_one_lineage_panel(tmp_path), output=out, query="query",
                     plot_format="png", methods=("hmm", "3seq"), min_methods=2,
                     informative_sites=False, cluster_lineages=False),
        logger,
    )
    gate = json.loads((out / "run_provenance.json").read_text())["run"][
        "min methods (agreement gate)"]
    assert gate == "2"


def test_merge_gap_reaches_the_coverage_diagnostic(tmp_path, logger, monkeypatch) -> None:
    seen: list[int] = []
    real = run_module.call_coverage_gaps

    def spy(result, window_size, params):
        seen.append(params.merge_gap)
        return real(result, window_size, params)

    monkeypatch.setattr(run_module, "call_coverage_gaps", spy)
    msa = _one_lineage_panel(tmp_path)
    for label, merge_gap in (("default", None), ("set", 250)):
        run_recomb(RecombParams(msa=msa, output=tmp_path / label, query="query",
                                plot_format="png", methods=("3seq",),
                                merge_gap=merge_gap), logger)
    assert seen == [1000, 250]  # the window size by default, the option when given
```

`test_a_testable_hmm_still_counts` is a guard (passes before and after).

- [ ] **Step 2: Run them, expected failure**

Run: `pytest tests/unit/test_untestable_callers.py -q`
Expected (observed on plan G's candidate code without this task): 2 failed -- `AssertionError: 3SEQ's region is reported once the gate is lowered to the callers that ran` / `assert []`, and `assert [1000, 1000] == [1000, 250]`.

- [ ] **Step 3: Implement**

(a) `src/tessera/recomb/run.py`, directly after `    n_ran = len(params.methods) - len(not_run)`:

```python
    # A one-state HMM: lineage clustering pooled every reference into one lineage, so the
    # HMM had no second state and could not call a region whatever the data (plan F
    # reports it as "not run"). Like an unrunnable barcode caller it does not count toward
    # the agreement gate: a gate it cannot help meet would only suppress what the callers
    # that did test found.
    hmm_one_state = "hmm" in params.methods and (
        len(result.similarities) - sum(len(c) - 1 for c in reference_clusters) < 2
    )
    if hmm_one_state:
        n_ran -= 1
        logger.warning(
            "The hmm caller was left a single state by lineage clustering and could not "
            "test; it does not count toward --min-methods."
        )
```

`reference_clusters` holds only multi-member clusters (Task 11), so `len(result.similarities) - sum(len(c) - 1 ...)` is the number of HMM states after pooling -- the condition plan F Task 3 tests on the pooled view.

(b) Replace

```python
    gate_lowered = bool(not_run) and min_agree < params.min_methods
```

with

```python
    gate_lowered = (bool(not_run) or hmm_one_state) and min_agree < params.min_methods
```

(c) In the `CoverageParams.with_defaults(` call, after `rel_drop=params.coverage_rel_drop,` add

```python
        # --merge-gap merges coverage gaps too (None -> the window size, as before).
        merge_gap=params.merge_gap,
```

(d) `src/tessera/cli/cmd_recomb.py`, the `--merge-gap` help becomes

```python
        help="Merge same-parent regions (heuristic caller) and coverage gaps separated by "
        "at most this many bp. Default: the window size.",
```

(e) `docs/detection-methods.md`:
- in the `--min-methods` paragraph, change "The threshold is clamped to the number of callers actually run" to "The threshold is clamped to the number of callers that could test -- a barcode caller without typed references and an HMM left a single state by lineage clustering do not count --";
- at the end of the lineage-clustering limitation paragraph (after "run with `--no-cluster-lineages` to keep the HMM's vote."), add "If clustering does leave the HMM a single state, it does not count toward `--min-methods`.";
- in the output table's `coverage_gaps.tsv` row, after "with a `kind`:", insert "gaps closer than `--merge-gap` (default one window size) are merged;".

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_untestable_callers.py tests/unit/test_barcode.py tests/unit/test_coverage.py tests/unit/test_run.py -q`
Expected: all pass.

- [ ] **Step 5: Full checks and harness gate**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q` (expected `799 passed`), then:

```bash
validation/data/audit-g/gate.sh task16
python validation/run_hybrids.py > validation/data/audit-g/task16/hybrids.txt 2>&1
```

Acceptance: every file identical to `task14`'s (Task 15 runs no gate; measured identical at 3 replicates on plan G's candidate code: no harness passes `--merge-gap`, and no harness panel leaves the HMM one state -- the specificity tips are 2.4 % apart, the near-identical panels are compared on informative columns, and the redundant panels keep three other lineages). Hybrids: identical verdicts. Any difference means the change reached a run that sets neither condition: revert and record.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/run.py src/tessera/cli/cmd_recomb.py docs/detection-methods.md tests/unit/test_untestable_callers.py
git commit -m "$(cat <<'MSG'
Do not count a one-state HMM toward --min-methods; pass --merge-gap to coverage

When lineage clustering pools every reference into one lineage the HMM has
a single state and cannot call a region. It still counted toward the
agreement gate, so --method hmm,3seq --min-methods 2 could only be met by
a caller that tested nothing: on five references 0.25 % apart with a real
6 kb tract every 3SEQ region was suppressed and the run read as negative.
Like an unrunnable barcode caller, it no longer counts; the gate is
clamped, logged and recorded in the provenance as for the barcode. With
the default four callers and --min-methods 2 nothing changes.

--merge-gap was never passed to the coverage diagnostic, which always
merged gaps within one window size. It is now; the default is unchanged,
so only runs that set the option change, and their coverage gaps and
donor-absent regions follow it.

Both items were left by plan F because they change reported regions.
Harness: unchanged (no harness sets either condition).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

---

### Task 17: Changelog, final gate, pull request  (all items)

**Files:**
- Modify: `CHANGELOG.md` (`## [Unreleased]`, lines 7-218)

**Interfaces:**
- Consumes: the gate outputs under `validation/data/audit-g/` from every earlier task.

- [ ] **Step 1: Changelog**

In `CHANGELOG.md`, under `## [Unreleased]`, add a paragraph after the opening paragraph and the bullets below to the existing sections. **Delete the bullet of any task whose gate failed and was reverted**, and replace every number with the one you measured at 10 replicates.

Opening paragraph (append after the existing one):

```markdown
Caller-statistics fixes from the fourth audit
(`docs/superpowers/specs/2026-10-12-audit-swot-gap.md`, items A1-A14, B2, B4). The
lineage-cluster changes (A1, A2), the descent floor (A3, A6), the last-window fix and the
Bootscan null can change which regions are reported on panels with near-duplicate
references or at the genome ends; the others change reported values only.
```

Under `### Fixed`, a new sub-heading `Callers (audit 4):` with:

```markdown
- **Near-duplicate references could silence the HMM.** A lineage cluster's per-window
  evidence was summed over its members, so three copies of the backbone outweighed a
  single-genome donor and the HMM lost a 1 kb tract (HMM among the callers 3/10 runs, 8/10
  without the copies). The evidence is now the member mean (<n>/10). (A1)
- **A redundant donor lineage was split into two regions and lost Bootscan.** Clusters are
  now computed once and shared: the merge treats two members of one cluster as one donor,
  and Bootscan counts support per cluster (one row <n>/10, Bootscan among the callers
  <n>/10, against 7/10 and 1/10). (A2)
- **The HMM sign test returned p = 0 above 1000 discordant sites.** It is now exact at every
  size, in log space, and never below the smallest normal double. (A9)
- **3SEQ discarded short tracts before testing them, and both triplet callers pruned the
  Benjamini-Hochberg family.** 3SEQ tests every donor with at least 20 discriminating sites
  (its exact DP is 10-40 times faster); MaxChi counts shallow donors in its family with
  p = 1. (A3, A6)
- **A merged region reported one member's identity over the union span**, and its
  `n_windows` was a maximum over different units. Both are restated over the merged span;
  `n_windows` counts base-pair windows for every caller. (A11, A8)
- **HMM `n_windows` and `posterior_support` described the untrimmed segment**; they now
  describe the reported span. Segment ends are padded by the median window step. (B4, A10)
- **Bootscan's permutation null dropped wins** when the window count was not a multiple of
  the block, giving anti-conservative p-values. (A5)
- **The last `(width - window) % step` columns were in no window.** A final window flush
  with the end is added.
- **A tract at the first discriminating site had a zero-width breakpoint interval** (3SEQ,
  MaxChi, GENECONV). (A14)
- **Memory:** the PHI gamete matrices, the Hudson-Kaplan pair list and GENECONV's
  permutations no longer exceed the scaling caps (136 MB to 25 MB, 381 MB to 0.25 MB, 680
  MB to 250 MB on the test inputs); MaxChi's null no longer overflows int64 at bacterial
  scale. (A13, A7)
```

Under `### Changed`:

```markdown
- **`--min-region`, when given, bounds every reported region**, whichever caller found it.
  Unset, only the HMM applies a floor, to its segments. The help text said otherwise. (B2, A4)
- **An HMM left a single state by lineage clustering no longer counts toward
  `--min-methods`**, as an unrunnable barcode caller does not: a gate it cannot help meet
  only suppressed the regions of the callers that did test. (B7)
- **`--merge-gap` also sets how far apart coverage gaps may be and still be merged**, and
  so shapes donor-absent regions; the default (one window size) is unchanged. (B6)
```

Under `### Documentation`:

```markdown
- 3SEQ and MaxChi are no longer described as independent votes: they score the same tract
  and their p-values correlate at 0.95 under the clonal null. The open `--min-methods`
  decision record gains these numbers and the options for the barcode caller, which votes
  without a significance test. (A12, A14)
```

Under `### Added`:

```markdown
- `validation/run_regimes.py` simulates panels with a redundant backbone or donor lineage.
```

- [ ] **Step 2: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: all green, `799 passed`.

- [ ] **Step 3: Final gate**

```bash
validation/data/audit-g/gate.sh final
python validation/run_hybrids.py > validation/data/audit-g/final/hybrids.txt 2>&1
diff validation/data/audit-g/baseline/validation_regions.tsv validation/data/audit-g/final/validation_regions.tsv
```

- [ ] **Step 4: Commit the changelog**

```bash
git add CHANGELOG.md
git commit -m "$(cat <<'MSG'
Record the audit-4 caller-statistics changes in the changelog

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
MSG
)"
```

- [ ] **Step 5: Pull request (only when the maintainer asks for it)**

Push and open the PR against `fix-audit-caller-behaviour` (or `main` if that branch has merged). The description:

```markdown
## Summary
<one line per task kept; "reverted: <task>, gate result <numbers>" per task reverted>

## Harness results, before and after
<baseline and final: spec_mm1 / spec_mm2 TOTAL and positive-control lines; every
regime line at both gates; run_validation verdicts and the validation_regions diff;
hybrids summary>

## Gate results per task
<the acceptance lines of each task's gate, with the measured values>

## Skipped
<every gate that could not run and why -- e.g. "run_hybrids.py: no network"; Tasks 10,
11, 12, 14 and 16 must not merge before it has run>

## Found, not fixed
<Self-review notes 3, 4 and 5 of the plan, with their numbers>

## Decision needed
<the --min-methods decision record: the A12/A14 addendum; the half-window measurement>

<the pull-request attribution line required by the session's instructions>
```

---

## Self-review notes

Every task's tests were written in a scratch worktree of `4b35e39`, run against the unchanged source to confirm the failure recorded in its Step 2 (27 tests fail there for the stated reasons; the guards named in each task pass), then run against the cumulative candidate: `796 passed, 5 deselected` (799 with Task 16), `ruff` and `mypy` clean.

1. **A4 is handled through B2, not separately.** Keeping the sign test on the untrimmed segment is plan C's recorded decision (its deviation 1). What was wrong was the option's contract: Task 13 makes `--min-region`, when given, bound every reported region, HMM rows included, and states the default behaviour in the help and the documentation.

2. **A14, first half (the barcode caller votes without a significance test), is a maintainer decision.** Changing what counts toward `--min-methods` changes what the agreement gate means, and nothing in the harnesses exercises the barcode caller (opt-in, needs a typed panel). Task 15 records the options in the open decision record; no task changes the count. The second half of A14 (zero-width terminal breakpoint intervals) is fixed in Task 9, which also covers `maxchi.py:170`, the same line the spec did not list.

3. **A2 interacts with plan C's C1 (commit `a53ea53`), and part of the redundancy case remains.** Under informative-site windowing, clustering compares a pair on informative columns only, in windows of 40 sites, at the 98.5 % floor; one differing informative site in a window (97.5 %) keeps a pair apart. Measured on a near-identical panel (`simulate_near_identical`, 30 kb, 0.4 % divergence) with a copy of the donor 0.03 % away from it: `lineage_clusters` returned no cluster, and the event was reported as two rows naming the two copies (`hmm,maxchi,bootscan` and `3seq`), before and after Task 11. An exact copy is pooled and gives one row (Review Focus 3). Remedying this means choosing C1's window size or floor for informative columns, which needs the near-identical regimes and `run_hybrids.py` (mpox, VZV) as its gate; it is not attempted here. A1 does not interact with C1 beyond this: under informative-site windowing the clusters it averages are the ones C1 forms.

4. **The literal reading of A10 (pad by half the window width) changes the HMM's significance, not only its coordinates.** Measured in "What was measured": false regions at `--min-methods 1` fall from 7/12 to 1/12 on the specificity panels, but the HMM loses its calls on the near-identical panels (3/3 to 0/3) and the positive-control breakpoint error rises from 13 bp to 208 bp. It belongs to option B of the open decision (a scan-aware HMM gate) and is reported there, not implemented. Task 4 replaces only the single sampled gap by the median step.

5. **GENECONV has the same pre-test floor and family pruning as A3/A6** (`geneconv._MIN_RUN = 4`, candidates dropped before BH). It is opt-in and the spec does not name it, so it is left alone; the fix would mirror Task 12's MaxChi rule.

6. **MaxChi stays blind to tracts with a descent shallower than four sites** (Task 12). Removing its floor and permuting every donor made a 16-genome scan ten times slower (4.3 s to 45 s) for donors whose exact 3SEQ p was at least 0.21 on every clonal panel measured (144 donors). Counting them in MaxChi's family with p = 1 fixes A6 for MaxChi; A3 is covered by 3SEQ, which scores the same tract (A12).

7. **A6's reachability, measured:** on the 12 clonal specificity panels the floor removed 144 donors from the families, all with p >= 0.21, and no region changed at either gate. The inflation the spec describes is real arithmetic (Task 12's test) but did not change a call on any simulated panel.

8. **A7 is SUSPECTED.** Task 7 is a widening with a guard test at 120 000 sites, as allowed; the overflow was not reached through `_exceedances`.

9. **The spec's A5 example could not be reproduced as stated.** With 12 wins at positions 40-51 and block 10, the old code gives p 0.0920 at `n_win` 100 and 0.0970 at 103 (it rises, the spec reports a fall from 0.0930 to 0.0855); the spec does not say where its wins were. The defect itself is confirmed by Task 6's test (all wins, `n_win` 11, block 10: p 0.49 where it must be 1).

10. **A10's informative-site numbers (33 columns against a median of 237) and A13's 892 MB at the cap were not re-measured on the audit's panels;** the mechanisms were confirmed with the inputs in Tasks 4 and 8 (136 MB at z = 2000 scales to about 850 MB at z = 5000).

11. **Clustering now runs whenever `--cluster-lineages` is on (Task 11), also when the HMM is not among `--method`.** The cost is the existing all-pairs comparison, bounded by `cluster_max_refs = 200`.

12. **Plan F handed over two items (Task 16).** A one-state HMM now leaves the `--min-methods` count, as an unrunnable barcode caller does; `--merge-gap` now reaches the coverage diagnostic, with the default unchanged. Both were verified on the plan's candidate code in a scratch worktree (the change moved by `git diff` / `git apply`): the two failing tests fail as recorded, then pass; the fast suite is otherwise unchanged; the harnesses are identical (no harness sets either condition). Plans F and G edit the same `--min-region` and `--merge-gap` help strings, `RecombParams` comments and one `docs/detection-methods.md` paragraph (F Task 5 against G Tasks 4, 13 and 16), and F Task 3's test asserts the old gate count; each of those G tasks says what to change in F's text whichever plan lands first.

13. **Not covered by this plan's tests:** a real panel with near-duplicate references (A1 and A2 were measured on simulations only), and any change to `run_validation.py` or `run_hybrids.py` results; both are executor gates. B4's report wording, B1, B3 and B5-B9 belong to the reporting plan.
