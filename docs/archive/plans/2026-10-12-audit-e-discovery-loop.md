# Audit Plan E: Panel Discovery Loop Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the reference-recruitment loop (`detect`, `build-panel`, `fill-references`, `find-references`) download what its documentation says, stop on the gaps it searches, keep the lineage labels it already has, and say plainly when a requested step did not run or an option cannot act.

**Architecture:** Every change is local to `src/tessera/discover/` and the three panel commands' CLI modules. The stop test moves from the worst gap to a coverage shortfall over all gaps (Task 1); the download keeps one reference per gap and, inside the fill loop, excludes what curation removes so the next candidate gets a turn (Task 2); the remaining tasks are independent fixes to seeding, typing, pool selection, caching and option handling. All tests mock BLAST, efetch, the `datasets` CLI, skani and Nextclade; nothing in this plan touches the network or an aligner.

**Tech Stack:** Python 3.11+, numpy, pytest, ruff, mypy. No new dependency.

**Spec:** `docs/superpowers/specs/2026-10-12-audit-swot-gap.md` (section 1.4 items D1-D10; from section 1.5 item E5 the `--candidate-pool`/`--nextclade` precedence, the inert `--pool-consensus`, and the candidate-pool reader's extension list; the "Discover side" test gaps of section 3).

## Global Constraints

- No new runtime dependency.
- Modest scientific language in code, docs and messages; reported numbers must be faithful (state what passes, what fails, what was skipped).
- "Could not test" must never be reported as "tested, found nothing".
- A new behaviour needs a test that fails without the change; if a test is a guard that passes before and after, the plan says so.
- Branch: create `fix-audit4-discovery` from `fix-audit-caller-behaviour` (or from `main` once that branch is merged). Commit messages explain the problem, the evidence and the trade-off, and end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Append the aligner env to PATH, never prepend (prepending swaps the Python interpreter): `export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"`.
- The dev Mac has 18 GB RAM: run one harness at a time.
- Each task ends with `ruff check src tests validation`, `mypy`, and `pytest -m "not requires_binary" -q` green before its commit.
- No Unicode in code or docs beyond what files already contain; ASCII punctuation.
- **Plan C Task 6 is held in `git stash` and in `validation/data/audit-c/task6/task6.patch`.** It edits the same lines of `_type_panel` and `_select_from` as Task 5 here. Do not pop or drop that stash. If it lands before this branch, Task 5's edit becomes "add `datasets_rows=datasets_rows,` beside the `pool_labels=` argument it adds"; the two do not conflict in meaning (pool labels and NCBI datasets lineages are different sources).
- **CHANGELOG.** Every task adds its entry to `CHANGELOG.md` under `## [Unreleased]`. Fixes go in a group headed `Panel recruitment (fourth audit):` at the end of `### Fixed` (after the `Callers:` group); Task 1 creates the group. Behaviour changes a user will notice go as bullets in `### Changed`. If another fourth-audit plan has already added a group for its own items, put this group after it.

### Before Task 1

```bash
git switch fix-audit-caller-behaviour
git switch -c fix-audit4-discovery
pip install -e ".[dev]"
ruff check src tests validation && mypy && pytest -m "not requires_binary" -q
```

Expected: `All checks passed!`, `Success: no issues found in 79 source files`, `757 passed, 5 deselected`. If any differs, stop: the baseline is not the one this plan was verified against (commit `4b35e39`).

### What gates these changes

These tasks change what a recruited panel contains, so the question "which harness shows it" has to be answered per task:

- `validation/run_hybrids.py` builds each case's pool from the Nextclade tree and calls `select_regional` directly (`validation/run_hybrids.py:1093`), with `drop_siblings=False` and a selection window of at most a quarter of the genome (`window_params`, `validation/run_hybrids.py:202-211`). It never calls `fill_references`, `_download`, `_seed_windowed`, `_fetch_diverse`, `_type_panel` or the CLI. Task 7 changes only the `drop_siblings=True` branch and Task 8 only queries shorter than two selection windows, so neither is reachable from the harness either. **The hybrid harness cannot gate any task in this plan, and none of them changes its output.** Running it is not required by `CONTRIBUTING.md` for these tasks (no caller, default or region-calling change).
- `validation/run_specificity.py`, `run_validation.py` and `run_regimes.py` start from an alignment or a fixed collection; none recruits a panel.
- **The gate is the unit and integration tests named in each task**, all offline. They are listed per task under Step 4.
- Tasks 1 and 2 change what `detect` and `build-panel` recruit on real data. No offline harness measures that. Task 2 Step 7 asks for one paired real-data run before merge (network and aligner needed); if it cannot be run, the PR says so under "Skipped" rather than claiming the change is neutral.

## Review Focus

Inputs the spec implies but its examples do not exercise, most likely to bite first. Each has a test in the task that owns the code.

1. **A gap whose best candidate is a sibling of the query, with curation on** (the `detect` preset always curates). One reference per gap then downloads the sibling, curation removes it, and the gap gets nothing unless the removed genome is excluded and the next round searches again. *Task 2, `test_a_download_removed_by_curation_is_excluded_and_searched_past`.*
2. **A BLAST hit that is a short HSP at 100 % identity.** Ranking by `pct_identity` alone would download a genome that matched 60 bp of a 400 bp sub-interval. *Task 2, the `SHORT` candidate in `test_download_takes_one_reference_per_gap`.*
3. **`fill-references --candidate-pool dir` without `--seed-source local`.** The pool was read only under `local`; the default BLAST seeding ran and the pool was ignored without a word. *Task 11, `test_fill_references_refuses_a_candidate_pool_under_another_seed_source`.*
4. **A short query whose halves skani cannot match.** skani is a genome-scale tool and may report nothing for a 750-1450 bp chunk; splitting must not turn a working single-window selection into an empty panel. *Task 8, `test_select_regional_falls_back_to_the_whole_query_when_no_half_matches`.*
5. **`--lineage-map` together with the NCBI lineage sidecar on the default typing path.** Once the sidecar is read there, the user's map must still win. *Task 5, `test_a_user_lineage_map_still_overrides_the_sidecar_on_the_default_path`.*

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `src/tessera/discover/iterate.py` | modify | stop test (T1), curation exclusion and unchanged-panel round (T2), saturation (T3), empty fetch (T4), typing (T5), curation note (T6), cache key (T9), seed-option checks (T11) |
| `src/tessera/discover/run.py` | modify | one reference per gap (T2) |
| `src/tessera/discover/pool.py` | modify | empty package (T4), all-sibling warning (T7), short-query windows (T8), pool reader (T11) |
| `src/tessera/discover/lineage_assign.py` | modify | non-clade markers (T10) |
| `src/tessera/cli/cmd_fill_references.py` | modify | `--min-improvement` help and range (T1) |
| `src/tessera/cli/cmd_build_panel.py` | modify | shared seed-source rule (T11), `--candidate-pool` help |
| `src/tessera/cli/cmd_detect.py` | modify | `--candidate-pool` and `--pool-consensus` help (T11) |
| `tests/unit/test_iterate.py` | modify | T1-T6, T9 tests |
| `tests/unit/test_discover.py` | modify | T2 tests; two existing tests change (T2) |
| `tests/unit/test_pool.py` | modify | T7, T8, T11 tests |
| `tests/unit/test_lineage_assign.py` | modify | T10 test |
| `tests/unit/test_cli_input_checks.py` | modify | T1, T11 CLI tests |
| `tests/unit/test_cli_nextclade.py` | modify | T11 guard |
| `docs/reference-panels.md`, `CHANGELOG.md` | modify | each task documents its own change |

---

### Task 1: Stop on the coverage shortfall, not the worst gap  (audit items: D2)

The loop searches only the `--top-gaps` longest gaps (`iterate.py:445`) but stopped when the minimum best-similarity over **all** gaps failed to rise (`iterate.py:434`). A short gap outside the searched set pins that minimum, so the loop stopped after round 2 while the searched gaps were improving.

**Choice of measure.** The lead offered two candidates: the under-covered base count already in `RoundResult.undercovered_bp`, or the mean over the searched gaps. Neither is used:
- The base count does not move while a gap improves but stays below the threshold, because the gap keeps its length. The spec's own reproduction (three long gaps from about 0.70 to about 0.90 under a 0.94 threshold) is that case; the test below shows the base count constant at 6200 bp for three rounds of improvement.
- The mean over the searched gaps compares different gap sets: when the longest gap closes, the next-longest enters the set with its own, often lower, similarity and reads as a regression.
- The **coverage shortfall**, `sum(length_bp * (threshold - mean_best))` over all gaps, falls both when a gap closes and when it improves without closing; an unsearched gap adds a constant and cannot hold it still. It needs only fields the gaps already carry.

`--min-improvement` keeps its default value (0.01) but becomes a fraction of the previous round's shortfall. That is a change of meaning for a user-facing option, so the help, the docs and `### Changed` say so, and the CLI now rejects values outside 0-1.

**Files:**
- Modify: `src/tessera/discover/iterate.py` (`FillParams.min_improvement` line 106; import line 33; `RoundResult` lines 180-187; `_grow_collection` lines 387 and 422-440; `_progress_section` caption lines 1024-1029; `_write_trace` lines 1036-1044)
- Modify: `src/tessera/cli/cmd_fill_references.py` (help lines 105-108; validation after line 213)
- Modify: `docs/reference-panels.md` (lines 104-107), `CHANGELOG.md`
- Test: `tests/unit/test_iterate.py`, `tests/unit/test_cli_input_checks.py`

**Interfaces:**
- Produces: `coverage_shortfall(gaps: list[CoverageGap], threshold: float) -> float` in `tessera.discover.iterate`.
- Produces: `RoundResult.shortfall: float = 0.0` (new field between `threshold` and `added`; every existing construction is positional up to `threshold` or uses keywords, so nothing else changes).
- Produces: `fill_summary.tsv` gains a last column `coverage_shortfall` (appended, so existing column positions are unchanged).
- Consumes: `CoverageGap.length_bp`, `CoverageGap.mean_best`, and the `threshold` returned by `call_coverage_gaps`. Every gap's windows are below the threshold by construction (`coverage.py`), so each term is positive; `max(0.0, ...)` guards rounding only.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_iterate.py`:

```python


# --- audit 4: the stop test follows the gaps the loop searches (D2) -----------------

def _sized_gap(start: int, length: int, mean_best: float) -> CoverageGap:
    return CoverageGap(
        msa_start=start, msa_end=start + length, query_start=start, query_end=start + length,
        length_bp=length, n_windows=max(1, length // 100), best_label="refA",
        mean_best=mean_best, kind="divergent",
    )


def test_stop_test_is_not_pinned_by_a_gap_the_search_never_reaches(
    monkeypatch, tmp_path, logger
):
    """Three long gaps are searched every round and improve (0.70, 0.80, 0.90) without
    clearing the 0.94 threshold, so their length does not change. A short fourth gap is
    below --top-gaps, never searched, and stays at 0.55. The old test compared the worst
    gap (always 0.55) and stopped after round 2; an under-covered base count would also
    stop, because the long gaps keep their length while they improve."""
    query, coll, out = _setup(tmp_path)

    def gaps(long_best: float) -> list[CoverageGap]:
        return [_sized_gap(0, 2000, long_best), _sized_gap(3000, 2000, long_best),
                _sized_gap(6000, 2000, long_best), _sized_gap(9000, 200, 0.55)]

    _common_mocks(monkeypatch, [
        (gaps(0.70), 0.94), (gaps(0.80), 0.94), (gaps(0.90), 0.94), ([], 0.95),
    ])
    _one_new_hit_per_round(monkeypatch)

    trace = fill_references(
        FillParams(query=query, collection=coll, output=out, max_rounds=5, top_gaps=3),
        logger,
    )

    assert [r.round for r in trace] == [1, 2, 3, 4]
    assert trace[-1].n_gaps == 0
    assert [r.undercovered_bp for r in trace[:3]] == [6200, 6200, 6200]
    assert [r.shortfall for r in trace[:3]] == [1518.0, 918.0, 318.0]


def test_fill_summary_records_the_coverage_shortfall(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    _common_mocks(monkeypatch, [([_gap(0.84)], 0.94), ([], 0.95)])
    _one_new_hit_per_round(monkeypatch)

    fill_references(FillParams(query=query, collection=coll, output=out), logger)

    header, first, _ = (out / "fill_summary.tsv").read_text().splitlines()
    row = dict(zip(header.split("\t"), first.split("\t"), strict=True))
    assert row["coverage_shortfall"] == "50.0"  # 500 bp x (0.94 - 0.84)
```

(`_one_new_hit_per_round` is defined further up the same module; it is resolved at call time.)

Append to `tests/unit/test_cli_input_checks.py`:

```python


@pytest.mark.parametrize("value", ["-0.1", "1.5"])
def test_fill_references_rejects_a_min_improvement_outside_zero_one(
    tmp_path: Path, no_pipeline, value: str
) -> None:
    """--min-improvement is a fraction of the previous round's coverage shortfall."""
    result = runner.invoke(app, [
        "fill-references", "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        "--min-improvement", value,
    ])
    _rejected(result, "Invalid --min-improvement")
```

- [ ] **Step 2: Run them, expected failure**

```bash
pytest tests/unit/test_iterate.py -k "pinned or shortfall" -q
pytest tests/unit/test_cli_input_checks.py -k min_improvement -q
```

Observed at `4b35e39`:
- `test_stop_test_is_not_pinned_by_a_gap_the_search_never_reaches`: `assert [1, 2] == [1, 2, 3, 4]` (the loop stopped at round 2).
- `test_fill_summary_records_the_coverage_shortfall`: `KeyError: 'coverage_shortfall'`.
- both CLI cases: `assert 'Unexpected error' not in ...`, output `[ERROR] Unexpected error: the command reached the pipeline; it should have been rejected`.

- [ ] **Step 3: Implement**

`src/tessera/discover/iterate.py`, the import:

```diff
-from ..recomb.coverage import CoverageParams, call_coverage_gaps
+from ..recomb.coverage import CoverageGap, CoverageParams, call_coverage_gaps
```

`FillParams`:

```diff
-    min_improvement: float = 0.01  # stop if the worst gap's best-sim gains less than this
+    # Stop once a round reduces the coverage shortfall by less than this fraction of the
+    # previous round's shortfall (see coverage_shortfall).
+    min_improvement: float = 0.01
```

`RoundResult`, and the new function directly after it:

```diff
     worst_best_sim: float
     threshold: float
+    shortfall: float = 0.0  # coverage_shortfall of this round's gaps
     added: list[str] = field(default_factory=list)
 
 
+def coverage_shortfall(gaps: list[CoverageGap], threshold: float) -> float:
+    """Under-covered length weighted by how far each gap falls below the threshold.
+
+    The sum over gaps of ``length_bp * (threshold - mean_best)``. It falls when a gap
+    closes and also when a gap's best match improves without yet clearing the threshold,
+    which an under-covered base count does not register. A gap that no search reaches
+    adds the same amount every round, so it cannot hold the measure still the way it
+    held the worst-gap similarity still.
+    """
+    return sum(g.length_bp * max(0.0, threshold - g.mean_best) for g in gaps)
+
+
 def fill_references(params: FillParams, logger: logging.Logger) -> list[RoundResult]:
```

`_grow_collection`:

```diff
     built_from: set[str] = set()  # labels the last MSA was built from
-    prev_worst: float | None = None
+    prev_shortfall: float | None = None
```

```diff
         undercovered = sum(max(0, g.query_end - g.query_start) for g in gaps)
         worst = min((g.mean_best for g in gaps), default=1.0)
-        rr = RoundResult(rnd, len(gaps), undercovered, round(worst, 4), round(threshold, 4))
+        shortfall = coverage_shortfall(gaps, threshold)
+        rr = RoundResult(
+            rnd, len(gaps), undercovered, round(worst, 4), round(threshold, 4),
+            round(shortfall, 1),
+        )
         trace.append(rr)
         logger.info(
-            "Round %d: %d gap(s), %d bp under-covered, worst best-similarity %.3f.",
-            rnd, len(gaps), undercovered, worst,
+            "Round %d: %d gap(s), %d bp under-covered, coverage shortfall %.1f, "
+            "worst best-similarity %.3f.",
+            rnd, len(gaps), undercovered, shortfall, worst,
         )
 
         if not gaps:
             logger.info("Converged: no coverage gaps remain.")
             break
-        if prev_worst is not None and (worst - prev_worst) < params.min_improvement:
+        # Judged on the shortfall, not on the worst gap: only the --top-gaps longest gaps
+        # are searched, so a short gap outside that set held the worst-gap similarity
+        # still and stopped the loop while the searched gaps were improving.
+        if (
+            prev_shortfall is not None
+            and (prev_shortfall - shortfall) < params.min_improvement * prev_shortfall
+        ):
             logger.info(
-                "Stopping: coverage stopped improving (worst best-sim %.3f vs %.3f last round).",
-                worst, prev_worst,
+                "Stopping: coverage stopped improving (shortfall %.1f against %.1f last "
+                "round; a reduction below %g of it).",
+                shortfall, prev_shortfall, params.min_improvement,
             )
             break
-        prev_worst = worst
+        prev_shortfall = shortfall
```

`_progress_section` caption:

```diff
         f'ending with {final_size} references. The run stops when the gaps close, no new '
-        'reference is found, or a round no longer improves the worst gap.</p>'
+        'reference is found, or a round no longer reduces the coverage shortfall (the '
+        'under-covered length weighted by how far each gap falls below the '
+        'threshold).</p>'
```

`_write_trace`:

```diff
     header = ["round", "n_gaps", "undercovered_bp", "worst_best_similarity",
-              "threshold", "references_added"]
+              "threshold", "references_added", "coverage_shortfall"]
 ...
                 r.round, r.n_gaps, r.undercovered_bp, r.worst_best_sim,
-                r.threshold, ",".join(r.added) or "-",
+                r.threshold, ",".join(r.added) or "-", r.shortfall,
             ])) + "\n")
```

`src/tessera/cli/cmd_fill_references.py`:

```diff
     min_improvement: float = typer.Option(
         0.01, "--min-improvement",
-        help="Stop once a round improves the worst gap's best similarity by less than this.",
+        help="Stop once a round reduces the coverage shortfall (under-covered length "
+        "weighted by how far each gap falls below the threshold) by less than this "
+        "fraction of the previous round's shortfall (0-1).",
     ),
```

```diff
         _require_range(max_rounds, "--max-rounds", lo=1)
+        _require_range(min_improvement, "--min-improvement", lo=0.0, hi=1.0)
         _require_scan_windows(window_size, window_step)
```

`docs/reference-panels.md`, the paragraph after the `fill-references` example (lines 104-107):

```diff
-It stops when the gaps close, when no new reference can be found, or when a round no
-longer improves the worst gap (so a genuinely hypervariable region is reported, not
-chased forever). Each round is recorded in `filled/fill_summary.tsv`, and the query's
-own record is auto-excluded from its FASTA header. This needs an aligner and Entrez
+It stops when the gaps close, when no new reference can be found, or when a round no
+longer reduces the coverage shortfall -- the under-covered length weighted by how far
+each gap falls below the threshold -- by at least `--min-improvement` (default 0.01, a
+fraction of the previous round's shortfall), so a genuinely hypervariable region is
+reported, not chased forever. Only the `--top-gaps` longest gaps are searched each round;
+a shorter gap outside that set adds the same amount every round and does not stop the
+loop. Each round, with its `coverage_shortfall`, is recorded in
+`filled/fill_summary.tsv`, and the query's own record is auto-excluded from its FASTA
+header. This needs an aligner and Entrez
```

`CHANGELOG.md`, at the end of `### Fixed` (create the group):

```markdown
Panel recruitment (fourth audit):

- `fill-references` stopped after round 2 whenever there were more coverage gaps than
  `--top-gaps`: the stop test watched the worst gap, while the search reached only the
  longest ones, so an unsearched short gap held the measure still.
```

and in `### Changed`:

```markdown
- **The fill loop stops on the coverage shortfall.** A round must reduce the under-covered
  length, weighted by how far each gap falls below the threshold, by at least
  `--min-improvement` of the previous round's value. **`--min-improvement` is now a
  fraction (0-1, default unchanged at 0.01)**, where it was an absolute gain in the worst
  gap's best similarity. `fill_summary.tsv` gains a last column, `coverage_shortfall`.
```

- [ ] **Step 4: Run, expected pass**

```bash
pytest tests/unit/test_iterate.py tests/unit/test_cli_input_checks.py -q
```

The gate for this task is these tests plus the existing `test_loop_stops_when_coverage_stalls` (a guard: an unchanged gap still stops the loop), `test_last_round_downloads_are_aligned` and `test_loop_converges_when_gaps_close`, all in `tests/unit/test_iterate.py`.

- [ ] **Step 5: Full checks**

```bash
ruff check src tests validation && mypy && pytest -m "not requires_binary" -q
```

Expected: `761 passed, 5 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/discover/iterate.py src/tessera/cli/cmd_fill_references.py \
        tests/unit/test_iterate.py tests/unit/test_cli_input_checks.py \
        docs/reference-panels.md CHANGELOG.md
git commit -F - <<'EOF'
Stop the fill loop on the coverage shortfall, not the worst gap

The loop searches only the --top-gaps longest gaps but stopped when the
minimum best-similarity over all gaps failed to rise. A short gap outside
the searched set pins that minimum, so a run with more gaps than
--top-gaps stopped after round 2 while the searched gaps improved
(fourth audit, item D2).

The stop test now compares the coverage shortfall: the sum over gaps of
length times the distance below the threshold. It falls when a gap closes
and when a gap improves without closing; an unsearched gap adds a
constant. The under-covered base count was considered and not used: it
stays constant while a gap improves below the threshold, which is the
audit's own reproduction (6200 bp for three improving rounds in the
test).

Trade-off: --min-improvement changes meaning, from an absolute gain in
the worst gap's similarity to a fraction of the previous shortfall. The
default value (0.01) is kept, values outside 0-1 are refused, and the
help, docs and changelog say so. fill_summary.tsv gains a last column.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 2: One new reference per gap  (audit items: D1; the recorded item "curated-away siblings are re-downloaded each round")

`_download` (`run.py:272-292`) fetched every deduplicated new accession, although its docstring, the `fill-references` module docstring (`iterate.py:4-5`), `find-references --download`'s help (`cmd_find_references.py:58-61`) and `docs/reference-panels.md:66-70` all say one per gap. Sub-interval searching (`_tiles`, `run.py:125-139`) multiplies the candidates: at the defaults (`--max-hits 5`, `--subtile 400`) a 1000 bp gap is searched in 4 sub-intervals and can add up to 20 genomes, a 6 kb gap in 29 and up to 145, per gap per round.

**Which commands change.** `_download` is shared, so `find-references --download`, `fill-references`, `detect` and `build-panel` all change. `find-references`' documentation promises "the best new one", not a bulk download, so it is not kept bulk. `candidate_references.tsv` still lists every candidate.

**How the best is chosen.** Highest identity over the sub-interval the hit was found in: `pct_identity * query_coverage / 100`, then `pct_identity`, then BLAST order. The lead's suggestion was the highest-identity hit; that was not taken as is, because `pct_identity` is measured over the HSP only and a 60 bp HSP at 100 % would outrank a hit spanning the sub-interval at 95 % (Review Focus item 2). If the best fetch fails, the gap's next candidate is tried. A gap whose best candidate was already fetched for an earlier gap is served by it.

**The interaction this creates, and the second half of the task.** With one per gap, the best candidate is often a sibling of the query (it matches everywhere). Under `--curate`, which the `detect` and `build-panel` presets always set, curation then removes it. Before this change the gap's other candidates arrived in the same round; now nothing would, and worse, the removed genome was never excluded (`exclude` is not extended), so every later round would pick it again, and the unchanged panel would read as a stall and stop the loop. Two changes in `_grow_collection` close that:
- every genome curation drops is added to `exclude`, so later searches move on to the gap's next candidate;
- a round whose previous round's downloads were all curated away is not judged by the stall test (the panel is the one the last round saw, so an unchanged shortfall says nothing about whether a search could help); it searches again.

Reviewers should approve or reject these together with the cap: the cap alone would regress `detect` on exactly the masking case curation exists for.

**Files:**
- Modify: `src/tessera/discover/run.py` (`_download`, lines 272-292)
- Modify: `src/tessera/discover/iterate.py` (`_grow_collection`: state and `curate_if_pending` at lines 389-401; the stall test from Task 1)
- Modify: `docs/reference-panels.md` (lines 66-70, 85-88, 297-301), `CHANGELOG.md`
- Test: `tests/unit/test_discover.py`, `tests/unit/test_iterate.py`
- **Existing tests that change:** `tests/unit/test_discover.py::test_curate_after_download_keeps_the_users_own_genomes` (lines 153-207) and `::test_curate_after_download_into_a_new_directory` (lines 210-253) pin the bulk download: each gives the gap two hits, `NEW123` (90 % over 95 %) and `NEWSIB` (97 % over 96 %), and expected both downloaded with `NEWSIB` then curated away. With one per gap only `NEWSIB` is downloaded. Their intent (the user's genomes are kept, a downloaded sibling is dropped) is unchanged; only the expected file lists change (Step 3).

**Interfaces:**
- Produces in `tessera.discover.run`: `_subinterval_identity(candidate: Candidate) -> tuple[float, float]`, `_ranked_new_by_gap(candidates: list[Candidate]) -> list[list[Candidate]]`. `_download(candidates, dest, logger) -> list[Candidate]` keeps its signature.
- Consumes: `Candidate.gap` (`msa_start`, `msa_end` identify a gap), `Hit.pct_identity`, `Hit.query_coverage`; curation rows from `curate_collection_dir` with `role` ending in `-dropped` (`panel.py:244-252`); `_base_accession` from `run.py`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_discover.py`:

```python


# --- audit 4: one new reference per gap (D1) ---------------------------------------

def _gap_at(start: int):
    from tessera.recomb.coverage import CoverageGap

    return CoverageGap(
        msa_start=start, msa_end=start + 1000, query_start=start, query_end=start + 1000,
        length_bp=1000, n_windows=10, best_label="ref", mean_best=0.8, kind="divergent",
    )


def _cand(gap, accession: str, identity: float, coverage: float, *, present: bool = False):
    from tessera.discover.run import Candidate

    return Candidate(gap, Hit(accession, accession, identity, coverage, 1e-9), present)


def _recording_efetch(monkeypatch, fetched: list[str], *, failing: tuple[str, ...] = ()):
    def fake_efetch(accession, dest, logger):
        fetched.append(accession)
        if accession in failing:
            raise RuntimeError("efetch timed out")
        path = dest / f"{accession}.fasta"
        path.write_text(f">{accession}\nACGT\n")
        return path

    monkeypatch.setattr(discover_run, "efetch_fasta", fake_efetch)


def test_download_takes_one_reference_per_gap(monkeypatch, tmp_path, logger):
    """Sub-interval searches return up to --max-hits each; 40 candidates from one gap
    used to become 40 downloads. The best new candidate per gap is fetched: identity
    over the searched sub-interval, so a 100 % HSP over 15 % of it does not win."""
    gap1, gap2 = _gap_at(0), _gap_at(2000)
    candidates = [_cand(gap1, f"G1_{i:02d}", 85.0 + i * 0.1, 95.0) for i in range(40)]
    candidates += [
        _cand(gap2, "OLD", 99.0, 99.0, present=True),  # already in the collection
        _cand(gap2, "SHORT", 100.0, 15.0),  # a short, perfect HSP
        _cand(gap2, "G2_A", 91.0, 95.0),  # 86.5 over the sub-interval
        _cand(gap2, "G2_B", 93.0, 90.0),  # 83.7 over the sub-interval
    ]
    fetched: list[str] = []
    _recording_efetch(monkeypatch, fetched)

    downloaded = discover_run._download(candidates, tmp_path, logger)

    assert fetched == ["G1_39", "G2_A"]
    assert [c.hit.accession for c in downloaded] == ["G1_39", "G2_A"]
    assert [c.gap for c in downloaded] == [gap1, gap2]


def test_download_tries_the_next_candidate_when_a_fetch_fails(monkeypatch, tmp_path, logger):
    gap = _gap_at(0)
    fetched: list[str] = []
    _recording_efetch(monkeypatch, fetched, failing=("BEST",))

    downloaded = discover_run._download(
        [_cand(gap, "BEST", 95.0, 95.0), _cand(gap, "NEXT", 90.0, 95.0)], tmp_path, logger
    )

    assert fetched == ["BEST", "NEXT"]
    assert [c.hit.accession for c in downloaded] == ["NEXT"]


def test_a_reference_best_for_two_gaps_is_fetched_once(monkeypatch, tmp_path, logger):
    """The second gap is served by the reference the first one fetched; it is not given
    a further, worse reference of its own."""
    gap1, gap2 = _gap_at(0), _gap_at(2000)
    fetched: list[str] = []
    _recording_efetch(monkeypatch, fetched)

    discover_run._download(
        [_cand(gap1, "SHARED", 95.0, 95.0), _cand(gap2, "SHARED", 94.0, 95.0),
         _cand(gap2, "WORSE", 80.0, 95.0)],
        tmp_path, logger,
    )

    assert fetched == ["SHARED"]
```

Append to `tests/unit/test_iterate.py`:

```python


# --- audit 4: how many genomes one round may add (D1) ------------------------------

def test_a_round_adds_at_most_one_reference_per_searched_gap(monkeypatch, tmp_path, logger):
    """Two 1000 bp gaps, each searched in four sub-intervals of 400 bp returning five
    hits each: the real collect_candidates and _download, with BLAST and efetch mocked."""
    from tessera.discover import run as discover_run

    query, coll, out = _setup(tmp_path)
    _common_mocks(monkeypatch, [([_sized_gap(0, 1000, 0.8), _sized_gap(2000, 1000, 0.8)],
                                 0.94), ([], 0.95)])
    monkeypatch.setattr(iterate, "read_fasta", lambda path: [("q", "ACGT" * 1000)])
    calls = iter(range(1000))

    def fake_blast(seq, *, max_hits, logger, email=None, cache_dir=None, **kw):
        n = next(calls)
        return [Hit(f"H{n}_{k}", "t", 90.0 - k, 95.0, 1e-9) for k in range(max_hits)]

    fetched: list[str] = []

    def fake_efetch(accession, dest, logger):
        fetched.append(accession)
        (dest / f"{accession}.fasta").write_text(f">{accession}\nACGT\n")
        return dest / f"{accession}.fasta"

    monkeypatch.setattr(discover_run, "blast_subsequence", fake_blast)
    monkeypatch.setattr(discover_run, "efetch_fasta", fake_efetch)

    trace = fill_references(
        FillParams(query=query, collection=coll, output=out, top_gaps=2, max_hits=5),
        logger,
    )

    assert next(calls) == 8  # 2 gaps x 4 sub-intervals were searched
    assert len(trace[0].added) == 2
    assert fetched == trace[0].added


def test_a_download_removed_by_curation_is_excluded_and_searched_past(
    monkeypatch, tmp_path, logger
):
    """With one reference per gap, the best hit for a gap is often a sibling of the query,
    which curation then removes. The removed genome is excluded from later searches, and
    a round built from an unchanged panel is not read as a stall: the next round
    searches again and recruits the gap's next candidate."""
    from tessera.discover import run as discover_run

    query, coll, out = _setup(tmp_path)
    _common_mocks(monkeypatch, [([_gap(0.80)], 0.94), ([_gap(0.80)], 0.94)])
    excluded_per_round: list[set[str]] = []

    def collect(targets, query_row, existing, *, exclude, **kw):
        excluded_per_round.append(set(exclude))
        hits = [] if "SIB" in exclude else [Hit("SIB", "sibling", 99.0, 98.0, 0.0)]
        return [Candidate(targets[0], h, False)
                for h in hits + [Hit("DONOR", "donor", 88.0, 30.0, 1e-9)]]

    monkeypatch.setattr(iterate, "collect_candidates", collect)

    def fake_efetch(accession, dest, logger):
        (dest / f"{accession}.fasta").write_text(f">{accession}\nACGT\n")
        return dest / f"{accession}.fasta"

    monkeypatch.setattr(discover_run, "efetch_fasta", fake_efetch)
    _fake_skani(monkeypatch, {"refA": (92.0, 95.0), "SIB": (99.0, 98.0),
                              "DONOR": (88.0, 30.0)})

    trace = fill_references(
        FillParams(query=query, collection=coll, output=out, curate=True, max_rounds=2),
        logger,
    )

    assert [r.added for r in trace] == [[], ["DONOR"]]
    assert "SIB" not in excluded_per_round[0] and "SIB" in excluded_per_round[1]
    assert sorted(p.name for p in (out / "collection").iterdir()) == [
        "DONOR.fasta", "refA.fasta",
    ]
```

- [ ] **Step 2: Run them, expected failure**

```bash
pytest tests/unit/test_discover.py tests/unit/test_iterate.py -q
```

Observed (on `4b35e39`; Task 1 does not change these outcomes):
- `test_download_takes_one_reference_per_gap`: `assert ['G1_00', 'G1... 'G1_05', ...] == ['G1_39', 'G2_A']`, "Left contains 41 more items".
- `test_a_reference_best_for_two_gaps_is_fetched_once`: `assert ['SHARED', 'WORSE'] == ['SHARED']`.
- `test_a_round_adds_at_most_one_reference_per_searched_gap`: `assert 40 == 2`.
- `test_a_download_removed_by_curation_is_excluded_and_searched_past`: `assert [['DONOR'], []] == [[], ['DONOR']]` (both arrived in round 1 under the bulk download).
- `test_download_tries_the_next_candidate_when_a_fetch_fails` **passes before and after**: the bulk download also tried both. It is a guard for the new per-gap loop.

- [ ] **Step 3: Implement**

`src/tessera/discover/run.py`, replace `_download` (lines 272-292) with:

```python
def _subinterval_identity(candidate: Candidate) -> tuple[float, float]:
    """Rank key: identity over the searched sub-interval, then the HSP's own identity.

    ``pct_identity`` is measured over the aligned part only, so a 60 bp HSP at 100 %
    would outrank a hit that spans the sub-interval at 95 %. Weighting by query coverage
    counts the unaligned remainder as unmatched.
    """
    hit = candidate.hit
    return (hit.pct_identity * hit.query_coverage / 100.0, hit.pct_identity)


def _ranked_new_by_gap(candidates: list[Candidate]) -> list[list[Candidate]]:
    """The new (not yet held) candidates of each gap, best first, gaps in search order."""
    by_gap: dict[tuple[int, int], list[Candidate]] = {}
    for c in candidates:
        if not c.in_collection:
            by_gap.setdefault((c.gap.msa_start, c.gap.msa_end), []).append(c)
    # sorted() is stable, so equal keys keep BLAST's own order.
    return [sorted(group, key=_subinterval_identity, reverse=True) for group in by_gap.values()]


def _download(
    candidates: list[Candidate], dest: Path, logger: logging.Logger
) -> list[Candidate]:
    """Fetch the best new reference per gap. Returns what was saved.

    A gap searched in sub-intervals contributes up to ``--max-hits`` candidates per
    sub-interval; all of them are reported in ``candidate_references.tsv``, but only the
    best is downloaded (see :func:`_subinterval_identity`). When that fetch fails the
    next candidate is tried. A gap whose best candidate was already fetched for an
    earlier gap is served by it and gets no second reference.
    """
    ranked = _ranked_new_by_gap(candidates)
    if not ranked:
        logger.info("No new references to download (all hits already in the collection).")
        return []
    logger.info("Downloading the best new reference for each of %d gap(s) to %s",
                len(ranked), dest)
    downloaded: list[Candidate] = []
    fetched: set[str] = set()
    failed: set[str] = set()
    for options in ranked:
        for candidate in options:
            accession = candidate.hit.accession
            if accession in fetched:
                break  # this gap is served by a reference another gap fetched
            if accession in failed:
                continue
            try:
                out = efetch_fasta(accession, dest, logger)
            except Exception as exc:  # noqa: BLE001 - report and try the gap's next candidate
                logger.warning("  ! failed to download %s: %s", accession, exc)
                failed.add(accession)
                continue
            logger.info("  + %s -> %s (gap query %d-%d, %.1f%% identity)",
                        accession, out.name, candidate.gap.query_start,
                        candidate.gap.query_end, candidate.hit.pct_identity)
            fetched.add(accession)
            downloaded.append(candidate)
            break
    return downloaded
```

`src/tessera/discover/iterate.py`, `_grow_collection` state and `curate_if_pending`:

```diff
     backbone: Path | None = None
     pending_round: RoundResult | None = None  # the round whose downloads await curation
+    # Set when curation removed every reference the previous round downloaded: the next
+    # build then sees the panel the last one saw, and an unchanged shortfall says nothing
+    # about whether a further search could help.
+    all_added_curated_away = False
 
     def curate_if_pending() -> None:
-        nonlocal curate_pending, pending_round
+        nonlocal curate_pending, pending_round, all_added_curated_away
         if not curate_pending:
             return
         _curate_round(params, collection, backbone, panel_rows, logger)
+        # A genome curation removed would otherwise be found and downloaded again by
+        # every later search, in place of the gap's next candidate.
+        exclude.update(
+            _base_accession(row["genome"]) for row in panel_rows.values()
+            if row["role"].endswith("-dropped")
+        )
         if pending_round is not None:
             kept = _labels(collection)
-            pending_round.added = [a for a in pending_round.added if a in kept]
+            downloaded = pending_round.added
+            pending_round.added = [a for a in downloaded if a in kept]
+            all_added_curated_away = bool(downloaded) and not pending_round.added
         curate_pending, pending_round = False, None
```

The stall test (as left by Task 1):

```diff
         # Judged on the shortfall, not on the worst gap: only the --top-gaps longest gaps
         # are searched, so a short gap outside that set held the worst-gap similarity
         # still and stopped the loop while the searched gaps were improving.
-        if (
+        if all_added_curated_away:
+            logger.info(
+                "Curation removed every reference the last round added; searching again "
+                "with them excluded."
+            )
+            all_added_curated_away = False
+        elif (
             prev_shortfall is not None
             and (prev_shortfall - shortfall) < params.min_improvement * prev_shortfall
         ):
```

(`prev_shortfall = shortfall` after the block runs on both paths, so the next round compares against this one.)

`tests/unit/test_discover.py`, the two existing tests:

```diff
@@ test_curate_after_download_keeps_the_users_own_genomes
-    assert sorted(p.name for p in coll.iterdir()) == [
-        "NEW123.fasta", "refA.fasta", "refB.fasta", "userC.fasta",
-    ]
+    # One reference per gap: NEWSIB (97 % over 96 % of the gap) outranks NEW123 and is the
+    # only download; curation then removes it and leaves the user's genomes alone.
+    assert sorted(p.name for p in coll.iterdir()) == [
+        "refA.fasta", "refB.fasta", "userC.fasta",
+    ]
@@ test_curate_after_download_into_a_new_directory
-    assert sorted(p.name for p in fresh.iterdir()) == ["NEW123.fasta"]
+    # NEWSIB is the gap's one download, and curation is free to drop it.
+    assert sorted(p.name for p in fresh.iterdir()) == []
```

The remaining assertions of both tests (`rows["userC"] == "sibling-kept"`, `rows["NEWSIB"] == "sibling-dropped"`, the collection left untouched) stay as they are and still pass.

`docs/reference-panels.md`, the `find-references` paragraph (lines 66-70):

```diff
 When a gap is found, `tessera find-references` searches NCBI for the missing genome:
 it BLASTs the under-covered query subsequence against `nt`, reports candidate
 references (accession, identity, whether already in your collection), and -- with
 `--download` -- fetches the best new one into a collection directory to re-run with.
+Every candidate is listed in `candidate_references.tsv`; one is downloaded per gap,
+the one with the highest identity over the sub-interval it was found in (identity times
+query coverage, so a short perfect match does not outrank a hit that spans it).
```

The sub-interval paragraph (lines 85-87):

```diff
 flanks' donor and hides the tract's. (On a norovirus capsid gap this is the difference
-between recruiting only the flank reference and also recruiting the divergent core's
-distinct donor.) Resolution has two levers: `--window-size` controls whether a short
+between finding only the flank reference and also finding the divergent core's distinct
+donor. Both are listed as candidates; one reference is downloaded per gap per run, and
+`fill-references` can take another in a later round, when what remains of the gap is
+searched again.) Resolution has two levers: `--window-size` controls whether a short
```

The `find-references --curate` paragraph (lines 297-301), appended:

```diff
 have dropped them. Both need skani (and skDER for dereplication):
 `conda install -c bioconda skani skder`.
+
+Because one reference is downloaded per gap, a gap whose best candidate is a sibling of
+the query is left without a new reference once curation removes it; its other
+candidates are in `candidate_references.tsv`. Inside `fill-references --curate` (and
+`detect`, which curates) a genome curation removes is excluded from later searches,
+and the next round takes the gap's next candidate.
```

`CHANGELOG.md`, in the `Panel recruitment (fourth audit):` group:

```markdown
- The fill loop re-downloaded, every round, the genomes curation had just removed, because
  they were never added to the exclusion list. They now are.
```

and in `### Changed`:

```markdown
- **One new reference is downloaded per coverage gap**, as documented. Every sub-interval
  search contributed all of its hits, so at the defaults one round could add up to 20
  genomes for a 1000 bp gap and 145 for a 6 kb gap, all of which went into the next
  alignment. The best candidate is the one with the highest identity over the searched
  sub-interval; `candidate_references.tsv` still lists them all. This applies to
  `find-references --download`, `fill-references`, `detect` and `build-panel`. With
  `--curate`, a gap whose best candidate is a sibling of the query gets its next
  candidate in the following round rather than in the same one; `find-references`, which
  runs once, leaves such a gap without a download.
```

- [ ] **Step 4: Run, expected pass**

```bash
pytest tests/unit/test_discover.py tests/unit/test_iterate.py -q
```

The gate for this task: the four new tests above, the two changed tests, and the existing curation-loop tests in `tests/unit/test_iterate.py` (`test_last_round_downloads_are_curated_before_the_final_build`, `test_curate_never_removes_the_user_reference`, `test_curate_runs_on_a_supplied_collection_before_the_first_build`), all of which still pass unchanged.

- [ ] **Step 5: Full checks**

```bash
ruff check src tests validation && mypy && pytest -m "not requires_binary" -q
```

Expected: `766 passed, 5 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/discover/run.py src/tessera/discover/iterate.py \
        tests/unit/test_discover.py tests/unit/test_iterate.py \
        docs/reference-panels.md CHANGELOG.md
git commit -F - <<'EOF'
Download one reference per gap, and search past what curation removes

_download fetched every new accession from every sub-interval search,
although its docstring, the find-references help and the docs say one
per gap. At the defaults one round could add up to 20 genomes for a
1000 bp gap and 145 for a 6 kb gap, all of them aligned in the next
round (fourth audit, item D1; the spec's 60 and 435 assume 15 hits per
sub-interval, not the default 5).

The best candidate per gap is now the one with the highest identity over
the sub-interval it was found in (identity times query coverage), so a
short perfect HSP does not win. A failed fetch falls through to the
gap's next candidate. Every candidate is still listed in
candidate_references.tsv. find-references, fill-references, detect and
build-panel all change; find-references' docs promised one per gap.

Trade-off, and the second change: the best candidate is often a sibling
of the query, which --curate (always on in detect and build-panel) then
removes. Under the bulk download the gap's other candidates arrived in
the same round; now they would not, and because removed genomes were
never excluded, every later round would pick the sibling again and the
unchanged panel would read as a stall. Curation-removed genomes are now
added to the exclusion list, and a round whose predecessor's downloads
were all removed skips the stall test and searches again. A gap whose
best hit is a sibling therefore costs one round. find-references, which
runs once, leaves such a gap without a download; its docs say so.

Two find-references curation tests pinned the bulk download and now
expect only the sibling to be fetched; their intent is unchanged.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

- [ ] **Step 7: Paired real-data check before merge (not a gate; record the outcome)**

Needs network, Entrez Direct, `mafft` and skani. Fetch the validation data first if absent (`python validation/fetch.py`).

```bash
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
q="$PWD/validation/data/norovirus_gii/query/GII.P16-GII.1.fasta"
git worktree add --detach /tmp/planE_base fix-audit-caller-behaviour
(cd /tmp/planE_base && PYTHONPATH=src python -m tessera detect -q "$q" -o /tmp/planE_before --max-rounds 2)
python -m tessera detect -q "$q" -o /tmp/planE_after --max-rounds 2
for d in /tmp/planE_before /tmp/planE_after; do
  echo "$d"; cat "$d/fill_summary.tsv"; ls "$d/collection" | wc -l
  cut -f1-8 "$d/recombination_regions.tsv"
done
git worktree remove --force /tmp/planE_base
```

Use a worktree, not `git stash`, for the baseline: the repository's stash holds plan C Task 6 (see Global Constraints), and a second entry would shift it. Record both summaries in the PR: panel size, rounds run, and whether the reported regions and donors differ. If this cannot be run, write "paired detect check: not run (reason)" under "Skipped" in the PR. Do not describe the change as neutral on real data without it.

---

### Task 3: Saturation needs a kept hit  (audit items: D4)

`_seed_windowed` (`iterate.py:935-1003`) decided saturation from `any_hits` (hits before filtering, line 963) and then tested only the kept hits (line 1000-1002). When every hit was the query's own record or an `--exclude` accession, the kept list was empty, `not any(...)` was vacuously true, and the run reported saturation: a negative-lineage sweep and, failing that, a whole-taxon NCBI Virus download, logged with a false explanation. `_seed_whole` (line 923) already uses `bool(kept) and all(...)`. This task makes the two seeders follow one rule.

**Files:**
- Modify: `src/tessera/discover/iterate.py` (`_seed_windowed`, lines 957-964 and 997-1002)
- Modify: `CHANGELOG.md`
- Test: `tests/unit/test_iterate.py`

**Interfaces:**
- `_seed_windowed(query_seq, params, exclude, logger, *, drop_siblings) -> tuple[list[str], bool]` keeps its signature; only the boolean changes for an all-filtered input.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_iterate.py`:

```python


# --- audit 4: saturation needs a kept hit (D4) -------------------------------------

@pytest.mark.parametrize(
    ("hit", "exclude"),
    [(Hit("OWN", "the query's own record", 99.9, 99.0, 0.0), set()),  # a self-hit
     (Hit("BAD1", "excluded by the user", 97.0, 98.0, 0.0), {"BAD1"})],  # --exclude
)
def test_windowed_seeding_is_not_saturated_when_every_hit_was_filtered(
    monkeypatch, tmp_path, logger, hit, exclude
):
    """Every search returned something, and the filters removed all of it. Nothing was
    kept, so nothing says the query's lineage saturates NCBI; reporting saturation would
    start a negative-lineage sweep and a whole-taxon download on a false explanation."""
    monkeypatch.setattr(iterate, "blast_subsequence", lambda seq, **k: [hit])
    params = FillParams(query=tmp_path / "q.fasta", collection=None, output=tmp_path / "o",
                        seed_window=100)

    windowed = iterate._seed_windowed("A" * 300, params, exclude, logger, drop_siblings=True)
    whole = iterate._seed_whole("A" * 300, params, exclude, logger)

    assert windowed == ([], False)
    assert whole == ([], False)  # the whole-query seeder already got this right


def test_windowed_seeding_is_saturated_when_every_kept_hit_is_a_sibling(
    monkeypatch, tmp_path, logger
):
    """Guard (passes before and after): the true saturation case still reports it."""
    monkeypatch.setattr(iterate, "blast_subsequence",
                        lambda seq, **k: [Hit("SIB", "sibling", 98.0, 99.0, 0.0)])
    params = FillParams(query=tmp_path / "q.fasta", collection=None, output=tmp_path / "o",
                        seed_window=100)

    _, saturated = iterate._seed_windowed("A" * 300, params, set(), logger, drop_siblings=True)

    assert saturated is True
```

- [ ] **Step 2: Run them, expected failure**

```bash
pytest tests/unit/test_iterate.py -k saturated -q
```

Observed: both parameter cases of `test_windowed_seeding_is_not_saturated_when_every_hit_was_filtered` fail with `assert ([], True) == ([], False)`. `test_windowed_seeding_is_saturated_when_every_kept_hit_is_a_sibling` passes before and after (guard).

- [ ] **Step 3: Implement**

`src/tessera/discover/iterate.py`, `_seed_windowed`:

```diff
-    per_window_hits: list[list] = []
-    any_hits = False  # did any search return anything, before filtering?
+    per_window_hits: list[list] = []  # each window's hits after the self/--exclude filter
     for chunk in chunks:
         if len(chunk) < MIN_SUBSEQ:
             continue
         hits = _blast_or_none(chunk, params, logger)
-        any_hits = any_hits or bool(hits)
         per_window_hits.append([h for h in hits if _keep(h, exclude, params)])
```

```diff
-    # "Only siblings were found" needs something to have been found. When every
-    # search failed or came back empty there is no evidence about what NCBI holds, and
-    # calling that saturation would switch the seed source under a false explanation.
-    saturated = any_hits and not any(
-        not _is_sibling_hit(h) for kept in per_window_hits for h in kept
-    )
+    # "Only siblings were found" needs something to have been kept. When every search
+    # failed, came back empty, or returned only the query's own record or --exclude
+    # accessions, there is no evidence about what NCBI holds, and calling that saturation
+    # would switch the seed source under a false explanation. Same rule as _seed_whole.
+    kept_hits = [h for kept in per_window_hits for h in kept]
+    saturated = bool(kept_hits) and all(_is_sibling_hit(h) for h in kept_hits)
     return selected[:SEED_TOTAL_CAP], saturated
```

`CHANGELOG.md`, `Panel recruitment (fourth audit):` group:

```markdown
- Windowed BLAST seeding reported that the query's lineage saturates NCBI when every hit
  had been filtered out (the query's own record, or `--exclude` accessions), and then ran
  the negative-lineage search and the NCBI Virus download on that false explanation. It
  now needs a kept hit, as whole-query seeding already did.
```

- [ ] **Step 4: Run, expected pass**

```bash
pytest tests/unit/test_iterate.py -q
```

Gate: the three tests above plus the existing seeding tests (`test_parents_mode_falls_back_when_only_siblings`, `test_saturation_auto_switches_to_ncbi_virus_diversity`, `test_parents_mode_suppresses_siblings_and_seeds_parents`), unchanged.

- [ ] **Step 5: Full checks**

```bash
ruff check src tests validation && mypy && pytest -m "not requires_binary" -q
```

Expected: `769 passed, 5 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/discover/iterate.py tests/unit/test_iterate.py CHANGELOG.md
git commit -F - <<'EOF'
Report seeding saturation only when a kept hit says so

_seed_windowed computed "saturated" from whether any search returned
anything before filtering, then tested only the kept hits. When every
hit was the query's own record or an --exclude accession, the kept list
was empty and the all-siblings test held vacuously, so the run started
a negative-lineage sweep and a whole-taxon NCBI Virus download, logged
as saturation (fourth audit, item D4). _seed_whole already required a
kept hit; both seeders now use that rule.

No trade-off: an all-filtered window still falls back as before, and the
true all-sibling case still reports saturation (guard test).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 4: Zero RefSeq genomes broaden instead of aborting  (audit items: D5)

`_fetch_diverse` requests RefSeq first and broadens to complete genomes below three (`iterate.py:775-788`), but `fetch_ncbi_virus` raised when the package held no genomes at all (`pool.py:170-174`), so the broadening was unreachable for exactly zero. A taxon with no RefSeq genome and 50 complete ones aborted; with two RefSeq genomes it broadened and returned all 50.

`fetch_ncbi_virus` has one caller (`_fetch_diverse`), so it now returns `[]` for an empty package and the caller decides; `_fetch_diverse` raises the same message when the final set, after broadening, is empty. The raise happens inside `atomic_cache_dir`, which discards the staging tree on any exception (`core/cache.py:158-162`), so no empty panel is cached.

**Files:**
- Modify: `src/tessera/discover/pool.py` (`fetch_ncbi_virus` docstring lines 142-150, empty-package branch lines 170-174)
- Modify: `src/tessera/discover/iterate.py` (`_fetch_diverse`, before line 793)
- Modify: `CHANGELOG.md`
- Test: `tests/unit/test_iterate.py`

**Interfaces:**
- `fetch_ncbi_virus(...) -> list[Path]` now returns `[]` when the package has no `genomic.fna`, where it raised `UserInputError`.
- `_fetch_diverse(params, logger) -> list[Path]` raises `UserInputError("NCBI Virus returned no genomes for '<taxon>'. ...")` when nothing is left after broadening.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_iterate.py`:

```python


# --- audit 4: zero RefSeq genomes broaden (D5) -------------------------------------

def _fake_datasets(monkeypatch, *, refseq: int, complete: int) -> list[list[str]]:
    """Stand in for the datasets CLI: write a package holding ``refseq`` or ``complete``
    genomes. A package with none carries no genomic.fna, as datasets writes it."""
    import zipfile

    from tessera.discover import pool as pool_mod

    commands: list[list[str]] = []

    def fake_run_tool(tool, command, **kwargs):
        commands.append(command)
        n = refseq if "--refseq" in command else complete
        with zipfile.ZipFile(command[command.index("--filename") + 1], "w") as zf:
            zf.writestr("ncbi_dataset/data/data_report.jsonl", "")
            if n:
                zf.writestr("ncbi_dataset/data/genomic.fna",
                            "".join(f">G{i}.1 some virus\nACGT\n" for i in range(n)))

    monkeypatch.setattr(pool_mod, "datasets_available", lambda: True)
    monkeypatch.setattr(pool_mod, "run_tool", fake_run_tool)
    monkeypatch.setattr(pool_mod.eutils_throttle, "wait", lambda logger: None)
    return commands


def test_a_taxon_with_no_refseq_genome_broadens_to_complete_genomes(
    monkeypatch, tmp_path, logger
):
    commands = _fake_datasets(monkeypatch, refseq=0, complete=50)
    out = tmp_path / "o"
    out.mkdir()
    params = FillParams(query=tmp_path / "q.fasta", collection=None, output=out,
                        taxon="Some virus", cache_dir=tmp_path / "cache")

    genomes = iterate._fetch_diverse(params, logger)

    assert len(genomes) == 50
    assert ["--refseq" in c for c in commands] == [True, False]


def test_a_taxon_with_no_genome_at_all_is_a_clear_error_and_caches_nothing(
    monkeypatch, tmp_path, logger
):
    _fake_datasets(monkeypatch, refseq=0, complete=0)
    out = tmp_path / "o"
    out.mkdir()
    params = FillParams(query=tmp_path / "q.fasta", collection=None, output=out,
                        taxon="Some virus", cache_dir=tmp_path / "cache")

    with pytest.raises(UserInputError, match="NCBI Virus returned no genomes"):
        iterate._fetch_diverse(params, logger)

    cache_root = tmp_path / "cache" / "ncbi_virus"
    assert not cache_root.exists() or not any(p.is_dir() for p in cache_root.iterdir())
```

- [ ] **Step 2: Run them, expected failure**

```bash
pytest tests/unit/test_iterate.py -k taxon -q
```

Observed: `test_a_taxon_with_no_refseq_genome_broadens_to_complete_genomes` fails with `tessera.core.errors.UserInputError: NCBI Virus returned no genomes for 'Some virus'. Try --taxon or a broader scope.` raised from the RefSeq request. `test_a_taxon_with_no_genome_at_all_is_a_clear_error_and_caches_nothing` **passes before and after** (guard: the error message and the empty cache must survive the move of the raise).

- [ ] **Step 3: Implement**

`src/tessera/discover/pool.py`, `fetch_ncbi_virus`:

```diff
     downstream split/align work; a heavily sequenced taxon still downloads in full,
     so a curated ``--candidate-pool`` is preferable there.
+
+    Returns an empty list when the package holds no genomes for this scope. Whether that
+    is an error is the caller's decision: an empty RefSeq set is the case for broadening
+    to complete genomes, not for stopping the run.
     """
```

```diff
         fna = data_dir / "genomic.fna"
         if not fna.exists():
-            raise UserInputError(
-                f"NCBI Virus returned no genomes for '{taxon}'. Try --taxon or a broader scope."
-            )
+            logger.info("NCBI Virus holds no %s genomes for '%s'.", scope, taxon)
+            return []
```

`src/tessera/discover/iterate.py`, `_fetch_diverse`, inside the `with atomic_cache_dir(cache) as staging:` block, immediately before `if len(fetched) > params.fetch_limit:`:

```diff
+        if not fetched:
+            # Raised inside the staging context, so no empty panel is cached.
+            raise UserInputError(
+                f"NCBI Virus returned no genomes for '{taxon}'. Try --taxon or a broader "
+                "scope."
+            )
         if len(fetched) > params.fetch_limit:
```

`CHANGELOG.md`, `Panel recruitment (fourth audit):` group:

```markdown
- A taxon with no RefSeq genome stopped the run with "NCBI Virus returned no genomes"
  instead of broadening to complete genomes, which it did for a taxon with one or two
  RefSeq genomes. It now broadens; the error remains for a taxon with no genome at all.
```

- [ ] **Step 4: Run, expected pass**

```bash
pytest tests/unit/test_iterate.py tests/unit/test_pool.py -q
```

Gate: the two tests above plus `test_fetch_diverse_broadens_downloads_all_and_caches`, `test_interrupted_fetch_leaves_no_reusable_cache` and `test_fetch_scope_does_not_share_a_cache_slot`, unchanged.

- [ ] **Step 5: Full checks**

```bash
ruff check src tests validation && mypy && pytest -m "not requires_binary" -q
```

Expected: `771 passed, 5 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/discover/pool.py src/tessera/discover/iterate.py \
        tests/unit/test_iterate.py CHANGELOG.md
git commit -F - <<'EOF'
Broaden to complete genomes when a taxon has no RefSeq genome

_fetch_diverse broadens from RefSeq to complete genomes below three
RefSeq genomes, but fetch_ncbi_virus raised when the package held none,
so the broadening was unreachable for exactly zero: such a taxon
aborted, while one with two RefSeq genomes broadened (fourth audit,
item D5).

fetch_ncbi_virus, whose only caller is _fetch_diverse, now returns an
empty list for an empty package, and _fetch_diverse raises the same
message when nothing is left after broadening. The raise is inside the
atomic staging context, so no empty panel is cached (guard test).

Not verified: whether the datasets CLI writes an empty package or exits
non-zero for a scope with no genomes. The test follows the audit's
reproduction (no genomic.fna). If datasets exits non-zero, run_tool
raises ToolExecutionError and the run still stops; that path is unchanged.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 5: The default typing path reads the lineage sidecar  (audit items: D3; test gap "the default and deep typing branches compared on one panel")

Recruitment writes each NCBI Virus genome's exact lineage to `<output>/ncbi_lineages.tsv` (`_record_ncbi_lineages`, `iterate.py:485-498`). Only the deep branch of `_type_panel` read it (`iterate.py:526-535`); the default branch (`iterate.py:546-551`) re-mined headers with `genotype_from_title`, which needs a digit-bearing token of four characters or more, so `B.1`, `A` and `XBB` were lost and those genomes skipped lineage-aware handling. `build_lineage_map` already takes `datasets_rows` with the right priority (user map > datasets > mined title, `recomb/typing.py:167-192`); the fix passes it.

**Files:**
- Modify: `src/tessera/discover/iterate.py` (`_type_panel`, lines 525-551)
- Modify: `docs/reference-panels.md` (the last sentence of "Type the recruited panel automatically", lines 321-323), `CHANGELOG.md`
- Test: `tests/unit/test_iterate.py`

**Interfaces:**
- Consumes: `_read_ncbi_lineages(output) -> list[tuple[str, str]] | None`, `build_lineage_map(*, user_tsv, datasets_rows, title_by_label, organism)`.
- `_type_panel` keeps its signature; `lineages.tsv` rows from the sidecar carry source `ncbi-datasets` on both branches.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_iterate.py`:

```python


# --- audit 4: the default typing path reads the lineage sidecar (D3) ---------------

def test_default_and_deep_typing_agree_on_a_sidecar_typed_panel(tmp_path, logger):
    """Recruitment records each NCBI Virus genome's exact lineage in
    <output>/ncbi_lineages.tsv. Only the deep branch read it; the default branch mined
    the headers again with a rule that needs a digit-bearing token of four characters or
    more, so B.1, A and XBB were lost and those genomes skipped lineage-aware handling."""
    coll = tmp_path / "collection"
    coll.mkdir()
    lineages = {"OQ1.1": "B.1", "OQ2.1": "A", "OQ3.1": "XBB"}
    for acc, lineage in lineages.items():
        (coll / f"{acc}.fasta").write_text(f">{acc} {lineage}\nACGT\n")
    out = tmp_path / "out"
    out.mkdir()
    (out / "ncbi_lineages.tsv").write_text(
        "".join(f"{acc}\t{lin}\n" for acc, lin in lineages.items())
        + "NOTINPANEL.1\tBA.2\n"  # fetched, not selected: must not reach lineages.tsv
    )
    query = tmp_path / "q.fasta"
    query.write_text(">q\nACGT\n")

    def typed(deep: bool) -> dict[str, str]:
        params = FillParams(query=query, collection=None, output=out, deep_typing=deep)
        lineage_map, _ = iterate._type_panel(params, coll, "q", logger)
        return lineage_map

    default = typed(False)
    assert default == lineages
    assert typed(True) == default
    assert "NOTINPANEL" not in (out / "lineages.tsv").read_text()


def test_a_user_lineage_map_still_overrides_the_sidecar_on_the_default_path(tmp_path, logger):
    """Priority stays user map > NCBI datasets lineage > mined header."""
    coll = tmp_path / "collection"
    coll.mkdir()
    (coll / "OQ1.1.fasta").write_text(">OQ1.1 B.1\nACGT\n")
    out = tmp_path / "out"
    out.mkdir()
    (out / "ncbi_lineages.tsv").write_text("OQ1.1\tB.1\n")
    user_map = tmp_path / "mine.tsv"
    user_map.write_text("OQ1.1\tB.1.617.2\n")
    query = tmp_path / "q.fasta"
    query.write_text(">q\nACGT\n")

    params = FillParams(query=query, collection=None, output=out, lineage_map=user_map)
    lineage_map, _ = iterate._type_panel(params, coll, "q", logger)

    assert lineage_map == {"OQ1.1": "B.1.617.2"}
```

The deep branch in the first test reaches no network: every genome is typed by the sidecar, so `assign_lineages` has no untyped genome to send to Nextclade or skani.

- [ ] **Step 2: Run them, expected failure**

```bash
pytest tests/unit/test_iterate.py -k "sidecar" -q
```

Observed: `test_default_and_deep_typing_agree_on_a_sidecar_typed_panel` fails with `assert {} == {'OQ1.1': 'B.1', 'OQ2.1': 'A', 'OQ3.1': 'XBB'}`. `test_a_user_lineage_map_still_overrides_the_sidecar_on_the_default_path` **passes before and after** (guard for Review Focus item 5: before the change the sidecar was not read at all; after it, the user map must still win).

- [ ] **Step 3: Implement**

`src/tessera/discover/iterate.py`, `_type_panel`:

```diff
     coll_files = collection_genomes(collection)
+    # The consolidated sidecar records the full fetched NCBI-Virus set; restrict the
+    # datasets rows to accessions that survived selection so lineages.tsv describes the
+    # panel, not the download. Both branches read it: mining the headers again loses
+    # every lineage name without a digit-bearing token of four characters or more.
+    panel_accessions = {strip_sequence_extension(p.name) for p in coll_files}
+    recorded = _read_ncbi_lineages(params.output)
+    datasets_rows = (
+        [(acc, lin) for acc, lin in recorded if acc in panel_accessions]
+        if recorded else None
+    )
     if params.deep_typing:
-        # The consolidated sidecar records the full fetched NCBI-Virus set; restrict the
-        # datasets rows to accessions that survived selection so lineages.tsv describes the
-        # panel, not the download.
-        panel_accessions = {strip_sequence_extension(p.name) for p in coll_files}
-        recorded = _read_ncbi_lineages(params.output)
-        datasets_rows = (
-            [(acc, lin) for acc, lin in recorded if acc in panel_accessions]
-            if recorded else None
-        )
         lineage_rows = assign_lineages(
```

```diff
     else:
         lineage_rows = build_lineage_map(
             user_tsv=params.lineage_map,
+            datasets_rows=datasets_rows,
             title_by_label=titles_from_collection(coll_files),
             organism=params.taxon,
         )
```

`docs/reference-panels.md`, end of "Type the recruited panel automatically (`--deep-typing`)":

```diff
-the end of the run; it needs skani (the command errors early if skani is absent). References
-fetched from NCBI Virus are additionally tagged with the source `ncbi-datasets` -- their
-structured datasets lineage, captured during download in an `ncbi_lineages.tsv` sidecar,
-rather than a token mined back out of the header.
+the end of the run; it needs skani (the command errors early if skani is absent). With or
+without `--deep-typing`, references fetched from NCBI Virus are typed from their
+structured datasets lineage, captured during download in an `ncbi_lineages.tsv` sidecar,
+and tagged with the source `ncbi-datasets`, rather than from a token mined back out of the
+header (which loses short names such as `B.1` or `XBB`). A `--lineage-map` still takes
+precedence.
```

`CHANGELOG.md`, `Panel recruitment (fourth audit):` group:

```markdown
- Without `--deep-typing`, the panel was typed by mining headers again, ignoring the exact
  NCBI Virus lineages recorded during recruitment in `ncbi_lineages.tsv`; names without a
  digit-bearing token of four or more characters (`B.1`, `A`, `XBB`) were lost and those
  genomes skipped lineage-aware handling. Both typing paths now read the sidecar.
```

- [ ] **Step 4: Run, expected pass**

```bash
pytest tests/unit/test_iterate.py tests/unit/test_deep_typing.py tests/unit/test_ncbi_lineage_sidecar.py tests/unit/test_ncbi_lineage_consolidate.py -q
```

Gate: the two tests above and the existing typing tests named in the command.

- [ ] **Step 5: Full checks**

```bash
ruff check src tests validation && mypy && pytest -m "not requires_binary" -q
```

Expected: `773 passed, 5 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/discover/iterate.py tests/unit/test_iterate.py \
        docs/reference-panels.md CHANGELOG.md
git commit -F - <<'EOF'
Read the NCBI lineage sidecar on the default typing path

Recruitment records each NCBI Virus genome's exact lineage in
<output>/ncbi_lineages.tsv, but only --deep-typing read it. The default
path re-mined the headers with a rule that needs a digit-bearing token
of four characters or more, so B.1, A and XBB were lost; the audit
measured five of seven realistic headers untyped (fourth audit, item
D3). Untyped genomes skip lineage-aware selection and the recombinant-
lineage exclusion.

The sidecar rows, restricted to the panel as on the deep path, now go to
build_lineage_map on both paths. Its priority is unchanged: a user
--lineage-map still wins (guard test). A new test types one panel on
both paths and requires the same result.

Overlap: plan C Task 6 (held in a stash) edits the same call to add pool
labels; the two are independent sources and combine as two arguments.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 6: Say when a requested curation never ran  (audit items: D6)

`curate_pending` starts false for a freshly seeded collection by design (`iterate.py:389`; seeding applies its own sibling filter, and `test_curate_does_not_run_on_a_freshly_seeded_collection` pins that), and is set only after a round downloads (`iterate.py:461`). A seeded panel with no coverage gap breaks out of round 1 (`iterate.py:431-433`), so curation never runs. The report gates its panel section on the curation rows (`iterate.py:575-578`), so the run then wrote no panel table, no `panel_lineages.tsv` and no log line, which reads as "curated, nothing to drop". `detect` and `build-panel` always request curation, so this is their common path for a sibling-only seed.

Whether curation *should* run on a seeded panel before round 1 is a behaviour decision that the existing test records as deliberate; this task only reports. It is listed under Self-review notes for the maintainer.

**Files:**
- Modify: `src/tessera/discover/iterate.py` (new constant before `_report_sections`, line 564; `_report_sections` lines 575-578)
- Modify: `docs/reference-panels.md` (curation paragraph, lines 283-288), `CHANGELOG.md`
- Test: `tests/unit/test_iterate.py`

**Interfaces:**
- Produces: `CURATION_NOT_RUN: str` in `tessera.discover.iterate` (the log and report text).
- `_report_sections(...)` keeps its signature; with `params.curate` and no curation rows it logs a warning and adds a "Reference panel" section holding that text.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_iterate.py`:

```python


# --- audit 4: a requested curation that never ran is reported (D6) -----------------

def test_a_requested_curation_that_never_ran_is_reported(monkeypatch, tmp_path, caplog):
    """A freshly seeded panel is not curated before round 1, and a panel with no
    coverage gap downloads nothing, so curation never runs. The report used to drop its
    panel section and the log said nothing, which reads as "curated, nothing to drop"."""
    import logging

    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "ACGT" * 100 + "\n")
    out = tmp_path / "out"
    _common_mocks(monkeypatch, [([], 0.95)])
    monkeypatch.setattr(iterate, "skani_available", lambda: True)

    def seed(params, collection, query_records, exclude, logger):
        for name in ("S1", "S2"):
            (collection / f"{name}.fasta").write_text(f">{name}\nACGT\n")

    monkeypatch.setattr(iterate, "_seed_collection", seed)
    sections: dict = {}
    monkeypatch.setattr(iterate, "run_recomb", lambda p, logger, **kw: sections.update(kw))
    log = logging.getLogger("fill_loop_test")
    with caplog.at_level(logging.WARNING, logger="fill_loop_test"):
        fill_references(FillParams(query=query, collection=None, output=out, curate=True), log)

    assert "Curation (--curate) was requested but did not run" in caplog.text
    body = dict(sections["extra_sections"])["Reference panel"]
    assert "was requested but did not run" in body
    assert not (out / "panel_lineages.tsv").exists()


def test_no_curation_note_when_curation_was_not_requested(monkeypatch, tmp_path, logger):
    """Guard (passes before and after): without --curate there is no panel section."""
    query, coll, out = _setup(tmp_path)
    _common_mocks(monkeypatch, [([], 0.95)])
    sections: dict = {}
    monkeypatch.setattr(iterate, "run_recomb", lambda p, logger, **kw: sections.update(kw))

    fill_references(FillParams(query=query, collection=coll, output=out), logger)

    assert "Reference panel" not in dict(sections["extra_sections"])
```

(The logger is not under `tessera`, because the CLI tests switch that logger's propagation off and `caplog` sees only records reaching the root, as in `test_curate_warns_when_it_leaves_too_few_references`.)

- [ ] **Step 2: Run them, expected failure**

```bash
pytest tests/unit/test_iterate.py -k "curation_that_never or no_curation_note" -q
```

Observed: `test_a_requested_curation_that_never_ran_is_reported` fails with `assert 'Curation (--curate) was requested but did not run' in ''`. `test_no_curation_note_when_curation_was_not_requested` passes before and after (guard).

- [ ] **Step 3: Implement**

`src/tessera/discover/iterate.py`, immediately before `def _report_sections(`:

```python
CURATION_NOT_RUN = (
    "Curation (--curate) was requested but did not run on this panel. A freshly seeded "
    "collection is not curated before the first build (seeding applies its own sibling "
    "filter), and curation then runs only on references a round downloads; this run "
    "downloaded none that reached it, or found no whole-genome backbone to anchor the "
    "sibling test. Close relatives of the query that seeding kept may still be in the "
    "panel, and no panel_lineages.tsv was written."
)
```

`_report_sections`:

```diff
     if params.curate and panel_rows:
         table = list(panel_rows.values())
         write_panel_tsv(params.output / "panel_lineages.tsv", table, logger, lineage_map)
         extra_sections.append(("Reference panel", panel_table_html(table, lineage_map)))
+    elif params.curate:
+        # A curation that did not run must not read as one that ran and kept everything.
+        logger.warning(CURATION_NOT_RUN)
+        extra_sections.append(
+            ("Reference panel", f'<p class="cap">{html.escape(CURATION_NOT_RUN)}</p>')
+        )
```

`_report_sections` runs whether or not detection follows (`fill_references` calls it before checking `params.report`), so `build-panel`, which sets `report=False`, also logs the warning.

`docs/reference-panels.md`, curation paragraph (lines 283-288):

```diff
 supplied (before round 1), and each round's downloads (before the next build). A panel
 seeded from scratch is not curated before round 1, because seeding applies its own sibling
-filter. With `--reference`, that genome is never removed by curation; the sibling test
-stays anchored on the query's closest whole-genome relative, which the log names.
+filter; if no round then downloads anything (a seeded panel with no coverage gap),
+curation does not run at all. The log and the report's "Reference panel" section then say
+so, and no `panel_lineages.tsv` is written. With `--reference`, that genome is never
+removed by curation; the sibling test stays anchored on the query's closest whole-genome
+relative, which the log names.
```

`CHANGELOG.md`, `Panel recruitment (fourth audit):` group:

```markdown
- When `--curate` was requested (always, for `detect` and `build-panel`) but never ran --
  a seeded panel with no coverage gap downloads nothing for it to curate -- the run wrote
  no panel table and no log line. The log and the report's "Reference panel" section now
  state that curation did not run.
```

- [ ] **Step 4: Run, expected pass**

```bash
pytest tests/unit/test_iterate.py -q
```

Gate: the two tests above, plus `test_curate_does_not_run_on_a_freshly_seeded_collection` and `test_curate_runs_on_a_supplied_collection_before_the_first_build`, unchanged.

- [ ] **Step 5: Full checks**

```bash
ruff check src tests validation && mypy && pytest -m "not requires_binary" -q
```

Expected: `775 passed, 5 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/discover/iterate.py tests/unit/test_iterate.py \
        docs/reference-panels.md CHANGELOG.md
git commit -F - <<'EOF'
Report a requested curation that never ran

A freshly seeded collection is deliberately not curated before round 1,
and curation is owed only after a round downloads. A seeded panel with
no coverage gap therefore never reaches curation. The report gated its
panel section on the curation rows, so the run wrote no table, no
panel_lineages.tsv and no log line, which reads as a curation that ran
and kept everything (fourth audit, item D6). detect and build-panel
always request curation.

The log now carries a warning and the report a "Reference panel"
section saying curation did not run and why. Whether a seeded panel
should be curated before round 1 is left unchanged: an existing test
records it as deliberate, and changing it alters what detect recruits.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 7: Warn when every pool genome is a sibling  (audit items: D7)

`select_regional` drops the query's siblings only when some genome survives (`pool.py:328`): `siblings and len(siblings) < len(pool)`. When the whole pool is siblings the drop is skipped, which is right (an empty pool selects nothing), but it was silent, and the panel looked like any other. The seeding path logs "Every window's hits were siblings of the query; seeding the best hits instead (no distinct parental lineage found in NCBI)" (`iterate.py:985-988`); this task logs the pool counterpart, as a warning because it can mask a recombination.

**Files:**
- Modify: `src/tessera/discover/pool.py` (`select_regional`, lines 323-333)
- Modify: `docs/reference-panels.md` (sibling-drop paragraph, lines 236-240), `CHANGELOG.md`
- Test: `tests/unit/test_pool.py`

**Interfaces:**
- `select_regional(...) -> PoolSelection` unchanged in signature and output; it logs one warning when every genome in the reduced pool is a sibling and `drop_siblings` is on.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_pool.py`:

```python


def _patch_selection(monkeypatch, whole, window):
    monkeypatch.setattr(pool, "skani_available", lambda: True)
    monkeypatch.setattr(pool, "dereplicate", lambda g, ani, logger: (list(g), []))
    monkeypatch.setattr(pool, "skani_query_ani", _fake_ani(whole, window))


# --- audit 4: an all-sibling pool is reported (D7) ---------------------------------

def test_select_regional_warns_when_every_genome_is_a_sibling(monkeypatch, tmp_path, caplog):
    """Four genomes at 99 % ANI over 99 % coverage: the sibling drop would empty the pool,
    so it is skipped. That used to happen in silence and yield a normal-looking panel."""
    import logging

    names = ["S1", "S2", "S3", "S4"]
    genomes = _genomes(tmp_path, names)
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "A" * 3000 + "\n")
    _patch_selection(
        monkeypatch,
        {n: (99.0, 99.0) for n in names},
        {n: {"win_0": (99.0, 80.0), "win_1": (99.0, 80.0)} for n in names},
    )
    log = logging.getLogger("pool_test")
    with caplog.at_level(logging.WARNING, logger="pool_test"):
        sel = pool.select_regional(query, genomes, window=1500, per_window=2, logger=log)

    assert sel.selected  # still a panel: the siblings are all there is
    assert "Every pool genome is a sibling of the query" in caplog.text
    assert "no distinct parental lineage found in the pool" in caplog.text
```

(`_patch_selection` is also used by Task 8.)

- [ ] **Step 2: Run it, expected failure**

```bash
pytest tests/unit/test_pool.py -k every_genome -q
```

Observed: `assert 'Every pool genome is a sibling of the query' in ''`.

- [ ] **Step 3: Implement**

`src/tessera/discover/pool.py`, `select_regional`:

```diff
         if siblings and len(siblings) < len(pool):
             logger.info(
                 "Dropping %d sibling(s) of the query (>= %.0f%% ANI over near-full coverage).",
                 len(siblings), sibling_identity,
             )
             pool = [g for g in pool if g not in siblings]
+        elif siblings:
+            # Dropping them would empty the pool, so they stay; say so, as BLAST seeding
+            # does, rather than return a panel that looks like any other.
+            logger.warning(
+                "Every pool genome is a sibling of the query (>= %.0f%% ANI over near-full "
+                "coverage); selecting among them instead (no distinct parental lineage "
+                "found in the pool). Siblings match the query everywhere and can mask a "
+                "recombination.",
+                sibling_identity,
+            )
```

`docs/reference-panels.md`, sibling-drop paragraph:

```diff
 The sibling drop uses an absolute identity cutoff, which suits divergent-lineage
 recombinants (HIV subtypes ~12 % apart). For closely related parents that have no
 masking sibling (e.g. SARS-CoV-2 sublineages), add `--seed-keep-siblings`. NCBI Virus
 / `datasets` needs `conda install -c conda-forge ncbi-datasets-cli`.
+When every genome in the pool is a sibling, none is dropped (the pool would be empty),
+and the run logs a warning that no distinct parental lineage was found in the pool.
```

`CHANGELOG.md`, `Panel recruitment (fourth audit):` group:

```markdown
- Pool selection kept every genome without comment when all of them were siblings of the
  query, and returned a normal-looking panel. It now logs a warning, as BLAST seeding
  already did for the same situation.
```

- [ ] **Step 4: Run, expected pass**

```bash
pytest tests/unit/test_pool.py tests/unit/test_select_from_lineage.py tests/unit/test_lineage_select.py -q
```

Gate: the test above and the existing `test_select_regional_drops_siblings_and_picks_regional_parents`, `test_select_regional_keeps_siblings_when_disabled`, unchanged. The hybrid harness calls `select_regional` with `drop_siblings=False` and does not reach this branch.

- [ ] **Step 5: Full checks**

```bash
ruff check src tests validation && mypy && pytest -m "not requires_binary" -q
```

Expected: `776 passed, 5 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/discover/pool.py tests/unit/test_pool.py docs/reference-panels.md CHANGELOG.md
git commit -F - <<'EOF'
Warn when every genome in a candidate pool is a sibling of the query

select_regional skips the sibling drop when it would empty the pool,
which is right, but it did so silently, and the panel it returned looked
like any other (fourth audit, item D7: four genomes at 99 % ANI over
99 % coverage, none dropped, no message). BLAST seeding logs the same
situation. Pool selection now logs a warning in matching words; the
selection itself is unchanged.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 8: Match a short query region by region  (audit items: D9; test gap "the window count for a short query")

`select_regional` sets `n_windows = max(1, len(query_seq) // window)` (`pool.py:338`) and gives the last window the remainder. A query between one and two windows long therefore gets a single window over the whole query: a 2.9 kb genome at the default 1500 bp window is selected genome-wide, so a genome that matches only one half cannot be recruited for it. At the default window this affects queries of 1.5 to 3 kb, such as hepatitis delta (about 1.7 kb); hepatitis B (about 3.2 kb) already gets two windows.

**Rule chosen.** Only queries longer than one window and shorter than two change: they are cut into two equal halves. Longer queries keep their window grid exactly, so nothing changes for any genome of at least two windows (including every hybrid-harness case, whose selection window is at most a quarter of the genome). Rounding the window count instead would also move the last window of long genomes; that was not taken.

**Fallback.** skani is a genome-scale tool and has not been verified here on 750-1450 bp chunks. If neither half selects any genome, the whole-query ranking already computed (`whole`) supplies the panel, which is what the single window gave before (Review Focus item 4). The pre-change failure of that test arises because the mock also reports nothing for the single whole-query window; the test pins the new fallback, not a defect.

**Files:**
- Modify: `src/tessera/discover/pool.py` (`select_regional`, lines 337-338 and after line 363)
- Modify: `docs/reference-panels.md` (`local` bullet, lines 193-198), `CHANGELOG.md`
- Test: `tests/unit/test_pool.py`

**Interfaces:**
- `select_regional(...) -> PoolSelection` unchanged in signature. For `window < len(query) < 2 * window`, `PoolSelection.table` rows carry two windows, `(0, ceil(L/2))` and `(ceil(L/2), L)`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_pool.py`:

```python


# --- audit 4: a short query is still matched region by region (D9) -----------------

def test_select_regional_cuts_a_short_query_into_two_windows(monkeypatch, tmp_path, logger):
    """A 2.9 kb query at the default 1500 bp window was one window over the whole query,
    so a genome that matches only one half could not be recruited for that half."""
    par_a, par_b = _genomes(tmp_path, ["PAR_A", "PAR_B"])
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "A" * 2900 + "\n")
    _patch_selection(
        monkeypatch,
        {"PAR_A": (90.0, 85.0), "PAR_B": (88.0, 84.0)},
        {"PAR_A": {"win_0": (96.0, 80.0), "win_1": (85.0, 80.0)},
         "PAR_B": {"win_0": (85.0, 80.0), "win_1": (96.0, 80.0)}},
    )

    sel = pool.select_regional(query, [par_a, par_b], window=1500, per_window=1, logger=logger)

    assert {p.stem for p in sel.selected} == {"PAR_A", "PAR_B"}
    assert [(r["query_start"], r["query_end"]) for r in sel.table] == [(0, 1450), (1450, 2900)]


def test_select_regional_keeps_one_window_for_a_query_no_longer_than_it(
    monkeypatch, tmp_path, logger
):
    """Guard (passes before and after): nothing to split at or below one window."""
    par_a, par_b = _genomes(tmp_path, ["PAR_A", "PAR_B"])
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "A" * 1500 + "\n")
    _patch_selection(
        monkeypatch,
        {"PAR_A": (90.0, 85.0), "PAR_B": (88.0, 84.0)},
        {"PAR_A": {"win_0": (96.0, 80.0)}, "PAR_B": {"win_0": (85.0, 80.0)}},
    )

    sel = pool.select_regional(query, [par_a, par_b], window=1500, per_window=1, logger=logger)

    assert [(r["query_start"], r["query_end"]) for r in sel.table] == [(0, 1500)]


def test_select_regional_falls_back_to_the_whole_query_when_no_half_matches(
    monkeypatch, tmp_path, logger
):
    """skani may report nothing for a sub-window chunk (it is a genome-scale tool). The
    split must not then leave a short query with an empty panel where the single
    whole-query window used to select one: the whole-genome ranking is used instead."""
    par_a, par_b = _genomes(tmp_path, ["PAR_A", "PAR_B"])
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "A" * 2900 + "\n")
    _patch_selection(
        monkeypatch,
        {"PAR_A": (90.0, 85.0), "PAR_B": (88.0, 84.0)},
        {"PAR_A": {"win_0": (0.0, 0.0), "win_1": (0.0, 0.0)},
         "PAR_B": {"win_0": (0.0, 0.0), "win_1": (0.0, 0.0)}},
    )

    sel = pool.select_regional(query, [par_a, par_b], window=1500, per_window=1, logger=logger)

    assert [p.stem for p in sel.selected] == ["PAR_A"]
```

The existing `test_select_regional_drops_siblings_and_picks_regional_parents` (3000 bp, two windows) is the guard that a query of exactly two windows is unchanged.

- [ ] **Step 2: Run them, expected failure**

```bash
pytest tests/unit/test_pool.py -k "short_query or no_longer or falls_back" -q
```

Observed:
- `test_select_regional_cuts_a_short_query_into_two_windows`: `assert {'PAR_A'} == {'PAR_A', 'PAR_B'}`.
- `test_select_regional_falls_back_to_the_whole_query_when_no_half_matches`: `assert [] == ['PAR_A']` (new behaviour; see above).
- `test_select_regional_keeps_one_window_for_a_query_no_longer_than_it`: passes before and after (guard).

- [ ] **Step 3: Implement**

`src/tessera/discover/pool.py`, `select_regional`:

```diff
     query_seq = records[0][1].replace("-", "")
     n_windows = max(1, len(query_seq) // window)
+    # Integer division gives a query between one and two windows long a single window
+    # over the whole of it, which is whole-genome selection under another name (a 2.9 kb
+    # genome at the default 1500 bp). Cut such a query into two equal halves instead.
+    split_short = n_windows == 1 and len(query_seq) > window
+    if split_short:
+        n_windows = 2
+        window = -(-len(query_seq) // 2)  # ceiling: the two halves cover the query
     selected: list[Path] = []
```

After the `with tempfile.TemporaryDirectory() as tmp:` loop, before the final `logger.info("Regional selection: ...")`:

```diff
+    if not selected and split_short:
+        # skani is a genome-scale tool and may report nothing for a half-length chunk.
+        # Fall back to the whole-query ranking the single window used to give, rather
+        # than return an empty panel for a short genome.
+        logger.info(
+            "No pool genome matched either half of the %d bp query; selecting on the "
+            "whole query instead.", len(query_seq),
+        )
+        ranked = sorted(
+            (g for g in related if whole[g][0] > 0), key=lambda g: whole[g][0], reverse=True
+        )
+        selected = ranked[:per_window]
     logger.info(
         "Regional selection: %d genome(s) cover the query across %d window(s).",
```

`docs/reference-panels.md`, the `local` bullet:

```diff
   the divergent parents are no longer out-ranked. The query is then matched region by
-  region and the best-covering genomes are kept -- a focused parental panel, offline
-  and reproducible.
+  region (in windows of `--seed-window`, 1500 bp by default; a query between one and two
+  windows long is cut into two halves) and the best-covering genomes are kept -- a
+  focused parental panel, offline and reproducible.
```

`CHANGELOG.md`, `Panel recruitment (fourth audit):` group:

```markdown
- Regional pool selection treated a query between one and two selection windows long
  (1.5-3 kb at the default) as a single window, which is whole-genome selection. Such a
  query is now cut into two halves, with the whole-genome ranking as a fallback if neither
  half matches. Longer queries are unaffected.
```

- [ ] **Step 4: Run, expected pass**

```bash
pytest tests/unit/test_pool.py tests/unit/test_select_from_lineage.py tests/unit/test_lineage_select.py -q
```

Gate: the three tests above and the two existing `select_regional` tests. No harness reaches this branch (see "What gates these changes").

- [ ] **Step 5: Full checks**

```bash
ruff check src tests validation && mypy && pytest -m "not requires_binary" -q
```

Expected: `779 passed, 5 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/discover/pool.py tests/unit/test_pool.py docs/reference-panels.md CHANGELOG.md
git commit -F - <<'EOF'
Cut a query shorter than two selection windows into two halves

select_regional divides the query length by the window with integer
division, so a query between one and two windows long (a 2.9 kb genome
at the default 1500 bp) got one window over the whole of it: regional
selection degenerated to whole-genome selection for small genomes
(fourth audit, item D9).

Such a query is now cut into two equal halves. Longer queries keep
their window grid exactly, so no genome of two or more windows, and no
hybrid-harness case, is affected. Rounding the window count was not
used because it also moves the last window of long genomes.

Not verified: skani's behaviour on 750-1450 bp chunks (no binary here).
If neither half selects a genome, the whole-genome ranking already
computed supplies the panel, which is what the single window gave.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 9: Leave `--fetch-limit` out of the cache key  (audit items: D8)

`_fetch_diverse` keyed the NCBI Virus cache on `fetch_limit` (`iterate.py:762`) although the fetch never passes it (`iterate.py:776-792`); its only effect is the notice at line 793. Changing the flag therefore forced a full re-download of an identical genome set.

**Which way.** The CLI help (`cmd_fill_references.py:75-79`: "Notify (and dereplicate locally) when a broad NCBI Virus fetch returns more than this many complete genomes; the full set is downloaded either way"), the docs (`docs/reference-panels.md:246-249`: "A broad fetch is not truncated") and `_fetch_diverse`'s docstring ("rather than truncated to a biased sample") all describe a notice threshold, and `test_fetch_diverse_broadens_downloads_all_and_caches` pins that the limit is not passed. Applying it to the fetch would contradict all four and keep only the first N records in file order, a biased sample. So the option stays a notice threshold and leaves the cache key.

**Files:**
- Modify: `src/tessera/discover/iterate.py` (`FillParams.fetch_limit` comment line 92; `scope` line 762)
- Modify: `CHANGELOG.md`
- Test: `tests/unit/test_iterate.py`

**Interfaces:**
- The cache scope passed to `ncbi_virus_cache` and `write_cache_manifest` becomes `{"source_refseq": ...}`.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_iterate.py`:

```python


# --- audit 4: --fetch-limit is not part of the cache key (D8) ----------------------

def test_fetch_limit_does_not_split_the_cache(monkeypatch, tmp_path, logger):
    """--fetch-limit only sets when a notice is logged; the full set is downloaded either
    way (docs/reference-panels.md, and the option's help). As a cache-key field it forced a
    full re-download of an identical genome set whenever the value changed."""
    from tessera.discover import pool as pool_mod

    calls = {"fetch": 0}

    def fake_fetch(taxon, d, *, refseq=True, complete_only=False, released_after=None,
                   limit=None, logger):
        calls["fetch"] += 1
        written = []
        for i in range(4):
            p = d / f"G{i}.fasta"
            p.write_text(">x\nA\n")
            written.append(p)
        return written

    monkeypatch.setattr(pool_mod, "fetch_ncbi_virus", fake_fetch)
    out = tmp_path / "o"
    out.mkdir()
    common = {"query": tmp_path / "q.fasta", "collection": None, "output": out,
              "taxon": "SARS-CoV-2", "cache_dir": tmp_path / "cache"}

    first = iterate._fetch_diverse(FillParams(fetch_limit=2000, **common), logger)
    second = iterate._fetch_diverse(FillParams(fetch_limit=5, **common), logger)

    assert calls["fetch"] == 1
    assert [p.name for p in second] == [p.name for p in first]
```

- [ ] **Step 2: Run it, expected failure**

```bash
pytest tests/unit/test_iterate.py -k fetch_limit -q
```

Observed: `assert 2 == 1`.

- [ ] **Step 3: Implement**

`src/tessera/discover/iterate.py`:

```diff
-    fetch_limit: int = 2000  # cap a broad NCBI Virus fetch (complete genomes)
+    fetch_limit: int = 2000  # log a notice when a broad NCBI Virus fetch exceeds this
```

```diff
     # Keyed by what was asked for, not by what the fetch settled on: the RefSeq ->
     # complete broadening below is how a request is satisfied, not a separate request.
-    scope = {"source_refseq": params.source_refseq, "fetch_limit": params.fetch_limit}
+    # --fetch-limit is not part of the request: it only sets when a notice is logged, and
+    # the full set is downloaded whatever its value.
+    scope = {"source_refseq": params.source_refseq}
```

`CHANGELOG.md`, `Panel recruitment (fourth audit):` group:

```markdown
- `--fetch-limit` was part of the NCBI Virus cache key although it does not change what
  is fetched (it only sets when a notice is logged), so changing it re-downloaded an
  identical genome set. It no longer is. A cache written by an earlier version is keyed
  differently and is fetched once more on first use.
```

- [ ] **Step 4: Run, expected pass**

```bash
pytest tests/unit/test_iterate.py tests/unit/test_cache.py -q
```

Gate: the test above and `test_fetch_scope_does_not_share_a_cache_slot` (RefSeq and complete requests still get separate slots), unchanged.

- [ ] **Step 5: Full checks**

```bash
ruff check src tests validation && mypy && pytest -m "not requires_binary" -q
```

Expected: `780 passed, 5 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/discover/iterate.py tests/unit/test_iterate.py CHANGELOG.md
git commit -F - <<'EOF'
Leave --fetch-limit out of the NCBI Virus cache key

--fetch-limit was in the cache scope but never reaches the fetch; its
only effect is a log notice. Changing it forced a full re-download of
an identical genome set (fourth audit, item D8).

The option's help, the docs, the _fetch_diverse docstring and an
existing test all describe it as a notice threshold over a full,
untruncated download, so it is removed from the key rather than applied
to the fetch (which would keep the first N records in file order, a
biased sample).

Trade-off: caches written by earlier versions carry the old key and are
fetched once more on first use.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 10: Use the shared non-clade marker set  (audit items: D10)

`nextclade.py:50` defines `NON_CLADE_MARKERS = frozenset({"", "example", "NA", "?"})` and says consumers should treat all four as "no clade"; `lineage_assign._reference_tips` (`lineage_assign.py:121`) hardcoded `("NA", "example")`, so a tip labelled `?` became a clade named `?` in nearest-reference typing.

**Files:**
- Modify: `src/tessera/discover/lineage_assign.py` (import line 21; line 121)
- Modify: `CHANGELOG.md`
- Test: `tests/unit/test_lineage_assign.py`

**Interfaces:**
- Consumes: `tessera.discover.nextclade.NON_CLADE_MARKERS`. `lineage_assign` already imports from `.nextclade`, so no new import cycle.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_lineage_assign.py`:

```python


def test_reference_tips_treat_every_non_clade_marker_as_no_clade(monkeypatch, tmp_path):
    """nextclade.NON_CLADE_MARKERS lists four tokens that mean "no clade"; the tip parser
    hardcoded two of them, so a tip labelled "?" became a clade named "?"."""
    from tessera.discover.nextclade import NON_CLADE_MARKERS

    labels = ["B.1", *sorted(NON_CLADE_MARKERS)]
    tips = [tmp_path / f"t{i}.fasta" for i in range(len(labels))]
    headers = {
        t: f"acc{i} {label}".strip()
        for i, (t, label) in enumerate(zip(tips, labels, strict=True))
    }
    for p in tips:
        p.write_text(">x\nACGT\n")

    class _DS:
        path = "nextstrain/x"
        tag = "2020"

    monkeypatch.setattr(la, "resolve_dataset", lambda q, o, *, email, logger: _DS())
    monkeypatch.setattr(la, "build_pool", lambda ds, *, cache_dir, logger: tips)
    monkeypatch.setattr(la, "first_header", lambda p: headers[p])

    reps = la._reference_tips(query=_g("q"), nextclade_dataset=None, email=None,
                              cache_dir=tmp_path, logger=logging.getLogger("t"))

    assert list(reps.values()) == ["B.1"]
```

- [ ] **Step 2: Run it, expected failure**

```bash
pytest tests/unit/test_lineage_assign.py -k non_clade -q
```

Observed: `assert ['B.1', '?'] == ['B.1']`.

- [ ] **Step 3: Implement**

`src/tessera/discover/lineage_assign.py`:

```diff
-from .nextclade import build_pool, resolve_dataset
+from .nextclade import NON_CLADE_MARKERS, build_pool, resolve_dataset
```

```diff
         clade = parts[1].strip() if len(parts) > 1 else ""
-        if clade and clade not in ("NA", "example"):
+        if clade not in NON_CLADE_MARKERS:
             clade_of[tip] = clade
```

(`""` is in the set, so the empty-clade case the old `if clade` handled is still excluded.)

`CHANGELOG.md`, `Panel recruitment (fourth audit):` group:

```markdown
- Nearest-reference typing took a Nextclade tip labelled `?` for a clade named `?`; it
  now uses the shared set of "no clade" markers.
```

- [ ] **Step 4: Run, expected pass**

```bash
pytest tests/unit/test_lineage_assign.py -q
```

Gate: the test above and `test_reference_tips_parses_and_picks_one_per_clade`, unchanged.

- [ ] **Step 5: Full checks**

```bash
ruff check src tests validation && mypy && pytest -m "not requires_binary" -q
```

Expected: `781 passed, 5 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/discover/lineage_assign.py tests/unit/test_lineage_assign.py CHANGELOG.md
git commit -F - <<'EOF'
Use the shared non-clade marker set in nearest-reference typing

nextclade.NON_CLADE_MARKERS names four tokens that mean "no clade" and
says consumers should treat all four that way; _reference_tips
hardcoded two, so a tip labelled "?" became a clade named "?" (fourth
audit, item D10). It now uses the shared set.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

### Task 11: Refuse seed options that cannot act; read a pool like a collection  (audit items: E5 -- `--candidate-pool` with `--nextclade`, inert `--pool-consensus`, the candidate-pool reader's extension list)

Three related cases where an option did nothing and the run did not say so:
- `detect` and `build-panel` took `--candidate-pool` over `--nextclade`/`--nextclade-dataset` in silence (`iterate.py:157-162`, `cmd_build_panel.py:28-34`). The precedence was stated only in code docstrings; neither command's help nor `docs/reference-panels.md` defines one, and `--nextclade-dataset`'s help says it "implies --nextclade". **Refused** with a `UserInputError`, as contradictory inputs.
- `--pool-consensus` acts only in `_fetch_nextclade` (`iterate.py:731-732`); with any other seed source it did nothing. **Refused** rather than warned: the user asked for a consensus pool and would silently get a different panel, the contradiction is known before any network work, and the remedy (drop the flag or add `--nextclade`) is one word.
- In `fill-references`, `--candidate-pool` is read only under `--seed-source local` (`iterate.py:708-711`); with the default `blast` the pool was ignored (Review Focus item 3). **Refused** by the same check.

Both refusals live in `FillParams.__post_init__`, so they cover all three commands and run inside each command's `stage_errors` block before any work. A new `seed_source_for` holds the `detect`/`build-panel` rule in one place (it was duplicated).

The pool reader (`iter_pool_genomes`, `pool.py:68-78`) accepted `.fasta`, `.fa`, `.fna`, `.gz` and dropped anything else, while every collection reader takes each non-hidden file (`collection_genomes`, `core/io.py:118-128`) and rejects a non-FASTA file by content (`_require_fasta`, `core/io.py:200-215`); `type-lineages` already reads that way (`cmd_type_lineages.py:58-63`). The pool reader now does the same.

Not covered: with `--collection` given, seeding is skipped and every seed option is inert; that is not refused here (see Self-review notes).

**Files:**
- Modify: `src/tessera/discover/iterate.py` (`FillParams` after line 120; `for_detection` lines 148-162; new `seed_source_for` before `RoundResult`)
- Modify: `src/tessera/cli/cmd_build_panel.py` (remove `_seed_source`, lines 28-34; import at line 91; call at line 110; `--candidate-pool` help lines 60-63)
- Modify: `src/tessera/cli/cmd_detect.py` (`--candidate-pool` help lines 44-47; `--pool-consensus` help lines 79-83)
- Modify: `src/tessera/discover/pool.py` (imports line 29; `iter_pool_genomes` lines 68-78)
- Modify: `docs/reference-panels.md` (after the `fill-references`/`detect` examples, line 230), `CHANGELOG.md`
- Test: `tests/unit/test_cli_input_checks.py`, `tests/unit/test_cli_nextclade.py`, `tests/unit/test_pool.py`

**Interfaces:**
- Produces: `seed_source_for(candidate_pool: Path | None, nextclade: bool, nextclade_dataset: str | None) -> str` in `tessera.discover.iterate`; raises `UserInputError` when a pool and a Nextclade request are both given.
- Produces: `FillParams.__post_init__(self) -> None`, raising `UserInputError` for `pool_consensus` without `seed_source == "nextclade"` and for `candidate_pool` without `seed_source == "local"`. Every existing `FillParams(...)` construction in `tests/` and `validation/` was checked and none combines these.
- Changes: `iter_pool_genomes(pool: Path) -> list[Path]` returns `collection_genomes(pool)` and raises `UserInputError` for a file that is not FASTA. `cmd_build_panel._seed_source` is removed (no test refers to it).

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_cli_input_checks.py`:

```python


# --- seed-source options that cannot act together (audit 4, E5) ------------------

@pytest.mark.parametrize("command", ["detect", "build-panel"])
@pytest.mark.parametrize(
    "nextclade_args",
    [["--nextclade"], ["--nextclade-dataset", "nextstrain/sars-cov-2/XBB"]],
)
def test_a_candidate_pool_and_a_nextclade_request_are_refused_together(
    tmp_path: Path, no_pipeline, command: str, nextclade_args: list[str]
) -> None:
    """Both choose where the panel is recruited from. The pool used to win in silence,
    so a run asked to use Nextclade clades used the pool and did not say so."""
    pool_dir = tmp_path / "pool"
    pool_dir.mkdir()
    result = runner.invoke(app, [
        command, "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        "--candidate-pool", str(pool_dir), *nextclade_args,
    ])
    _rejected(result, "--candidate-pool and --nextclade/--nextclade-dataset")


@pytest.mark.parametrize(
    ("command", "args"),
    [("detect", []),
     ("detect", ["--candidate-pool", "POOL"]),
     ("fill-references", ["--seed-source", "local", "--candidate-pool", "POOL"]),
     ("fill-references", ["--seed-source", "ncbi-virus"])],
)
def test_pool_consensus_without_a_nextclade_pool_is_refused(
    tmp_path: Path, no_pipeline, command: str, args: list[str]
) -> None:
    """--pool-consensus acts only on a Nextclade pool; under any other seed source it
    did nothing and the run looked like one that had used it."""
    pool_dir = tmp_path / "pool"
    pool_dir.mkdir()
    args = [str(pool_dir) if a == "POOL" else a for a in args]
    result = runner.invoke(app, [
        command, "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        "--pool-consensus", *args,
    ])
    _rejected(result, "--pool-consensus")


def test_fill_references_refuses_a_candidate_pool_under_another_seed_source(
    tmp_path: Path, no_pipeline
) -> None:
    """fill-references reads --candidate-pool only with --seed-source local; without it
    the default BLAST seeding ran and the pool was ignored."""
    pool_dir = tmp_path / "pool"
    pool_dir.mkdir()
    result = runner.invoke(app, [
        "fill-references", "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        "--candidate-pool", str(pool_dir),
    ])
    _rejected(result, "--seed-source local")
```

Append to `tests/unit/test_cli_nextclade.py`:

```python


def test_detect_accepts_pool_consensus_with_nextclade(monkeypatch, tmp_path):
    """Guard (passes before and after): the one combination --pool-consensus is for."""
    captured = _capture(monkeypatch)
    q = tmp_path / "q.fasta"
    q.write_text(">q\nACGT\n")
    result = runner.invoke(app, ["detect", "-q", str(q), "-o", str(tmp_path / "o"),
                                 "--nextclade", "--pool-consensus"])
    assert result.exit_code == 0, result.output
    assert captured["params"].pool_consensus is True
```

Append to `tests/unit/test_pool.py`:

```python


# --- audit 4: a candidate pool is read like a collection (E5) ----------------------

def test_candidate_pool_reads_the_same_files_as_a_collection(tmp_path):
    """The pool reader accepted four extensions where every collection reader takes each
    non-hidden file, so a .fas or .fsa genome was dropped from the pool in silence."""
    import gzip

    from tessera.core.io import collection_genomes

    for name in ("a.fasta", "b.fas", "c.fsa", "d"):
        (tmp_path / name).write_text(f">{name}\nACGT\n")
    with gzip.open(tmp_path / "e.fa.gz", "wt") as fo:
        fo.write(">e\nACGT\n")
    (tmp_path / ".DS_Store").write_bytes(b"\x00\x01")

    genomes = pool.iter_pool_genomes(tmp_path)

    assert genomes == collection_genomes(tmp_path)
    assert [p.name for p in genomes] == ["a.fasta", "b.fas", "c.fsa", "d", "e.fa.gz"]


def test_candidate_pool_refuses_a_file_that_is_not_fasta(tmp_path):
    import pytest

    from tessera.core.errors import UserInputError

    (tmp_path / "a.fasta").write_text(">a\nACGT\n")
    (tmp_path / "README.txt").write_text("pool assembled 2026-09\n")

    with pytest.raises(UserInputError, match="does not look like a FASTA file"):
        pool.iter_pool_genomes(tmp_path)
```

- [ ] **Step 2: Run them, expected failure**

```bash
pytest tests/unit/test_cli_input_checks.py tests/unit/test_cli_nextclade.py tests/unit/test_pool.py -q
```

Observed:
- all nine new CLI cases: `assert 'Unexpected error' not in ...`, output `[ERROR] Unexpected error: the command reached the pipeline; it should have been rejected`.
- `test_candidate_pool_reads_the_same_files_as_a_collection`: list mismatch, `Right contains 3 more items, first extra item: .../c.fsa` (only `a.fasta` and `e.fa.gz` were read).
- `test_candidate_pool_refuses_a_file_that_is_not_fasta`: `Failed: DID NOT RAISE <class 'tessera.core.errors.UserInputError'>`.
- `test_detect_accepts_pool_consensus_with_nextclade`: passes before and after (guard).

- [ ] **Step 3: Implement**

`src/tessera/discover/iterate.py`, inside `FillParams`, after the last field (`deep_typing`) and before `for_detection`:

```python
    def __post_init__(self) -> None:
        # Options that act through one seed source are refused under another rather than
        # ignored: a run that dropped a requested pool or consensus in silence reads like
        # one that used it.
        if self.pool_consensus and self.seed_source != "nextclade":
            raise UserInputError(
                "--pool-consensus builds one consensus genome per clade of a Nextclade "
                f"pool, but the seed source here is '{self.seed_source}'. Add --nextclade "
                "(detect, build-panel) or --seed-source nextclade (fill-references), or "
                "leave out --pool-consensus."
            )
        if self.candidate_pool is not None and self.seed_source != "local":
            raise UserInputError(
                "--candidate-pool is read only with --seed-source local, but the seed "
                f"source here is '{self.seed_source}'. Add --seed-source local, or leave "
                "out --candidate-pool."
            )
```

`for_detection`:

```diff
         Recruit the parental lineages (drop the query's siblings, curate the
         panel) so a recombinant query whose own lineage is common in NCBI is
-        not masked. The seed source follows the supplied input: a local pool
-        wins, then a Nextclade dataset, then the default BLAST recruitment.
-        This is the single home for the preset so ``tessera detect`` does not
-        duplicate (and drift from) the field list.
+        not masked. The seed source follows the supplied input
+        (:func:`seed_source_for`). This is the single home for the preset so
+        ``tessera detect`` does not duplicate (and drift from) the field list.
         """
-        if candidate_pool:
-            seed_source = "local"
-        elif nextclade or nextclade_dataset:
-            seed_source = "nextclade"
-        else:
-            seed_source = "blast"
+        seed_source = seed_source_for(candidate_pool, nextclade, nextclade_dataset)
         return cls(
```

New module-level function, between the `FillParams` class and `@dataclass class RoundResult`:

```python
def seed_source_for(
    candidate_pool: Path | None, nextclade: bool, nextclade_dataset: str | None
) -> str:
    """The seed source ``detect`` and ``build-panel`` take from their inputs.

    A local pool gives ``local``, a Nextclade flag or dataset ``nextclade``, and neither
    the default BLAST recruitment. A pool together with a Nextclade request names two
    seed sources; the pool used to win in silence, so it is refused instead.
    """
    if candidate_pool and (nextclade or nextclade_dataset):
        raise UserInputError(
            "--candidate-pool and --nextclade/--nextclade-dataset each choose where the "
            "panel is recruited from; give one of them."
        )
    if candidate_pool:
        return "local"
    if nextclade or nextclade_dataset:
        return "nextclade"
    return "blast"
```

`src/tessera/cli/cmd_build_panel.py`:

```diff
-def _seed_source(candidate_pool, nextclade: bool, nextclade_dataset: str | None) -> str:
-    """Local pool wins, then Nextclade, then the default BLAST recruitment."""
-    if candidate_pool:
-        return "local"
-    if nextclade or nextclade_dataset:
-        return "nextclade"
-    return "blast"
-
-
 @app.command(name="build-panel")
```

```diff
     candidate_pool: Path | None = typer.Option(
         None, "--candidate-pool",
-        help="Use a local genome pool instead of fetching from NCBI (for mega-taxa).",
+        help="Use a local genome pool instead of fetching from NCBI (for mega-taxa). "
+        "Cannot be combined with --nextclade or --nextclade-dataset.",
     ),
```

```diff
-    from ..discover.iterate import FillParams, fill_references
+    from ..discover.iterate import FillParams, fill_references, seed_source_for
```

```diff
-            seed_source=_seed_source(candidate_pool, nextclade, nextclade_dataset),
+            seed_source=seed_source_for(candidate_pool, nextclade, nextclade_dataset),
```

`src/tessera/cli/cmd_detect.py`:

```diff
     candidate_pool: Path | None = typer.Option(
         None, "--candidate-pool",
-        help="Use a local genome pool instead of fetching from NCBI (for mega-taxa).",
+        help="Use a local genome pool instead of fetching from NCBI (for mega-taxa). "
+        "Cannot be combined with --nextclade or --nextclade-dataset.",
     ),
```

```diff
     pool_consensus: bool = typer.Option(
         False, "--pool-consensus/--no-pool-consensus",
-        help="With a Nextclade pool, use one denoised consensus genome per clade (a "
-        "stable per-lineage reference) instead of every tree tip.",
+        help="With a Nextclade pool (--nextclade or --nextclade-dataset), use one "
+        "denoised consensus genome per clade (a stable per-lineage reference) instead of "
+        "every tree tip. Refused without a Nextclade pool.",
     ),
```

`src/tessera/discover/pool.py`:

```diff
-from ..core.io import read_fasta, strip_sequence_extension, write_fasta_record
+from ..core.io import (
+    _require_fasta,
+    collection_genomes,
+    read_fasta,
+    strip_sequence_extension,
+    write_fasta_record,
+)
```

```diff
 def iter_pool_genomes(pool: Path) -> list[Path]:
-    """List the genome FASTA files in a pool directory."""
+    """List the genome files in a pool directory, read the way a collection is read.
+
+    Every non-hidden file is a genome, whatever its extension and gzip-compressed or
+    not, and a file that is not FASTA is an error rather than something to skip. A list
+    of accepted extensions dropped a ``.fas`` or ``.fsa`` genome in silence.
+    """
     if not pool.is_dir():
         raise UserInputError(f"Candidate pool {pool} is not a directory.")
-    genomes = sorted(
-        p for p in pool.iterdir()
-        if p.is_file() and p.suffix.lower() in {".fasta", ".fa", ".fna", ".gz"}
-    )
+    genomes = collection_genomes(pool)
     if not genomes:
-        raise UserInputError(f"Candidate pool {pool} has no FASTA genome files.")
+        raise UserInputError(f"Candidate pool {pool} has no genome files.")
+    for genome in genomes:
+        _require_fasta(genome)
     return genomes
```

`docs/reference-panels.md`, after the second example block of the seed-source section (the `detect --nextclade` / `fill-references --seed-source nextclade` block, line 230):

```markdown
`--candidate-pool` and `--nextclade`/`--nextclade-dataset` each choose where `detect` and
`build-panel` recruit from, so giving both is refused. `--pool-consensus` applies only to
a Nextclade pool and is refused with any other seed source, and in `fill-references`
`--candidate-pool` needs `--seed-source local`. A candidate pool is read like a
collection: every non-hidden file is a genome, whatever its extension, and a file that is
not FASTA is an error.
```

`CHANGELOG.md`, `Panel recruitment (fourth audit):` group:

```markdown
- The candidate-pool reader accepted four file extensions, so a `.fas` or `.fsa` genome
  was left out of the pool without notice. It now reads a pool as every other command
  reads a collection, and rejects a file that is not FASTA.
```

and in `### Changed`:

```markdown
- **Seed options that cannot act are refused.** `detect` and `build-panel` used
  `--candidate-pool` and ignored `--nextclade`/`--nextclade-dataset` when both were given;
  `--pool-consensus` did nothing without a Nextclade pool; and `fill-references` ignored
  `--candidate-pool` without `--seed-source local`. Each is now a usage error.
```

- [ ] **Step 4: Run, expected pass**

```bash
pytest tests/unit/test_cli_input_checks.py tests/unit/test_cli_nextclade.py tests/unit/test_pool.py tests/unit/test_iterate.py tests/unit/test_cli_keep_recombinant.py tests/unit/test_cli_deep_typing.py -q
```

Gate: the tests above, plus the existing `test_cli_nextclade.py` routing tests and the `_seed_from_pool` tests in `test_iterate.py` (`test_pool_seeding_skips_the_query_own_record_and_excluded_accessions` and the three `test_pool_exclusion_*`), which construct `FillParams(seed_source="local", candidate_pool=...)` and read pools through `iter_pool_genomes`; all unchanged.

- [ ] **Step 5: Full checks**

```bash
ruff check src tests validation && mypy && pytest -m "not requires_binary" -q
```

Expected: `793 passed, 5 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/discover/iterate.py src/tessera/discover/pool.py \
        src/tessera/cli/cmd_build_panel.py src/tessera/cli/cmd_detect.py \
        tests/unit/test_cli_input_checks.py tests/unit/test_cli_nextclade.py \
        tests/unit/test_pool.py docs/reference-panels.md CHANGELOG.md
git commit -F - <<'EOF'
Refuse seed options that cannot act, and read a pool like a collection

Three options did nothing in some combinations and the run did not say
so (fourth audit, item E5): detect and build-panel took --candidate-pool
over --nextclade/--nextclade-dataset; --pool-consensus was inert under
any seed source but nextclade; and fill-references ignored
--candidate-pool without --seed-source local. The pool precedence was
stated only in code docstrings, and --nextclade-dataset's help says it
implies --nextclade, so the combination is treated as contradictory.

All three are refused before any work, from FillParams.__post_init__
and a shared seed_source_for (which replaces build-panel's copy of the
rule). Refusing rather than warning: the remedy is one flag, and a
warning in a long log is easy to miss when the panel silently differs.

The candidate-pool reader accepted four extensions and dropped the rest,
while every collection reader takes each non-hidden file and rejects
non-FASTA content. It now reads the same way (collection_genomes and
_require_fasta, as type-lineages does), so a .fas genome is no longer
lost and a stray README is an error, not a skipped file.

Not changed: with --collection given, seeding is skipped and the seed
options are inert; that case is not refused here.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
```

---

## Self-review notes

**Verification.** Every test in this plan was run in a detached worktree of `4b35e39` with `PYTHONPATH=src`: each failed for the reason recorded in its Step 2 before the candidate implementation, and passed after it. With all eleven tasks applied the fast suite gave `793 passed, 5 deselected` (757 before), `ruff check src tests validation` was clean and `mypy` reported no issues in 79 source files. The per-task counts in Step 5 follow from the number of tests each task adds. No network call, aligner or harness was run.

**Existing tests that change.** Only two, both in Task 2: `tests/unit/test_discover.py::test_curate_after_download_keeps_the_users_own_genomes` and `::test_curate_after_download_into_a_new_directory`. They pinned the bulk download; their expected file lists change and their intent is kept.

**Deviations from the lead's guidance, with reasons.**
- D2 uses neither the under-covered base count nor the mean over searched gaps (Task 1 explains both failures; the base count stays at 6200 bp through three improving rounds in the test).
- D1 ranks by identity times query coverage, not by identity alone (Review Focus item 2).
- D1 is paired with excluding curation-removed genomes and an unchanged-panel exemption from the stall test. Without them the cap regresses `detect` whenever a gap's best hit is a sibling.

**Items not planned, with evidence.**
- *Suspected: the negative-lineage seeder searches only the first 1500 bases.* Confirmed as a code/docstring mismatch: `_seed_negative_lineage` says "from a whole-query BLAST" (`iterate.py:867`) and BLASTs `query_seq[:1500]` (`iterate.py:872`). The consequence for a recombinant query (the dominant token being one parent's lineage) is not reproduced. A whole-genome BLAST is much slower, so whether to search the whole query, sample several windows, or correct the docstring is the maintainer's choice.
- *Suspected: reference typing resolves its Nextclade dataset from whichever untyped genome sorts first* (`lineage_assign.py:153`, `query=untyped[0]`). Confirmed in code; the cost is a wasted call and an untyped panel only when an outgroup sorts first. Not reproduced, not planned.
- *Curating a seeded panel before round 1* (D6's alternative remedy). `test_curate_does_not_run_on_a_freshly_seeded_collection` records the current behaviour as deliberate, and changing it alters what `detect` recruits; Task 6 reports instead.
- *Seed options inert when `--collection` is given* (seeding is skipped). Not refused in Task 11; refusing would break a workflow where a user keeps one command line and adds `--collection`. A maintainer decision.
- *Discover-side test gap "a pool build over an existing cache".* The Nextclade pool cache install (`nextclade.build_pool`, its own temp-dir install) belongs with the core item about pool installs under concurrency and stray files (spec G7); it is left to the plan that owns `core/cache.py`.

**Unverified assumptions.**
- Task 4 follows the audit's reproduction of an empty package (no `genomic.fna`). If the `datasets` CLI instead exits non-zero for a scope with no genomes, `run_tool` raises and the run still stops; that path is unchanged and untested.
- Task 8's fallback exists because skani's output on 750-1450 bp chunks was not measured. With a binary present, one run on a hepatitis B or delta genome would show whether the halves select anything.
- Tasks 1 and 2 change real-data recruitment in ways no offline harness measures (Task 2 Step 7).

**Corrections to the spec.**
- D1's magnitudes (60 genomes for a 1000 bp gap, 435 for 6 kb) correspond to 15 hits per sub-interval. At the default `--max-hits 5` the bounds are 20 and 145 per gap per round (4 and 29 sub-intervals). The defect stands; the numbers are about three times the default case.
- The brief says the hybrid harness "does not call the CLI selection path". It does call `select_regional` directly (`validation/run_hybrids.py:1093`), but with `drop_siblings=False` and a selection window of at most a quarter of the genome, so neither changed branch (Tasks 7, 8) is reachable from it.
- D9 names hepatitis B as affected. At about 3.2 kb and the default 1500 bp window it gets two windows already; genomes of 1.5-3 kb (hepatitis delta at about 1.7 kb, for example) are the affected range.
