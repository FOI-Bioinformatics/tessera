# Plurality best-match major parent Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Derive the HMM caller's major parent as the plurality best-match genome (most windows won) instead of the most-occupied Viterbi state, so a uniformly-mediocre sibling can no longer be chosen as the backbone.

**Architecture:** Add a `plurality_major(result, labels)` helper in `recomb/hmm.py`; call it from `segment_query` in place of `_major_state`; remove the now-unused `_major_state`. Everything downstream (region calling, `reconcile_major`, the verdict) is unchanged. Then re-run the attribution harness as a mandatory gate.

**Tech Stack:** Python 3 (numpy), pytest, ruff. Design spec: `docs/superpowers/specs/2026-07-01-plurality-major-parent-design.md`.

## Global Constraints

- Only the HMM caller's major derivation changes. No change to `select_regional`/panel selection, the HMM segmentation model, sibling exclusion, the 3SEQ/MaxChi/Bootscan/barcode callers, or `reconcile_major`.
- `plurality_major` = the label that is the strict best match (highest similarity) in the most informative windows; ties broken by highest mean similarity over informative windows, then deterministic label order; `None` when no window is informative.
- No new runtime dependency (numpy is already imported in `hmm.py`); do not import private helpers across modules -- compute the per-window winner inline with numpy.
- Modest scientific language in code/docstrings. Ruff clean (line length 100, rules E/F/I/UP/B); CI coverage floor (80%) holds.
- Commit messages end with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.
- Branch `plurality-major-spec` (already checked out off main); do not commit on main.

---

### Task 1: `plurality_major` helper + `segment_query` swap

**Files:**
- Modify: `src/tessera/recomb/hmm.py` (add `plurality_major`; edit `segment_query` line 86; remove `_major_state` lines 205-207)
- Test: `tests/unit/test_hmm.py` (add tests; keep the existing ones)

**Interfaces:**
- Consumes: `WindowSimilarity` (`recomb/similarity.py`) with `.similarities: dict[str, list[float]]` and `.positions: list`.
- Produces: `plurality_major(result: WindowSimilarity, labels: list[str]) -> str | None`; `segment_query` returns this as its second element (same signature as today).

- [ ] **Step 1: Write the failing tests**

Add to `tests/unit/test_hmm.py` (keep the existing imports and tests; add `plurality_major` to the import):

```python
from tessera.recomb.hmm import plurality_major, segment_query


def _wsim(similarities: dict[str, list[float]]) -> WindowSimilarity:
    # plurality_major reads only .similarities; the other fields satisfy the
    # dataclass's required arguments.
    n = len(next(iter(similarities.values())))
    return WindowSimilarity(
        positions=list(range(n)), query_positions=list(range(n)),
        similarities=similarities, query="q", width=n,
        query_cumulative=np.arange(n + 1),
    )


def test_plurality_major_prefers_window_winner_over_uniform_genome():
    # X is the strict best in 7 of 10 windows; Y is uniformly second-best everywhere
    # (it would occupy the most HMM path but wins no window).
    x = [1.0] * 7 + [0.80] * 3
    y = [0.92] * 10
    z = [0.70] * 7 + [1.0] * 3  # wins the last 3
    assert plurality_major(_wsim({"X": x, "Y": y, "Z": z}), ["X", "Y", "Z"]) == "X"


def test_plurality_major_tiebreak_by_mean_similarity():
    # Both win 5 of 10 windows; A has the higher mean over informative windows.
    a = [1.0] * 5 + [0.0] * 5
    b = [0.6] * 5 + [0.61] * 5  # wins windows 5-9 (0.61 > 0.0), mean lower than A's wins
    # A wins 0-4 (1.0 > 0.6); B wins 5-9 (0.61 > 0.0). 5-5 tie -> mean similarity breaks it.
    assert plurality_major(_wsim({"A": a, "B": b}), ["A", "B"]) in {"A", "B"}
    # Make the tiebreak deterministic and checkable: A's mean (0.5) vs B's (0.605) -> B.
    assert plurality_major(_wsim({"A": a, "B": b}), ["A", "B"]) == "B"


def test_plurality_major_none_when_no_informative_window():
    nan = float("nan")
    assert plurality_major(_wsim({"A": [nan, nan], "B": [nan, nan]}), ["A", "B"]) is None


def test_plurality_major_single_label():
    assert plurality_major(_wsim({"A": [0.9, 0.8]}), ["A"]) == "A"
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/unit/test_hmm.py -k plurality_major -v`
Expected: FAIL with `ImportError: cannot import name 'plurality_major'`.

- [ ] **Step 3: Add the `plurality_major` helper**

In `src/tessera/recomb/hmm.py`, add this function (place it near `segment_query`, above it):

```python
def plurality_major(result: WindowSimilarity, labels: list[str]) -> str | None:
    """The label that is the strict best match in the most informative windows.

    This is the major parent as defined for the query -- the reference it most
    resembles across the genome -- rather than the state an HMM path happens to occupy
    most (which a uniformly-mediocre genome, best in no window, can win). Ties are broken
    by the highest mean similarity over informative windows, then by label order.
    Returns ``None`` when no window is informative.
    """
    if not labels:
        return None
    matrix = np.array([result.similarities[label] for label in labels], dtype=float)  # (S, W)
    wins = np.zeros(len(labels), dtype=int)
    for i in range(matrix.shape[1]):
        col = matrix[:, i]
        if np.isfinite(col).any():
            wins[int(np.nanargmax(col))] += 1
    if not wins.any():
        return None
    means = [
        float(np.nanmean(matrix[j])) if np.isfinite(matrix[j]).any() else float("-inf")
        for j in range(len(labels))
    ]
    return labels[max(range(len(labels)), key=lambda j: (int(wins[j]), means[j]))]
```

- [ ] **Step 4: Use it in `segment_query`; remove `_major_state`**

In `segment_query` (`recomb/hmm.py:86`), replace:

```python
    major = _major_state(path, labels)
```

with:

```python
    major = plurality_major(result, labels)
```

Then delete the now-unused `_major_state` function (lines 205-207):

```python
def _major_state(path: np.ndarray, labels: list[str]) -> str:
    counts = np.bincount(path, minlength=len(labels))
    return labels[int(np.argmax(counts))]
```

- [ ] **Step 5: Run the new tests to verify they pass**

Run: `python -m pytest tests/unit/test_hmm.py -k plurality_major -v`
Expected: PASS (4 passed).

- [ ] **Step 6: Run the full HMM + regions suites to confirm no regression**

Run: `python -m pytest tests/unit/test_hmm.py tests/unit/test_regions.py -v`
Expected: PASS -- including the existing `segment_query` tests that assert `major == "refA"` (refA is the dominant window-winner, so plurality agrees with the prior occupancy result).

- [ ] **Step 7: Full fast suite + ruff**

Run: `python -m pytest -m "not requires_binary" -q && ruff check src/tessera/recomb/hmm.py tests/unit/test_hmm.py`
Expected: all pass; ruff clean. (`np.bincount` was the only user of `path` in `_major_state`; `path` is still used by `_path_to_segments`, so no unused-variable warning.)

- [ ] **Step 8: Commit**

```bash
git add src/tessera/recomb/hmm.py tests/unit/test_hmm.py
git commit -m "$(cat <<'EOF'
Derive the HMM major parent by plurality best-match

segment_query now names the major parent as the genome that is the strict best
match in the most windows (tie-broken by mean similarity), matching the module
definition and the heuristic caller, instead of the most-occupied Viterbi state.
A uniformly-mediocre genome (best in no window) can no longer become the
backbone. Removes the now-unused _major_state.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Harness gate + record the result

**Files:**
- Modify: `validation/attribution-results.md` (append a section)

**Interfaces:** none (measurement + documentation).

- [ ] **Step 1: Run the attribution harness comparison (opt-in, needs the aligner env)**

Run:
```bash
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"   # skani/mafft/minimap2/skDER; keep the project python
python validation/run_hybrids.py --compare
```
Expected outcome to check (not assert programmatically): `rsv_a` flips to **PASS** (major parent now A.1, donor A.D.1.8), the other 15 runnable cases hold, `flu_h3n2_ha` still FAILs on its separate short-segment issue -> a clean **17/18**. This is a measurement gate, not a unit test; it needs the aligner env on PATH and is not part of CI.

- [ ] **Step 2: Record the result faithfully**

Append a section to `validation/attribution-results.md` stating the before (16/18) and after counts, that `rsv_a`'s major parent moved from the sibling A.D.5.2 to the true A.1 (donor A.D.1.8 unchanged), and any case that shifted. If any case regresses, record it plainly and flag that the merge decision is revisited -- do not present a regressed run as a clean pass. Use modest language; state what passed, what failed, and that the harness path (not CI) produced it.

- [ ] **Step 3: Commit**

```bash
git add validation/attribution-results.md
git commit -m "$(cat <<'EOF'
Record the plurality-major harness result

Attribution --compare re-run after the HMM major-parent change: rsv_a's major
moves from the sibling A.D.5.2 to the true A.1 (donor A.D.1.8 unchanged). Records
the measured pass count and any case that shifted.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Self-Review

**Spec coverage:**
- `plurality_major` helper (strict best-match plurality, mean-similarity tiebreak, `None` on no informative window) -> Task 1 Step 3.
- `segment_query` swap + `_major_state` removal -> Task 1 Step 4.
- CI unit tests (plurality winner over a uniform genome, tiebreak, none, single label) + existing HMM/region tests hold -> Task 1 Steps 1, 6.
- Mandatory harness gate (clean 17/18, faithful record) -> Task 2.
- No panel/segmentation/other-caller change; no new dependency -> Global Constraints, Task 1.

**Placeholder scan:** every code step shows complete code; the only prose-only steps are Task 2 (a measurement + a specified results write-up) -- no TODO/TBD.

**Type consistency:** `plurality_major(result: WindowSimilarity, labels: list[str]) -> str | None` matches its call in `segment_query` (which passes `result` and `labels`, both in scope there) and returns the same `str | None` `segment_query` already returns as its second element. `_major_state` is removed and has no remaining callers (only `segment_query` used it).

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-07-01-plurality-major-parent.md`. Two execution options:

1. **Subagent-Driven (recommended)** - fresh subagent per task, review between tasks, fast iteration.
2. **Inline Execution** - execute tasks in this session with checkpoints.

Which approach?
