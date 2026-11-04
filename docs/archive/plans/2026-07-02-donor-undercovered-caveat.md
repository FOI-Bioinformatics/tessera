# Donor-undercovered caveat Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** When a coverage gap overlaps a confident donor-present region, caveat that region (`donor_undercovered=True`) instead of silently dropping the gap, so the "true donor may be absent" signal survives alongside the closest-available attribution.

**Architecture:** Factor the gap/region reconciliation out of `run_recomb` into a pure `reconcile_gaps` helper in `recomb/coverage.py` (unit-testable without an aligner), change its rule from drop-on-overlap to caveat-on-overlap, and teach the harness scorer to accept the caveat. Reuses the existing `Region.donor_undercovered` field (already serialized to the TSV and shown in reports).

**Tech Stack:** Python 3 (numpy), pytest, ruff. Aligner env `recomfi-aln` for the opt-in gate. Design spec: `docs/superpowers/specs/2026-07-02-donor-undercovered-caveat-design.md`.

## Global Constraints

- Only the gap/region reconciliation and the harness scorer change. No change to the coverage threshold, the callers, `flag_undercovered_regions`, or the report layout.
- The 18 positive harness cases must not regress (true donor present -> no coverage gap over the donor region -> no new caveats). Verified on the gate (Task 3).
- No new runtime dependency. Modest language; ruff clean (line length 100, E/F/I/UP/B); CI coverage floor (80%) holds.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. Branch `donor-undercovered-caveat` (checked out); do not touch main.

---

### Task 1: `reconcile_gaps` helper + wire into `run_recomb`

**Files:** Modify `src/tessera/recomb/coverage.py` (add helper), `src/tessera/recomb/run.py` (lines 293-297, use it). Test: `tests/unit/test_coverage_reconcile.py` (create).

**Interfaces:**
- Consumes: `Region` (`recomb/regions.py`, fields incl. `msa_start`, `msa_end`, `donor_undercovered: bool = False`); `gaps_as_regions` (already used in run.py).
- Produces: `reconcile_gaps(regions: list[Region], gap_regions: list[Region]) -> list[Region]`.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_coverage_reconcile.py`:

```python
"""reconcile_gaps: a coverage gap caveats an overlapping called region; a
non-overlapping gap stays a donor-absent region."""

from __future__ import annotations

from tessera.recomb.coverage import reconcile_gaps
from tessera.recomb.regions import Region


def _region(msa_start: int, msa_end: int, *, donor_absent: bool = False) -> Region:
    return Region(
        minor_parent="m", major_parent="M", msa_start=msa_start, msa_end=msa_end,
        query_start=msa_start, query_end=msa_end, length_bp=msa_end - msa_start,
        n_windows=1, mean_sim_minor=0.98, mean_sim_major=0.93, margin=0.05,
        donor_absent=donor_absent,
    )


def test_overlapping_gap_caveats_region_and_is_not_emitted():
    called = _region(5000, 9000)
    gap = _region(5200, 5800, donor_absent=True)
    absent = reconcile_gaps([called], [gap])
    assert called.donor_undercovered is True   # the region is caveated
    assert absent == []                         # the gap is not double-emitted


def test_non_overlapping_gap_stays_donor_absent():
    called = _region(1000, 4000)
    gap = _region(5200, 5800, donor_absent=True)
    absent = reconcile_gaps([called], [gap])
    assert called.donor_undercovered is False
    assert absent == [gap]                       # genuine donor-absent region survives


def test_gap_overlapping_several_regions_caveats_all():
    a, b = _region(5000, 5500), _region(5400, 6000)
    absent = reconcile_gaps([a, b], [_region(5300, 5600, donor_absent=True)])
    assert a.donor_undercovered and b.donor_undercovered and absent == []
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/unit/test_coverage_reconcile.py -v`
Expected: FAIL with `ImportError: cannot import name 'reconcile_gaps'`.

- [ ] **Step 3: Add `reconcile_gaps` to `coverage.py`**

Add near `gaps_as_regions` in `src/tessera/recomb/coverage.py` (import `Region` if not already imported there):

```python
def reconcile_gaps(regions: list[Region], gap_regions: list[Region]) -> list[Region]:
    """Reconcile coverage-gap regions with the called regions.

    A gap that overlaps a called region caveats that region (``donor_undercovered``):
    the true donor for part of the region may be absent, so the reported donor is the
    closest available reference, not necessarily the real one. A gap that overlaps no
    called region is a genuine donor-absent region and is returned to be added to the
    output.
    """
    absent: list[Region] = []
    for gap in gap_regions:
        overlapping = [
            p for p in regions if gap.msa_start < p.msa_end and p.msa_start < gap.msa_end
        ]
        if overlapping:
            for p in overlapping:
                p.donor_undercovered = True
        else:
            absent.append(gap)
    return absent
```

- [ ] **Step 4: Use it in `run_recomb`**

In `src/tessera/recomb/run.py`, replace the reconciliation block (lines 293-297):

```python
    absent = [
        r for r in gaps_as_regions(coverage_gaps, result, major_parent)
        if not any(r.msa_start < p.msa_end and p.msa_start < r.msa_end for p in regions)
    ]
    regions = sorted(regions + absent, key=lambda r: r.msa_start)
```

with:

```python
    absent = reconcile_gaps(
        regions, gaps_as_regions(coverage_gaps, result, major_parent)
    )
    regions = sorted(regions + absent, key=lambda r: r.msa_start)
```

Add `reconcile_gaps` to the existing `from .coverage import (...)` block in run.py (alongside `gaps_as_regions`). The comment above the block should note the new caveat behaviour ("a gap overlapping a confident region caveats it as donor_undercovered; a gap covering no region is reported as donor-absent").

- [ ] **Step 5: Run the new tests + the existing recomb tests**

Run: `python -m pytest tests/unit/test_coverage_reconcile.py tests/unit -k "coverage or region or recomb or run" -v`
Expected: PASS (the new 3 + existing coverage/region tests unchanged).

- [ ] **Step 6: Ruff**

Run: `ruff check src/tessera/recomb/coverage.py src/tessera/recomb/run.py tests/unit/test_coverage_reconcile.py`
Expected: clean.

- [ ] **Step 7: Commit**

```bash
git add src/tessera/recomb/coverage.py src/tessera/recomb/run.py tests/unit/test_coverage_reconcile.py
git commit -m "$(cat <<'EOF'
Caveat a confident region when a coverage gap overlaps it

reconcile_gaps replaces the drop-on-overlap rule: a coverage gap that overlaps a
called region now flags that region donor_undercovered (the true donor for part
of the span may be absent; the reported donor is the closest available) instead
of being discarded. A gap covering no called region is still reported as a
donor-absent region.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Harness scorer accepts the caveat

**Files:** Modify `validation/run_hybrids.py` (`_score_panel_donor_absent`). Test: extend `tests/unit/test_harness_scoring.py`.

**Interfaces:**
- Consumes: the regions TSV `donor_undercovered` column (already written by `report_text.py`).
- Produces: updated `_score_panel_donor_absent` verdict.

- [ ] **Step 1: Write the failing tests**

Add to `tests/unit/test_harness_scoring.py`:

```python
def test_donor_absent_pass_on_undercovered_caveat(tmp_path):
    # No donor_absent region, but the span-overlapping region is caveated
    # donor_undercovered -> the tool acknowledged the missing donor -> PASS.
    _write_regions(tmp_path, [
        {"minor_parent": "gC", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "hmm,3seq", "donor_absent": "no",
         "donor_undercovered": "yes"}])
    setup = _setup(out=tmp_path, case_type="panel_donor_absent", clade_a="A", clade_b="B")
    res = rh._score_regions(tmp_path, {"gA": "A", "gC": "C"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is True


def test_donor_absent_fail_on_clean_confident_misattribution(tmp_path):
    # Confident, cross-clade, NOT undercovered -> a silent mis-attribution -> FAIL.
    _write_regions(tmp_path, [
        {"minor_parent": "gC", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "hmm,3seq", "donor_absent": "no",
         "donor_undercovered": "no"}])
    setup = _setup(out=tmp_path, case_type="panel_donor_absent", clade_a="A", clade_b="B")
    res = rh._score_regions(tmp_path, {"gA": "A", "gC": "C"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is False
```

Note: `_write_regions` in this file writes a fixed column set. Add `donor_undercovered` to its column list (default `"no"`) so these rows serialize; existing tests that omit it get `"no"`, unchanged.

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/unit/test_harness_scoring.py -k donor_absent -v`
Expected: the undercovered-caveat test FAILs (current scorer ignores `donor_undercovered`).

- [ ] **Step 3: Update `_score_panel_donor_absent`**

In `validation/run_hybrids.py`, change the verdict so a caveat counts as the honest signal and an undercovered region is not treated as a silent mis-attribution:

```python
    absent_hit = any(r.get("donor_absent") == "yes" and overlaps(r) for r in regions)
    present = [r for r in regions if r.get("donor_absent") != "yes"]

    def cross_clade(r):
        top = base_clade(clade_of(r["minor_parent"]) or "").split(".")[0]
        return top != setup.clade_a.split(".")[0]

    # A confident, cross-clade region that is NOT caveated is a silent mis-attribution.
    misattr = any(
        overlaps(r) and "," in (r.get("methods") or "")
        and r.get("donor_undercovered") != "yes" and cross_clade(r)
        for r in present
    )
    # The caller acknowledged the missing donor if the span is flagged donor-absent or a
    # span-overlapping region is caveated donor_undercovered.
    caveated = any(overlaps(r) and r.get("donor_undercovered") == "yes" for r in present)
    passed = (absent_hit or caveated) and not misattr
    return _base(setup, mode, runtime, present, n_refs=n_refs,
                 absent_hit=absent_hit, caveated=caveated, misattributed=misattr,
                 **{"pass": passed})
```

- [ ] **Step 4: Run the tests + ruff**

Run: `python -m pytest tests/unit/test_harness_scoring.py -v && ruff check validation/run_hybrids.py tests/unit/test_harness_scoring.py`
Expected: PASS, ruff clean.

- [ ] **Step 5: Commit**

```bash
git add validation/run_hybrids.py tests/unit/test_harness_scoring.py
git commit -m "$(cat <<'EOF'
Accept the donor_undercovered caveat in the donor-absent scorer

_score_panel_donor_absent now PASSes when the span is flagged donor-absent OR a
span-overlapping region is caveated donor_undercovered (the caller acknowledged
the missing donor), and treats only a clean, uncaveated confident cross-clade
region as a silent mis-attribution.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Harness gate + record

**Files:** Modify `validation/attribution-results.md`.

- [ ] **Step 1: Confirm the fix on the case (opt-in, aligner env)**

Run:
```bash
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
python validation/run_hybrids.py donorabsent_rsv
```
Expected: `donorabsent_rsv` now **PASS** (the span-overlapping region is caveated `donor_undercovered`, so the caller no longer silently mis-attributes). If it still FAILs, inspect `recombination_regions.tsv` for the `donor_undercovered` column and the region's overlap with the true span before adjusting.

- [ ] **Step 2: Full run -- confirm no regression**

Run: `python validation/run_hybrids.py`
Expected: the 18 positives + `lowdiv_rsv` + `neg_measles` hold; `donorabsent_rsv` flips to PASS -> sensitivity 20/20, specificity 1/1. Record the before/after (16/... -> ...) faithfully; a regression on any positive is a finding to report, not hide.

- [ ] **Step 3: Record**

Update `validation/attribution-results.md`: the donor-absent finding is now closed by the caveat behaviour (the caller keeps the closest-available attribution but flags `donor_undercovered` so the missing-donor signal is not lost). State the measured pass count and that no positive regressed.

- [ ] **Step 4: Commit**

```bash
git add validation/attribution-results.md
git commit -m "$(cat <<'EOF'
Record: donor-absent finding closed by the undercovered caveat

Harness re-run after the reconcile_gaps change: donorabsent_rsv now PASSes (the
donor region is kept but flagged donor_undercovered), with no regression to the
positives. Records the measured pass count.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Self-Review

**Spec coverage:** the caveat-on-overlap reconciliation (`reconcile_gaps` + run.py) -> Task 1; the harness scorer accepting the caveat -> Task 2; the harness gate + faithful record -> Task 3. The `donor_undercovered` field is reused (already serialized), so no report/format change.

**Placeholder scan:** every code step shows complete code; Task 3 is a measurement + specified write-up (a genuine FAIL-on-regression would be reported, not a TODO).

**Type consistency:** `reconcile_gaps(regions: list[Region], gap_regions: list[Region]) -> list[Region]` matches its call in run.py (`gaps_as_regions(...)` returns `list[Region]`) and mutates `Region.donor_undercovered` (a real `bool` field). The harness reads the TSV column `donor_undercovered` (values `"yes"`/`"no"`) written by `report_text.py`.

## Execution Handoff

Two options: **1. Subagent-Driven (recommended)** or **2. Inline Execution**. Which approach?
