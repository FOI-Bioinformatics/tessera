# Harder harness Phase 2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add hard-topology cases (multi-breakpoint, short tract, asymmetric, terminal) and a masking-sibling attribution case to the hybrid harness, via a `make_mosaic`/`true_spans` generalisation done additively so the Phase-1 cases and scorers are untouched.

**Architecture:** Extend `validation/run_hybrids.py`. Add pure `make_mosaic` and `pick_parents_n`; add `CaseSetup.true_spans` (default `[]`) populated for every case (single-span for the existing types); add a `mosaic` case type with a `pattern` knob and a `_score_mosaic` scorer, and a `mask_sibling` case type with a `_score_mask_sibling` scorer. `make_hybrid` and `q_start`/`q_end` stay; the Phase-1 scorers are unchanged.

**Tech Stack:** Python 3, pytest, ruff. Aligner env `recomfi-aln` for the opt-in gate. Design spec: `docs/superpowers/specs/2026-07-02-harder-harness-phase2-design.md`.

## Global Constraints

- Additive: the 24+ existing cases (all non-mosaic) behave byte-identically; `make_hybrid`, `q_start`, `q_end`, and every Phase-1 scorer are unchanged. `true_spans` is new data they don't have to use.
- `true_spans` holds only **non-backbone** donor spans as `(q_start, q_end, donor_clade)`.
- `AB_short` is **detection-gated only** (a miss is reported, never a FAIL). `mask_sibling` requires **exact** donor attribution.
- Harness-only; no product/pipeline/data-layer change; no new runtime dependency. Modest language; ruff clean (line length 100); coverage floor 80%.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. Branch `harder-harness-phase2` (checked out); do not touch main.

---

### Task 1: `make_mosaic` + `pick_parents_n` + `CaseSetup.true_spans`

**Files:** Modify `validation/run_hybrids.py` (add two functions; `CaseSetup` ~line 508; `_prepare_case` return ~line 641 and each construction branch). Test: `tests/unit/test_mosaic.py` (create).

**Interfaces:**
- Produces: `make_mosaic(reference, segments) -> (str, list[tuple[int,int,str]])` where `segments` is `list[(muts, frac_start, frac_end, donor_clade|None)]`; `pick_parents_n(tips, reference, n, *, floor) -> list[tuple[str,str]]`; `CaseSetup.true_spans: list[tuple[int,int,str]]`.
- Consumes: `reconstruct_gapped`, `pct_identity`, `clade_representative`, `MIN_MEMBERS`.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_mosaic.py`:

```python
"""make_mosaic span math and pick_parents_n selection (no aligner)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "run_hybrids", Path(__file__).resolve().parents[2] / "validation" / "run_hybrids.py")
rh = importlib.util.module_from_spec(_SPEC)
sys.modules["run_hybrids"] = rh
_SPEC.loader.exec_module(rh)


def test_make_mosaic_emits_non_backbone_spans_in_order():
    ref = "A" * 100
    # No mutations -> each segment reconstructs to the reference slice; boundaries at
    # 0/30/45/70/100. Non-backbone segments (clade set) become true_spans.
    segs = [([], 0.0, 0.30, None), ([], 0.30, 0.45, "B"),
            ([], 0.45, 0.70, None), ([], 0.70, 1.0, "C")]
    query, spans = rh.make_mosaic(ref, segs)
    assert len(query) == 100
    assert [c for _s, _e, c in spans] == ["B", "C"]
    assert spans[0] == (30, 45, "B") and spans[1] == (70, 100, "C")


def test_make_mosaic_single_insert_matches_make_hybrid_span():
    ref = "ACGT" * 25  # 100 bp
    segs = [([], 0.0, 0.35, None), ([], 0.35, 0.65, "B"), ([], 0.65, 1.0, None)]
    _q, spans = rh.make_mosaic(ref, segs)
    assert spans == [(35, 65, "B")]
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/unit/test_mosaic.py -v`
Expected: FAIL (`module 'run_hybrids' has no attribute 'make_mosaic'`).

- [ ] **Step 3: Add `make_mosaic` and `pick_parents_n`**

Add near `make_hybrid` in `validation/run_hybrids.py`:

```python
def make_mosaic(reference, segments):
    """Splice a mosaic from ``segments`` = list of (muts, frac_start, frac_end, donor_clade).
    A segment with ``donor_clade`` None is backbone (not a true span). Returns
    ``(query, true_spans)`` with ``true_spans`` = (q_start, q_end, donor_clade) in query
    coordinates for the non-backbone segments, in order."""
    length = len(reference)
    mosaic = ""
    true_spans = []
    for muts, f0, f1, clade in segments:
        b0, b1 = int(length * f0), int(length * f1)
        seg = reconstruct_gapped(reference, muts)[b0:b1]
        q0 = len(mosaic.replace("-", ""))
        mosaic += seg
        q1 = len(mosaic.replace("-", ""))
        if clade is not None:
            true_spans.append((q0, q1, clade))
    return mosaic.replace("-", "").upper(), true_spans


def pick_parents_n(tips, reference, n, *, floor=0.0):
    """Greedily pick the ``n`` most mutually-divergent clades (>= ``MIN_MEMBERS`` genomes,
    pairwise divergence >= ``floor``). Returns ``[(clade, source_genome), ...]`` with the
    largest clade first (backbone). Raises ``CaseSkipped`` when fewer than ``n`` qualify."""
    # Reuse the existing module-level clade_members (as pick_parents does): prefer
    # non-recombinant clades, allow recombinant ones only if too few remain.
    eligible = {c: m for c, m in clade_members(tips, exclude_recombinant=True).items()
                if len(m) >= MIN_MEMBERS}
    if len(eligible) < n:
        eligible = {c: m for c, m in clade_members(tips, exclude_recombinant=False).items()
                    if len(m) >= MIN_MEMBERS}
    if len(eligible) < n:
        raise CaseSkipped(f"fewer than {n} clades with >= {MIN_MEMBERS} genomes")
    reps = {c: clade_representative(m, tips) for c, m in eligible.items()}
    seq = {c: reconstruct_gapped(reference, tips[reps[c]][1]) for c in reps}

    def div(a, b):
        return 100.0 - pct_identity(seq[a], seq[b])

    clades = sorted(reps)
    first = max(((div(a, b), a, b) for i, a in enumerate(clades) for b in clades[i + 1:]),
                default=None)
    if first is None or first[0] < floor:
        raise CaseSkipped("no clade pair meets the divergence floor")
    chosen = [first[1], first[2]]
    while len(chosen) < n:
        rest = [c for c in clades if c not in chosen]
        nxt = max(rest, key=lambda c: min(div(c, x) for x in chosen), default=None)
        if nxt is None or min(div(nxt, x) for x in chosen) < floor:
            raise CaseSkipped(f"fewer than {n} clades meet the divergence floor")
        chosen.append(nxt)
    chosen.sort(key=lambda c: len(eligible[c]), reverse=True)  # backbone = largest
    return [(c, reps[c]) for c in chosen]
```

`clade_members(tips, *, exclude_recombinant)` is an existing module-level helper (used inside `pick_parents`'s `eligible_of` closure); `pick_parents_n` reuses it directly, so no refactor of `pick_parents` is needed.

- [ ] **Step 4: Add `CaseSetup.true_spans` and populate it**

Add to `CaseSetup` (after `q_end`): `true_spans: list[tuple[int, int, str]] = field(default_factory=list)`.

In `_prepare_case`, build `true_spans` per branch and pass it to the `CaseSetup(...)` return:
- `neg_pure` / `neg_within`: `true_spans = []`.
- `single_insert` / `low_div` / `panel_*` / `mask_sibling`: `true_spans = [(q_start, q_end, clade_b)]`.
- (mosaic sets it in Task 2.)
Add `true_spans=true_spans` to the `CaseSetup(...)` call. `q_start`/`q_end` stay exactly as today.

- [ ] **Step 5: Run tests + existing harness tests + ruff**

Run: `python -m pytest tests/unit/test_mosaic.py tests/unit/test_harness_scoring.py -v && ruff check validation/run_hybrids.py tests/unit/test_mosaic.py`
Expected: PASS (existing harness-scoring tests unaffected), ruff clean.

- [ ] **Step 6: Commit** ("Add make_mosaic, pick_parents_n, and CaseSetup.true_spans"; trailer).

---

### Task 2: `mosaic` case type + `_score_mosaic` + patterns

**Files:** Modify `validation/run_hybrids.py` (`_prepare_case` construction; new scorer + dispatch). Test: extend `tests/unit/test_harness_scoring.py`.

**Interfaces:** Consumes `make_mosaic`, `pick_parents_n`, `donor_match`, `clade_match`, `window_params`. Produces `_score_mosaic` (registered).

- [ ] **Step 1: Write the failing tests**

Add to `tests/unit/test_harness_scoring.py` (uses the existing `_setup`/`_write_regions`; note `_setup` needs a `true_spans` kwarg once CaseSetup has the field):

```python
def test_mosaic_pass_all_spans_recovered(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 300,
         "query_end": 450, "methods": "hmm", "donor_absent": "no"},
        {"minor_parent": "gC", "major_parent": "gA", "query_start": 700,
         "query_end": 1000, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="mosaic", clade_a="A",
                   true_spans=[(300, 450, "B"), (700, 1000, "C")], pattern="ABAC")
    res = rh._score_regions(tmp_path, {"gA": "A", "gB": "B", "gC": "C"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is True and res["spans_hit"] == 2


def test_mosaic_fail_when_a_span_missed(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 300,
         "query_end": 450, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="mosaic", clade_a="A",
                   true_spans=[(300, 450, "B"), (700, 1000, "C")], pattern="ABAC")
    res = rh._score_regions(tmp_path, {"gA": "A", "gB": "B", "gC": "C"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is False and res["spans_hit"] == 1


def test_mosaic_short_is_detection_gated(tmp_path):
    _write_regions(tmp_path, [  # a region exists but does not overlap the true span
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 10,
         "query_end": 40, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="mosaic", clade_a="A",
                   true_spans=[(500, 560, "B")], pattern="AB_short")
    res = rh._score_regions(tmp_path, {"gA": "A", "gB": "B"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is True  # detection-gated: a call exists, span recovery only reported
```

`_setup` must accept `pattern` and `true_spans`; add them (default `pattern=""`, `true_spans=None -> []`) and set `case.pattern` where the scorer reads it (store `pattern` on `CaseSetup` -- add a field `pattern: str = ""`).

- [ ] **Step 2: Run to verify failure** (`-k mosaic`; FAIL, no scorer/field).

- [ ] **Step 3: Add `pattern` to `CaseSetup`, `_score_mosaic`, register it**

Add `pattern: str = ""` to `CaseSetup` (populate from `case.get("pattern","")` in `_prepare_case`). Add and register `"mosaic": _score_mosaic`:

```python
def _score_mosaic(out_dir, clade_of, setup, n_refs, mode, runtime) -> dict:
    """Multi-span mosaic: each non-backbone span must be recovered with the right donor
    clade, backbone must be clade_a. AB_short is detection-gated (span recovery reported,
    not required); AB_terminal additionally requires the recovered region near query 0."""
    regions = parse_regions(out_dir / "recombination_regions.tsv")
    present = [r for r in regions if r.get("donor_absent") != "yes"]
    major_clade = clade_of(regions[0]["major_parent"]) if regions else "?"
    detected = len(present) >= 1
    backbone_ok = clade_match(major_clade, setup.clade_a)

    def hit(q0, q1, donor):
        return any(
            int(r["query_start"]) <= q1 and int(r["query_end"]) >= q0
            and donor_match(clade_of(r["minor_parent"]), donor, setup.clade_a)
            for r in present)

    spans_hit = sum(hit(q0, q1, d) for q0, q1, d in setup.true_spans)
    spans_total = len(setup.true_spans)
    if setup.pattern == "AB_short":
        passed = detected                      # a miss is reported, not failed
    elif setup.pattern == "AB_terminal":
        near0 = any(int(r["query_start"]) <= setup.window for r in present)
        passed = detected and backbone_ok and spans_hit == spans_total and near0
    else:                                       # ABAC, AB_9010
        passed = detected and backbone_ok and spans_hit == spans_total
    return _base(setup, mode, runtime, present, n_refs=n_refs, major_clade=major_clade,
                 backbone_ok=backbone_ok, spans_hit=spans_hit, spans_total=spans_total,
                 **{"pass": passed})
```

- [ ] **Step 4: Add the `mosaic` construction branch in `_prepare_case`**

Add an `elif case_type == "mosaic":` branch that reads `pattern = case["pattern"]`, chooses parents (`pick_parents_n` for `ABAC`, else `pick_parents`), builds `segments` per pattern (see the design's fractions), calls `make_mosaic`, and sets `query_seq` + `true_spans` + `clade_a`/`clade_b`/`divergence` (divergence = max pairwise among chosen). For `AB_short`, compute the donor fraction from `window_params(len(reference))` so the tract is sub-window. Set `q_start`/`q_end` to the envelope `min/max` of `true_spans`.

- [ ] **Step 5: Run tests + ruff; Step 6: Commit** ("Add mosaic case type with multi-breakpoint/short/terminal patterns"; trailer).

---

### Task 3: `mask_sibling` case type + `_score_mask_sibling`

**Files:** Modify `validation/run_hybrids.py`. Test: extend `tests/unit/test_harness_scoring.py`.

- [ ] **Step 1: Write the failing tests**

```python
def test_mask_sibling_pass_on_exact_donor(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="mask_sibling", clade_a="A.1", clade_b="A.D.1.8",
                   true_spans=[(100, 200, "A.D.1.8")])
    res = rh._score_regions(tmp_path, {"gA": "A.1", "gB": "A.D.1.8"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is True


def test_mask_sibling_fail_on_sibling_donor(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gS", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="mask_sibling", clade_a="A.1", clade_b="A.D.1.8",
                   true_spans=[(100, 200, "A.D.1.8")])
    res = rh._score_regions(tmp_path, {"gA": "A.1", "gS": "A.D.5.2"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is False  # a sibling sub-clade, not exact -> FAIL
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Add `_score_mask_sibling` + register + construction branch**

```python
def _score_mask_sibling(out_dir, clade_of, setup, n_refs, mode, runtime) -> dict:
    """Donor with a sibling sub-clade present: PASS requires the donor attributed EXACTLY
    (attribution_tier == 'exact'), not merely to a sibling. Regression guard for the
    plurality-major fix."""
    regions = parse_regions(out_dir / "recombination_regions.tsv")
    present = [r for r in regions if r.get("donor_absent") != "yes"]
    major_clade = clade_of(regions[0]["major_parent"]) if regions else "?"
    q0, q1, donor = setup.true_spans[0]

    def overlaps(r):
        return int(r["query_start"]) <= q1 and int(r["query_end"]) >= q0

    exact = any(overlaps(r) and attribution_tier(clade_of(r["minor_parent"]), donor) == "exact"
                for r in present)
    passed = len(present) >= 1 and clade_match(major_clade, setup.clade_a) and exact
    return _base(setup, mode, runtime, present, n_refs=n_refs, major_clade=major_clade,
                 donor_exact=exact, **{"pass": passed})
```

Construction: `elif case_type == "mask_sibling":` -- pin `clades: [A, B]`, standard splice from B (`make_hybrid`), then `CaseSkipped` unless a sibling sub-clade of B (same top-level, different sub-clade) is present in `members_by_clade`. `true_spans = [(q_start, q_end, clade_b)]`.

- [ ] **Step 4: Run tests + ruff; Step 5: Commit** ("Add mask_sibling case type requiring exact donor attribution"; trailer).

---

### Task 4: Register cases, harness gate, record

**Files:** `validation/run_hybrids.py` (HYBRIDS), `validation/attribution-results.md`, `validation/README.md`.

- [ ] **Step 1: Add case entries** reusing existing datasets, e.g.:
  - `{"name": "mosaic_dengue", "dataset": "nextstrain/dengue/all", "case_type": "mosaic", "pattern": "ABAC"}`
  - `{"name": "asym_measles", "dataset": "nextstrain/measles/genome/WHO-2012", "case_type": "mosaic", "pattern": "AB_9010"}`
  - `{"name": "short_wnv", "dataset": "nextstrain/wnv/all-lineages", "case_type": "mosaic", "pattern": "AB_short"}`
  - `{"name": "terminal_mumps", "dataset": "nextstrain/mumps/genome", "case_type": "mosaic", "pattern": "AB_terminal"}`
  - `{"name": "masksib_rsv", "dataset": "nextstrain/rsv/a/EPI_ISL_412866", "case_type": "mask_sibling", "clades": ["A.1", "A.D.1.8"]}`
  (Pins/patterns may be adjusted after a probe run; datasets with too few clades `CaseSkipped`.)
  Ensure `_run_default`'s per-row printer tolerates the mosaic/mask rows (they carry `major_clade`/`backbone_ok`; `spans_hit` is extra -- printed only if you add a column, otherwise ignored).

- [ ] **Step 2: Harness gate (opt-in, aligner env)**

```bash
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
python validation/run_hybrids.py
```
Confirm the Phase-1 set holds (single-insert byte-identical) and the new topology/mask cases score. A `AB_short` miss is reported, not a FAIL; a genuine FAIL on `ABAC`/`mask_sibling` is a real finding -- record it honestly (as with `donorabsent`), do not tune it away. Adjust taxa/pins only for *validity* (enough clades, in-band divergence), not to force green.

- [ ] **Step 3: Record** the outcome in `validation/attribution-results.md` and document the new patterns + `mask_sibling` in `validation/README.md` (the "Harder cases" subsection).

- [ ] **Step 4: Commit** ("Register Phase-2 topology + mask-sibling cases; record harness result"; trailer).

---

## Self-Review

**Spec coverage:** make_mosaic/true_spans + pick_parents_n + CaseSetup.true_spans (Task 1); mosaic patterns ABAC/9010/short/terminal + _score_mosaic (Task 2); mask_sibling + exact-donor scorer (Task 3); cases + gate + record (Task 4). Additive migration keeps Phase-1 scorers and single_insert byte-identical (Task 1 populates single-span true_spans; q_start/q_end unchanged).

**Placeholder scan:** complete code for the pure functions and scorers; the construction branches (Task 2 Step 4, Task 3 Step 3) are described with exact fractions/guards referencing the design; Task 4 is a measurement + specified write-up.

**Type consistency:** `make_mosaic(...) -> (str, list[tuple[int,int,str]])`; `CaseSetup.true_spans: list[tuple[int,int,str]]`; scorers loop it as `(q0,q1,donor)`; every scorer signature `(out_dir, clade_of, setup, n_refs, mode, runtime) -> dict` and is registered in `_score_regions`. `pick_parents_n(...) -> list[(clade, src)]` matches its mosaic-branch use.

## Execution Handoff

Two options: **1. Subagent-Driven (recommended)** or **2. Inline Execution**. Which approach?
