# Harder harness Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the sensitivity-only 18/18 hybrid harness into a precision/recall harness by adding a `case_type` discriminator and five must-pass hard case types (negatives, low-divergence, donor-absent, equidistant), scored by branch, with a sensitivity/specificity summary split.

**Architecture:** Extend `validation/run_hybrids.py` in place. Add `case_type` (and a `tier` key, used in Phase 3) to each `HYBRIDS` dict, default so all 24 existing cases are unchanged. Refactor `_score_regions` to dispatch by `case_type`; extract today's logic verbatim into `_score_single_insert`. Add per-type query construction in `_prepare_case`, per-type scoring functions, and a specificity line in `_run_default`.

**Tech Stack:** Python 3, pytest, ruff. Aligner env `recomfi-aln` for the opt-in run. Design spec: `docs/superpowers/specs/2026-07-01-harder-harness-phase1-design.md`.

## Global Constraints

- All 24 existing `HYBRIDS` cases default to `case_type="single_insert"` / `tier="must_pass"` and must behave byte-identically (the existing 18/18 must hold).
- No data-layer or pipeline changes; harness-only. No new runtime dependency.
- Negatives score `passed = not-called`; low_div requires backbone at same top-level clade (no <4% free pass); donor-absent must hit a `donor_absent` region; equidistant must attribute to B not C.
- Unsuitable data -> `CaseSkipped` (never a FAIL).
- Reuse `pick_parents`, `make_hybrid`, `reconstruct_gapped`, `clade_match`/`donor_match`/`attribution_tier`/`shared_clade_depth`, `coverage`'s `donor_absent` regions. Modest language; ruff clean (line length 100); coverage floor 80%.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. Branch `harder-harness-phase1` (checked out); do not touch main.

---

### Task 1: `case_type` plumbing + scoring dispatch (behavior-preserving)

**Files:** Modify `validation/run_hybrids.py` (`CaseSetup` :448-469; `_prepare_case` return :531-535; `_score_regions` :643-686). Test: `tests/unit/test_harness_scoring.py` (create).

**Interfaces:**
- Produces: `CaseSetup.case_type: str`; `_score_single_insert(out_dir, clade_of, setup, n_refs, mode, runtime) -> dict` (today's `_score_regions` body); `_score_regions` dispatches on `setup.case_type`.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_harness_scoring.py`:

```python
"""Unit tests for run_hybrids scoring branches (no aligner; synthetic regions TSV)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "run_hybrids",
    Path(__file__).resolve().parents[2] / "validation" / "run_hybrids.py",
)
rh = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(rh)


def _setup(**kw):
    d = dict(name="t", out=Path("/tmp"), clade_a="A", clade_b="B", divergence=10.0,
             query=Path("/tmp/q"), query_label="q", q_start=100, q_end=200,
             window=300, step=30, sel_window=1500, aligner="mafft", reference="",
             tips={}, pool=[], members_by_clade={}, case_type="single_insert")
    d.update(kw)
    return rh.CaseSetup(**d)


def _write_regions(tmp_path, rows):
    """rows: list of dicts with query_start/query_end/major_parent/minor_parent/methods/
    donor_absent. Writes recombination_regions.tsv + an empty profile."""
    cols = ["minor_parent", "major_parent", "query_start", "query_end", "methods",
            "donor_absent"]
    lines = ["\t".join(cols)]
    for r in rows:
        lines.append("\t".join(str(r.get(c, "")) for c in cols))
    (tmp_path / "recombination_regions.tsv").write_text("\n".join(lines) + "\n")
    (tmp_path / "recombination_profile.tsv").write_text("phi_p\trmin\n-\t0\n")


def test_single_insert_pass(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "hmm,3seq", "donor_absent": "no"},
    ])
    clade_of = {"gA": "A", "gB": "B"}.get
    setup = _setup(out=tmp_path)
    res = rh._score_regions(tmp_path, clade_of, setup, 5, "tip", 1.0)
    assert res["pass"] is True and res["backbone_ok"] and res["donor_ok"]
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/unit/test_harness_scoring.py -v`
Expected: FAIL (`CaseSetup.__init__` has no `case_type`).

- [ ] **Step 3: Add `case_type` to `CaseSetup`**

In `validation/run_hybrids.py`, add to `CaseSetup` (after `members_by_clade`):

```python
    case_type: str = "single_insert"
```

- [ ] **Step 4: Populate it in `_prepare_case`**

In the `_prepare_case` return (`CaseSetup(...)` at :531), add:

```python
        case_type=case.get("case_type", "single_insert"),
```

- [ ] **Step 5: Extract `_score_single_insert` and dispatch**

Rename the current `def _score_regions(` body to `def _score_single_insert(` (identical body). Then add a new dispatcher above it:

```python
def _score_regions(
    out_dir: Path, clade_of, setup: CaseSetup, n_refs: int, mode: str, runtime: float,
) -> dict:
    """Score a completed run by the case's type (default: single-insert hybrid)."""
    scorer = {
        "single_insert": _score_single_insert,
    }.get(setup.case_type, _score_single_insert)
    return scorer(out_dir, clade_of, setup, n_refs, mode, runtime)
```

(Later tasks register more scorers in this dict.)

- [ ] **Step 6: Run the test + the existing harness-related tests + ruff**

Run: `python -m pytest tests/unit/test_harness_scoring.py -v && ruff check validation/run_hybrids.py tests/unit/test_harness_scoring.py`
Expected: PASS, ruff clean.

- [ ] **Step 7: Commit**

```bash
git add validation/run_hybrids.py tests/unit/test_harness_scoring.py
git commit -m "$(cat <<'EOF'
Add case_type plumbing and scoring dispatch to the hybrid harness

CaseSetup gains case_type (default single_insert); _score_regions dispatches to
per-type scorers, with today's logic extracted verbatim into
_score_single_insert. Behaviour-preserving scaffolding for the hard case types.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Negative controls (`neg_pure`, `neg_within`) + specificity summary

**Files:** Modify `validation/run_hybrids.py` (`_prepare_case` query build; new scorers + `_score_regions` dispatch; `_run_default` summary :778-782). Test: extend `tests/unit/test_harness_scoring.py`.

**Interfaces:**
- Consumes: `make_hybrid`, `reconstruct_gapped`, `pick_parents`, `clade_of`.
- Produces: `_pick_within_clade(tips, reference, logger) -> (clade, src_a, src_b)`; `_score_neg_pure`/`_score_neg_within` (registered in the dispatch dict); a `specificity` line in `_run_default`.

- [ ] **Step 1: Write the failing tests**

Add to `tests/unit/test_harness_scoring.py`:

```python
def test_neg_pure_pass_when_no_regions(tmp_path):
    _write_regions(tmp_path, [])  # no calls
    setup = _setup(out=tmp_path, case_type="neg_pure", clade_b="", q_start=0, q_end=0)
    res = rh._score_regions(tmp_path, (lambda x: "A"), setup, 5, "tip", 1.0)
    assert res["pass"] is True and res["n_false_regions"] == 0


def test_neg_pure_fail_on_false_call(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 10,
         "query_end": 90, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="neg_pure", clade_b="", q_start=0, q_end=0)
    res = rh._score_regions(tmp_path, {"gA": "A", "gB": "B"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is False and res["n_false_regions"] == 1


def test_neg_within_passes_same_clade_call(tmp_path):
    # A region whose minor and major are the SAME top-level clade is not a cross-clade FP.
    _write_regions(tmp_path, [
        {"minor_parent": "gA2", "major_parent": "gA1", "query_start": 10,
         "query_end": 90, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="neg_within", clade_a="A", clade_b="A")
    res = rh._score_regions(tmp_path, {"gA1": "A.1", "gA2": "A.2"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is True


def test_neg_within_fails_cross_clade_call(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA1", "query_start": 10,
         "query_end": 90, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="neg_within", clade_a="A", clade_b="A")
    res = rh._score_regions(tmp_path, {"gA1": "A.1", "gB": "B"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is False
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/unit/test_harness_scoring.py -k neg -v`
Expected: FAIL (no `neg_pure`/`neg_within` scorers; `n_false_regions` absent).

- [ ] **Step 3: Add the two scorers and register them**

Add these functions near `_score_single_insert`, and add both to the dispatch dict in `_score_regions` (`"neg_pure": _score_neg_pure, "neg_within": _score_neg_within`):

```python
def _base(setup: CaseSetup, mode: str, runtime: float, present: list, **extra) -> dict:
    """The reporting fields shared by every scorer."""
    row = {
        "name": setup.name, "clade_a": setup.clade_a, "clade_b": setup.clade_b,
        "divergence": setup.divergence, "mode": mode, "runtime": runtime,
        "detected": len(present) >= 1, "n_regions": len(present),
        "case_type": setup.case_type,
    }
    row.update(extra)
    return row


def _score_neg_pure(out_dir, clade_of, setup, n_refs, mode, runtime) -> dict:
    """Non-recombinant query: PASS iff no (donor-present) region is called."""
    regions = parse_regions(out_dir / "recombination_regions.tsv")
    present = [r for r in regions if r.get("donor_absent") != "yes"]
    return _base(setup, mode, runtime, present, n_refs=n_refs,
                 n_false_regions=len(present), pass=(len(present) == 0))


def _score_neg_within(out_dir, clade_of, setup, n_refs, mode, runtime) -> dict:
    """Within-clade splice: PASS iff no region attributes a CROSS-top-level-clade donor."""
    regions = parse_regions(out_dir / "recombination_regions.tsv")
    present = [r for r in regions if r.get("donor_absent") != "yes"]
    cross = [r for r in present
             if base_clade(clade_of(r["minor_parent"])).split(".")[0]
             != base_clade(clade_of(r["major_parent"])).split(".")[0]]
    return _base(setup, mode, runtime, present, n_refs=n_refs,
                 n_false_regions=len(cross), pass=(len(cross) == 0))
```

(`base_clade` and `parse_regions` are existing module functions. Python keyword `pass` cannot be a dict key via `pass=`; use `**{"pass": ...}`:)

Replace the `pass=(...)` kwargs above with dict-literal returns to avoid the `pass` keyword:

```python
    row = _base(setup, mode, runtime, present, n_refs=n_refs, n_false_regions=len(present))
    row["pass"] = len(present) == 0
    return row
```

(and analogously for `_score_neg_within` with `len(cross)`). Fix `_base` not to receive `pass`.

- [ ] **Step 4: Query construction for negatives in `_prepare_case`**

Add `_pick_within_clade` near `pick_parents`:

```python
def _pick_within_clade(tips, reference, logger) -> tuple[str, str, str]:
    """Largest non-recombinant clade and its two most-divergent members (for neg_within)."""
    by_clade: dict[str, list[str]] = {}
    for acc, (clade, _muts) in tips.items():
        if clade and clade not in ("?", "NA") and not is_recombinant_clade(clade):
            by_clade.setdefault(clade, []).append(acc)
    ranked = sorted(by_clade.items(), key=lambda kv: len(kv[1]), reverse=True)
    for clade, members in ranked:
        if len(members) < 2:
            continue
        seqs = {a: reconstruct_gapped(reference, tips[a][1]) for a in members}
        best = max(((100.0 - pct_identity(seqs[a], seqs[b]), a, b)
                    for i, a in enumerate(members) for b in members[i + 1:]),
                   default=None)
        if best:
            return clade, best[1], best[2]
    raise CaseSkipped("no clade has two divergent members for a within-clade splice")
```

In `_prepare_case`, branch the query/pool construction on `case.get("case_type")` *before* the current pair-splice block (around :488). For `neg_pure`: pick clade-A via `pick_parents` (ignore B), write the unspliced clade-A central genome as the query, `q_start=q_end=0`, `clade_b=""`. For `neg_within`: use `_pick_within_clade` to get `(clade, src_a, src_b)`, set `clade_a=clade_b=clade`, splice with `make_hybrid`. Both then build the source-removed pool exactly as today (dropping their two sources). Keep the existing block for `single_insert`. (Implement as an `if/elif/else` that assigns `clade_a, clade_b, src_a, src_b, query_seq, q_start, q_end` before the shared pool build.)

- [ ] **Step 5: Specificity split in `_run_default`**

After the existing `passed`/`ran` counters (:778-782), add a negative/positive split and print a specificity line. Positive = `case_type in {"single_insert","low_div",...}`; negative = `case_type in {"neg_pure","neg_within"}`. Compute over `results` (each scored dict carries `case_type`; skipped/errored excluded):

```python
    scored = [r for r in results if r.get("pass") is not None]
    negs = [r for r in scored if r.get("case_type") in ("neg_pure", "neg_within")]
    pos = [r for r in scored if r not in negs]
    sens_pass = sum(1 for r in pos if r["pass"])
    spec_pass = sum(1 for r in negs if r["pass"])
    false_calls = sum(r.get("n_false_regions", 0) for r in negs)
    print(f"\nsensitivity {sens_pass}/{len(pos)}  "
          f"specificity {spec_pass}/{len(negs)} ({false_calls} false call(s))  "
          f"({skipped} skipped, {errored} error)")
```

Keep the existing `passed/ran` line too (or replace it with the above -- the reviewer decides; the design requires both numbers visible).

- [ ] **Step 6: Run tests + ruff**

Run: `python -m pytest tests/unit/test_harness_scoring.py -v && ruff check validation/run_hybrids.py tests/unit/test_harness_scoring.py`
Expected: PASS, ruff clean.

- [ ] **Step 7: Commit** (message: "Add negative-control case types and specificity summary to the harness"; trailer.)

---

### Task 3: Low-divergence attribution (`low_div`)

**Files:** Modify `validation/run_hybrids.py` (`pick_parents` :372-377; `_prepare_case` :488-491; new `_score_low_div`). Test: extend `tests/unit/test_harness_scoring.py`.

**Interfaces:**
- Produces: `pick_parents(..., objective="max")` param; `divergence_band` handling; `_score_low_div` (registered).

- [ ] **Step 1: Write the failing tests**

```python
def test_low_div_requires_backbone_top_level(tmp_path):
    # donor correct top-level (B), backbone WRONG top-level (C) -> FAIL even at low div.
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gC", "query_start": 120,
         "query_end": 180, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="low_div", divergence=2.0,
                   clade_a="A", clade_b="B")
    res = rh._score_regions(tmp_path, {"gB": "B.1", "gC": "C.1"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is False  # backbone C != A, not relaxed away


def test_low_div_pass_top_level(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="low_div", divergence=2.0,
                   clade_a="A", clade_b="B")
    res = rh._score_regions(tmp_path, {"gA": "A.2", "gB": "B.1"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is True
```

- [ ] **Step 2: Run to verify failure** (`-k low_div`; FAIL, no scorer).

- [ ] **Step 3: Add `_score_low_div` and register it**

```python
def _score_low_div(out_dir, clade_of, setup, n_refs, mode, runtime) -> dict:
    """Low-divergence attribution: require detection + donor top-level + backbone
    top-level clade (no <4% free pass). Exact-vs-sibling is reported, not gating."""
    regions = parse_regions(out_dir / "recombination_regions.tsv")
    present = [r for r in regions if r.get("donor_absent") != "yes"]
    major_clade = clade_of(regions[0]["major_parent"]) if regions else "?"
    span = [r for r in present
            if int(r["query_start"]) <= setup.q_end and int(r["query_end"]) >= setup.q_start]
    donor_hits = [r for r in span
                  if attribution_tier(clade_of(r["minor_parent"]), setup.clade_b) != "mismatch"]
    detected = len(present) >= 1
    backbone_ok = attribution_tier(major_clade, setup.clade_a) != "mismatch"
    passed = detected and len(donor_hits) >= 1 and backbone_ok
    donor_obs = max((clade_of(r["minor_parent"]) for r in span),
                    key=lambda c: shared_clade_depth(c, setup.clade_b), default="?")
    return _base(setup, mode, runtime, present, n_refs=n_refs, major_clade=major_clade,
                 backbone_ok=backbone_ok, donor_ok=len(donor_hits) >= 1,
                 backbone_tier=attribution_tier(major_clade, setup.clade_a),
                 donor_tier=attribution_tier(donor_obs, setup.clade_b),
                 **{"pass": passed})
```

(`_base` must accept `**extra` and NOT define `pass`; `**{"pass": passed}` sets it.)

- [ ] **Step 4: `pick_parents` objective + `divergence_band`**

In `pick_parents`, add `objective: str = "max"` param; when `"min"`, select the least-divergent pair with `div >= case min_divergence`. Thread `case.get("pair_objective", "max")` from `_prepare_case`'s call. After computing `divergence` (:488), if `case.get("divergence_band")` is set and `divergence` is outside `[lo, hi]`, `raise CaseSkipped(...)`.

- [ ] **Step 5: Run tests + ruff.** Expected PASS, clean.

- [ ] **Step 6: Commit** ("Add low-divergence attribution case type with an honest top-level floor"; trailer).

---

### Task 4: Panel-adversarial (`panel_donor_absent`, `panel_equidistant`)

**Files:** Modify `validation/run_hybrids.py` (`_prepare_case` pool build :518-528; new scorers). Test: extend `tests/unit/test_harness_scoring.py`.

- [ ] **Step 1: Write the failing tests**

```python
def test_donor_absent_pass_on_gap(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "coverage", "donor_absent": "yes"}])
    setup = _setup(out=tmp_path, case_type="panel_donor_absent", clade_a="A", clade_b="B")
    res = rh._score_regions(tmp_path, {"gA": "A"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is True


def test_donor_absent_fail_on_misattribution(tmp_path):
    # No donor_absent region; a confident cross-clade present region mis-attributes the span.
    _write_regions(tmp_path, [
        {"minor_parent": "gC", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "hmm,3seq", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="panel_donor_absent", clade_a="A", clade_b="B")
    res = rh._score_regions(tmp_path, {"gA": "A", "gC": "C"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is False


def test_equidistant_pass_when_B_wins(tmp_path):
    _write_regions(tmp_path, [
        {"minor_parent": "gB", "major_parent": "gA", "query_start": 120,
         "query_end": 180, "methods": "hmm", "donor_absent": "no"}])
    setup = _setup(out=tmp_path, case_type="panel_equidistant", clade_a="A", clade_b="B")
    res = rh._score_regions(tmp_path, {"gA": "A", "gB": "B", "gC": "C"}.get, setup, 5, "tip", 1.0)
    assert res["pass"] is True
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Add the two scorers + register**

```python
def _score_panel_donor_absent(out_dir, clade_of, setup, n_refs, mode, runtime) -> dict:
    """True donor removed from the panel: PASS iff the span is flagged donor-absent and
    not confidently mis-attributed to a present cross-clade donor."""
    regions = parse_regions(out_dir / "recombination_regions.tsv")
    def overlaps(r):
        return int(r["query_start"]) <= setup.q_end and int(r["query_end"]) >= setup.q_start
    absent_hit = any(r.get("donor_absent") == "yes" and overlaps(r) for r in regions)
    present = [r for r in regions if r.get("donor_absent") != "yes"]
    misattr = any(overlaps(r) and "," in (r.get("methods") or "")
                  and base_clade(clade_of(r["minor_parent"])).split(".")[0] != setup.clade_a.split(".")[0]
                  for r in present)
    passed = absent_hit and not misattr
    return _base(setup, mode, runtime, present, n_refs=n_refs,
                 absent_hit=absent_hit, misattributed=misattr, **{"pass": passed})


def _score_panel_equidistant(out_dir, clade_of, setup, n_refs, mode, runtime) -> dict:
    """Two equidistant candidate donors (B pinned as truth): PASS iff donor attributed to
    B, not the decoy C. Guards the plurality_major tie-break."""
    res = _score_single_insert(out_dir, clade_of, setup, n_refs, mode, runtime)
    # _score_single_insert already requires donor_match(observed, clade_b, clade_a); a C
    # attribution fails donor_ok. Reuse its verdict directly.
    return res
```

- [ ] **Step 4: Pool manipulation in `_prepare_case`**

For `panel_donor_absent`: after building `pool`/`members_by_clade`, drop every genome/member whose clade `clade_match(clade, clade_b)`, and `CaseSkipped` if that empties the pool. For `panel_equidistant`: the case pins `clades: [A, B]` and a `decoy_clade` C via a new `case` field; assert both B and C remain in `members_by_clade` (else `CaseSkipped`); construction is otherwise the standard single splice from B.

- [ ] **Step 5: Run tests + ruff. Step 6: Commit** ("Add donor-absent and equidistant panel-adversarial case types"; trailer).

---

### Task 5: Register concrete cases, docs, and the harness measurement

**Files:** Modify `validation/run_hybrids.py` (`HYBRIDS` list); `validation/README.md`; `validation/attribution-results.md`.

- [ ] **Step 1: Add case entries** to `HYBRIDS` -- a small, curated set reusing existing datasets, e.g.:
  - `{"name": "neg_measles", "dataset": "nextstrain/measles/genome/WHO-2012", "case_type": "neg_pure"}`
  - `{"name": "neg_within_hiv", "dataset": "community/neherlab/hiv-1/hxb2", "case_type": "neg_within"}`
  - `{"name": "lowdiv_dengue", "dataset": "nextstrain/dengue/all", "case_type": "low_div", "pair_objective": "min", "divergence_band": [1.0, 4.0], "min_divergence": 1.0}`
  - `{"name": "donorabsent_rsv", "dataset": "nextstrain/rsv/a/EPI_ISL_412866", "case_type": "panel_donor_absent"}`
  - `{"name": "equidist_dengue", "dataset": "nextstrain/dengue/all", "case_type": "panel_equidistant", "clades": ["<A>", "<B>"], "decoy_clade": "<C>"}` (fill clades after a probe run).

- [ ] **Step 2: Run the harness (opt-in, aligner env), tune, and record**

```bash
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
python validation/run_hybrids.py            # default run: sensitivity + specificity lines
```
Confirm the existing 18/18 positives hold and the new must-pass hard cases behave (a genuine FAIL is a finding to record, not hide). Pin `clades`/`decoy_clade` for the equidistant case from the logged clade list. Record the outcome (sensitivity + specificity + per-new-case verdicts) in `validation/attribution-results.md`, and document the new case types + the precision/recall framing in `validation/README.md`.

- [ ] **Step 3: Commit** ("Register Phase-1 hard cases; record harness precision/recall result"; trailer).

---

## Self-Review

**Spec coverage:** case_type/tier plumbing + dispatch (Task 1); neg_pure/neg_within + specificity split (Task 2); low_div + honest floor (Task 3); donor-absent + equidistant (Task 4); concrete cases + docs + measurement (Task 5). Existing 18/18 preserved by the `single_insert` default (Task 1) and verified in Task 5.

**Placeholder scan:** code steps carry complete code; the `pass`-keyword pitfall is called out explicitly (use `**{"pass": ...}`); Task 5's clade pins are intentionally filled from a probe run (documented), not a TODO.

**Type consistency:** every scorer has signature `(out_dir, clade_of, setup, n_refs, mode, runtime) -> dict` and is registered in the `_score_regions` dispatch dict; all return a dict carrying `case_type`, `pass`, and (for negatives) `n_false_regions`, which `_run_default` reads. `CaseSetup.case_type` is set in `_prepare_case` and read in dispatch.

## Execution Handoff

Two options: **1. Subagent-Driven (recommended)** or **2. Inline Execution**. Which approach?
