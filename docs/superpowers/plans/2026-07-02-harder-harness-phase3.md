# Harder harness Phase 3 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add frontier cases (inter-species, reassortment) behind a `--frontier` flag, scored XFAIL/XPASS/KNOWN-LIMIT in a separate block that can never count as a must-pass regression.

**Architecture:** Extend `validation/run_hybrids.py`. A `tier` gate hides frontier cases unless `--frontier`; a `_run_frontier` runner prints a separate table and always returns 0. Inter-species splices two fully-reconstructed genomes by fraction (`make_cross_hybrid`) over a combined two-species panel; reassortment concatenates per-clade segment consensuses. A `_score_frontier` measures the identity envelope / junction.

**Tech Stack:** Python 3, pytest, ruff. Aligner env `recomfi-aln` for the opt-in gate. Design spec: `docs/superpowers/specs/2026-07-02-harder-harness-phase3-design.md`.

## Global Constraints

- Faithful reporting is load-bearing: without `--frontier`, only `tier == "must_pass"` runs and the headline is unchanged; frontier verdicts are XPASS/XFAIL/KNOWN-LIMIT/SKIP in a separate table and never touch `passed`/`errored` or the return code.
- Harness-only (no `src/` change); no new runtime dependency. Modest language; ruff clean (line length 100); coverage floor 80%.
- Cross-species tips carry mutations relative to *different* references -- do NOT use `make_hybrid` across datasets; reconstruct each genome then splice by fraction.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. Branch `harder-harness-phase3` (checked out); do not touch main.

---

### Task 1: `tier`/`--frontier` gate + frontier runner scaffolding

**Files:** Modify `validation/run_hybrids.py` (`main`, add `_run_frontier`). Test: `tests/unit/test_frontier_gate.py` (create).

**Interfaces:** Produces the `tier` case key (default `"must_pass"`); `main` filters by tier + `--frontier`; `_run_frontier(cases, logger) -> int` (always returns 0).

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_frontier_gate.py`:

```python
"""The tier/--frontier gate: frontier cases run only under --frontier."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "run_hybrids", Path(__file__).resolve().parents[2] / "validation" / "run_hybrids.py")
rh = importlib.util.module_from_spec(_SPEC)
sys.modules["run_hybrids"] = rh
_SPEC.loader.exec_module(rh)


def test_select_cases_hides_frontier_without_flag():
    cases = [{"name": "a", "tier": "must_pass"}, {"name": "f", "tier": "frontier"}]
    assert [c["name"] for c in rh._select_cases(cases, names=[], frontier=False)] == ["a"]
    assert [c["name"] for c in rh._select_cases(cases, names=[], frontier=True)] == ["f"]


def test_select_cases_default_tier_is_must_pass():
    cases = [{"name": "a"}]  # no tier -> must_pass
    assert [c["name"] for c in rh._select_cases(cases, names=[], frontier=False)] == ["a"]
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/unit/test_frontier_gate.py -v`
Expected: FAIL (`module 'run_hybrids' has no attribute '_select_cases'`).

- [ ] **Step 3: Add `_select_cases` and wire `--frontier` into `main`**

Add a helper and use it in `main`:

```python
def _select_cases(cases, *, names, frontier):
    """Filter HYBRIDS by name and by tier: frontier cases run only when `frontier` is set,
    must_pass cases only when it is not."""
    want = "frontier" if frontier else "must_pass"
    return [c for c in cases
            if (not names or c["name"] in names)
            and c.get("tier", "must_pass") == want]
```

In `main`, parse the flag and route:

```python
def main(argv: list[str]) -> int:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    logger = logging.getLogger("tessera")
    compare = "--compare" in argv
    frontier = "--frontier" in argv
    names = [a for a in argv if not a.startswith("-")]
    cases = _select_cases(HYBRIDS, names=names, frontier=frontier)
    if frontier:
        return _run_frontier(cases, logger)
    if compare:
        return _run_compare(cases, logger)
    return _run_default(cases, logger)
```

- [ ] **Step 4: Add `_run_frontier`**

```python
def _run_frontier(cases: list[dict], logger: logging.Logger) -> int:
    """Run the frontier (known-limitation) cases and print a separate XPASS/XFAIL/
    KNOWN-LIMIT/SKIP table. Always returns 0: a frontier outcome is never a regression."""
    print(f"\nTessera FRONTIER probes -- {len(cases)} case(s) (known limitations)\n" + "=" * 72)
    for case in cases:
        try:
            rec = _run_frontier_case(case, logger)
            print(f"[{rec['verdict']:11}] {case['name']:16} {rec.get('detail', '')}")
        except CaseSkipped as exc:
            print(f"[SKIP       ] {case['name']:16} {exc}")
        except Exception as exc:  # noqa: BLE001 - report, never fail the batch
            logger.exception("[%s] frontier ERROR", case["name"])
            print(f"[ERROR      ] {case['name']:16} {type(exc).__name__}: {exc}")
    print("\n(frontier probes measure the envelope; they never count as a regression)")
    return 0
```

`_run_frontier_case(case, logger) -> dict` (with a `verdict` and `detail`) is added in Tasks 2-3; for Task 1 add a stub that raises `CaseSkipped("no frontier case types yet")` so the runner is testable, then replace it.

- [ ] **Step 5: Run tests + ruff**

Run: `python -m pytest tests/unit/test_frontier_gate.py tests/unit/test_harness_scoring.py -v && ruff check validation/run_hybrids.py tests/unit/test_frontier_gate.py`
Expected: PASS (must-pass path unchanged), ruff clean.

- [ ] **Step 6: Commit** ("Add tier/--frontier gate and the frontier runner"; trailer).

---

### Task 2: `inter_species` frontier case

**Files:** Modify `validation/run_hybrids.py`. Test: extend `tests/unit/test_frontier_gate.py`.

**Interfaces:** Produces `make_cross_hybrid(genome_a, genome_b, insert=INSERT) -> (query, q_start, q_end)`; `_prepare_frontier_species` / the inter_species branch of `_run_frontier_case`; `_score_frontier_envelope(identity, detected, backbone_ok, donor_ok) -> (verdict, detail)`.

- [ ] **Step 1: Write the failing tests**

```python
def test_make_cross_hybrid_splices_by_fraction():
    a = "A" * 100
    b = "B" * 100
    q, s, e = rh.make_cross_hybrid(a, b, insert=(0.35, 0.65))
    assert q == "A" * 35 + "B" * 30 + "A" * 35 and (s, e) == (35, 65)


def test_envelope_known_limit_below_floor():
    v, _d = rh._score_frontier_envelope(0.74, detected=True, backbone_ok=False, donor_ok=False)
    assert v == "KNOWN-LIMIT"


def test_envelope_xpass_in_envelope_detected_attributed():
    v, _d = rh._score_frontier_envelope(0.83, detected=True, backbone_ok=True, donor_ok=True)
    assert v == "XPASS"


def test_envelope_xfail_in_envelope_not_attributed():
    v, _d = rh._score_frontier_envelope(0.83, detected=True, backbone_ok=False, donor_ok=False)
    assert v == "XFAIL"
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Add `make_cross_hybrid` and `_score_frontier_envelope`**

```python
def make_cross_hybrid(genome_a, genome_b, insert=INSERT):
    """Splice two already-reconstructed, ungapped genomes by length fraction (their
    coordinate systems differ, so a reference-relative splice is invalid). Returns
    (query, q_start, q_end) in query coordinates."""
    la, lb = len(genome_a), len(genome_b)
    a1, a2 = int(la * insert[0]), int(la * insert[1])
    b1, b2 = int(lb * insert[0]), int(lb * insert[1])
    query = genome_a[:a1] + genome_b[b1:b2] + genome_a[a2:]
    return query.upper(), a1, a1 + (b2 - b1)


def _score_frontier_envelope(identity, detected, backbone_ok, donor_ok):
    """Envelope verdict for an inter-species case. Below the HMM identity floor (0.80)
    attribution is not meaningful -> KNOWN-LIMIT; in-envelope -> XPASS when detected and
    attributed, else XFAIL. Returns (verdict, detail)."""
    if identity < 0.80:
        return "KNOWN-LIMIT", f"identity {identity:.2f} below the 0.80 HMM floor"
    if detected and backbone_ok and donor_ok:
        return "XPASS", f"identity {identity:.2f}, detected + attributed"
    return "XFAIL", f"identity {identity:.2f}, in-envelope but not attributed"
```

- [ ] **Step 4: Add the inter_species prepare + run branch**

Add `_prepare_frontier_species(case, logger) -> (setup, identity)`: resolve `case["dataset"]` and `case["second_dataset"]`; `build_pool` both; `collect_tips` both (prefix each clade with `"A:"`/`"B:"`); reconstruct the two source genomes (backbone clade's central genome from A, donor clade's from B) and `make_cross_hybrid`; build a combined source-removed pool + `members_by_clade` (prefixed clades). Compute `identity = pct_identity(genome_a, genome_b)/100`. Assemble a `CaseSetup` (reuse the existing dataclass; `true_spans=[(q_start,q_end,"B:"+donor_clade)]`).

In `_run_frontier_case`, dispatch on `case["case_type"]`: for `"inter_species"`, call the prepare, then `_build_and_score(setup, "tip", methods, out, logger)` to get the scored dict, then
`verdict, detail = _score_frontier_envelope(identity, rec["detected"], rec["backbone_ok"], rec["donor_ok"])`; return `{"verdict": verdict, "detail": detail}`.

(If `_build_and_score`'s representation invariant SKIPs the cross-species panel, let `CaseSkipped` propagate to the runner -- reported as SKIP.)

- [ ] **Step 5: Run tests + ruff; Step 6: Commit** ("Add inter_species frontier case with envelope scoring"; trailer).

---

### Task 3: `reassortant` frontier case (detection probe)

**Files:** Modify `validation/run_hybrids.py`. Test: extend `tests/unit/test_frontier_gate.py`.

- [ ] **Step 1: Write the failing test** (junction-overlap check):

```python
def test_reassortant_junction_hit(tmp_path):
    # A region overlapping the HA|NA junction at query 1500 -> XPASS (detection-gated).
    cols = ["minor_parent", "major_parent", "query_start", "query_end", "methods",
            "donor_absent", "donor_undercovered"]
    (tmp_path / "recombination_regions.tsv").write_text(
        "\t".join(cols) + "\n" + "\t".join(["m", "M", "1400", "1700", "hmm", "no", "no"]) + "\n")
    (tmp_path / "recombination_profile.tsv").write_text("phi_p\trmin\n-\t0\n")
    v, _d = rh._score_frontier_reassortant(tmp_path, junction=1500)
    assert v == "XPASS"


def test_reassortant_junction_miss(tmp_path):
    cols = ["minor_parent", "major_parent", "query_start", "query_end", "methods",
            "donor_absent", "donor_undercovered"]
    (tmp_path / "recombination_regions.tsv").write_text(
        "\t".join(cols) + "\n" + "\t".join(["m", "M", "100", "300", "hmm", "no", "no"]) + "\n")
    (tmp_path / "recombination_profile.tsv").write_text("phi_p\trmin\n-\t0\n")
    v, _d = rh._score_frontier_reassortant(tmp_path, junction=1500)
    assert v == "XFAIL"
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Add `_score_frontier_reassortant`**

```python
def _score_frontier_reassortant(out_dir, junction):
    """Detection-gated: XPASS iff a called region overlaps the HA|NA junction. The
    whole-segment-swap labeling is a documented KNOWN-LIMIT (not scored)."""
    regions = parse_regions(out_dir / "recombination_regions.tsv")
    present = [r for r in regions if r.get("donor_absent") != "yes"]
    hit = any(int(r["query_start"]) <= junction <= int(r["query_end"]) for r in present)
    if hit:
        return "XPASS", f"breakpoint called at the HA|NA junction (~{junction} bp)"
    return "XFAIL", f"no breakpoint at the HA|NA junction (~{junction} bp)"
```

- [ ] **Step 4: Add the reassortant prepare + run branch**

`_prepare_frontier_reassortant(case, logger)`: resolve the HA dataset and `case["second_dataset"]` (NA); for each segment build per-clade consensus genomes (`consensus_sequence` over `clade_members`); pick a backbone clade X and donor clade Y per segment (two most-divergent within that segment); the panel references are `HA_consensus[c] ++ NA_consensus[c]` written per clade (compound label); the query is `HA_consensus[X] ++ NA_consensus[Y]`; `junction = len(HA_consensus[X])`. Build a collection from the concatenated references, align, `run_recomb`, and return `(out_dir, junction)`. In `_run_frontier_case` dispatch `"reassortant"` -> this prepare + `_score_frontier_reassortant`. `CaseSkipped` when a segment has < 2 clades or a consensus is empty.

(This is a probe: if the concatenated-panel alignment/scan does not run, let `CaseSkipped`/the runner's `ERROR` path report it -- do not force a result.)

- [ ] **Step 5: Run tests + ruff; Step 6: Commit** ("Add reassortant frontier detection probe"; trailer).

---

### Task 4: Register frontier cases, gate, record

**Files:** `validation/run_hybrids.py` (HYBRIDS), `validation/attribution-results.md`, `validation/README.md`.

- [ ] **Step 1: Register the two frontier cases**

```python
    {"name": "interspecies_rsv", "tier": "frontier", "case_type": "inter_species",
     "dataset": "nextstrain/rsv/a/EPI_ISL_412866",
     "second_dataset": "nextstrain/rsv/b/EPI_ISL_1653999"},
    {"name": "reassort_flu", "tier": "frontier", "case_type": "reassortant",
     "dataset": "nextstrain/flu/h3n2/ha/EPI1857216",
     "second_dataset": "nextstrain/flu/h3n2/na/EPI1857215"},
```

- [ ] **Step 2: Gate (opt-in, aligner env)**

```bash
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
python validation/run_hybrids.py --frontier       # the two probes
python validation/run_hybrids.py                   # must-pass headline unchanged
```
Record whatever the probes measure -- the RSV-A x RSV-B identity + verdict (expected KNOWN-LIMIT or XFAIL near the 0.80 floor) and the flu reassortant junction outcome. A SKIP/ERROR from an assembly that will not run is an honest probe result; do not tune to force an XPASS.

- [ ] **Step 3: Record** in `attribution-results.md` (a Phase-3 section: the measured envelope, framed as known-limitations) and document `--frontier` in `README.md`.

- [ ] **Step 4: Commit** ("Register frontier cases; record the measured envelope"; trailer).

---

## Self-Review

**Spec coverage:** tier/--frontier gate + separate runner (Task 1); inter_species cross-splice + envelope scoring (Task 2); reassortant concatenation + junction probe (Task 3); registration + gate + record (Task 4). Faithful reporting: `_run_frontier` returns 0 and prints a separate table; must-pass path untouched.

**Placeholder scan:** complete code for the gate, `make_cross_hybrid`, and both `_score_frontier_*` functions; the two-dataset prepare functions are described precisely (they are I/O-heavy and exercised by the opt-in gate); Task 4 is a measurement + specified write-up.

**Type consistency:** `_select_cases(cases, *, names, frontier)`; `make_cross_hybrid(...) -> (str, int, int)`; `_score_frontier_envelope(...) -> (verdict, detail)` and `_score_frontier_reassortant(...) -> (verdict, detail)` both return a 2-tuple consumed by `_run_frontier_case`; frontier verdicts are strings printed by `_run_frontier`, never folded into must-pass counters.

## Execution Handoff

Two options: **1. Subagent-Driven (recommended)** or **2. Inline Execution**. Which approach?
