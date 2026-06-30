# Donor Re-attribution Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an opt-in post-merge step that re-labels a called recombinant region's donor to the clade whose denoised consensus best matches the query over that region, closing the `rsv_a` attribution FAIL without touching detection or the backbone.

**Architecture:** A new pure module `recomb/reattribute.py` builds one majority-consensus base array per clade from the typed panel MSA and, per region, re-attributes the donor to the best non-backbone clade when it beats the current donor by a similarity margin. `run_recomb` calls it once after `consensus_regions(...)`, only when the flag is on and the panel is typed. Surfaced as `--reattribute-donors` (default off); validated on the `--compare` harness before any default change.

**Tech Stack:** Python 3.11+, numpy (no new dependency); Typer CLI; pytest.

## Global Constraints

- Dependency-free: no new runtime dependency; reimplement in numpy (reuse `CANONICAL_BASES`, the `barcode`/`consensus` majority-base primitive).
- Attribution fixes are **post-hoc only** — never change what competes for detection or the backbone (the reverted-competition lesson). This step changes only `Region.minor_parent` and `Region.mean_sim_minor`.
- Opt-in this cycle (default off). Promotion to default-when-typed is a separate later cycle gated on a clean harness.
- Modest scientific language; faithful reporting.
- Commit messages end with the project's co-author trailer. Branch before committing on the default branch.

---

### Task 1: `recomb/reattribute.py` — the re-attribution function

**Files:**
- Create: `src/tessera/recomb/reattribute.py`
- Test: `tests/unit/test_reattribute.py`

**Interfaces:**
- Consumes: `Region` (`recomb/regions.py`, dataclass; required fields `minor_parent, major_parent, msa_start, msa_end, query_start, query_end, length_bp, n_windows, mean_sim_minor, mean_sim_major, margin`); `CANONICAL_BASES = frozenset(b"ACGT")` (`recomb/similarity.py`); `lineage_of(label, lineage_map) -> str | None` (`recomb/typing.py`); `WindowSimilarity` (only `.rows: dict[str, np.ndarray]` and `.query: str` are used).
- Produces: `reattribute_donors(regions, result, lineage_map, major_clade, *, margin=0.03, min_sites=20, logger=None) -> list[Region]`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/unit/test_reattribute.py
"""Pure donor re-attribution (no aligner/network)."""

from __future__ import annotations

import numpy as np

from tessera.recomb.regions import Region
from tessera.recomb.reattribute import reattribute_donors


class _Result:
    """Minimal stand-in for WindowSimilarity: reattribute_donors uses only rows/query."""

    def __init__(self, rows: dict[str, np.ndarray], query: str):
        self.rows = rows
        self.query = query


def _enc(seq: str) -> np.ndarray:
    return np.frombuffer(seq.encode(), dtype=np.uint8)


def _region(minor: str, major: str, lo: int, hi: int) -> Region:
    return Region(minor_parent=minor, major_parent=major, msa_start=lo, msa_end=hi,
                  query_start=lo, query_end=hi, length_bp=hi - lo, n_windows=1,
                  mean_sim_minor=0.5, mean_sim_major=0.9, margin=0.0)


# Width-20 panel. Over span [5:15) clade C matches the query fully; the current donor
# (clade B) matches 6/10; clade A is the backbone (matches everywhere).
_QUERY = "A" * 20
_A = "A" * 20                                  # backbone member
_B = "AAAAA" + "TTTT" + "A" * 11               # 4 mismatches in the span -> 6/10
_C = "A" * 20                                  # full match in the span


def _panel() -> tuple[_Result, dict[str, str]]:
    rows = {"q": _enc(_QUERY), "a1": _enc(_A), "b1": _enc(_B), "c1": _enc(_C)}
    lineage_map = {"a1": "A", "b1": "B", "c1": "C"}
    return _Result(rows, "q"), lineage_map


def test_reattributes_to_better_clade_above_margin():
    result, lm = _panel()
    regions = [_region("b1", "a1", 5, 15)]
    out = reattribute_donors(regions, result, lm, "A", margin=0.1, min_sites=4)
    assert out[0].minor_parent == "c1"               # re-labelled to clade C's rep
    assert out[0].major_parent == "a1"               # backbone untouched
    assert out[0].mean_sim_minor == 1.0


def test_unchanged_below_margin():
    result, lm = _panel()
    regions = [_region("b1", "a1", 5, 15)]
    out = reattribute_donors(regions, result, lm, "A", margin=0.5, min_sites=4)
    assert out[0].minor_parent == "b1"               # 0.4 gain < 0.5 margin


def test_never_reattributes_to_the_backbone_clade():
    # Only clades A (backbone) and B; A scores best in the span but is excluded, so the
    # current donor B is kept rather than re-labelled to A.
    rows = {"q": _enc(_QUERY), "a1": _enc(_A), "b1": _enc(_B)}
    lm = {"a1": "A", "b1": "B"}
    out = reattribute_donors([_region("b1", "a1", 5, 15)], _Result(rows, "q"), lm,
                             "A", margin=0.1, min_sites=4)
    assert out[0].minor_parent == "b1"


def test_noop_without_lineage_map():
    result, _ = _panel()
    regions = [_region("b1", "a1", 5, 15)]
    assert reattribute_donors(regions, result, None, "A") is regions


def test_noop_with_fewer_than_two_clades():
    rows = {"q": _enc(_QUERY), "a1": _enc(_A)}
    out = reattribute_donors([_region("a1", "a1", 5, 15)], _Result(rows, "q"),
                             {"a1": "A"}, "A", margin=0.1, min_sites=4)
    assert out[0].minor_parent == "a1"


def test_skips_span_with_too_few_comparable_sites():
    result, lm = _panel()
    regions = [_region("b1", "a1", 5, 8)]            # 3 columns < min_sites
    out = reattribute_donors(regions, result, lm, "A", margin=0.1, min_sites=4)
    assert out[0].minor_parent == "b1"
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_reattribute.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'tessera.recomb.reattribute'`.

- [ ] **Step 3: Write the implementation**

```python
# src/tessera/recomb/reattribute.py
"""Post-hoc donor re-attribution.

The other callers attribute a recombinant region's donor by competing individual genomes;
at low divergence a near-identical adjacent-clade genome can win by chance. This step
refines a *called* region's donor -- after detection and the backbone are fixed -- to the
clade whose denoised consensus best matches the query over that region, guarded by a
similarity margin. It never changes detection or the backbone (attribution fixes are
post-hoc), and is a no-op on an untyped panel.
"""

from __future__ import annotations

import logging
from dataclasses import replace

import numpy as np

from .regions import Region
from .similarity import CANONICAL_BASES, WindowSimilarity
from .typing import LineageMap, lineage_of

_BASES = np.array(sorted(CANONICAL_BASES), dtype=np.uint8)  # A C G T byte codes


def _clade_consensus(
    rows: dict[str, np.ndarray], query: str, lineage_map: LineageMap
) -> tuple[dict[str, np.ndarray], dict[str, str]]:
    """Per-clade majority canonical base per column (0 where no canonical base), and a
    representative member label per clade."""
    by_clade: dict[str, list[str]] = {}
    for label in rows:
        if label == query:
            continue
        clade = lineage_of(label, lineage_map)
        if clade:
            by_clade.setdefault(clade, []).append(label)
    cons: dict[str, np.ndarray] = {}
    rep: dict[str, str] = {}
    for clade, members in by_clade.items():
        stack = np.array([rows[m] for m in members])  # (n_members, width)
        counts = np.stack([(stack == b).sum(axis=0) for b in _BASES])  # (4, width)
        majority = _BASES[np.argmax(counts, axis=0)]
        cons[clade] = np.where(counts.sum(axis=0) > 0, majority, 0).astype(np.uint8)
        rep[clade] = members[0]
    return cons, rep


def reattribute_donors(
    regions: list[Region],
    result: WindowSimilarity,
    lineage_map: LineageMap | None,
    major_clade: str | None,
    *,
    margin: float = 0.03,
    min_sites: int = 20,
    logger: logging.Logger | None = None,
) -> list[Region]:
    """Refine each region's donor to the best-matching non-backbone clade consensus.

    For each region, score every clade's consensus against the query over the region's MSA
    columns (fraction of matching canonical positions). Re-label the donor to the best
    non-``major_clade`` clade when it beats the region's current donor clade by ``margin``;
    otherwise leave the region unchanged. A no-op without a lineage map or with fewer than
    two scorable clades. Returns a new list; input regions are not mutated.
    """
    if not lineage_map:
        return regions
    cons, rep = _clade_consensus(result.rows, result.query, lineage_map)
    if len(cons) < 2:
        return regions
    query = result.rows[result.query]
    canon = np.isin(query, _BASES)

    def score(clade: str, lo: int, hi: int) -> float | None:
        c = cons[clade]
        comparable = canon[lo:hi] & (c[lo:hi] != 0)
        n = int(comparable.sum())
        if n < min_sites:
            return None
        match = int((query[lo:hi][comparable] == c[lo:hi][comparable]).sum())
        return match / n

    out: list[Region] = []
    for region in regions:
        lo, hi = region.msa_start, region.msa_end
        current = lineage_of(region.minor_parent, lineage_map)
        cur = score(current, lo, hi) if current in cons else None
        best, best_score = None, -1.0
        for clade in cons:
            if clade == major_clade:
                continue
            s = score(clade, lo, hi)
            if s is not None and s > best_score:
                best, best_score = clade, s
        if best is not None and best != current and best_score - (cur or 0.0) >= margin:
            if logger is not None:
                logger.info(
                    "Re-attributed donor %s -> %s over query %d-%d (consensus sim %.3f "
                    "vs %.3f).", region.minor_parent, rep[best], region.query_start,
                    region.query_end, best_score, cur or 0.0,
                )
            out.append(replace(region, minor_parent=rep[best],
                               mean_sim_minor=round(best_score, 4)))
        else:
            out.append(region)
    return out
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_reattribute.py -q`
Expected: PASS (6 passed).

- [ ] **Step 5: Lint and commit**

```bash
ruff check src/tessera/recomb/reattribute.py tests/unit/test_reattribute.py
git add src/tessera/recomb/reattribute.py tests/unit/test_reattribute.py
git commit -m "Add post-hoc donor re-attribution function

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Wire re-attribution into `run_recomb`

**Files:**
- Modify: `src/tessera/recomb/run.py` (`RecombParams`, and the merge block at `run.py:242-247`)
- Test: covered by Task 1 unit tests + Task 4 harness; no new unit test (this is integration into a binary-dependent pipeline).

**Interfaces:**
- Consumes: `reattribute_donors(...)` from Task 1; existing `regions`, `result`, `lineage_map`, `major_parent` already in `run_recomb` scope.
- Produces: `RecombParams.reattribute_donors: bool = False`, `RecombParams.reattribute_margin: float = 0.03`.

- [ ] **Step 1: Add the two fields to `RecombParams`**

In `src/tessera/recomb/run.py`, in the `@dataclass class RecombParams` (after the `organism` field, around `run.py:87`), add:

```python
    # Opt-in post-hoc donor re-attribution: re-label a region's donor to the clade whose
    # consensus best matches the query over that region (typed panels only). Off by
    # default; the backbone and detection are untouched. See recomb/reattribute.py.
    reattribute_donors: bool = False
    reattribute_margin: float = 0.03
```

- [ ] **Step 2: Import the function**

In `src/tessera/recomb/run.py`, beside the existing `from .ensemble import consensus_regions, reconcile_major` (`run.py:21`), add:

```python
from .reattribute import reattribute_donors
```

- [ ] **Step 3: Call it after the merge**

In `src/tessera/recomb/run.py`, immediately after the `consensus_regions(...)` call that assigns `regions, method_breakdown` (`run.py:243-247`), insert:

```python
    if params.reattribute_donors and lineage_map:
        regions = reattribute_donors(
            regions, result, lineage_map, major_parent,
            margin=params.reattribute_margin, logger=logger,
        )
```

- [ ] **Step 4: Verify the suite still passes**

Run: `pytest -m "not requires_binary" -q`
Expected: PASS (existing count + the 6 from Task 1; no regressions).

- [ ] **Step 5: Lint and commit**

```bash
ruff check src/tessera/recomb/run.py
git add src/tessera/recomb/run.py
git commit -m "Call donor re-attribution after the ensemble merge (opt-in)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: CLI flags on recomb / detect / fill-references

**Files:**
- Modify: `src/tessera/cli/cmd_recomb.py` (option + `RecombParams(...)` construction)
- Modify: `src/tessera/cli/cmd_detect.py`, `src/tessera/cli/cmd_fill_references.py` (option), `src/tessera/discover/iterate.py` (`FillParams` field + thread into the `RecombParams` it builds)
- Test: none new (thin Typer wiring; covered by the harness in Task 4).

**Interfaces:**
- Consumes: `RecombParams.reattribute_donors`, `RecombParams.reattribute_margin` (Task 2); `FillParams` (`discover/iterate.py`).
- Produces: `--reattribute-donors/--no-reattribute-donors` and `--reattribute-margin` on the three commands.

- [ ] **Step 1: Add the options to `cmd_recomb.py`**

In `src/tessera/cli/cmd_recomb.py`, add two options to the `recomb(...)` signature (after the existing `phi_window` option):

```python
    reattribute_donors: bool = typer.Option(
        False, "--reattribute-donors/--no-reattribute-donors",
        help="Refine each region's donor to the clade whose consensus best matches the "
        "query over that region (needs a typed panel; backbone unchanged). Off by default.",
    ),
    reattribute_margin: float = typer.Option(
        0.03, "--reattribute-margin",
        help="A donor is re-attributed only when a clade beats the current donor's "
        "consensus similarity by at least this much.",
    ),
```

- [ ] **Step 2: Pass them into `RecombParams` in `cmd_recomb.py`**

In the `RecombParams(...)` construction in `src/tessera/cli/cmd_recomb.py`, add:

```python
            reattribute_donors=reattribute_donors,
            reattribute_margin=reattribute_margin,
```

- [ ] **Step 3: Add the `FillParams` field**

In `src/tessera/discover/iterate.py`, in `@dataclass class FillParams` (beside `lineage_map`), add:

```python
    reattribute_donors: bool = False  # opt-in donor re-attribution in the detection step
```

In the same file, where `FillParams` builds its `RecombParams` for the detection run (`run_recomb(RecombParams(... methods=params.methods, lineage_map=lineage_map or None, organism=...)`), add:

```python
                reattribute_donors=params.reattribute_donors,
```

Also add `reattribute_donors=reattribute_donors` to the `FillParams.for_detection(...)` classmethod's constructor call and signature (so `detect` forwards it).

- [ ] **Step 4: Add the option to `cmd_detect.py` and `cmd_fill_references.py`**

In both `src/tessera/cli/cmd_detect.py` and `src/tessera/cli/cmd_fill_references.py`, add the option to the command signature:

```python
    reattribute_donors: bool = typer.Option(
        False, "--reattribute-donors/--no-reattribute-donors",
        help="Refine each region's donor to the best-matching clade consensus "
        "(typed panel; backbone unchanged). Off by default.",
    ),
```

and pass `reattribute_donors=reattribute_donors` into the `FillParams.for_detection(...)` call (detect) and the `FillParams(...)` construction (fill-references).

- [ ] **Step 5: Verify and commit**

Run: `pytest -m "not requires_binary" -q && ruff check src tests`
Expected: PASS; clean.

```bash
git add src/tessera/cli/cmd_recomb.py src/tessera/cli/cmd_detect.py \
        src/tessera/cli/cmd_fill_references.py src/tessera/discover/iterate.py
git commit -m "Add --reattribute-donors to recomb/detect/fill-references

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Harness `--compare` config + docs (the acceptance gate)

**Files:**
- Modify: `validation/run_hybrids.py` (`COMPARE_CONFIGS`, `_build_and_score` to pass the flag)
- Modify: `validation/README.md` (document the `reattribute` config)
- Modify: `docs/gap-analysis.md` and `validation/attribution-results.md` (record the outcome once measured)

**Interfaces:**
- Consumes: `RecombParams.reattribute_donors` (Task 2); the existing `_build_and_score` / `COMPARE_CONFIGS` machinery.
- Produces: a `reattribute` column in `--compare` and `attribution_compare.tsv`.

- [ ] **Step 1: Let `_build_and_score` set the flag**

In `validation/run_hybrids.py`, give `_build_and_score` an optional parameter and pass it into `RecombParams`:

```python
def _build_and_score(
    setup: CaseSetup, panel_mode: str, methods: tuple[str, ...], out_dir: Path,
    logger: logging.Logger, *, reattribute: bool = False,
) -> dict:
    ...
    windowing = run_recomb(RecombParams(msa=msa, output=out_dir, query=setup.query_label,
                                        window_size=setup.window, window_step=setup.step,
                                        organism=setup.name, methods=methods,
                                        lineage_map=lineage_map,
                                        reattribute_donors=reattribute), logger)
```

- [ ] **Step 2: Add a config tuple shape + the `reattribute` config**

In `validation/run_hybrids.py`, extend each `COMPARE_CONFIGS` entry to carry a `reattribute` flag and add the new config. Change the tuples to `(label, panel_mode, methods, reattribute)` and update the two loops that unpack them (`compare_case` and `_run_compare`) accordingly:

```python
COMPARE_CONFIGS = [
    ("baseline", "tip", ("hmm", "3seq", "maxchi", "bootscan"), False),
    ("+barcode", "tip", ("hmm", "3seq", "maxchi", "bootscan", "barcode"), False),
    ("consensus", "consensus", ("hmm", "3seq", "maxchi", "bootscan"), False),
    ("consensus+barcode", "consensus", ("hmm", "3seq", "maxchi", "bootscan", "barcode"), False),
    ("mixed", "mixed", ("hmm", "3seq", "maxchi", "bootscan"), False),
    # post-hoc donor re-attribution on the tip panel (backbone unchanged)
    ("reattribute", "tip", ("hmm", "3seq", "maxchi", "bootscan"), True),
]
```

In `compare_case`, unpack and pass the flag:

```python
    for label, mode, names, reattribute in COMPARE_CONFIGS:
        methods = parse_methods(",".join(names))
        sub = setup.out / ("cmp_" + re.sub(r"[^\w]+", "_", label))
        try:
            records[label] = _build_and_score(setup, mode, methods, sub, logger,
                                               reattribute=reattribute)
        except CaseSkipped:
            ...
```

In `_run_compare`, update the two `for label, _, _ in COMPARE_CONFIGS:` loops to `for label, _, _, _ in COMPARE_CONFIGS:`.

- [ ] **Step 3: Verify the harness still imports and the unit suite passes**

Run: `pytest -m "not requires_binary" -q && ruff check validation/run_hybrids.py`
Expected: PASS; clean. (The harness end-to-end run is opt-in and not part of CI.)

- [ ] **Step 4: Document the config in `validation/README.md`**

In the `--compare` subsection of `validation/README.md`, add a line noting the `reattribute` column: post-hoc donor re-attribution on the tip panel, expected to recover `rsv_a`'s donor (backbone unchanged) where the consensus panel could not.

- [ ] **Step 5: Commit**

```bash
git add validation/run_hybrids.py validation/README.md
git commit -m "Add the reattribute config to the --compare attribution harness

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 6: Run the gate (needs the aligner env; not CI)**

```bash
export PATH="$HOME/miniforge3/envs/recomfi-aln/bin:$PATH"
python validation/run_hybrids.py --compare
```

Expected acceptance: the `reattribute` column flips `rsv_a` to PASS (donor `exact`, backbone still `exact`) with **zero regressions** on the 16 baseline PASSes and no new donor `mismatch`. Record the measured outcome in `validation/attribution-results.md` and, if clean, note in `docs/gap-analysis.md` that G1 has a working caller-level fix (promotion to default-when-typed is a separate later cycle). If a regression appears, raise `--reattribute-margin` and re-measure before promoting.

---

## Notes for the implementer

- The harness end-to-end run (Task 4, Step 6) needs `mafft`/`skani` on PATH and contacts the Nextclade server; it is **not** part of CI. Everything else (Tasks 1-4 Steps 1-5) is CI-verifiable.
- Do not change detection, the backbone (`reconcile_major`), support, q-values, breakpoints, or method agreement anywhere. This feature only re-labels `minor_parent` and recomputes `mean_sim_minor`.
- `flu_h3n2_ha` is out of scope (a short-segment backbone mismatch, a different gap).
