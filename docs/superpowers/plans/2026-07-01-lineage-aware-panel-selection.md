# Lineage-aware panel selection + recombinant toggle Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** When references carry lineage labels, build the panel by lineage -- one query-closest representative per non-recombinant lineage, recombinant (CRF/URF/X) lineages excluded by default with a toggle to keep them -- instead of clade-blind ANI dereplication.

**Architecture:** Add a `_lineage_select` reduction step inside `select_regional` (discover/pool.py) that fires when a `lineage_of` map is supplied; untyped genomes fall back to the existing ANI `dereplicate`. Thread the map and a `keep_recombinant` toggle through the product path (`_select_from` / `FillParams` / the `detect` + `fill-references` CLI) and the validation harness. With no map, behaviour is exactly today's global path.

**Tech Stack:** Python 3 (numpy/pandas/typer), skani/skDER external tools, pytest, ruff. Design spec: `docs/superpowers/specs/2026-07-01-lineage-aware-panel-selection-design.md`.

## Global Constraints

- No new runtime dependency; reuse `is_recombinant_lineage` (`recomb/typing.py`), `skani_query_ani`, `dereplicate` (`discover/panel.py`), `build_lineage_map` / `lineage_map_from_rows` / `titles_from_collection` (`recomb/typing.py`).
- Backward compatible: `select_regional(..., lineage_of=None)` must reproduce today's exact behaviour (same skani calls on the same genome set).
- Toggle default is **exclude** recombinant lineages: `--keep-recombinant-lineages/--no-keep-recombinant-lineages` default off; `FillParams.keep_recombinant: bool = False`.
- Use modest scientific language in code and docstrings. No Unicode in nextflow files (none touched here).
- Ruff clean (`ruff check src tests validation`, line length 100); CI coverage floor (80%) holds.
- Commit messages end with the trailer `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.
- Branch `lineage-select-spec` (already checked out off main); do not commit on main.

---

### Task 1: `_lineage_select` reduction step in `select_regional`

**Files:**
- Modify: `src/tessera/discover/pool.py` (import block ~line 33; `select_regional` line 199-235)
- Test: `tests/unit/test_lineage_select.py` (create)

**Interfaces:**
- Consumes: `dereplicate(genomes, *, ani, logger) -> (reps, redundant)` and `skani_query_ani(query, refs, logger) -> dict[Path, (ani, af)]` from `discover/panel.py`; `is_recombinant_lineage(name) -> bool` and `strip_sequence_extension(name) -> str`.
- Produces: `_lineage_select(genomes, lineage_of, query_ani, *, keep_recombinant, derep_ani, logger) -> list[Path]`; and `select_regional(..., lineage_of: dict[str, str] | None = None, keep_recombinant: bool = False)`.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_lineage_select.py`:

```python
"""Unit tests for the lineage-aware panel reduction in discover/pool.py.

``_lineage_select`` is exercised directly with a stubbed ANI map and a stubbed
``dereplicate`` so no skani/skDER binary is needed.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from tessera.discover import pool as pool_mod

_LOG = logging.getLogger("tessera.test")


@pytest.fixture
def stub_dereplicate(monkeypatch):
    """Stub ``dereplicate`` to keep the first genome of the input (a stand-in for
    ANI collapse) and record what it was called with."""
    calls = {}

    def fake(genomes, *, ani, logger):
        calls["genomes"] = list(genomes)
        calls["ani"] = ani
        return (genomes[:1], genomes[1:])

    monkeypatch.setattr(pool_mod, "dereplicate", fake)
    return calls


def _p(name: str) -> Path:
    return Path(f"/pool/{name}.fasta")


def test_keeps_one_query_closest_rep_per_lineage(stub_dereplicate):
    a1, a2, b1 = _p("a1"), _p("a2"), _p("b1")
    lineage_of = {"a1": "A.1", "a2": "A.1", "b1": "B.1"}
    query_ani = {a1: (95.0, 90.0), a2: (98.0, 90.0), b1: (80.0, 90.0)}
    kept = pool_mod._lineage_select(
        [a1, a2, b1], lineage_of, query_ani,
        keep_recombinant=False, derep_ani=99.0, logger=_LOG,
    )
    assert set(kept) == {a2, b1}  # a2 beats a1 within A.1; b1 is B.1's only member


def test_excludes_recombinant_lineages_by_default(stub_dereplicate):
    a1, x1 = _p("a1"), _p("x1")
    lineage_of = {"a1": "A.1", "x1": "CRF01_AE"}
    query_ani = {a1: (95.0, 90.0), x1: (99.0, 90.0)}
    kept = pool_mod._lineage_select(
        [a1, x1], lineage_of, query_ani,
        keep_recombinant=False, derep_ani=99.0, logger=_LOG,
    )
    assert kept == [a1]  # the CRF is dropped even though it is the closest match


def test_keeps_recombinant_lineages_when_toggled(stub_dereplicate):
    a1, x1 = _p("a1"), _p("x1")
    lineage_of = {"a1": "A.1", "x1": "CRF01_AE"}
    query_ani = {a1: (95.0, 90.0), x1: (99.0, 90.0)}
    kept = pool_mod._lineage_select(
        [a1, x1], lineage_of, query_ani,
        keep_recombinant=True, derep_ani=99.0, logger=_LOG,
    )
    assert set(kept) == {a1, x1}  # the CRF competes as its own lineage


def test_untyped_genomes_route_to_ani_fallback(stub_dereplicate):
    a1, u1, u2, u3 = _p("a1"), _p("u1"), _p("u2"), _p("u3")
    lineage_of = {"a1": "A.1"}  # u1..u3 untyped
    query_ani = {a1: (95.0, 90.0)}
    kept = pool_mod._lineage_select(
        [a1, u1, u2, u3], lineage_of, query_ani,
        keep_recombinant=False, derep_ani=99.0, logger=_LOG,
    )
    # a1 kept as its lineage rep; the 3 untyped go through dereplicate (stub keeps u1).
    assert set(kept) == {a1, u1}
    assert stub_dereplicate["genomes"] == [u1, u2, u3]
    assert stub_dereplicate["ani"] == 99.0


def test_empty_map_sends_everything_to_fallback(stub_dereplicate):
    u1, u2, u3 = _p("u1"), _p("u2"), _p("u3")
    kept = pool_mod._lineage_select(
        [u1, u2, u3], {}, {},
        keep_recombinant=False, derep_ani=99.0, logger=_LOG,
    )
    assert kept == [u1]  # all untyped -> the ANI fallback (stub keeps the first)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_lineage_select.py -v`
Expected: FAIL with `AttributeError: module 'tessera.discover.pool' has no attribute '_lineage_select'`.

- [ ] **Step 3: Add the import**

In `src/tessera/discover/pool.py`, add after the existing `from .panel import ...` line (line 33):

```python
from ..recomb.typing import is_recombinant_lineage
```

(Place it with the other `from ..` imports, keeping ruff's import order -- `..core` / `..recomb` groups sort before `.blast` / `.panel`. Move it above the `.blast` import if ruff `I` complains.)

- [ ] **Step 4: Add the `_lineage_select` helper**

In `src/tessera/discover/pool.py`, add this function immediately before `def select_regional` (line 199):

```python
def _lineage_select(
    genomes: list[Path],
    lineage_of: dict[str, str],
    query_ani: dict[Path, tuple[float, float]],
    *,
    keep_recombinant: bool,
    derep_ani: float,
    logger: logging.Logger,
) -> list[Path]:
    """Reduce a pool by lineage: keep one query-closest representative per lineage.

    Genomes with a lineage label are grouped by lineage and the member with the
    highest whole-genome ANI to the query is kept; recombinant (CRF/URF/X) lineages
    are dropped unless ``keep_recombinant`` is set, since they carry both parents'
    segments and mask the true parents. Genomes with no label fall back to ANI
    dereplication (the pre-lineage behaviour). Returns the union.
    """
    typed_by_lineage: dict[str, list[Path]] = {}
    untyped: list[Path] = []
    for g in genomes:
        lineage = lineage_of.get(strip_sequence_extension(g.name))
        if not lineage:
            untyped.append(g)
        elif keep_recombinant or not is_recombinant_lineage(lineage):
            typed_by_lineage.setdefault(lineage, []).append(g)
    reps = [
        max(members, key=lambda g: query_ani.get(g, (0.0, 0.0))[0])
        for members in typed_by_lineage.values()
    ]
    if len(untyped) > 2:
        kept_untyped, _ = dereplicate(untyped, ani=derep_ani, logger=logger)
    else:
        kept_untyped = untyped
    logger.info(
        "Lineage selection: %d representative(s) from %d lineage(s); "
        "%d untyped genome(s) -> %d after ANI dereplication.",
        len(reps), len(typed_by_lineage), len(untyped), len(kept_untyped),
    )
    return reps + kept_untyped
```

- [ ] **Step 5: Wire it into `select_regional`**

In `src/tessera/discover/pool.py`, add the two params to the `select_regional` signature. Change lines 207-211 from:

```python
    dereplicate_pool: bool = True,
    drop_siblings: bool = True,
    sibling_identity: float = DEFAULT_SIBLING_IDENTITY,
    sibling_coverage: float = DEFAULT_SIBLING_COVERAGE,
    logger: logging.Logger,
```

to:

```python
    dereplicate_pool: bool = True,
    lineage_of: dict[str, str] | None = None,
    keep_recombinant: bool = False,
    drop_siblings: bool = True,
    sibling_identity: float = DEFAULT_SIBLING_IDENTITY,
    sibling_coverage: float = DEFAULT_SIBLING_COVERAGE,
    logger: logging.Logger,
```

Then replace the reduction block (lines 230-237), which currently reads:

```python
    pool = genomes
    if dereplicate_pool and len(genomes) > 2:
        pool, _ = dereplicate(genomes, ani=derep_ani, logger=logger)
        logger.info(
            "Dereplicated pool: %d representative(s) from %d genome(s).", len(pool), len(genomes)
        )

    whole = skani_query_ani(query_fasta, pool, logger)
```

with:

```python
    if lineage_of:
        # Typed references: reduce by lineage, not by clade-blind ANI. Compute the
        # whole-genome query ANI on the full pool first, so a lineage's representative
        # is its query-closest member and the sibling/related filters below can reuse it.
        whole = skani_query_ani(query_fasta, genomes, logger)
        pool = _lineage_select(
            genomes, lineage_of, whole,
            keep_recombinant=keep_recombinant, derep_ani=derep_ani, logger=logger,
        )
    else:
        pool = genomes
        if dereplicate_pool and len(genomes) > 2:
            pool, _ = dereplicate(genomes, ani=derep_ani, logger=logger)
            logger.info(
                "Dereplicated pool: %d representative(s) from %d genome(s).",
                len(pool), len(genomes),
            )
        whole = skani_query_ani(query_fasta, pool, logger)
```

Also extend the `select_regional` docstring (line 219, before the closing `"""`) with a sentence:

```python
    When ``lineage_of`` labels the genomes, the pool is reduced by lineage -- one
    query-closest representative per lineage, recombinant lineages excluded unless
    ``keep_recombinant`` -- instead of by ANI; untyped genomes still dereplicate.
```

- [ ] **Step 6: Run the new tests to verify they pass**

Run: `pytest tests/unit/test_lineage_select.py -v`
Expected: PASS (5 passed).

- [ ] **Step 7: Run the existing pool tests to confirm no regression**

Run: `pytest tests/unit -k "pool or regional or select" -v`
Expected: PASS (the `lineage_of=None` path is unchanged).

- [ ] **Step 8: Ruff**

Run: `ruff check src/tessera/discover/pool.py tests/unit/test_lineage_select.py`
Expected: clean (no output).

- [ ] **Step 9: Commit**

```bash
git add src/tessera/discover/pool.py tests/unit/test_lineage_select.py
git commit -m "$(cat <<'EOF'
Add lineage-aware reduction to select_regional

When select_regional is given a lineage_of map it keeps one query-closest
representative per lineage and drops recombinant (CRF/URF/X) lineages unless
keep_recombinant is set; untyped genomes fall back to ANI dereplication. With
no map the existing global path is unchanged.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Product wiring -- `_select_from` and `FillParams`

**Files:**
- Modify: `src/tessera/discover/iterate.py` (`FillParams` dataclass ~line 60; `FillParams.for_detection` ~line 110-160; `_select_from` line 567-574)
- Test: `tests/unit/test_select_from_lineage.py` (create)

**Interfaces:**
- Consumes: `select_regional(..., lineage_of=..., keep_recombinant=...)` from Task 1; `build_lineage_map`, `lineage_map_from_rows`, `titles_from_collection` (already imported in `iterate.py`, lines 29-38).
- Produces: `FillParams.keep_recombinant: bool = False`, forwarded by `FillParams.for_detection`; `_select_from` now passes a `lineage_of` map derived from the pool headers.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_select_from_lineage.py`:

```python
"""_select_from should type the pool from its headers and forward lineage
selection into select_regional."""

from __future__ import annotations

import logging
from pathlib import Path

from tessera.discover import iterate as it
from tessera.discover.iterate import FillParams

_LOG = logging.getLogger("tessera.test")


def _write(dir_: Path, acc: str, note: str) -> Path:
    p = dir_ / f"{acc}.fasta"
    p.write_text(f">{acc} {note}\nACGT\n")
    return p


def test_select_from_passes_lineage_map(monkeypatch, tmp_path):
    pool = tmp_path / "pool"
    pool.mkdir()
    g1 = _write(pool, "ACC1", "genotype A")
    g2 = _write(pool, "ACC2", "genotype B")
    captured = {}

    def fake_select_regional(query, genomes, **kwargs):
        captured.update(kwargs)
        return it.PoolSelection(selected=list(genomes), table=[])

    monkeypatch.setattr("tessera.discover.pool.select_regional", fake_select_regional)

    params = FillParams(
        query=tmp_path / "q.fasta", collection=None, output=tmp_path / "out",
        taxon="Human immunodeficiency virus 1", keep_recombinant=False,
    )
    it._select_from(params, [g1, g2], _LOG)

    assert captured["keep_recombinant"] is False
    assert captured["lineage_of"]  # a non-empty map was built and forwarded
    assert set(captured["lineage_of"]) == {"ACC1", "ACC2"}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `pytest tests/unit/test_select_from_lineage.py -v`
Expected: FAIL with `TypeError: __init__() got an unexpected keyword argument 'keep_recombinant'`.

- [ ] **Step 3: Add the `FillParams` field**

In `src/tessera/discover/iterate.py`, add after `reattribute_donors: bool = False` (the last field of `FillParams`, ~line 106):

```python
    keep_recombinant: bool = False  # keep recombinant (CRF/URF/X) lineages in a typed panel
```

- [ ] **Step 4: Forward it through `for_detection`**

In `FillParams.for_detection`, add a parameter after `reattribute_donors: bool = False,` (~line 129):

```python
        keep_recombinant: bool = False,
```

and add to the returned `cls(...)` call, after `reattribute_donors=reattribute_donors,` (~line 158):

```python
            keep_recombinant=keep_recombinant,
```

- [ ] **Step 5: Type the pool in `_select_from` and forward selection**

In `src/tessera/discover/iterate.py`, replace `_select_from` (lines 567-574) with:

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
    lineage_of = lineage_map_from_rows(rows)
    return select_regional(
        params.query, genomes, window=params.seed_window, per_window=SEED_PER_WINDOW,
        derep_ani=params.derep_diverse_ani, drop_siblings=not params.seed_keep_siblings,
        lineage_of=lineage_of or None, keep_recombinant=params.keep_recombinant,
        logger=logger,
    )
```

- [ ] **Step 6: Run the new test to verify it passes**

Run: `pytest tests/unit/test_select_from_lineage.py -v`
Expected: PASS (1 passed).

- [ ] **Step 7: Run the iterate/fill tests to confirm no regression**

Run: `pytest tests/unit -k "iterate or fill or select_from" -v`
Expected: PASS.

- [ ] **Step 8: Ruff**

Run: `ruff check src/tessera/discover/iterate.py tests/unit/test_select_from_lineage.py`
Expected: clean.

- [ ] **Step 9: Commit**

```bash
git add src/tessera/discover/iterate.py tests/unit/test_select_from_lineage.py
git commit -m "$(cat <<'EOF'
Type the pool and forward lineage selection in _select_from

_select_from now builds a lineage map from the pool headers (reusing
build_lineage_map) and forwards it plus keep_recombinant into select_regional,
so a typed pool is reduced by lineage. FillParams gains keep_recombinant
(default False), forwarded through for_detection.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: CLI toggle on `detect` and `fill-references`

**Files:**
- Modify: `src/tessera/cli/cmd_detect.py` (option block ~line 71-75; `FillParams.for_detection(...)` call ~line 94-105)
- Modify: `src/tessera/cli/cmd_fill_references.py` (option block ~line 155; `FillParams(...)` call ~line 177-194)
- Test: `tests/unit/test_cli_keep_recombinant.py` (create)

**Interfaces:**
- Consumes: `FillParams.keep_recombinant` and `FillParams.for_detection(keep_recombinant=...)` from Task 2.
- Produces: a `--keep-recombinant-lineages/--no-keep-recombinant-lineages` option (default off) on both commands.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_cli_keep_recombinant.py`:

```python
"""The keep-recombinant toggle reaches FillParams from both CLI commands."""

from __future__ import annotations

from tessera.discover.iterate import FillParams


def test_default_excludes_recombinant():
    params = FillParams(query=None, collection=None, output=None)
    assert params.keep_recombinant is False


def test_for_detection_forwards_keep_recombinant():
    params = FillParams.for_detection(
        query=None, output=None, keep_recombinant=True,
    )
    assert params.keep_recombinant is True
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `pytest tests/unit/test_cli_keep_recombinant.py -v`
Expected: FAIL (`for_detection` has no `keep_recombinant` yet -- this guards Task 2's wiring at the boundary the CLI uses).

Note: if Task 2 is already merged this test passes at Step 2; that is acceptable -- proceed to add the CLI options below, which the test does not itself exercise but which Steps 6-7 verify via `--help`.

- [ ] **Step 3: Add the option to `detect`**

In `src/tessera/cli/cmd_detect.py`, add after the `reattribute_donors` option (after line 75, before the closing `) -> None:`):

```python
    keep_recombinant_lineages: bool = typer.Option(
        False, "--keep-recombinant-lineages/--no-keep-recombinant-lineages",
        help="Keep recombinant (CRF/URF/X) lineages in the panel. Off by default: "
        "they carry both parents' segments and mask the true parents. Turn on to "
        "recruit the best-matching genome even if it is a known hybrid.",
    ),
```

Then pass it into the `FillParams.for_detection(...)` call, after `reattribute_donors=reattribute_donors,` (line 104):

```python
            keep_recombinant=keep_recombinant_lineages,
```

- [ ] **Step 4: Add the option to `fill-references`**

In `src/tessera/cli/cmd_fill_references.py`, add after the `reattribute_donors` option (after line ~158, before the closing `) -> None:`), mirroring the detect help text:

```python
    keep_recombinant_lineages: bool = typer.Option(
        False, "--keep-recombinant-lineages/--no-keep-recombinant-lineages",
        help="Keep recombinant (CRF/URF/X) lineages in the panel. Off by default: "
        "they carry both parents' segments and mask the true parents. Turn on to "
        "recruit the best-matching genome even if it is a known hybrid.",
    ),
```

Then pass it into the `FillParams(...)` call, after `reattribute_donors=reattribute_donors,` (line 194):

```python
            keep_recombinant=keep_recombinant_lineages,
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `pytest tests/unit/test_cli_keep_recombinant.py -v`
Expected: PASS (2 passed).

- [ ] **Step 6: Verify the detect option is registered**

Run: `python -m tessera.cli.main detect --help | grep -A1 keep-recombinant`
Expected: shows `--keep-recombinant-lineages / --no-keep-recombinant-lineages`.

- [ ] **Step 7: Verify the fill-references option is registered**

Run: `python -m tessera.cli.main fill-references --help | grep -A1 keep-recombinant`
Expected: shows the same option.

- [ ] **Step 8: Ruff**

Run: `ruff check src/tessera/cli/cmd_detect.py src/tessera/cli/cmd_fill_references.py tests/unit/test_cli_keep_recombinant.py`
Expected: clean.

- [ ] **Step 9: Commit**

```bash
git add src/tessera/cli/cmd_detect.py src/tessera/cli/cmd_fill_references.py tests/unit/test_cli_keep_recombinant.py
git commit -m "$(cat <<'EOF'
Add --keep-recombinant-lineages toggle to detect and fill-references

Default off: recombinant (CRF/URF/X) lineages are excluded from a typed panel
because they mask the true parents. Turn on to recruit the best-matching genome
even when it is a known hybrid.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: Validation harness -- measure the gain

**Files:**
- Modify: `validation/run_hybrids.py` (the tip-branch `select_regional` call ~line 575)

**Interfaces:**
- Consumes: `select_regional(..., lineage_of=..., keep_recombinant=False)` from Task 1; `clade_of_label(label, tips)` (harness, line 416); `strip_sequence_extension`.
- Produces: no new public interface; the harness now reduces the tip pool by lineage.

- [ ] **Step 1: Build the lineage map and pass it in**

In `validation/run_hybrids.py`, in the `else` branch that calls `select_regional` (around line 575), build a `lineage_of` map from the tips and pass it. Replace:

```python
            selected = select_regional(setup.query, setup.pool, window=setup.sel_window,
                                       per_window=2, drop_siblings=False, logger=logger).selected
```

with:

```python
            # Type the pool by tree clade so the panel is reduced by lineage (one
            # query-closest representative per clade, recombinant clades excluded) --
            # this keeps a true parent that ANI would collapse (rsv_a's A.D.1.8) and
            # drops the masking CRF clades (hiv1).
            lineage_of = {
                strip_sequence_extension(g.name): clade_of_label(g.name, setup.tips)
                for g in setup.pool
            }
            lineage_of = {k: v for k, v in lineage_of.items() if v and v not in ("?", "NA")}
            selected = select_regional(
                setup.query, setup.pool, window=setup.sel_window, per_window=2,
                drop_siblings=False, lineage_of=lineage_of or None,
                keep_recombinant=False, logger=logger,
            ).selected
```

Confirm `strip_sequence_extension` is imported in `run_hybrids.py` (it is used already at lines 565-566); if not, add it to the existing `from tessera.core.io import ...` line.

- [ ] **Step 2: Ruff**

Run: `ruff check validation/run_hybrids.py`
Expected: clean.

- [ ] **Step 3: Run the harness comparison (opt-in, needs the aligner env)**

Run:
```bash
export PATH="$HOME/miniforge3/envs/recomfi-aln/bin:$PATH"
python validation/run_hybrids.py --compare
```
Expected: `rsv_a` and `hiv1` both PASS while the other 16 pathogens hold -- a clean **17/18**. (This step is a measurement gate, not a unit test; it needs mafft/minimap2/skani/skDER on PATH and is not part of CI.)

- [ ] **Step 4: Record the result**

Append a short section to `validation/attribution-results.md` stating the before (16/18 baseline) and after (17/18 with lineage selection) counts, which pathogens changed (`rsv_a`, `hiv1`), and the remaining failure (`flu_h3n2_ha`, backbone mismatch -- unchanged). Use modest language; state what passed, what failed, and that the harness path (not CI) produced it.

- [ ] **Step 5: Commit**

```bash
git add validation/run_hybrids.py validation/attribution-results.md
git commit -m "$(cat <<'EOF'
Reduce the harness tip pool by lineage and record 17/18

The hybrid harness now types its tip pool by tree clade and passes it to
select_regional, so the panel is built by lineage: rsv_a keeps its true parent
(A.D.1.8) that ANI collapsed and hiv1 drops the masking CRF clades. Records the
measured move from 16/18 to a clean 17/18.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Self-Review

**Spec coverage:**
- `select_regional` gains `lineage_of` + `keep_recombinant`, reduction by lineage with ANI fallback -> Task 1.
- Product wiring (`_select_from` builds the map from headers; `FillParams.keep_recombinant`, forwarded through `for_detection`) -> Task 2.
- CLI toggle on `detect` + `fill-references` (default exclude) -> Task 3.
- Harness measurement targeting a clean 17/18, recorded in `attribution-results.md` -> Task 4.
- Backward compatibility (`lineage_of=None` == today) -> Task 1 Step 5 (`else` branch preserves the exact prior sequence) + Step 7 regression run.
- Reuse of `is_recombinant_lineage` / `skani_query_ani` / `dereplicate` / `build_lineage_map` -> Tasks 1-2, no new dependency.

**Placeholder scan:** every code step shows complete code; no TODO/TBD; the only prose-only step (Task 4 Step 4) is a results write-up whose content is specified.

**Type consistency:** `_lineage_select(genomes, lineage_of, query_ani, *, keep_recombinant, derep_ani, logger)` and the two new `select_regional` params match between Task 1's definition and Task 2/4's call sites; `FillParams.keep_recombinant` and `for_detection(keep_recombinant=...)` match between Task 2 and Task 3; `lineage_of` is `dict[str, str] | None` at every call (`... or None` guards the empty-map case). `query_ani` values are `(ani, af)` tuples, consistent with `skani_query_ani`'s return.

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-07-01-lineage-aware-panel-selection.md`. Two execution options:

1. **Subagent-Driven (recommended)** - fresh subagent per task, review between tasks, fast iteration.
2. **Inline Execution** - execute tasks in this session with checkpoints.

Which approach?
