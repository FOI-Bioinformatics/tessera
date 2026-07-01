# Auto-typing in recruitment + ncbi-datasets source tag Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** An opt-in `--deep-typing` flag runs the full lineage-assignment ladder once at the end of `fill-references`/`detect` recruitment (instead of tier-1 header mining), and recruited NCBI-Virus genomes are tagged `ncbi-datasets` via a fetch sidecar.

**Architecture:** Two parts. Part B first (the plumbing the ladder consumes): `assign_lineages` gains a `datasets_rows` pass-through; `fetch_ncbi_virus` persists its lineage map as a sidecar; recruitment consolidates it to `<output>/ncbi_lineages.tsv`. Part A then makes `_type_panel` call the full ladder (with those rows) when `--deep-typing` is on, wired through `FillParams` and both CLI commands. Off = today's tier-1 path unchanged.

**Tech Stack:** Python 3 (typer), skani/Nextclade (already-optional external tools), pytest, ruff. Design spec: `docs/superpowers/specs/2026-07-01-auto-typing-in-recruitment-design.md`.

## Global Constraints

- No new runtime dependency; reuse `assign_lineages` (`discover/lineage_assign.py`), `build_lineage_map` / `write_lineage_map` / `titles_from_collection` / `genotype_from_title` / `first_header` / `LINEAGES_TSV` (`recomb/typing.py`), `skani_available` (`discover/panel.py`), `_lineage_labels` / `fetch_ncbi_virus` (`discover/pool.py`).
- `--deep-typing/--no-deep-typing` default **off** on both `detect` and `fill-references`; `FillParams.deep_typing: bool = False`.
- With `deep_typing` off, `_type_panel` must be byte-for-byte today's tier-1 behaviour.
- Sidecar file name is `ncbi_lineages.tsv`, two columns `accession<TAB>lineage`, distinct from the 3-column `lineages.tsv` (`LINEAGES_TSV`).
- Modest scientific language in code/docstrings. No Unicode in nextflow files (none touched). Ruff clean (line length 100, rules E/F/I/UP/B); CI coverage floor (80%) holds.
- Commit messages end with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.
- Branch `auto-typing-spec` (already checked out off main); do not commit on main.

---

### Task 1: `datasets_rows` pass-through in `assign_lineages`

**Files:**
- Modify: `src/tessera/discover/lineage_assign.py` (`assign_lineages` signature ~line 130-142)
- Test: `tests/unit/test_lineage_assign.py` (add one test; create the file if absent)

**Interfaces:**
- Consumes: `build_lineage_map(*, user_tsv, datasets_rows, title_by_label, organism)` (`recomb/typing.py:167`), whose priority is user > ncbi-datasets > title.
- Produces: `assign_lineages(..., datasets_rows: Iterable[tuple[str, str]] | None = None, ...)` forwarding `datasets_rows` into `build_lineage_map`.

- [ ] **Step 1: Write the failing test**

Add to `tests/unit/test_lineage_assign.py` (create with this header if the file does not exist):

```python
"""Unit tests for discover/lineage_assign.py (tool calls stubbed)."""

from __future__ import annotations

import logging
from pathlib import Path

from tessera.discover import lineage_assign as la

_LOG = logging.getLogger("tessera.test")


def test_datasets_rows_override_title_and_tag_source(tmp_path, monkeypatch):
    # One genome whose title mines to genotype "1" (source 'title'); a datasets row
    # for the same accession must win and be tagged 'ncbi-datasets'.
    g = tmp_path / "ACC1.fasta"
    g.write_text(">ACC1 genotype 1\nACGT\n")
    # No Nextclade dataset, no de-novo edges: keep the ladder on tier-1 only.
    monkeypatch.setattr(la, "_reference_tips", lambda **k: {})
    rows = la.assign_lineages(
        [g], taxon=None, datasets_rows=[("ACC1", "DENV1")], logger=_LOG,
    )
    by_label = {label: (lineage, source) for label, lineage, source in rows}
    assert by_label["ACC1"] == ("DENV1", "ncbi-datasets")
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/unit/test_lineage_assign.py::test_datasets_rows_override_title_and_tag_source -v`
Expected: FAIL with `TypeError: assign_lineages() got an unexpected keyword argument 'datasets_rows'`.

- [ ] **Step 3: Add the parameter and forward it**

In `src/tessera/discover/lineage_assign.py`, add `datasets_rows` to the `assign_lineages` signature (keyword-only, after `nextclade_dataset`):

```python
def assign_lineages(
    genomes: list[Path], *, user_lineage_map: Path | None = None,
    taxon: str | None = None, nextclade_dataset: str | None = None,
    datasets_rows: Iterable[tuple[str, str]] | None = None,
    ref_ani_floor: float = DEFAULT_REF_ANI_FLOOR, cluster_ani: float = DEFAULT_CLUSTER_ANI,
    email: str | None = None, cache_dir: Path | None = None, logger: logging.Logger,
) -> list[tuple[str, str, str]]:
```

Add the import if not present (top of file): `from collections.abc import Iterable`.

Then forward it in the `build_lineage_map(...)` call (~line 138):

```python
    rows = build_lineage_map(
        user_tsv=user_lineage_map,
        datasets_rows=datasets_rows,
        title_by_label=titles_from_collection(genomes),
        organism=taxon,
    )
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `python -m pytest tests/unit/test_lineage_assign.py -v`
Expected: PASS.

- [ ] **Step 5: Run the existing lineage-assign tests + ruff**

Run: `python -m pytest tests/unit -k "lineage_assign or assign_lineages" -v && ruff check src/tessera/discover/lineage_assign.py tests/unit/test_lineage_assign.py`
Expected: PASS, ruff clean.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/discover/lineage_assign.py tests/unit/test_lineage_assign.py
git commit -m "$(cat <<'EOF'
Add datasets_rows pass-through to assign_lineages

assign_lineages forwards an optional datasets_rows list into build_lineage_map,
so a structured NCBI-datasets lineage overrides a mined title token for the same
accession and is tagged 'ncbi-datasets'.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: `fetch_ncbi_virus` writes the `ncbi_lineages.tsv` sidecar

**Files:**
- Modify: `src/tessera/discover/pool.py` (near `fetch_ncbi_virus` line 151 and module constants ~line 40)
- Test: `tests/unit/test_ncbi_lineage_sidecar.py` (create)

**Interfaces:**
- Consumes: the `labels: dict[str, str]` (accession -> lineage) already built at `fetch_ncbi_virus` line 151 by `_lineage_labels`.
- Produces: module constant `NCBI_LINEAGES_TSV = "ncbi_lineages.tsv"`; helper `write_ncbi_lineages(dest: Path, labels: dict[str, str]) -> Path`; a call to it inside `fetch_ncbi_virus` so `<dest>/ncbi_lineages.tsv` exists after a fetch.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_ncbi_lineage_sidecar.py`:

```python
"""The NCBI-datasets lineage sidecar written during fetch."""

from __future__ import annotations

from pathlib import Path

from tessera.discover.pool import NCBI_LINEAGES_TSV, write_ncbi_lineages


def test_write_ncbi_lineages_round_trip(tmp_path: Path):
    dest = tmp_path / "cache"
    dest.mkdir()
    out = write_ncbi_lineages(dest, {"ACC1": "DENV1", "ACC2": "DENV4"})
    assert out == dest / NCBI_LINEAGES_TSV
    lines = sorted(out.read_text().splitlines())
    assert lines == ["ACC1\tDENV1", "ACC2\tDENV4"]


def test_write_ncbi_lineages_empty_writes_nothing(tmp_path: Path):
    dest = tmp_path / "cache"
    dest.mkdir()
    out = write_ncbi_lineages(dest, {})
    assert out is None
    assert not (dest / NCBI_LINEAGES_TSV).exists()
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/unit/test_ncbi_lineage_sidecar.py -v`
Expected: FAIL with `ImportError: cannot import name 'NCBI_LINEAGES_TSV'`.

- [ ] **Step 3: Add the constant and helper**

In `src/tessera/discover/pool.py`, add near the other module constants (after line 43, `DEFAULT_DEREP_ANI`):

```python
NCBI_LINEAGES_TSV = "ncbi_lineages.tsv"  # sidecar: accession<TAB>datasets-lineage
```

Add the helper (place it just above `fetch_ncbi_virus`, ~line 110):

```python
def write_ncbi_lineages(dest: Path, labels: dict[str, str]) -> Path | None:
    """Persist an accession -> NCBI-datasets-lineage map as a 2-column sidecar in
    ``dest`` (``accession<TAB>lineage``). Returns the path, or ``None`` when there is
    nothing to write. Lets a later typing pass tag these accessions ``ncbi-datasets``
    after the source ``data_report.jsonl`` has been discarded."""
    if not labels:
        return None
    out = dest / NCBI_LINEAGES_TSV
    with open(out, "w") as fo:
        for acc, lineage in sorted(labels.items()):
            fo.write(f"{acc}\t{lineage}\n")
    return out
```

- [ ] **Step 4: Call it inside `fetch_ncbi_virus`**

In `src/tessera/discover/pool.py`, `fetch_ncbi_virus`, right after the existing
`labels = _lineage_labels(data_report.jsonl)` line (line 151) and before `_split_fasta`:

```python
        labels = _lineage_labels(data_report.jsonl)
        write_ncbi_lineages(dest, labels)
        written = _split_fasta(fna, dest, labels, limit=limit)
```

(Only the middle line is added; the surrounding two lines already exist.)

- [ ] **Step 5: Run the tests + ruff**

Run: `python -m pytest tests/unit/test_ncbi_lineage_sidecar.py -v && ruff check src/tessera/discover/pool.py tests/unit/test_ncbi_lineage_sidecar.py`
Expected: PASS, ruff clean.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/discover/pool.py tests/unit/test_ncbi_lineage_sidecar.py
git commit -m "$(cat <<'EOF'
Persist NCBI-datasets lineages as a fetch sidecar

fetch_ncbi_virus writes an accession<TAB>lineage sidecar (ncbi_lineages.tsv)
into its dest from the data_report.jsonl labels, before the report is discarded,
so recruitment can later tag those accessions 'ncbi-datasets'.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Consolidate and read the sidecar in recruitment

**Files:**
- Modify: `src/tessera/discover/iterate.py` (import block ~line 40; `_fetch_diverse` line 525-567; add two helpers near `_type_panel` line 361)
- Test: `tests/unit/test_ncbi_lineage_consolidate.py` (create)

**Interfaces:**
- Consumes: `NCBI_LINEAGES_TSV` from `discover/pool.py` (Task 2).
- Produces: `_record_ncbi_lineages(output: Path, source_dir: Path) -> None` (merge/dedup `source_dir`'s sidecar into `<output>/ncbi_lineages.tsv`); `_read_ncbi_lineages(output: Path) -> list[tuple[str, str]] | None` (read `<output>/ncbi_lineages.tsv`, or `None` when absent). A call to `_record_ncbi_lineages(params.output, cache)` in `_fetch_diverse`.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_ncbi_lineage_consolidate.py`:

```python
"""Consolidating and reading the NCBI-datasets lineage sidecar in recruitment."""

from __future__ import annotations

from pathlib import Path

from tessera.discover.iterate import _read_ncbi_lineages, _record_ncbi_lineages
from tessera.discover.pool import NCBI_LINEAGES_TSV


def _sidecar(d: Path, rows: dict[str, str]) -> None:
    d.mkdir(parents=True, exist_ok=True)
    (d / NCBI_LINEAGES_TSV).write_text("".join(f"{a}\t{l}\n" for a, l in rows.items()))


def test_read_absent_returns_none(tmp_path: Path):
    assert _read_ncbi_lineages(tmp_path) is None


def test_record_then_read_merges_and_dedups(tmp_path: Path):
    out = tmp_path / "out"
    out.mkdir()
    src1 = tmp_path / "cache1"
    src2 = tmp_path / "cache2"
    _sidecar(src1, {"ACC1": "DENV1", "ACC2": "DENV4"})
    _sidecar(src2, {"ACC2": "DENV4", "ACC3": "DENV2"})  # ACC2 duplicate
    _record_ncbi_lineages(out, src1)
    _record_ncbi_lineages(out, src2)
    rows = dict(_read_ncbi_lineages(out))
    assert rows == {"ACC1": "DENV1", "ACC2": "DENV4", "ACC3": "DENV2"}


def test_record_missing_source_is_noop(tmp_path: Path):
    out = tmp_path / "out"
    out.mkdir()
    _record_ncbi_lineages(out, tmp_path / "does_not_exist")
    assert _read_ncbi_lineages(out) is None
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/unit/test_ncbi_lineage_consolidate.py -v`
Expected: FAIL with `ImportError: cannot import name '_read_ncbi_lineages'`.

- [ ] **Step 3: Add the import**

In `src/tessera/discover/iterate.py`, add `NCBI_LINEAGES_TSV` to the existing `from .pool import ...` (there is a local import of pool symbols; if pool is only imported inside functions, add a module-level `from .pool import NCBI_LINEAGES_TSV` with the other `from .` imports near line 40).

- [ ] **Step 4: Add the two helpers**

In `src/tessera/discover/iterate.py`, add near `_type_panel` (before it, ~line 360):

```python
def _record_ncbi_lineages(output: Path, source_dir: Path) -> None:
    """Merge ``source_dir``'s NCBI-datasets lineage sidecar into the run's consolidated
    ``<output>/ncbi_lineages.tsv`` (dedup by accession). No-op when the source has none."""
    src = source_dir / NCBI_LINEAGES_TSV
    if not src.exists():
        return
    merged: dict[str, str] = dict(_read_ncbi_lineages(output) or [])
    for line in src.read_text().splitlines():
        acc, _, lineage = line.partition("\t")
        if acc and lineage:
            merged[acc] = lineage
    with open(output / NCBI_LINEAGES_TSV, "w") as fo:
        for acc, lineage in sorted(merged.items()):
            fo.write(f"{acc}\t{lineage}\n")


def _read_ncbi_lineages(output: Path) -> list[tuple[str, str]] | None:
    """Read ``<output>/ncbi_lineages.tsv`` into ``(accession, lineage)`` rows, or
    ``None`` when the sidecar is absent."""
    path = output / NCBI_LINEAGES_TSV
    if not path.exists():
        return None
    rows: list[tuple[str, str]] = []
    for line in path.read_text().splitlines():
        acc, _, lineage = line.partition("\t")
        if acc and lineage:
            rows.append((acc, lineage))
    return rows
```

- [ ] **Step 5: Wire the consolidation into `_fetch_diverse`**

In `src/tessera/discover/iterate.py`, `_fetch_diverse`, record the cache's sidecar just before every `return` of genomes. The cleanest single point is to capture the return value; change the three return sites (cache-hit line 546, refseq-enough line 551, and the final line 567) so each records first. Concretely, wrap the returns:

At the cache-hit branch (line 541-546), before `return existing`:
```python
        _record_ncbi_lineages(params.output, cache)
        return existing
```
At the refseq-enough branch (line 550-551), before `return fetched`:
```python
        if len(fetched) >= SEED_MIN_DIVERSE:
            _record_ncbi_lineages(params.output, cache)
            return fetched
```
At the final return (line 567):
```python
    _record_ncbi_lineages(params.output, cache)
    return fetched
```

(`params.output` exists at every point; `cache` is the `ncbi_virus_cache(...)` dir assigned at line 539. The sidecar was written into `cache` by `fetch_ncbi_virus` on the fetching run and persists there across cached runs.)

- [ ] **Step 6: Run the tests + the iterate suite + ruff**

Run: `python -m pytest tests/unit/test_ncbi_lineage_consolidate.py tests/unit -k "iterate or fill or fetch" -v && ruff check src/tessera/discover/iterate.py tests/unit/test_ncbi_lineage_consolidate.py`
Expected: PASS, ruff clean.

- [ ] **Step 7: Commit**

```bash
git add src/tessera/discover/iterate.py tests/unit/test_ncbi_lineage_consolidate.py
git commit -m "$(cat <<'EOF'
Consolidate NCBI-datasets lineage sidecars in recruitment

_fetch_diverse merges each NCBI-Virus cache's ncbi_lineages.tsv into one
<output>/ncbi_lineages.tsv; _read_ncbi_lineages reads it back (None when absent).
This is the datasets_rows source the deep-typing ladder will consume.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: `--deep-typing` in `_type_panel` (FillParams + guard + ladder branch)

**Files:**
- Modify: `src/tessera/discover/iterate.py` (`FillParams` ~line 106; `FillParams.for_detection` ~line 115-160; `fill_references` guard ~line 181; `_type_panel` line 361-386; import `assign_lineages`)
- Test: `tests/unit/test_deep_typing.py` (create)

**Interfaces:**
- Consumes: `assign_lineages(genomes, *, user_lineage_map, taxon, nextclade_dataset, datasets_rows, email, cache_dir, logger)` (Task 1); `_read_ncbi_lineages` (Task 3); `skani_available` (`discover/panel.py`).
- Produces: `FillParams.deep_typing: bool = False` (forwarded by `for_detection`); `_type_panel` runs the full ladder when `deep_typing` is set, else the current tier-1 path.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_deep_typing.py`:

```python
"""_type_panel routing: tier-1 by default, the full ladder under --deep-typing."""

from __future__ import annotations

import logging
from pathlib import Path

from tessera.discover import iterate as it
from tessera.discover.iterate import FillParams, _type_panel
from tessera.recomb.typing import LINEAGES_TSV

_LOG = logging.getLogger("tessera.test")


def _panel(tmp_path: Path) -> tuple[FillParams, Path]:
    out = tmp_path / "out"
    coll = out / "collection"
    coll.mkdir(parents=True)
    (coll / "ACC1.fasta").write_text(">ACC1 genotype 1\nACGT\n")
    q = tmp_path / "q.fasta"
    q.write_text(">QRY genotype 2\nACGT\n")
    return FillParams(query=q, collection=None, output=out, taxon=None), out


def test_default_uses_tier1_not_ladder(tmp_path, monkeypatch):
    params, out = _panel(tmp_path)
    called = {"ladder": False}
    monkeypatch.setattr(it, "assign_lineages",
                        lambda *a, **k: called.__setitem__("ladder", True) or [])
    _type_panel(params, out / "collection", "QRY", _LOG)
    assert called["ladder"] is False
    assert (out / LINEAGES_TSV).exists()


def test_deep_typing_calls_ladder_with_datasets_rows(tmp_path, monkeypatch):
    params, out = _panel(tmp_path)
    params.deep_typing = True
    (out / "ncbi_lineages.tsv").write_text("ACC1\tDENV1\n")
    seen = {}

    def fake_ladder(genomes, **kwargs):
        seen.update(kwargs)
        return [("ACC1", "DENV1", "ncbi-datasets")]

    monkeypatch.setattr(it, "assign_lineages", fake_ladder)
    lineage_map, _ = _type_panel(params, out / "collection", "QRY", _LOG)
    assert seen["datasets_rows"] == [("ACC1", "DENV1")]
    assert lineage_map["ACC1"] == "DENV1"
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/unit/test_deep_typing.py -v`
Expected: FAIL (`FillParams` has no `deep_typing`, and `iterate` has no `assign_lineages` attribute to patch).

- [ ] **Step 3: Add the `FillParams` field**

In `src/tessera/discover/iterate.py`, add after `keep_recombinant: bool = False` (last field of `FillParams`):

```python
    deep_typing: bool = False  # type the recruited panel with the full ladder, not just headers
```

- [ ] **Step 4: Forward it through `for_detection`**

In `FillParams.for_detection`, add a parameter after `keep_recombinant: bool = False,`:

```python
        deep_typing: bool = False,
```

and add to the returned `cls(...)` call, after `keep_recombinant=keep_recombinant,`:

```python
            deep_typing=deep_typing,
```

- [ ] **Step 5: Add the import and the skani guard**

In `src/tessera/discover/iterate.py`, add a module-level import with the other `from .` imports:

```python
from .lineage_assign import assign_lineages
```

In `fill_references`, extend the existing `--curate` skani guard (line 181-184) with a deep-typing check right after it:

```python
    if params.deep_typing and not skani_available():
        raise UserInputError(
            "--deep-typing runs the lineage ladder (Nextclade + de-novo) and needs skani. "
            "Install with: conda install -c bioconda skani skder"
        )
```

(`skani_available` is already imported in this module for the `--curate` guard; reuse it.)

- [ ] **Step 6: Branch `_type_panel` on the flag**

In `src/tessera/discover/iterate.py`, replace the body of `_type_panel` (lines 371-376, the `coll_files` + `lineage_rows = build_lineage_map(...)` block) with a branch:

```python
    coll_files = [p for p in collection.iterdir() if p.is_file()]
    if params.deep_typing:
        lineage_rows = assign_lineages(
            coll_files,
            user_lineage_map=params.lineage_map,
            taxon=params.taxon,
            nextclade_dataset=params.nextclade_dataset,
            datasets_rows=_read_ncbi_lineages(params.output),
            email=params.email,
            cache_dir=params.cache_dir,
            logger=logger,
        )
    else:
        lineage_rows = build_lineage_map(
            user_tsv=params.lineage_map,
            title_by_label=titles_from_collection(coll_files),
            organism=params.taxon,
        )
```

The rest of `_type_panel` (the query row append, `write_lineage_map`, `lineage_map_from_rows`, the log line, the return) is unchanged.

- [ ] **Step 7: Run the tests to verify they pass**

Run: `python -m pytest tests/unit/test_deep_typing.py -v`
Expected: PASS (2 passed).

- [ ] **Step 8: Run the iterate suite + ruff**

Run: `python -m pytest tests/unit -k "iterate or fill or type_panel or deep_typing" -v && ruff check src/tessera/discover/iterate.py tests/unit/test_deep_typing.py`
Expected: PASS, ruff clean.

- [ ] **Step 9: Commit**

```bash
git add src/tessera/discover/iterate.py tests/unit/test_deep_typing.py
git commit -m "$(cat <<'EOF'
Run the full lineage ladder in recruitment under --deep-typing

FillParams.deep_typing (default off) makes _type_panel type the recruited panel
with assign_lineages (Nextclade nearest-neighbour + de-novo), passing the
consolidated ncbi_lineages.tsv as datasets_rows, instead of tier-1 header
mining. A skani guard mirrors --curate. Off is unchanged.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: `--deep-typing` CLI option on `detect` and `fill-references`

**Files:**
- Modify: `src/tessera/cli/cmd_detect.py` (option block after `keep_recombinant_lineages` ~line 76-81; `FillParams.for_detection(...)` call ~line 94-105)
- Modify: `src/tessera/cli/cmd_fill_references.py` (option block after `keep_recombinant_lineages` ~line 160-165; `FillParams(...)` call ~line 177-195)
- Test: `tests/unit/test_cli_deep_typing.py` (create)

**Interfaces:**
- Consumes: `FillParams.deep_typing` and `FillParams.for_detection(deep_typing=...)` (Task 4).
- Produces: a `--deep-typing/--no-deep-typing` option (default off) on both commands.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_cli_deep_typing.py`:

```python
"""The deep-typing toggle reaches FillParams from both CLI commands."""

from __future__ import annotations

from tessera.discover.iterate import FillParams


def test_default_deep_typing_off():
    assert FillParams(query=None, collection=None, output=None).deep_typing is False


def test_for_detection_forwards_deep_typing():
    params = FillParams.for_detection(query=None, output=None, deep_typing=True)
    assert params.deep_typing is True
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/unit/test_cli_deep_typing.py -v`
Expected: FAIL (`for_detection` has no `deep_typing`) unless Task 4 is already merged, in which case it PASSes at Step 2 -- acceptable; still add the CLI options below (Steps 5-6 verify them via `--help`).

- [ ] **Step 3: Add the option to `detect`**

In `src/tessera/cli/cmd_detect.py`, after the `keep_recombinant_lineages` option block, before the closing `) -> None:`:

```python
    deep_typing: bool = typer.Option(
        False, "--deep-typing/--no-deep-typing",
        help="Type the recruited panel with the full lineage ladder (Nextclade "
        "nearest-neighbour + de-novo ANI clustering) instead of header mining alone. "
        "Needs skani; off by default.",
    ),
```

Pass it into the `FillParams.for_detection(...)` call, after `keep_recombinant=keep_recombinant_lineages,`:

```python
            deep_typing=deep_typing,
```

- [ ] **Step 4: Add the option to `fill-references`**

In `src/tessera/cli/cmd_fill_references.py`, after the `keep_recombinant_lineages` option block, before the closing `) -> None:`, add the identical option:

```python
    deep_typing: bool = typer.Option(
        False, "--deep-typing/--no-deep-typing",
        help="Type the recruited panel with the full lineage ladder (Nextclade "
        "nearest-neighbour + de-novo ANI clustering) instead of header mining alone. "
        "Needs skani; off by default.",
    ),
```

Pass it into the `FillParams(...)` call, after `keep_recombinant=keep_recombinant_lineages,`:

```python
            deep_typing=deep_typing,
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `python -m pytest tests/unit/test_cli_deep_typing.py -v`
Expected: PASS (2 passed).

- [ ] **Step 6: Verify both options render**

Run: `python -m tessera.cli.main detect --help | grep -A1 deep-typing && python -m tessera.cli.main fill-references --help | grep -A1 deep-typing`
Expected: both show `--deep-typing / --no-deep-typing`.

- [ ] **Step 7: Ruff**

Run: `ruff check src/tessera/cli/cmd_detect.py src/tessera/cli/cmd_fill_references.py tests/unit/test_cli_deep_typing.py`
Expected: clean.

- [ ] **Step 8: Commit**

```bash
git add src/tessera/cli/cmd_detect.py src/tessera/cli/cmd_fill_references.py tests/unit/test_cli_deep_typing.py
git commit -m "$(cat <<'EOF'
Add --deep-typing CLI option to detect and fill-references

Default off: type the recruited panel with header mining only. On: run the full
lineage ladder (Nextclade nearest-neighbour + de-novo). Threaded to
FillParams.deep_typing (via for_detection for detect).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Documentation

**Files:**
- Modify: `docs/reference-panels.md` (the lineage/typing section)

**Interfaces:** none (docs only).

- [ ] **Step 1: Document `--deep-typing` and the source tag**

In `docs/reference-panels.md`, in the section that covers lineage typing (near where `type-lineages` / `lineages.tsv` is described), add a short subsection: `--deep-typing` runs the full `type-lineages` ladder on the recruited panel at the end of a run (needs skani; off by default), so recruited genomes get Nextclade nearest-neighbour / de-novo labels rather than header-mined tokens alone; recruited NCBI-Virus genomes are tagged with source `ncbi-datasets` via the `ncbi_lineages.tsv` sidecar. Modest language; state that it is opt-in and adds a Nextclade/skani pass at the end.

- [ ] **Step 2: Commit**

```bash
git add docs/reference-panels.md
git commit -m "$(cat <<'EOF'
Document --deep-typing and the ncbi-datasets source tag

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Self-Review

**Spec coverage:**
- Part A: `--deep-typing` flag, default off, ladder in `_type_panel`, skani guard, CLI on both commands -> Tasks 4, 5.
- Part B Piece 1: `assign_lineages` `datasets_rows` pass-through -> Task 1.
- Part B Piece 2: fetch sidecar + consolidate + read -> Tasks 2, 3, consumed in Task 4.
- Backward compatibility (`deep_typing` off == tier-1) -> Task 4 Step 6 `else` branch + `test_default_uses_tier1_not_ladder`.
- No new dependency; reuse enumerated helpers -> all tasks.
- Docs -> Task 6.

**Placeholder scan:** every code step shows complete code; the only prose-only steps are Task 6 (doc content specified) -- no TODO/TBD.

**Type consistency:** `datasets_rows` is `Iterable[tuple[str, str]] | None` in `assign_lineages` (Task 1) and the value passed by `_type_panel` is `_read_ncbi_lineages(...) -> list[tuple[str, str]] | None` (Task 3) -- compatible. `NCBI_LINEAGES_TSV` defined in `pool.py` (Task 2) and imported in `iterate.py` (Task 3). `FillParams.deep_typing: bool` matches `for_detection(deep_typing=...)` (Task 4) and the CLI `deep_typing` variable (Task 5). `_read_ncbi_lineages` / `_record_ncbi_lineages` signatures match between Task 3's definition and Task 4's / `_fetch_diverse`'s call sites.

**Task ordering:** Task 1 (param) -> Task 2 (sidecar write) -> Task 3 (consolidate/read) -> Task 4 (consume both) -> Task 5 (CLI) -> Task 6 (docs). Each task's deliverable is independently testable.

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-07-01-auto-typing-in-recruitment.md`. Two execution options:

1. **Subagent-Driven (recommended)** - fresh subagent per task, review between tasks, fast iteration.
2. **Inline Execution** - execute tasks in this session with checkpoints.

Which approach?
