# reassort v2 (intragenic per-segment scan) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an opt-in `--scan-segments` flag to `tessera reassort` that runs the existing `run_recomb` pipeline on each assigned segment against its per-clade-consensus Nextclade panel, reporting within-segment recombination alongside the reassortment verdict.

**Architecture:** A new module `reassort/scan.py` builds a per-clade-consensus panel for one segment, aligns the segment query to it, runs the full `run_recomb` scan into `out/<segment>/`, and summarizes the regions. `assign.py` threads the already-resolved dataset out of `_type_segment` and, when `scan_segments` is on, scans each assigned segment. The CLI adds `--scan-segments`/`--aligner`, writes a `segment_scan.tsv` rollup, and prints a summary. Everything reuses `build_pool`/`build_msa`/`run_recomb` unchanged.

**Tech Stack:** Python 3.11+, skani + Nextclade + an aligner (mafft/minimap2) external backends, typer CLI, pytest.

## Global Constraints

- No new runtime dependency in the package. The aligner is the same optional external backend `detect` already uses; it becomes required only when `--scan-segments` is passed.
- No change to `detect`/`recomb`/the reassortment call. The scan reuses `run_recomb` as-is.
- A missing/unknown aligner is an up-front error when `--scan-segments` is set (mirrors the `skani_available()` check). A per-segment scan failure is non-fatal (that segment recorded not-scanned; the run continues).
- Modest scientific language in code and docs; faithful reporting (every input segment gets a `segment_scan.tsv` row; the opt-in e2e is recorded as measured).
- `ReassortmentResult` keeps `segments`/`verdict`/`groups`/`pair_notes` (positions unchanged; the harness builds `ReassortmentResult(segments=..., verdict=...)`); `scans` is added with a default so existing consumers are unaffected.
- Ruff clean (`ruff check src tests validation`, line length 100, E/F/I/UP/B); CI coverage floor 80% holds.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

---

## File Structure

- `src/tessera/reassort/scan.py` (create) -- `SegmentScan`, `window_params`, `require_aligner`, `scan_segment`, `_summarize_regions`, `_clade_of_header`. No change to `constellation.py`/`detect`.
- `src/tessera/reassort/assign.py` (modify) -- `_type_segment` returns the dataset; `assign_segments` gains `scan_segments`/`aligner`/`output`; `ReassortmentResult.scans`.
- `src/tessera/reassort/__init__.py` (modify) -- export `SegmentScan`, `scan_segment`.
- `src/tessera/cli/cmd_reassort.py` (modify) -- `--scan-segments`, `--aligner`; `segment_scan.tsv`; rollup summary.
- `tests/unit/test_reassort_scan.py` (create) -- `_summarize_regions`, stubbed `scan_segment`.
- `tests/unit/test_reassort_assign.py` (modify) -- stubbed `assign_segments(scan_segments=True)`.
- `tests/unit/test_cli_reassort.py` (modify) -- `--scan-segments`/`--aligner` exposed.
- `docs/reference-panels.md` (modify) -- document `--scan-segments`.
- `validation/attribution-results.md` (modify) -- record the opt-in e2e.

---

## Task 1: `scan.py` -- per-segment intragenic scan

**Files:**
- Create: `src/tessera/reassort/scan.py`
- Create: `tests/unit/test_reassort_scan.py`

**Interfaces:**
- Consumes: `build_pool`, `NextcladeDataset` (`discover.nextclade`); `nextclade_cache` (`core.cache`); `MsaParams`, `build_msa` (`msa.build`); `RecombParams`, `run_recomb` (`recomb.run`); `DEFAULT_METHODS` (`recomb.regions`); `LINEAGES_TSV`, `write_lineage_map`, `lineage_map_from_rows`, `first_header` (`recomb.typing`); `strip_sequence_extension`, `write_fasta_record` (`core.io`); `UserInputError` (`core.errors`); `registry` (`aligners.base`).
- Produces:
  - `@dataclass SegmentScan(segment: str, scanned: bool, recombinant: bool, n_regions: int, note: str)`
  - `window_params(genome_len: int) -> tuple[int, int, int]`
  - `require_aligner(aligner: str) -> None` (raises on unknown name or missing binary)
  - `scan_segment(segment: str, seq: str, dataset, out_dir: Path, *, aligner: str, cache_dir: Path | None, logger) -> SegmentScan`
  - `_summarize_regions(path: Path) -> tuple[int, bool]`

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_reassort_scan.py`:

```python
"""Unit tests for the per-segment intragenic scan (stubbed; no aligner/binaries)."""

from __future__ import annotations

import logging
from pathlib import Path

from tessera.reassort import scan
from tessera.reassort.scan import SegmentScan, _summarize_regions, scan_segment, window_params

LOG = logging.getLogger("test")


class _DS:
    path = "flu/ha"
    tag = "tag"


def _regions_tsv(tmp_path, rows):
    p = tmp_path / "recombination_regions.tsv"
    header = "query_start\tquery_end\tdonor_absent"
    p.write_text(header + "\n" + "\n".join(rows) + ("\n" if rows else ""))
    return p


def test_summarize_counts_present_regions(tmp_path):
    p = _regions_tsv(tmp_path, ["100\t200\tno", "300\t400\tno"])
    assert _summarize_regions(p) == (2, True)


def test_summarize_excludes_donor_absent(tmp_path):
    p = _regions_tsv(tmp_path, ["100\t200\tno", "300\t400\tyes"])
    assert _summarize_regions(p) == (1, True)


def test_summarize_empty_file_is_not_recombinant(tmp_path):
    p = _regions_tsv(tmp_path, [])
    assert _summarize_regions(p) == (0, False)


def test_summarize_missing_file_is_not_recombinant(tmp_path):
    assert _summarize_regions(tmp_path / "nope.tsv") == (0, False)


def test_window_params_scales_to_short_segment():
    window, step, sel = window_params(1700)
    assert 120 <= window <= 500
    assert step == max(20, window // 10)
    assert sel >= window


def test_scan_segment_single_clade_panel_is_skipped(tmp_path, monkeypatch):
    # A one-clade panel has nothing to recombine between -> scanned=False, note explains why.
    one = tmp_path / "pool" / "A_consensus.fasta"
    one.parent.mkdir(parents=True)
    one.write_text(">A_consensus A\nACGT\n")
    monkeypatch.setattr(scan, "build_pool", lambda ds, *, cache_dir, logger, per_clade_consensus: [one])
    monkeypatch.setattr(scan, "nextclade_cache", lambda p, t, override=None: Path("/x"))
    result = scan_segment("HA", "ACGTACGT", _DS(), tmp_path / "out", aligner="mafft",
                          cache_dir=None, logger=LOG)
    assert result == SegmentScan("HA", False, False, 0, "single-clade panel")


def test_scan_segment_runs_recomb_and_summarizes(tmp_path, monkeypatch):
    # Two-clade panel: build_msa/run_recomb are stubbed; run_recomb writes a regions file that
    # _summarize_regions reads back -> scanned + recombinant.
    pool = tmp_path / "pool"
    pool.mkdir()
    tips = []
    for c in ("A", "B"):
        t = pool / f"{c}_consensus.fasta"
        t.write_text(f">{c}_consensus {c}\nACGTACGT\n")
        tips.append(t)
    monkeypatch.setattr(scan, "build_pool", lambda ds, *, cache_dir, logger, per_clade_consensus: tips)
    monkeypatch.setattr(scan, "nextclade_cache", lambda p, t, override=None: Path("/x"))
    monkeypatch.setattr(scan, "build_msa", lambda params, logger: params.output)

    def fake_run_recomb(params, logger):
        (params.output / "recombination_regions.tsv").write_text(
            "query_start\tquery_end\tdonor_absent\n100\t200\tno\n")
        return "bp"
    monkeypatch.setattr(scan, "run_recomb", fake_run_recomb)

    result = scan_segment("HA", "ACGTACGT", _DS(), tmp_path / "out", aligner="mafft",
                          cache_dir=None, logger=LOG)
    assert result.segment == "HA" and result.scanned and result.recombinant
    assert result.n_regions == 1
    assert (tmp_path / "out" / "HA" / "recombination_regions.tsv").exists()


def test_scan_segment_failure_is_non_fatal(tmp_path, monkeypatch):
    pool = tmp_path / "pool"
    pool.mkdir()
    tips = []
    for c in ("A", "B"):
        t = pool / f"{c}_consensus.fasta"
        t.write_text(f">{c}_consensus {c}\nACGT\n")
        tips.append(t)
    monkeypatch.setattr(scan, "build_pool", lambda ds, *, cache_dir, logger, per_clade_consensus: tips)
    monkeypatch.setattr(scan, "nextclade_cache", lambda p, t, override=None: Path("/x"))

    def boom(params, logger):
        raise RuntimeError("aligner blew up")
    monkeypatch.setattr(scan, "build_msa", boom)

    result = scan_segment("HA", "ACGT", _DS(), tmp_path / "out", aligner="mafft",
                          cache_dir=None, logger=LOG)
    assert result.scanned is False
    assert "scan failed" in result.note
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tests/unit/test_reassort_scan.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'tessera.reassort.scan'`.

- [ ] **Step 3: Implement `scan.py`**

Create `src/tessera/reassort/scan.py`:

```python
"""Intragenic per-segment recombination scan (the ``reassort --scan-segments`` flag).

For each assigned segment, build a per-clade-consensus panel from its Nextclade dataset,
align the segment query to it, and run the ordinary single-backbone recombination scan
(:func:`tessera.recomb.run.run_recomb`) inside that one segment. This is orthogonal to the
whole-segment reassortment call: reassortment asks which parent each segment came from; this
asks whether a single segment is itself a within-segment mosaic of two lineages.
"""

from __future__ import annotations

import logging
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

from ..core.cache import nextclade_cache
from ..core.errors import UserInputError
from ..core.io import strip_sequence_extension, write_fasta_record
from ..discover.nextclade import build_pool
from ..msa.build import MsaParams, build_msa
from ..recomb.regions import DEFAULT_METHODS
from ..recomb.run import RecombParams, run_recomb
from ..recomb.typing import (
    LINEAGES_TSV,
    first_header,
    lineage_map_from_rows,
    write_lineage_map,
)


@dataclass
class SegmentScan:
    segment: str
    scanned: bool          # False if it could not be scanned (single-clade panel or error)
    recombinant: bool      # True if run_recomb called at least one present region
    n_regions: int
    note: str


def window_params(genome_len: int) -> tuple[int, int, int]:
    """Adaptive (recomb window, step, selection window) scaled to a short segment length,
    so a ~1-2.5 kb flu segment is not given a window wider than its alignment."""
    window = max(120, min(500, genome_len // 12))
    step = max(20, window // 10)
    select_window = max(window, min(1500, genome_len // 4))
    return window, step, select_window


def require_aligner(aligner: str) -> None:
    """Fail up front if ``aligner`` is unknown or its binary is not on PATH."""
    from ..aligners.base import registry
    names = set(registry.names())
    if aligner not in names:
        raise UserInputError(
            f"Unknown aligner '{aligner}'. Available: {', '.join(sorted(names))}."
        )
    registry.create(aligner).preflight()  # raises MissingBinaryError if the binary is absent


def _clade_of_header(path: Path) -> str:
    """The clade token from a consensus genome's ``>{label} {clade}`` header."""
    parts = first_header(path).split(None, 1)
    return parts[1].strip() if len(parts) > 1 and parts[1].strip() else "?"


def _summarize_regions(path: Path) -> tuple[int, bool]:
    """Count present (non ``donor_absent``) regions in a recomb TSV -> ``(n, n > 0)``."""
    if not path.exists():
        return 0, False
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        return 0, False
    header = lines[0].split("\t")
    absent_idx = header.index("donor_absent") if "donor_absent" in header else None
    n = 0
    for ln in lines[1:]:
        fields = ln.split("\t")
        if absent_idx is not None and absent_idx < len(fields) and fields[absent_idx] == "yes":
            continue
        n += 1
    return n, n > 0


def scan_segment(
    segment: str, seq: str, dataset, out_dir: Path, *,
    aligner: str, cache_dir: Path | None, logger: logging.Logger,
) -> SegmentScan:
    """Scan one segment for intragenic recombination. Never raises: a failure is recorded as
    ``scanned=False`` so the caller can continue with the other segments."""
    seg_dir = out_dir / re.sub(r"[^\w.-]+", "_", segment)
    try:
        panel = build_pool(
            dataset,
            cache_dir=nextclade_cache(dataset.path, dataset.tag, override=cache_dir),
            logger=logger, per_clade_consensus=True,
        )
        if len(panel) < 2:
            return SegmentScan(segment, False, False, 0, "single-clade panel")

        seg_dir.mkdir(parents=True, exist_ok=True)
        collection = seg_dir / "collection"
        if collection.exists():
            shutil.rmtree(collection)
        collection.mkdir(parents=True)
        rows = []
        for p in panel:
            shutil.copy(p, collection / p.name)
            rows.append((strip_sequence_extension(p.name), _clade_of_header(p), "consensus"))
        query = seg_dir / "query.fasta"
        with open(query, "w") as fo:
            write_fasta_record(fo, segment, seq)
        rows.append((segment, "query", "query"))
        write_lineage_map(seg_dir / LINEAGES_TSV, rows)
        lineage_map = lineage_map_from_rows(rows)

        window, step, _sel = window_params(len(seq))
        msa = seg_dir / "panel.msa.fasta"
        build_msa(MsaParams(query=query, collection=collection, output=msa, aligner=aligner), logger)
        run_recomb(RecombParams(msa=msa, output=seg_dir, query=segment,
                                window_size=window, window_step=step, organism=segment,
                                methods=DEFAULT_METHODS, lineage_map=lineage_map), logger)
    except Exception as exc:  # noqa: BLE001 - a per-segment scan failure is non-fatal
        logger.info("[%s] intragenic scan failed (%s); not scanned.", segment, exc)
        return SegmentScan(segment, False, False, 0, f"scan failed: {exc}")

    n_regions, recombinant = _summarize_regions(seg_dir / "recombination_regions.tsv")
    return SegmentScan(segment, True, recombinant, n_regions,
                       f"{n_regions} region(s)" if recombinant else "none")
```

- [ ] **Step 4: Run to verify they pass**

Run: `python -m pytest tests/unit/test_reassort_scan.py -q`
Expected: PASS (8 passed).

- [ ] **Step 5: Lint**

Run: `ruff check src/tessera/reassort/scan.py tests/unit/test_reassort_scan.py`
Expected: `All checks passed!`

- [ ] **Step 6: Commit**

```bash
git add src/tessera/reassort/scan.py tests/unit/test_reassort_scan.py
git commit -m "Add the intragenic per-segment scan module

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Wire the scan into `assign_segments`

**Files:**
- Modify: `src/tessera/reassort/assign.py`
- Modify: `src/tessera/reassort/__init__.py`
- Modify: `tests/unit/test_reassort_assign.py`

**Interfaces:**
- Consumes: `SegmentScan`, `require_aligner`, `scan_segment` from Task 1.
- Produces:
  - `ReassortmentResult` gains `scans: list[SegmentScan] = field(default_factory=list)`.
  - `assign_segments(query, *, dataset_overrides=None, ani_floor=DEFAULT_ANI_FLOOR, margin=DEFAULT_MARGIN, email=None, cache_dir=None, output=None, scan_segments=False, aligner="mafft", logger) -> ReassortmentResult`.
  - `_type_segment(...)` returns `(SegmentAssignment, candidates, universe, dataset_or_none)`.

- [ ] **Step 1: Write the failing test**

Add to `tests/unit/test_reassort_assign.py` (reuses the module-level `_DS`, `_patch`, `_write_query`, `LOG`):

```python
def test_scan_segments_scans_assigned_only(tmp_path, monkeypatch):
    # HA and NA both assign; scan_segments=True -> each assigned segment gets a SegmentScan;
    # scan_segment is stubbed so no aligner is needed. require_aligner is stubbed to a no-op.
    ha_tip = tmp_path / "HA_pool" / "strainA.fasta"
    na_tip = tmp_path / "NA_pool" / "strainB.fasta"
    for t in (ha_tip, na_tip):
        t.parent.mkdir(parents=True, exist_ok=True)
        t.write_text(">x\nACGT\n")

    def resolve(fasta, override, *, email, logger):
        return _DS("HA_ds") if "HA" in fasta.read_text() else _DS("NA_ds")

    _patch(monkeypatch, resolve=resolve,
           tips_by_path={"HA_ds": [ha_tip], "NA_ds": [na_tip]},
           ani_by_path={"HA_pool": {ha_tip: (99.0, 0.99)},
                        "NA_pool": {na_tip: (99.0, 0.99)}})
    monkeypatch.setattr(assign, "require_aligner", lambda aligner: None)
    seen = []

    def fake_scan(segment, seq, dataset, out_dir, *, aligner, cache_dir, logger):
        from tessera.reassort.scan import SegmentScan
        seen.append(segment)
        return SegmentScan(segment, True, segment == "HA", 1 if segment == "HA" else 0,
                           "1 region(s)" if segment == "HA" else "none")
    monkeypatch.setattr(assign, "scan_segment", fake_scan)

    q = _write_query(tmp_path, [("HA", "HAxx"), ("NA", "NAyy")])
    result = assign_segments(q, output=tmp_path / "out", scan_segments=True, logger=LOG)
    assert seen == ["HA", "NA"]
    scans = {s.segment: s for s in result.scans}
    assert scans["HA"].recombinant is True and scans["NA"].recombinant is False


def test_scan_segments_marks_unassigned(tmp_path, monkeypatch):
    # HA unassigned (no eligible tip); scan_segments=True -> HA gets a not-scanned row, no scan call.
    ha_tip = tmp_path / "HA_pool" / "strainA.fasta"
    na_tip = tmp_path / "NA_pool" / "strainB.fasta"
    for t in (ha_tip, na_tip):
        t.parent.mkdir(parents=True, exist_ok=True)
        t.write_text(">x\nACGT\n")

    def resolve(fasta, override, *, email, logger):
        return _DS("HA_ds") if "HA" in fasta.read_text() else _DS("NA_ds")

    _patch(monkeypatch, resolve=resolve,
           tips_by_path={"HA_ds": [ha_tip], "NA_ds": [na_tip]},
           ani_by_path={"HA_pool": {ha_tip: (10.0, 0.99)},   # below ani_floor -> unassigned
                        "NA_pool": {na_tip: (99.0, 0.99)}})
    monkeypatch.setattr(assign, "require_aligner", lambda aligner: None)

    def fake_scan(segment, seq, dataset, out_dir, *, aligner, cache_dir, logger):
        from tessera.reassort.scan import SegmentScan
        return SegmentScan(segment, True, False, 0, "none")
    monkeypatch.setattr(assign, "scan_segment", fake_scan)

    q = _write_query(tmp_path, [("HA", "HAxx"), ("NA", "NAyy")])
    result = assign_segments(q, output=tmp_path / "out", scan_segments=True, logger=LOG)
    scans = {s.segment: s for s in result.scans}
    assert scans["HA"].scanned is False and scans["HA"].note == "unassigned"
    assert scans["NA"].scanned is True


def test_scan_segments_needs_output(tmp_path, monkeypatch):
    from tessera.core.errors import UserInputError
    monkeypatch.setattr(assign, "skani_available", lambda: True)
    monkeypatch.setattr(assign, "require_aligner", lambda aligner: None)
    q = _write_query(tmp_path, [("HA", "AAAA")])
    with pytest.raises(UserInputError, match="output"):
        assign_segments(q, scan_segments=True, logger=LOG)
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tests/unit/test_reassort_assign.py -q`
Expected: FAIL -- `assign_segments` has no `scan_segments`/`output` kwargs yet (`TypeError`).

- [ ] **Step 3: Update `_type_segment` to return the dataset**

In `src/tessera/reassort/assign.py`, change the two non-exception returns of `_type_segment` to include the dataset, and the exception return to include `None`. The full updated function:

```python
def _type_segment(seg, seq, overrides, ani_floor, email, cache_dir, tmp, logger):
    """Type one segment. Returns ``(SegmentAssignment, candidates, universe, dataset)`` where
    ``candidates`` is ``[(strain, ani)]`` best-first (empty if unassigned), ``universe`` is every
    strain in the dataset (empty if unassigned), and ``dataset`` is the resolved Nextclade dataset
    (``None`` if resolution/typing failed). Never raises for this segment."""
    seg_fasta = Path(tmp) / f"{strip_sequence_extension(seg)}.fasta"
    with open(seg_fasta, "w") as fo:
        write_fasta_record(fo, seg, seq)
    try:
        dataset = resolve_dataset(seg_fasta, overrides.get(seg), email=email, logger=logger)
        tips = build_pool(
            dataset, cache_dir=nextclade_cache(dataset.path, dataset.tag, override=cache_dir),
            logger=logger,
        )
        ani = skani_query_ani(seg_fasta, tips, logger)
    except Exception as exc:  # noqa: BLE001 - a per-segment failure is non-fatal
        logger.info("[%s] could not type segment (%s); unassigned.", seg, exc)
        return SegmentAssignment(seg, "?", None, None, 0.0, "unassigned"), [], set(), None

    eligible = [t for t in tips
                if ani.get(t, (0.0, 0.0))[0] >= ani_floor and ani.get(t, (0.0, 0.0))[1] >= MIN_AF]
    ranked = sorted(eligible, key=lambda t: ani[t][0], reverse=True)
    if not ranked:
        return SegmentAssignment(seg, dataset.path, None, None, 0.0, "unassigned"), [], set(), dataset

    candidates = [(strip_sequence_extension(t.name), ani[t][0]) for t in ranked[:TOP_K]]
    universe = {strip_sequence_extension(t.name) for t in tips}
    best = ranked[0]
    logger.info("[%s] %s -> %s (%s), ANI %.1f%%", seg, dataset.path,
                strip_sequence_extension(best.name), _clade_of_tip(best), ani[best][0])
    assignment = SegmentAssignment(seg, dataset.path, strip_sequence_extension(best.name),
                                   _clade_of_tip(best), ani[best][0], "assigned")
    return assignment, candidates, universe, dataset
```

- [ ] **Step 4: Add the scan wiring to `assign_segments`**

In `src/tessera/reassort/assign.py`, add the imports at the top (after the `constellation` import):

```python
from .scan import SegmentScan, require_aligner, scan_segment
```

Add `scans` to the result dataclass:

```python
@dataclass
class ReassortmentResult:
    segments: list[SegmentAssignment] = field(default_factory=list)
    verdict: str = "undetermined"
    groups: list[ParentGroup] = field(default_factory=list)
    pair_notes: list[str] = field(default_factory=list)
    scans: list[SegmentScan] = field(default_factory=list)
```

Replace the `assign_segments` signature and body from the `if not skani_available()` check through the `return result` with:

```python
def assign_segments(
    query: Path,
    *,
    dataset_overrides: dict[str, str] | None = None,
    ani_floor: float = DEFAULT_ANI_FLOOR,
    margin: float = DEFAULT_MARGIN,
    email: str | None = None,
    cache_dir: Path | None = None,
    output: Path | None = None,
    scan_segments: bool = False,
    aligner: str = "mafft",
    logger: logging.Logger,
) -> ReassortmentResult:
    """Type each segment record of ``query`` against its Nextclade dataset and call
    reassortment. Each record's first-token name is the segment name (e.g. ``HA``). With
    ``scan_segments`` each assigned segment is additionally scanned for intragenic
    recombination (needs an ``aligner`` and an ``output`` directory)."""
    if not skani_available():
        raise UserInputError(
            "reassort needs skani. Install with: conda install -c bioconda skani skder"
        )
    if scan_segments:
        if output is None:
            raise UserInputError("--scan-segments needs an output directory (-o/--output).")
        require_aligner(aligner)
    overrides = dataset_overrides or {}
    records = read_fasta(query)
    if not records:
        raise UserInputError(f"Query FASTA {query} has no sequence.")
    names = [name for name, _ in records if name]
    dupes = sorted({n for n in names if names.count(n) > 1})
    if dupes:
        raise UserInputError(
            f"Query has duplicate segment names ({', '.join(dupes)}); each record must be a "
            f"distinct segment with a unique name."
        )

    result = ReassortmentResult()
    candidates: dict[str, list[tuple[str, float]]] = {}
    universes: dict[str, set[str]] = {}
    to_scan: dict[str, tuple[str, object]] = {}  # segment -> (seq, dataset) for assigned segments
    with tempfile.TemporaryDirectory() as tmp:
        for name, seq in records:
            seg = name or f"segment_{len(result.segments) + 1}"
            assignment, cand, universe, dataset = _type_segment(
                seg, seq, overrides, ani_floor, email, cache_dir, tmp, logger)
            result.segments.append(assignment)
            if assignment.status == "assigned":
                candidates[seg] = cand
                universes[seg] = universe
                to_scan[seg] = (seq, dataset)

    call = call_constellation(candidates, universes, margin)
    result.verdict = call.verdict
    result.groups = call.groups
    result.pair_notes = call.pair_notes

    if scan_segments:
        for s in result.segments:
            if s.status == "assigned":
                seq, dataset = to_scan[s.segment]
                result.scans.append(scan_segment(
                    s.segment, seq, dataset, output,
                    aligner=aligner, cache_dir=cache_dir, logger=logger))
            else:
                result.scans.append(SegmentScan(s.segment, False, False, 0, "unassigned"))
    return result
```

- [ ] **Step 5: Update `__init__.py`**

Replace `src/tessera/reassort/__init__.py` with:

```python
"""Per-segment reassortment detection (the ``tessera reassort`` command)."""

from .assign import ReassortmentResult, SegmentAssignment, assign_segments
from .constellation import Constellation, ParentGroup, call_constellation
from .scan import SegmentScan, scan_segment

__all__ = [
    "Constellation",
    "ParentGroup",
    "ReassortmentResult",
    "SegmentAssignment",
    "SegmentScan",
    "assign_segments",
    "call_constellation",
    "scan_segment",
]
```

- [ ] **Step 6: Run the tests**

Run: `python -m pytest tests/unit/test_reassort_assign.py tests/unit/test_reassort_scan.py -q`
Expected: PASS -- the 5 existing assign tests plus the 3 new (8 total) and the 8 scan tests from Task 1.

- [ ] **Step 7: Lint + full fast suite**

Run: `ruff check src/tessera/reassort tests/unit/test_reassort_assign.py`
Expected: `All checks passed!`
Run: `python -m pytest -m "not requires_binary" -q`
Expected: PASS, no failures.

- [ ] **Step 8: Commit**

```bash
git add src/tessera/reassort/assign.py src/tessera/reassort/__init__.py \
        tests/unit/test_reassort_assign.py
git commit -m "Wire the intragenic scan into assign_segments (--scan-segments)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: CLI `--scan-segments`/`--aligner`, rollup output, docs, e2e

**Files:**
- Modify: `src/tessera/cli/cmd_reassort.py`
- Modify: `tests/unit/test_cli_reassort.py`
- Modify: `docs/reference-panels.md`
- Modify: `validation/attribution-results.md`

**Interfaces:**
- Consumes: `assign_segments(..., output=..., scan_segments=..., aligner=...)`, `ReassortmentResult.scans` (`list[SegmentScan]` with `.segment`/`.scanned`/`.recombinant`/`.n_regions`/`.note`).

- [ ] **Step 1: Write the failing CLI tests**

Add to `tests/unit/test_cli_reassort.py`:

```python
def test_reassort_exposes_scan_segments_option():
    names = {opt for p in _reassort_command().params for opt in getattr(p, "opts", [])}
    assert "--scan-segments" in names


def test_reassort_exposes_aligner_option():
    names = {opt for p in _reassort_command().params for opt in getattr(p, "opts", [])}
    assert "--aligner" in names
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tests/unit/test_cli_reassort.py -q`
Expected: FAIL (the two new options are not yet parameters).

- [ ] **Step 3: Add the options, pass-through, and rollup output**

In `src/tessera/cli/cmd_reassort.py`, add the two options after `margin` (before `) -> None:`):

```python
    scan_segments: bool = typer.Option(
        False, "--scan-segments",
        help="Also scan each assigned segment for intragenic recombination (needs an aligner).",
    ),
    aligner: str = typer.Option(
        "mafft", "--aligner", help="Aligner backend for --scan-segments."),
```

Pass them through the `assign_segments` call (replace that call):

```python
        result = assign_segments(
            query, dataset_overrides=overrides,
            email=email or os.environ.get("NCBI_EMAIL"),
            ani_floor=ani_floor, margin=margin, output=output,
            scan_segments=scan_segments, aligner=aligner, logger=logger,
        )
```

After the existing `logger.info("Wrote %s and %s", tsv, ctsv)` line, add the segment-scan rollup:

```python
        if result.scans:
            stsv = output / "segment_scan.tsv"
            with open(stsv, "w") as fo:
                fo.write("segment\tintragenic_recombination\tn_regions\tnote\n")
                for sc in result.scans:
                    flag = "yes" if sc.recombinant else ("no" if sc.scanned else "n/a")
                    fo.write(f"{sc.segment}\t{flag}\t{sc.n_regions}\t{sc.note}\n")
            rollup = " | ".join(f"{sc.segment}: {sc.note}" for sc in result.scans)
            logger.info("Intragenic scan: %s", rollup)
            logger.info("Wrote %s", stsv)
```

- [ ] **Step 4: Run the CLI tests**

Run: `python -m pytest tests/unit/test_cli_reassort.py -q`
Expected: PASS (5 passed).

- [ ] **Step 5: Full fast suite + lint**

Run: `python -m pytest -m "not requires_binary" -q`
Expected: PASS.
Run: `ruff check src tests validation`
Expected: `All checks passed!`

- [ ] **Step 6: Update `docs/reference-panels.md`**

In the `## Detect reassortment across segments (reassort)` section, after the sentence ending "...not intragenic recombination within a segment." append a new paragraph:

```markdown
Pass `--scan-segments` to additionally look *inside* each segment. After assignment, each assigned
segment is aligned to its own per-clade-consensus panel and run through the ordinary recombination
scan, writing the full `run_recomb` output to `out/<segment>/` and a rollup to `out/segment_scan.tsv`
(`segment`, `intragenic_recombination`, `n_regions`, `note`). This asks a different question from the
reassortment verdict: not which parent each whole segment came from, but whether a single segment is
itself a within-segment mosaic of two lineages. It needs an aligner (`--aligner`, default `mafft`);
a missing or unknown aligner is reported up front. A segment whose panel has fewer than two clades,
or that is unassigned, is reported as not scanned rather than dropped.
```

- [ ] **Step 7: Opt-in e2e and record**

Run (aligner env; not CI):
```bash
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
export NCBI_EMAIL="andreas.sjodin@gmail.com"
```
Build a two-record query where HA is an intragenic mosaic of two H3N2 HA clades and NA is a single
clade (reuse the harness's `_segment_consensus` to fetch two HA clade consensuses, splice them at
mid-length into one HA record, and take one NA clade consensus), then run:
```bash
python -m tessera reassort -q <mosaic.fasta> -o out/ \
    --dataset HA=nextstrain/flu/h3n2/ha/EPI1857216 --dataset NA=nextstrain/flu/h3n2/na/EPI1857215 \
    --scan-segments
```
Expected: `out/segment_scan.tsv` shows HA `intragenic_recombination=yes` (>= 1 region) and NA `no`.
Record the measured result in `validation/attribution-results.md` (a short paragraph after the
`reassort_flu` bullet): the command, the HA/NA scan outcome, and that it is an opt-in probe. State the
measured result faithfully; if HA is not detected, record what was measured and stop for review
rather than retuning.

- [ ] **Step 8: Commit**

```bash
git add src/tessera/cli/cmd_reassort.py tests/unit/test_cli_reassort.py \
        docs/reference-panels.md validation/attribution-results.md
git commit -m "reassort CLI: --scan-segments / --aligner, rollup, docs, e2e

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review

**Spec coverage:**
- `--scan-segments` flag on reassort -> Task 3 CLI; wired in Task 2.
- `scan.py` per-segment build_pool(consensus) + build_msa + run_recomb -> Task 1 `scan_segment`.
- Full run_recomb report per segment in `out/<segment>/` -> Task 1 (run_recomb `output=seg_dir`).
- Reuse resolved dataset (no second resolve) -> Task 2 `_type_segment` returns the dataset; `to_scan` stores `(seq, dataset)`.
- Aligner required only for the flag; up-front error -> Task 1 `require_aligner`; called in Task 2 before typing; CLI passes `--aligner`.
- Non-fatal per-segment scan -> Task 1 `scan_segment` try/except; test `test_scan_segment_failure_is_non_fatal`.
- Single-clade panel skip -> Task 1; test `test_scan_segment_single_clade_panel_is_skipped`.
- Cache separation for the consensus panel -> handled by `build_pool(per_clade_consensus=True)` itself (it appends `_consensus` to the cache dir); the plan relies on that, no extra work.
- `segment_scan.tsv` + rollup -> Task 3.
- Unassigned segment reported not-scanned -> Task 2; test `test_scan_segments_marks_unassigned`.
- Verdict vocabulary + `ReassortmentResult` back-compat -> Task 2 adds `scans` with a default; harness `_score_frontier_reassortant` and `test_frontier_gate.py` untouched.
- Docs + faithful e2e -> Task 3 Steps 6-7.

**Placeholder scan:** every code step shows complete code; the only prose-described step is the e2e mosaic construction (Task 3 Step 7), which is an opt-in manual run, not code to commit.

**Type consistency:** `SegmentScan(segment, scanned, recombinant, n_regions, note)` is identical across Task 1 (definition), Task 2 (construction of the unassigned row and the stub), and Task 3 (consumes `.segment/.scanned/.recombinant/.n_regions/.note`). `scan_segment(segment, seq, dataset, out_dir, *, aligner, cache_dir, logger)` matches between Task 1 and its Task 2 call site. `_type_segment` returns a 4-tuple in Task 2's definition and is unpacked as a 4-tuple at its one call site. `assign_segments` gains `output`/`scan_segments`/`aligner` consistently in Task 2 and Task 3.

**Note for the implementer:** the existing `tests/unit/test_reassort_assign.py` unpacks `assign_segments` results only (never `_type_segment` directly), so changing `_type_segment`'s arity is safe. `build_pool` is called with the keyword `per_clade_consensus=True`; the stub in Task 1's tests must accept that keyword (the provided stubs do).
```
