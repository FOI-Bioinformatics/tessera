# Post-1.2.0 Audit, Plan A: Data Safety and Alignment Correctness -- Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the audit's data-loss and wrong-alignment defects (spec items A1-A10) so that no command deletes a user's genomes, every genome that is downloaded or supplied is aligned in the right orientation and register, and the reassort alignment-fraction filter does what its name says.

**Architecture:** Nine independent fixes, each a test that fails on `main` followed by the smallest change that passes it. None touches a region caller, a default threshold of the scan, or the output schema of `recomb`; they change what reaches the scan (alignments, panels) and one reassort filter. Two fixes restructure existing code rather than patch it: the `fill-references` round loop (Task 6) and the hand-off from aligner adapters to the MAF converter (Task 7).

**Tech Stack:** Python 3.11+, numpy / pandas / biopython / typer (no new dependency), pytest, ruff, mypy. External binaries only in `requires_binary` tests: mafft 7.526, minimap2 2.31.

**Spec:** `docs/superpowers/specs/2026-10-01-post-1.2.0-audit-design.md` (section "Plan A"). Read it first; this plan argues from it and names two places where verification showed the spec's design had to be narrowed (Task 6, Task 9).

**How this plan was produced.** Every test and every change below was applied to a scratch worktree of `main` (`457bdfb`) one task at a time. For each task the tests were added first and run (the failures listed under "Expected" are the ones observed), then the change was applied and the task's tests, the whole fast suite, `ruff` and `mypy` were run. The final tree passed `621 passed, 5 deselected` on the fast suite and `5 passed` on the `requires_binary` tests (mafft, minimap2). The `diff` blocks are the exact per-task diffs from that worktree and apply in task order with `git apply`. Task 8 is the exception to "verified against the real tool"; it says so.

**Applying a code block.** New tests are given as Python to append, insert or create; separate an appended or inserted block from its neighbours with two blank lines. **Applying a `diff` block.** Either make the edit by hand, or save the block to a file and run `git apply <file>` from the repository root. Each block is relative to the state the previous task left, so apply tasks in order (Tasks 2-9 are otherwise independent of each other except where "Interfaces" says so).

## Global Constraints

- No new runtime dependency (CLAUDE.md). Nothing in this plan adds one.
- Modest scientific language in code, docstrings, messages and docs; reported numbers must be faithful.
- A new behaviour needs a test that fails without the change (CONTRIBUTING.md). Tests in this plan that pass before the change are labelled as guards.
- "Could not test" must never read as "tested, found nothing": a genome that cannot be placed is a visible all-gap row plus a warning, not a missing row.
- Tests must not touch the network. Anything that could reach NCBI, BLAST or Nextclade is stubbed in the test itself (`blast_subsequence`, `efetch_fasta`, `efetch_available`, `skani_query_ani`, `build_pool`, `resolve_dataset`).
- Aligner binaries: **append** the env to `PATH` (`PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"`), never prepend it -- prepending swaps the Python interpreter. Commands below that need a binary show the prefix.
- `caplog` only sees records that reach the root logger, and the CLI tests set `propagate = False` on the `tessera` logger. A test that asserts on log text must pass a logger outside the `tessera` namespace.
- Ruff: line length 100, rules `E/F/I/UP/B`. `mypy src` is blocking in CI.
- Branch before committing (never commit on `main`). Commit messages end with the trailer `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`.
- Claude cannot merge or close PRs in this repository; open the PR, get CI green, and hand the merge to the maintainer.
- This plan changes no caller, default or region-calling code, so the hybrid and specificity harnesses are not a gate here. Task 3 changes a reassort filter; Task 10 lists the optional reassort benchmark.

## Review Focus

Inputs the spec implies but does not spell out, most likely to bite first. Each has a test in the task that owns the code.

1. **`--curate` on a supplied panel as close to the query as its own lineage.** Every genome but the backbone is classed as a sibling and one reference is left, on which detection cannot run. Expected: a warning that says so and suggests running without `--curate`. Task 6, `test_curate_warns_when_it_leaves_too_few_references`.
2. **`--reference` given with its file extension (`MYREF.fasta`)**, which `tessera msa` accepts. Expected: it is still the curation backbone. Task 6, `test_curate_accepts_a_reference_given_with_its_extension`.
3. **`find-references --download` naming a directory that does not exist yet, with `--curate`.** Expected: works; a downloaded sibling is dropped; the collection is untouched. Task 5, `test_curate_after_download_into_a_new_directory`.
4. **Two segment records whose names sanitise alike (`seg/1`, `seg_1`) in one `reassort --scan-segments` run.** Expected: two scan directories, neither overwriting the other. Task 4, `test_scan_segments_gives_alike_named_segments_separate_directories` and `test_unique_scan_dirs_separates_names_that_sanitise_alike`.
5. **progressiveMauve output that is not pairwise** (the query missing from the XMFA). Expected: an error naming the projection file, not an MSA that is silently one row short. Task 8, `test_progressivemauve_rejects_a_projection_that_is_not_pairwise`.

## File Structure

No file is created under `src/`. Responsibilities after this plan:

| File | Change | Responsibility |
|---|---|---|
| `src/tessera/core/io.py` | modify | FASTA reading (whitespace-insensitive), genome staging (cleaned copies), `safe_filename_stem` |
| `src/tessera/aligners/mafft.py` | modify | MAFFT command line (`--adjustdirection`) |
| `src/tessera/aligners/sibeliaz.py` | modify | hands the backbone's contig order and genome list to the MAF converter; rejects nameless records |
| `src/tessera/aligners/cactus.py` | modify | hands the genome list to the MAF converter |
| `src/tessera/aligners/progressivemauve.py` | modify | rows named by staged label, taken by position |
| `src/tessera/converters/maf_to_fasta.py` | modify | backbone layout from the caller; all-gap rows for unplaced genomes; truncated-block error |
| `src/tessera/converters/xmfa_to_fasta.py`, `mafft_merge.py` | modify | errors that name the file |
| `src/tessera/discover/panel.py` | modify | `curate_collection_dir(protect=)` and the two "kept" roles |
| `src/tessera/discover/run.py` | modify | find-references curates only its own downloads |
| `src/tessera/discover/iterate.py` | modify | the fill loop's round order and final build |
| `src/tessera/reassort/assign.py`, `scan.py` | modify | `MIN_AF` in percent; unique, contained scan directories |
| `src/tessera/cli/cmd_fill_references.py`, `cmd_find_references.py` | modify | `--curate` help text |
| `tests/integration/test_mafft_strand.py` | create | real-mafft strand tests |
| `tests/integration/test_whitespace_genomes.py` | create | real mafft / minimap2 on a reference with trailing spaces |
| `tests/unit/*.py` (10 files) | modify | one or more tests per fix, listed per task |
| `docs/aligners.md`, `docs/reference-panels.md`, `CHANGELOG.md` | modify | Task 10 |

---

### Task 1: Branch and baseline

**Files:** none.

**Interfaces:**
- Consumes: `main` at or after `457bdfb`.
- Produces: branch `fix-audit-data-safety`.

- [ ] **Step 1: Create the branch from an up-to-date main**

```bash
git checkout main && git pull --ff-only
git checkout -b fix-audit-data-safety
```

- [ ] **Step 2: Record the baseline**

Run: `ruff check src tests validation && pytest -m "not requires_binary" -q && mypy src`
Expected: `All checks passed!`, `586 passed, 1 deselected`, `Success: no issues found in 79 source files`.

If the pass count differs, stop and find out why before changing anything: the expected counts in later tasks are relative to this baseline.

- [ ] **Step 3: Confirm the aligners used by the binary tests are reachable**

Run: `PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin" sh -c 'mafft --version; minimap2 --version'`
Expected: a mafft version (`v7.526` here) and a minimap2 version (`2.31` here). If either is missing, Tasks 2 and 9 can still be implemented, but their `requires_binary` tests will skip; say so in the PR rather than reporting them as passed.

---
### Task 2: A1 -- MAFFT backend reorients reverse-strand input

**Spec item(s):** A1

**Files:**
- Modify: `src/tessera/aligners/mafft.py` (module docstring; the `add_genome` command at lines 76-83)
- Test: `tests/unit/test_aligner_backends.py` (append), `tests/integration/test_mafft_strand.py` (create; `requires_binary`)

**Interfaces:**
- Consumes: nothing from other tasks.
- Produces: no new names. The MAFFT command line gains `--adjustdirection`, placed before `--addfragments`. `merge_added_fragments` is unchanged: it is positional, so the `_R_<name>` records MAFFT emits for reversed sequences are merged like any other.

The adapter adds each genome onto the backbone with `mafft --keeplength --addfragments`. MAFFT aligns a sequence in the orientation it is given unless told otherwise, so a genome (or one contig of a draft) on the opposite strand comes out as noise that still fills the row. The fix is one flag. The integration test needs the real binary and is the test that proves the behaviour; the unit test only pins the command line so CI (which has no aligner) still guards the flag.

- [ ] **Step 1: Write the failing tests**

Create `tests/integration/test_mafft_strand.py`:

```python
"""MAFFT backend on reverse-strand input (needs the mafft binary).

An assembly carries no promise about strand: a whole genome, or one contig of a draft,
may be the reverse complement of the backbone. The row must still be a real alignment.
"""

from __future__ import annotations

import logging
import random
import shutil
from pathlib import Path

import pytest

from tessera.aligners.base import AlignParams
from tessera.aligners.mafft import MafftAligner
from tessera.core.io import read_fasta

_LOG = logging.getLogger("tessera.test")

pytestmark = pytest.mark.requires_binary

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _genomes() -> tuple[str, str]:
    """A 6 kb random reference and a query 3 % diverged from it."""
    rng = random.Random(1)
    ref = "".join(rng.choice("ACGT") for _ in range(6000))
    qry = "".join(
        rng.choice([b for b in "ACGT" if b != c]) if rng.random() < 0.03 else c for c in ref
    )
    return ref, qry


def _identity_by_third(ref_row: str, qry_row: str) -> list[float]:
    out = []
    for lo in range(0, len(ref_row), 2000):
        pairs = [
            (a, b) for a, b in zip(ref_row[lo:lo + 2000].upper(), qry_row[lo:lo + 2000].upper(),
                                   strict=True)
            if a in "ACGT" and b in "ACGT"
        ]
        out.append(sum(a == b for a, b in pairs) / len(pairs) if pairs else 0.0)
    return out


def _align(tmp_path: Path, query_fasta: str) -> tuple[str, str]:
    ref, _ = _genomes()
    (tmp_path / "ref.fasta").write_text(f">ref\n{ref}\n")
    (tmp_path / "qry.fasta").write_text(query_fasta)
    result = MafftAligner().align(
        [tmp_path / "ref.fasta", tmp_path / "qry.fasta"], tmp_path / "ref.fasta",
        tmp_path / "out", AlignParams(threads=1), _LOG,
    )
    msa = dict(read_fasta(result.msa_fasta))
    return msa["ref"], msa["qry"]


@pytest.fixture(autouse=True)
def _needs_mafft() -> None:
    if shutil.which("mafft") is None:
        pytest.skip("mafft not installed")


def test_reverse_complemented_query_aligns(tmp_path: Path) -> None:
    _, qry = _genomes()
    ref_row, qry_row = _align(tmp_path, f">q\n{_revcomp(qry)}\n")
    assert len(ref_row) == len(qry_row) == 6000
    assert all(identity > 0.95 for identity in _identity_by_third(ref_row, qry_row))


def test_draft_query_with_one_reversed_contig_aligns(tmp_path: Path) -> None:
    _, qry = _genomes()
    contigs = f">c1\n{qry[:2000]}\n>c2\n{_revcomp(qry[2000:4000])}\n>c3\n{qry[4000:]}\n"
    ref_row, qry_row = _align(tmp_path, contigs)
    assert len(ref_row) == len(qry_row) == 6000
    assert all(identity > 0.95 for identity in _identity_by_third(ref_row, qry_row))
```

Append to the end of `tests/unit/test_aligner_backends.py`:

```python
def test_mafft_adjusts_direction_of_added_sequences(monkeypatch, tmp_path: Path) -> None:
    """A draft contig or a whole genome may be on the opposite strand to the backbone.
    Without --adjustdirection MAFFT aligns it as given and the row is noise."""
    captured: dict[str, list[str]] = {}

    def fake_run(caps, cmd, **kw):
        captured["cmd"] = [str(c) for c in cmd]
        raise RuntimeError("stop after capture")

    monkeypatch.setattr(mafft_mod, "run_tool", fake_run)
    ref, qry = _two_genomes(tmp_path)
    with pytest.raises(RuntimeError):
        mafft_mod.MafftAligner().align([ref, qry], ref, tmp_path / "out",
                                       AlignParams(threads=1), _LOG)

    cmd = captured["cmd"]
    assert "--adjustdirection" in cmd
    # The option must precede --addfragments, whose two arguments are positional.
    assert cmd.index("--adjustdirection") < cmd.index("--addfragments")


def test_merge_added_fragments_ignores_reversed_name_prefix(tmp_path: Path) -> None:
    # `mafft --adjustdirection` renames a record it reverse-complemented to `_R_<name>`.
    # The merge is positional (first record = reference, the rest = contigs), so the
    # renamed contig must still be merged.
    aligned = tmp_path / "a.fasta"
    aligned.write_text(
        ">ref\nACGTACGT\n"
        ">contig1\nAC------\n"
        ">_R_contig2\n----ACGT\n"
    )
    _, merged = merge_added_fragments(aligned)
    assert merged == "AC--ACGT"
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin" pytest tests/unit/test_aligner_backends.py tests/integration/test_mafft_strand.py -q`

Expected: `3 failed, 10 passed`. The failures are:

```
tests/unit/test_aligner_backends.py::test_mafft_adjusts_direction_of_added_sequences
tests/integration/test_mafft_strand.py::test_reverse_complemented_query_aligns
tests/integration/test_mafft_strand.py::test_draft_query_with_one_reversed_contig_aligns
```

They fail because the unit test because `--adjustdirection` is not in the command; the two integration tests because identity over the reversed part is about 0.4, not above 0.95.

- [ ] **Step 3: Implement**

`src/tessera/aligners/mafft.py`:

````diff
diff --git a/src/tessera/aligners/mafft.py b/src/tessera/aligners/mafft.py
index 28aa705..9311780 100644
--- a/src/tessera/aligners/mafft.py
+++ b/src/tessera/aligners/mafft.py
@@ -7,8 +7,10 @@ contract each genome is added onto the backbone with
 ``mafft --addfragments <genome> --keeplength <reference>``: ``--keeplength``
 keeps the output in reference coordinates (insertions relative to the backbone
 are dropped) and ``--addfragments`` is designed for fragmented assemblies, so a
-multi-contig query is handled cleanly. A genome's contigs are then merged into a
-single reference-anchored row.
+multi-contig query is handled cleanly. ``--adjustdirection`` lets MAFFT reverse-
+complement an added sequence that is on the opposite strand to the backbone (it renames
+such a record ``_R_<name>``; the merge below is positional and ignores names). A genome's
+contigs are then merged into a single reference-anchored row.
 """
 
 from __future__ import annotations
@@ -76,8 +78,8 @@ class MafftAligner(Aligner):
             aligned = out_dir / f"{genome.stem}.aln.fasta"
             run_tool(
                 self.capabilities,
-                ["mafft", "--thread", threads, "--keeplength", *tuning,
-                 "--addfragments", str(genome.resolve()), str(ref_fasta.resolve())],
+                ["mafft", "--thread", threads, "--keeplength", "--adjustdirection",
+                 *tuning, "--addfragments", str(genome.resolve()), str(ref_fasta.resolve())],
                 logger=logger,
                 log_prefix=f"mafft:{genome.stem}",
                 stdout_path=aligned,
````

- [ ] **Step 4: Run the tests to verify they pass**

Run: `PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin" pytest tests/unit/test_aligner_backends.py tests/integration/test_mafft_strand.py -q`

Expected: all pass (the Step 2 command now gives `13 passed`).

- [ ] **Step 5: Run the fast gate**

Run: `ruff check src tests validation && pytest -m "not requires_binary" -q && mypy src`

Expected: `All checks passed!`, `588 passed, 3 deselected`, `Success: no issues found in 79 source files`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/aligners/mafft.py tests/integration/test_mafft_strand.py tests/unit/test_aligner_backends.py
git commit -m "$(cat <<'EOF'
Reorient reverse-strand genomes in the MAFFT backend

A genome or contig on the opposite strand to the backbone was aligned as given and
came out at chance-level identity. Run MAFFT with --adjustdirection.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

---
### Task 3: A2 -- reassort alignment-fraction filter in skani's unit

**Spec item(s):** A2

**Files:**
- Modify: `src/tessera/reassort/assign.py:28` (`MIN_AF`)
- Test: `tests/unit/test_reassort_assign.py` (mock values move to percent; one test appended)

**Interfaces:**
- Consumes: `tessera.discover.panel.skani_query_ani(query_fasta, refs, logger) -> dict[Path, tuple[float, float]]`, which returns `(ANI, Align_fraction_query)` both in percent (0-100). `tests/unit/test_panel.py::test_skani_query_ani_parses_and_marks_absent` already pins that scale.
- Produces: `tessera.reassort.assign.MIN_AF: float = 50.0` (percent).

`MIN_AF` was written as a fraction (0.5) but is compared with a percentage, so it never excluded anything. The existing tests mock `skani_query_ani` with fractions (`0.10`, `0.99`), which is why they passed: they pin the mock's unit, not skani's. Move every mock in this file to percent first, then add a test built from real skani output (a full-length tip at ANI 98.76 / AF 100.00 and a tip covering a fifth of the segment at ANI 99.29 / AF 20.26).

`Align_fraction_query` is the share of the *query segment* that aligned, so a short but complete segment still scores near 100; the filter removes tips that match only part of the segment.

- [ ] **Step 1: Write the failing tests**

Edit `tests/unit/test_reassort_assign.py` (existing lines change, and one test is added at the end):

````diff
diff --git a/tests/unit/test_reassort_assign.py b/tests/unit/test_reassort_assign.py
index 4289bec..bc488be 100644
--- a/tests/unit/test_reassort_assign.py
+++ b/tests/unit/test_reassort_assign.py
@@ -51,8 +51,8 @@ def test_low_af_tip_is_dropped_from_ranking(tmp_path, monkeypatch):
     _patch(monkeypatch,
            resolve=resolve,
            tips_by_path={"HA_ds": [ha_tip], "NA_ds": [na_tip]},
-           ani_by_path={"HA_pool": {ha_tip: (99.0, 0.10)},   # AF 0.10 < MIN_AF -> dropped
-                        "NA_pool": {na_tip: (99.0, 0.99)}})
+           ani_by_path={"HA_pool": {ha_tip: (99.0, 10.0)},   # AF 10 % < MIN_AF -> dropped
+                        "NA_pool": {na_tip: (99.0, 99.0)}})
     q = _write_query(tmp_path, [("HA", "HAxx"), ("NA", "NAyy")])
     result = assign_segments(q, logger=LOG)
     status = {s.segment: s.status for s in result.segments}
@@ -116,7 +116,7 @@ def test_skani_short_segment_is_non_fatal(tmp_path, monkeypatch):
     def skani(q, refs, logger):
         if refs[0].parent.name == "HA_pool":
             raise ToolExecutionError(["skani", "dist"], 1, "sequence too short")
-        return {na_tip: (99.0, 0.99)}
+        return {na_tip: (99.0, 99.0)}
 
     def build_pool(ds, *, cache_dir, logger):
         return [ha_tip] if ds.path == "HA_ds" else [na_tip]
@@ -152,7 +152,7 @@ def test_resolve_dataset_failure_is_non_fatal(tmp_path, monkeypatch):
     monkeypatch.setattr(assign, "nextclade_cache", lambda p, t, override=None: Path("/x"))
     monkeypatch.setattr(assign, "build_pool", lambda ds, *, cache_dir, logger: [na_tip])
     monkeypatch.setattr(assign, "skani_query_ani",
-                        lambda q, refs, logger: {na_tip: (99.0, 0.99)})
+                        lambda q, refs, logger: {na_tip: (99.0, 99.0)})
     monkeypatch.setattr(assign, "_clade_of_tip", lambda tip: "cladeX")
 
     q = _write_query(tmp_path, [("HA", "HAxx"), ("NA", "NAyy")])
@@ -208,8 +208,8 @@ def test_scan_segments_scans_assigned_only(tmp_path, monkeypatch):
 
     _patch(monkeypatch, resolve=resolve,
            tips_by_path={"HA_ds": [ha_tip], "NA_ds": [na_tip]},
-           ani_by_path={"HA_pool": {ha_tip: (99.0, 0.99)},
-                        "NA_pool": {na_tip: (99.0, 0.99)}})
+           ani_by_path={"HA_pool": {ha_tip: (99.0, 99.0)},
+                        "NA_pool": {na_tip: (99.0, 99.0)}})
     monkeypatch.setattr(assign, "require_aligner", lambda aligner: None)
     seen = []
 
@@ -240,8 +240,8 @@ def test_scan_segments_marks_unassigned(tmp_path, monkeypatch):
 
     _patch(monkeypatch, resolve=resolve,
            tips_by_path={"HA_ds": [ha_tip], "NA_ds": [na_tip]},
-           ani_by_path={"HA_pool": {ha_tip: (10.0, 0.99)},   # below ani_floor -> unassigned
-                        "NA_pool": {na_tip: (99.0, 0.99)}})
+           ani_by_path={"HA_pool": {ha_tip: (10.0, 99.0)},   # below ani_floor -> unassigned
+                        "NA_pool": {na_tip: (99.0, 99.0)}})
     monkeypatch.setattr(assign, "require_aligner", lambda aligner: None)
 
     def fake_scan(segment, seq, dataset, out_dir, *, aligner, cache_dir, logger):
@@ -291,3 +291,31 @@ def test_cap_candidates_still_bounds_the_tail():
 
 def test_cap_candidates_empty():
     assert cap_candidates([], margin=0.5, top_k=25) == []
+
+
+def test_alignment_fraction_filter_uses_skani_percent_scale(tmp_path, monkeypatch):
+    # skani reports Align_fraction_query in percent (0-100). A tip that aligns over a fifth
+    # of the segment at a higher ANI must not outrank a full-length match: the values below
+    # are what skani 0.3 printed for a 30 kb query against a full-length tip and a tip
+    # covering only its first 6 kb.
+    full = tmp_path / "HA_pool" / "full.fasta"
+    partial = tmp_path / "HA_pool" / "partial.fasta"
+    na_tip = tmp_path / "NA_pool" / "strainB.fasta"
+    for t in (full, partial, na_tip):
+        t.parent.mkdir(parents=True, exist_ok=True)
+        t.write_text(">x\nACGT\n")
+
+    def resolve(fasta, override, *, email, logger):
+        return _DS("HA_ds") if "HA" in fasta.read_text() else _DS("NA_ds")
+
+    _patch(monkeypatch,
+           resolve=resolve,
+           tips_by_path={"HA_ds": [full, partial], "NA_ds": [na_tip]},
+           ani_by_path={"HA_pool": {full: (98.76, 100.0), partial: (99.29, 20.26)},
+                        "NA_pool": {na_tip: (99.0, 99.0)}})
+    q = _write_query(tmp_path, [("HA", "HAxx"), ("NA", "NAyy")])
+    result = assign_segments(q, logger=LOG)
+    ha = next(s for s in result.segments if s.segment == "HA")
+    assert ha.status == "assigned"
+    assert ha.strain == "full"
+    assert ha.ani == pytest.approx(98.76)
````

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_reassort_assign.py -q`

Expected: `2 failed, 13 passed`. The failures are:

```
tests/unit/test_reassort_assign.py::test_low_af_tip_is_dropped_from_ranking
tests/unit/test_reassort_assign.py::test_alignment_fraction_filter_uses_skani_percent_scale
```

They fail because with the mocks in percent and `MIN_AF` still 0.5, a tip at AF 10 is not dropped, and the partial tip (AF 20.26) outranks the full-length one.

- [ ] **Step 3: Implement**

`src/tessera/reassort/assign.py`:

````diff
diff --git a/src/tessera/reassort/assign.py b/src/tessera/reassort/assign.py
index 7a57912..eb77d2b 100644
--- a/src/tessera/reassort/assign.py
+++ b/src/tessera/reassort/assign.py
@@ -25,7 +25,9 @@ from .constellation import DEFAULT_MARGIN, ParentGroup, call_constellation
 from .scan import SegmentScan, require_aligner, scan_segment
 
 DEFAULT_ANI_FLOOR = 80.0  # a segment below this ANI to every tip is left unassigned
-MIN_AF = 0.5              # a tip aligning over less than this fraction of the segment is ignored
+# A tip aligning over less than this share of the segment is ignored. In PERCENT (0-100),
+# the unit skani reports Align_fraction_query in and `skani_query_ani` returns.
+MIN_AF = 50.0
 TOP_K = 25                # internal cap on candidate strains kept per segment
 
 
````

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_reassort_assign.py tests/unit/test_cli_reassort.py -q`

Expected: all pass (the Step 2 command now gives `15 passed`).

- [ ] **Step 5: Run the fast gate**

Run: `ruff check src tests validation && pytest -m "not requires_binary" -q && mypy src`

Expected: `All checks passed!`, `589 passed, 3 deselected`, `Success: no issues found in 79 source files`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/reassort/assign.py tests/unit/test_reassort_assign.py
git commit -m "$(cat <<'EOF'
Compare the reassort alignment-fraction filter in percent

MIN_AF was a fraction compared against skani's percentage, so the filter never
fired and a partially aligning tip could become a segment's nearest strain.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

---
### Task 4: A7 -- segment names cannot leave or share a scan directory

**Spec item(s):** A7

**Files:**
- Modify: `src/tessera/core/io.py` (new `safe_filename_stem`, `import re`)
- Modify: `src/tessera/reassort/scan.py` (`unique_scan_dirs`, `scan_segment(..., dir_name=)`)
- Modify: `src/tessera/reassort/assign.py` (use the shared helper; pass `dir_name`)
- Test: `tests/unit/test_reassort_scan.py` (append), `tests/unit/test_reassort_assign.py` (two stubs gain a keyword; one test appended)

**Interfaces:**
- Consumes: nothing from other tasks.
- Produces:
  - `tessera.core.io.safe_filename_stem(name: str, fallback: str = "sequence") -> str`
  - `tessera.reassort.scan.unique_scan_dirs(segments: list[str]) -> dict[str, str]` (segment name -> directory name, unique)
  - `tessera.reassort.scan.scan_segment(segment, seq, dataset, out_dir, *, aligner, cache_dir, logger, dir_name: str | None = None) -> SegmentScan`

`assign._type_segment` already sanitises a segment name with `.strip(".")`; `scan.scan_segment` uses the same regex without it, so a record named `..` makes `seg_dir` the parent of the output and the scan then removes and rebuilds `<parent>/collection`. Move the sanitiser to one helper in `core/io.py` and use it in both places. Separately, two names can sanitise to the same string; `assign_segments` knows all the names, so it assigns each a unique directory and passes it down.

The two existing stubs of `scan_segment` in `test_reassort_assign.py` take a fixed keyword list and must accept the new `dir_name`.

- [ ] **Step 1: Write the failing tests**

Edit `tests/unit/test_reassort_assign.py` (existing lines change, and one test is added at the end):

````diff
diff --git a/tests/unit/test_reassort_assign.py b/tests/unit/test_reassort_assign.py
index bc488be..94489be 100644
--- a/tests/unit/test_reassort_assign.py
+++ b/tests/unit/test_reassort_assign.py
@@ -213,7 +213,8 @@ def test_scan_segments_scans_assigned_only(tmp_path, monkeypatch):
     monkeypatch.setattr(assign, "require_aligner", lambda aligner: None)
     seen = []
 
-    def fake_scan(segment, seq, dataset, out_dir, *, aligner, cache_dir, logger):
+    def fake_scan(segment, seq, dataset, out_dir, *, aligner, cache_dir, logger,
+                  dir_name=None):
         from tessera.reassort.scan import SegmentScan
         seen.append(segment)
         return SegmentScan(segment, True, segment == "HA", 1 if segment == "HA" else 0,
@@ -244,7 +245,8 @@ def test_scan_segments_marks_unassigned(tmp_path, monkeypatch):
                         "NA_pool": {na_tip: (99.0, 99.0)}})
     monkeypatch.setattr(assign, "require_aligner", lambda aligner: None)
 
-    def fake_scan(segment, seq, dataset, out_dir, *, aligner, cache_dir, logger):
+    def fake_scan(segment, seq, dataset, out_dir, *, aligner, cache_dir, logger,
+                  dir_name=None):
         from tessera.reassort.scan import SegmentScan
         return SegmentScan(segment, True, False, 0, "none")
     monkeypatch.setattr(assign, "scan_segment", fake_scan)
@@ -319,3 +321,28 @@ def test_alignment_fraction_filter_uses_skani_percent_scale(tmp_path, monkeypatc
     assert ha.status == "assigned"
     assert ha.strain == "full"
     assert ha.ani == pytest.approx(98.76)
+
+
+def test_scan_segments_gives_alike_named_segments_separate_directories(tmp_path, monkeypatch):
+    # "seg/1" and "seg_1" sanitise to the same directory name; each scan must get its own.
+    tip = tmp_path / "pool" / "strainA.fasta"
+    tip.parent.mkdir(parents=True)
+    tip.write_text(">x\nACGT\n")
+    _patch(monkeypatch,
+           resolve=lambda fasta, override, *, email, logger: _DS("ds"),
+           tips_by_path={"ds": [tip]},
+           ani_by_path={"pool": {tip: (99.0, 99.0)}})
+    monkeypatch.setattr(assign, "require_aligner", lambda aligner: None)
+    dirs = {}
+
+    def fake_scan(segment, seq, dataset, out_dir, *, aligner, cache_dir, logger,
+                  dir_name=None):
+        from tessera.reassort.scan import SegmentScan
+        dirs[segment] = dir_name
+        return SegmentScan(segment, True, False, 0, "none")
+    monkeypatch.setattr(assign, "scan_segment", fake_scan)
+
+    q = _write_query(tmp_path, [("seg/1", "AAAA"), ("seg_1", "CCCC")])
+    assign_segments(q, output=tmp_path / "out", scan_segments=True, logger=LOG)
+
+    assert dirs == {"seg/1": "seg_1", "seg_1": "seg_1_2"}
````

Append to the end of `tests/unit/test_reassort_scan.py`:

```python
def _two_clade_pool(tmp_path, monkeypatch):
    pool = tmp_path / "pool"
    pool.mkdir()
    tips = []
    for c in ("A", "B"):
        t = pool / f"{c}_consensus.fasta"
        t.write_text(f">{c}_consensus {c}\nACGTACGT\n")
        tips.append(t)
    _stub_pool(monkeypatch, tips)
    monkeypatch.setattr(scan, "build_msa", lambda params, logger: params.output)
    monkeypatch.setattr(scan, "run_recomb", lambda params, logger: "bp")


def test_scan_segment_name_cannot_leave_the_output_directory(tmp_path, monkeypatch):
    # The segment name is a FASTA header. ".." used to resolve the scan directory to the
    # parent of the output, whose `collection/` was then removed and rebuilt.
    _two_clade_pool(tmp_path, monkeypatch)
    out = tmp_path / "parent" / "out"
    out.mkdir(parents=True)
    precious = tmp_path / "parent" / "collection" / "precious.fasta"
    precious.parent.mkdir()
    precious.write_text(">x\nACGT\n")

    for name in ("..", ".", "../../escaped", "A/California/07/2009|HA"):
        scan_segment(name, "ACGTACGT", _DS(), out, aligner="mafft", cache_dir=None, logger=LOG)

    assert precious.exists()
    outside = [p for p in (tmp_path / "parent").rglob("*")
               if p.is_file() and out not in p.parents and p != precious]
    assert outside == []
    # Nothing is written into the output root itself: every segment has its own directory.
    assert [p for p in out.iterdir() if p.is_file()] == []


def test_safe_filename_stem():
    from tessera.core.io import safe_filename_stem

    assert safe_filename_stem("HA") == "HA"
    assert safe_filename_stem("A/California/07/2009|HA") == "A_California_07_2009_HA"
    assert safe_filename_stem("..", fallback="segment") == "segment"
    assert safe_filename_stem(".", fallback="segment") == "segment"
    assert safe_filename_stem("../x") == "_x"
    assert safe_filename_stem("", fallback="segment") == "segment"


def test_unique_scan_dirs_separates_names_that_sanitise_alike():
    from tessera.reassort.scan import unique_scan_dirs

    assert unique_scan_dirs(["HA", "NA"]) == {"HA": "HA", "NA": "NA"}
    # "seg/1" and "seg_1" both sanitise to "seg_1"; they must not share a directory.
    dirs = unique_scan_dirs(["seg/1", "seg_1", "seg 1"])
    assert dirs == {"seg/1": "seg_1", "seg_1": "seg_1_2", "seg 1": "seg_1_3"}
    assert len(set(dirs.values())) == 3
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_reassort_scan.py tests/unit/test_reassort_assign.py tests/unit/test_input_safety.py -q`

Expected: `4 failed, 51 passed`. The failures are:

```
tests/unit/test_reassort_scan.py::test_scan_segment_name_cannot_leave_the_output_directory
tests/unit/test_reassort_scan.py::test_safe_filename_stem
tests/unit/test_reassort_scan.py::test_unique_scan_dirs_separates_names_that_sanitise_alike
tests/unit/test_reassort_assign.py::test_scan_segments_gives_alike_named_segments_separate_directories
```

They fail because `safe_filename_stem` and `unique_scan_dirs` do not exist yet, `..` escapes the output directory, and `assign_segments` passes no `dir_name`.

- [ ] **Step 3: Implement**

`src/tessera/core/io.py`:

````diff
diff --git a/src/tessera/core/io.py b/src/tessera/core/io.py
index 62d524d..c6dfa58 100644
--- a/src/tessera/core/io.py
+++ b/src/tessera/core/io.py
@@ -11,6 +11,7 @@ from __future__ import annotations
 
 import gzip
 import logging
+import re
 import shutil
 from collections.abc import Sequence
 from pathlib import Path
@@ -99,6 +100,16 @@ def strip_sequence_extension(name: str) -> str:
     return name
 
 
+def safe_filename_stem(name: str, fallback: str = "sequence") -> str:
+    """A file or directory name derived from an untrusted label (a FASTA header).
+
+    Path separators and other punctuation become ``_``, and leading/trailing dots are
+    removed so the result can never be ``.`` or ``..`` or climb out of the directory it
+    is joined to. ``fallback`` is returned when nothing usable is left.
+    """
+    return re.sub(r"[^\w.-]+", "_", name).strip(".") or fallback
+
+
 def collection_genomes(directory: Path) -> list[Path]:
     """The genome files of a collection directory, sorted by name.
 
````

`src/tessera/reassort/assign.py`:

````diff
diff --git a/src/tessera/reassort/assign.py b/src/tessera/reassort/assign.py
index eb77d2b..04ba6b0 100644
--- a/src/tessera/reassort/assign.py
+++ b/src/tessera/reassort/assign.py
@@ -10,19 +10,23 @@ model that the intragenic recombination scan uses.
 from __future__ import annotations
 
 import logging
-import re
 import tempfile
 from dataclasses import dataclass, field
 from pathlib import Path
 
 from ..core.cache import nextclade_cache
 from ..core.errors import ToolExecutionError, UserInputError
-from ..core.io import read_fasta, strip_sequence_extension, write_fasta_record
+from ..core.io import (
+    read_fasta,
+    safe_filename_stem,
+    strip_sequence_extension,
+    write_fasta_record,
+)
 from ..discover.nextclade import NON_CLADE_MARKERS, build_pool, resolve_dataset
 from ..discover.panel import skani_available, skani_query_ani
 from ..recomb.typing import first_header
 from .constellation import DEFAULT_MARGIN, ParentGroup, call_constellation
-from .scan import SegmentScan, require_aligner, scan_segment
+from .scan import SegmentScan, require_aligner, scan_segment, unique_scan_dirs
 
 DEFAULT_ANI_FLOOR = 80.0  # a segment below this ANI to every tip is left unassigned
 # A tip aligning over less than this share of the segment is ignored. In PERCENT (0-100),
@@ -87,7 +91,7 @@ def _type_segment(seg, seq, overrides, ani_floor, margin, email, cache_dir, tmp,
     error, or any unexpected error) propagates so it surfaces rather than reading as unassigned."""
     # The segment name comes from a FASTA header, which may hold path separators
     # ("A/California/07/2009|HA") or climb out of the temp directory ("../x").
-    safe = re.sub(r"[^\w.-]+", "_", strip_sequence_extension(seg)).strip(".") or "segment"
+    safe = safe_filename_stem(strip_sequence_extension(seg), fallback="segment")
     seg_fasta = Path(tmp) / f"{safe}.fasta"
     with open(seg_fasta, "w") as fo:
         write_fasta_record(fo, seg, seq)
@@ -186,13 +190,15 @@ def assign_segments(
     result.pair_notes = call.pair_notes
 
     if scan_segments:
+        scan_dirs = unique_scan_dirs([s.segment for s in result.segments])
         for s in result.segments:
             if s.status == "assigned":
                 seq, dataset = to_scan[s.segment]
                 assert output is not None  # guarded above when scan_segments is set
                 result.scans.append(scan_segment(
                     s.segment, seq, dataset, output,
-                    aligner=aligner, cache_dir=cache_dir, logger=logger))
+                    aligner=aligner, cache_dir=cache_dir, logger=logger,
+                    dir_name=scan_dirs[s.segment]))
             else:
                 result.scans.append(SegmentScan(s.segment, False, False, 0, "unassigned"))
     return result
````

`src/tessera/reassort/scan.py`:

````diff
diff --git a/src/tessera/reassort/scan.py b/src/tessera/reassort/scan.py
index 44f492c..d5138fd 100644
--- a/src/tessera/reassort/scan.py
+++ b/src/tessera/reassort/scan.py
@@ -10,14 +10,13 @@ asks whether a single segment is itself a within-segment mosaic of two lineages.
 from __future__ import annotations
 
 import logging
-import re
 import shutil
 from dataclasses import dataclass
 from pathlib import Path
 
 from ..core.cache import nextclade_cache
 from ..core.errors import UserInputError
-from ..core.io import strip_sequence_extension, write_fasta_record
+from ..core.io import safe_filename_stem, strip_sequence_extension, write_fasta_record
 from ..discover.nextclade import NON_CLADE_MARKERS, build_pool
 from ..msa.build import MsaParams, build_msa
 from ..recomb.regions import DEFAULT_METHODS
@@ -87,13 +86,40 @@ def _summarize_regions(path: Path) -> tuple[int, bool]:
     return n, n > 0
 
 
+def unique_scan_dirs(segments: list[str]) -> dict[str, str]:
+    """Map each segment name to its own scan-directory name.
+
+    Segment names are FASTA headers, so two different names can sanitise to the same
+    string (``seg/1`` and ``seg_1``); their scans would then overwrite each other. Later
+    duplicates get a numeric suffix, in input order.
+    """
+    out: dict[str, str] = {}
+    used: set[str] = set()
+    for segment in segments:
+        base = safe_filename_stem(segment, fallback="segment")
+        name, n = base, 1
+        while name in used:
+            n += 1
+            name = f"{base}_{n}"
+        used.add(name)
+        out[segment] = name
+    return out
+
+
 def scan_segment(
     segment: str, seq: str, dataset, out_dir: Path, *,
     aligner: str, cache_dir: Path | None, logger: logging.Logger,
+    dir_name: str | None = None,
 ) -> SegmentScan:
     """Scan one segment for intragenic recombination. Never raises: a failure is recorded as
-    ``scanned=False`` so the caller can continue with the other segments."""
-    seg_name = re.sub(r"[^\w.-]+", "_", segment)
+    ``scanned=False`` so the caller can continue with the other segments.
+
+    ``dir_name`` is the scan directory under ``out_dir`` (see :func:`unique_scan_dirs`);
+    by default it is derived from the segment name. Either way it is a sanitised single
+    path component: the segment name comes from a FASTA header and must not be able to
+    name ``..`` or the output root itself.
+    """
+    seg_name = safe_filename_stem(dir_name or segment, fallback="segment")
     seg_dir = out_dir / seg_name
     try:
         pool = build_pool(
````

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_reassort_scan.py tests/unit/test_reassort_assign.py tests/unit/test_input_safety.py tests/unit/test_cli_reassort.py -q`

Expected: all pass (the Step 2 command now gives `55 passed`).

- [ ] **Step 5: Run the fast gate**

Run: `ruff check src tests validation && pytest -m "not requires_binary" -q && mypy src`

Expected: `All checks passed!`, `593 passed, 3 deselected`, `Success: no issues found in 79 source files`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/core/io.py src/tessera/reassort/assign.py src/tessera/reassort/scan.py tests/unit/test_reassort_assign.py tests/unit/test_reassort_scan.py
git commit -m "$(cat <<'EOF'
Keep reassort segment scans inside their own directories

A segment record named '..' resolved the scan directory to the parent of the
output, whose collection/ was then replaced. Sanitise segment names with one shared
helper and give names that sanitise alike separate scan directories.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

---
### Task 5: A6 -- find-references curation removes only what it downloaded

**Spec item(s):** A6

**Files:**
- Modify: `src/tessera/discover/panel.py` (`curate_collection_dir(..., protect=)`, two new role labels)
- Modify: `src/tessera/discover/run.py` (`find_references`, `_curate_download`)
- Modify: `src/tessera/cli/cmd_find_references.py` (`--curate` help text)
- Test: `tests/unit/test_panel.py` (insert), `tests/unit/test_discover.py` (insert)

**Interfaces:**
- Consumes: nothing from other tasks.
- Produces:
  - `tessera.discover.panel.curate_collection_dir(query_fasta, collection, backbone, *, ani_margin=..., af_min=..., derep_ani=..., protect: Iterable[Path] = (), logger) -> CurationResult`. Protected files are compared but never deleted; in `CurationResult.table` their role reads `sibling-kept` / `redundant-kept`, they are appended to `.kept` and removed from `.siblings` / `.redundant`.
  - `tessera.discover.run._curate_download(params, query_label, query_row, preexisting: list[Path], logger) -> None`
  - `panel_lineages.tsv` may now contain the roles `sibling-kept` and `redundant-kept`.

`find-references -c C --download C --curate` is the invocation the docs give for growing a collection in place. `_curate_download` curates the whole download directory, so it deletes the user's own genomes when they look like siblings or near-duplicates. Record what the directory held before the download and pass it as `protect`.

Reproducing this exposed a second defect in the same function: the curation backbone is picked from `--collection` *after* the download, so when `--download` is the collection a downloaded sibling (closer to the query than anything else) is picked as the backbone. The first new test covers both: `NEWSIB` must be dropped and `userC` must survive.

The third test (`..._into_a_new_directory`) already passes before the change. It pins that the fix does not break the case where `--download` names a directory that does not exist yet.

- [ ] **Step 1: Write the failing tests**

Insert into `tests/unit/test_discover.py`, immediately before the line `def test_download_without_efetch_is_a_clear_error(monkeypatch, tmp_path, logger):`:

```python
def test_curate_after_download_keeps_the_users_own_genomes(monkeypatch, tmp_path, logger):
    """`find-references -c C --download C --curate` is the documented way to grow a
    collection in place. Curation used to prune C itself, deleting the user's files."""
    from tessera.discover import panel

    coll = tmp_path / "collection"
    coll.mkdir()
    for name in ("refA", "refB", "userC"):
        (coll / f"{name}.fasta").write_text(f">{name}\nACGT\n")

    def fake_blast(seq, *, max_hits, logger, email=None, cache_dir=None):
        return [
            Hit("NEW123", "novel donor virus", 90.0, 95.0, 1e-40),
            Hit("NEWSIB", "a relative of the query", 97.0, 96.0, 1e-60),
        ]

    def fake_efetch(accession, collection_dir, logger):
        path = collection_dir / f"{accession}.fasta"
        path.write_text(f">{accession}\nACGT\n")
        return path

    def fake_ani(query, refs, logger):
        # refA is the backbone; refB and userC are its whole-genome twins (would be
        # dropped as siblings); NEWSIB is a downloaded sibling; NEW123 a regional donor.
        table = {
            "refA": (97.0, 95.0), "refB": (96.0, 95.0), "userC": (96.8, 96.0),
            "NEW123": (88.0, 30.0), "NEWSIB": (99.0, 98.0),
        }
        return {r: table[r.name.split(".")[0]] for r in refs}

    monkeypatch.setattr(discover_run, "blast_subsequence", fake_blast)
    monkeypatch.setattr(discover_run, "efetch_available", lambda: True)
    monkeypatch.setattr(discover_run, "efetch_fasta", fake_efetch)
    monkeypatch.setattr(panel, "skani_available", lambda: True)
    monkeypatch.setattr(panel, "skder_available", lambda: False)
    monkeypatch.setattr(panel, "skani_query_ani", fake_ani)

    find_references(
        FindRefParams(
            msa=_msa(tmp_path), query="q", output=tmp_path / "out",
            window_size=60, window_step=30, top_gaps=1,
            collection=coll, download=coll, curate=True,
        ),
        logger,
    )

    assert sorted(p.name for p in coll.iterdir()) == [
        "NEW123.fasta", "refA.fasta", "refB.fasta", "userC.fasta",
    ]
    rows = dict(
        line.split("\t")[:2]
        for line in (tmp_path / "out" / "panel_lineages.tsv").read_text().splitlines()[1:]
    )
    assert rows["userC"] == "sibling-kept"
    assert rows["NEWSIB"] == "sibling-dropped"


def test_curate_after_download_into_a_new_directory(monkeypatch, tmp_path, logger):
    """--download may name a directory that does not exist yet. Nothing in it predates
    the run, so curation is free to drop a downloaded sibling."""
    from tessera.discover import panel

    coll = tmp_path / "collection"
    coll.mkdir()
    (coll / "refA.fasta").write_text(">refA\nACGT\n")
    fresh = tmp_path / "downloads"

    def fake_blast(seq, *, max_hits, logger, email=None, cache_dir=None):
        return [
            Hit("NEW123", "novel donor virus", 90.0, 95.0, 1e-40),
            Hit("NEWSIB", "a relative of the query", 97.0, 96.0, 1e-60),
        ]

    def fake_efetch(accession, collection_dir, logger):
        collection_dir.mkdir(parents=True, exist_ok=True)
        path = collection_dir / f"{accession}.fasta"
        path.write_text(f">{accession}\nACGT\n")
        return path

    def fake_ani(query, refs, logger):
        table = {"refA": (92.0, 95.0), "NEW123": (88.0, 30.0), "NEWSIB": (99.0, 98.0)}
        return {r: table[r.name.split(".")[0]] for r in refs}

    monkeypatch.setattr(discover_run, "blast_subsequence", fake_blast)
    monkeypatch.setattr(discover_run, "efetch_available", lambda: True)
    monkeypatch.setattr(discover_run, "efetch_fasta", fake_efetch)
    monkeypatch.setattr(panel, "skani_available", lambda: True)
    monkeypatch.setattr(panel, "skder_available", lambda: False)
    monkeypatch.setattr(panel, "skani_query_ani", fake_ani)

    find_references(
        FindRefParams(
            msa=_msa(tmp_path), query="q", output=tmp_path / "out",
            window_size=60, window_step=30, top_gaps=1,
            collection=coll, download=fresh, curate=True,
        ),
        logger,
    )

    assert sorted(p.name for p in fresh.iterdir()) == ["NEW123.fasta"]
    assert sorted(p.name for p in coll.iterdir()) == ["refA.fasta"]
```

Insert into `tests/unit/test_panel.py`, immediately before the line `# --- typed Lineage column (conditional) ----------------------------------------`:

```python
def test_curate_collection_dir_never_deletes_protected_files(monkeypatch, tmp_path, logger):
    """Files the caller marks as protected take part in the comparison but stay on disk."""
    coll = tmp_path / "coll"
    coll.mkdir()
    backbone, own_sibling, new_sibling, parent = _genomes(
        coll, ["A1", "user_rel", "downloaded_rel", "C"]
    )
    ani = {
        backbone: (92.0, 81.0), own_sibling: (97.0, 94.0),
        new_sibling: (97.5, 95.0), parent: (89.0, 94.0),
    }
    monkeypatch.setattr(panel, "skani_available", lambda: True)
    monkeypatch.setattr(panel, "skder_available", lambda: False)
    monkeypatch.setattr(panel, "skani_query_ani", lambda *a, **k: ani)

    result = panel.curate_collection_dir(
        tmp_path / "q.fasta", coll, backbone, protect=[own_sibling], logger=logger,
    )

    assert sorted(p.name for p in coll.iterdir()) == ["A1.fasta", "C.fasta", "user_rel.fasta"]
    roles = {r["genome"]: r["role"] for r in result.table}
    assert roles == {
        "A1": "backbone", "C": "representative",
        "user_rel": "sibling-kept", "downloaded_rel": "sibling-dropped",
    }
    assert own_sibling in result.kept
    assert result.siblings == [new_sibling]


def test_panel_html_labels_a_kept_sibling():
    html = panel.panel_table_html([_panel_row("user_rel", role="sibling-kept")])
    assert "kept (sibling, pre-existing)" in html
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_panel.py tests/unit/test_discover.py -q`

Expected: `3 failed, 20 passed`. The failures are:

```
tests/unit/test_panel.py::test_curate_collection_dir_never_deletes_protected_files
tests/unit/test_panel.py::test_panel_html_labels_a_kept_sibling
tests/unit/test_discover.py::test_curate_after_download_keeps_the_users_own_genomes
```

They fail because `curate_collection_dir` has no `protect` parameter, the role label is unknown, and the user's `refB`/`userC` are deleted.

- [ ] **Step 3: Implement**

`src/tessera/cli/cmd_find_references.py`:

````diff
diff --git a/src/tessera/cli/cmd_find_references.py b/src/tessera/cli/cmd_find_references.py
index 7896bc0..c38217a 100644
--- a/src/tessera/cli/cmd_find_references.py
+++ b/src/tessera/cli/cmd_find_references.py
@@ -61,8 +61,9 @@ def find_references(
     ),
     curate: bool = typer.Option(
         False, "--curate",
-        help="After download, drop the query's siblings and dereplicate (needs skani/skDER, "
-        "--collection as the backbone source).",
+        help="After download, drop the query's siblings and near-duplicates among the "
+        "new downloads (needs skani/skDER, --collection as the backbone source). Genomes "
+        "already in the download directory are never removed.",
     ),
     sibling_margin: float = typer.Option(
         3.0, "--sibling-margin",
````

`src/tessera/discover/panel.py`:

````diff
diff --git a/src/tessera/discover/panel.py b/src/tessera/discover/panel.py
index 3f87d10..822a1a5 100644
--- a/src/tessera/discover/panel.py
+++ b/src/tessera/discover/panel.py
@@ -26,6 +26,7 @@ import html
 import logging
 import shutil
 import tempfile
+from collections.abc import Iterable
 from dataclasses import dataclass, field
 from pathlib import Path
 
@@ -266,22 +267,48 @@ def curate_collection_dir(
     ani_margin: float = DEFAULT_SIBLING_MARGIN,
     af_min: float = DEFAULT_AF_MIN,
     derep_ani: float = DEFAULT_DEREP_ANI,
+    protect: Iterable[Path] = (),
     logger: logging.Logger,
 ) -> CurationResult:
     """Curate the genomes in ``collection`` in place: drop siblings/redundant files.
 
     Runs :func:`curate_panel`, then deletes the dropped genome files from disk so a
     subsequent MSA rebuild sees only the diverse, sibling-free panel.
+
+    ``protect`` lists files that must stay on disk whatever the comparison says -- the
+    genomes a directory held before this run added to it. They still take part in the
+    sibling and redundancy comparison (a new download that duplicates one is dropped),
+    and a protected genome that would have been removed is reported with the role
+    ``sibling-kept`` / ``redundant-kept`` instead, so the table says what is on disk.
     """
     genomes = collection_genomes(collection)
     result = curate_panel(
         query_fasta, genomes, backbone,
         ani_margin=ani_margin, af_min=af_min, derep_ani=derep_ani, logger=logger,
     )
+    protected = {Path(p).resolve() for p in protect}
     keep = {p.resolve() for p in result.kept} | {backbone.resolve()}
+    spared: list[Path] = []
     for g in genomes:
-        if g.resolve() not in keep:
+        if g.resolve() in keep:
+            continue
+        if g.resolve() in protected:
+            spared.append(g)
+        else:
             g.unlink()
+    if spared:
+        spared_labels = {strip_sequence_extension(g.name) for g in spared}
+        for row in result.table:
+            if row["genome"] in spared_labels:
+                row["role"] = row["role"].replace("-dropped", "-kept")
+        spared_set = {g.resolve() for g in spared}
+        result.kept = [*result.kept, *spared]
+        result.siblings = [g for g in result.siblings if g.resolve() not in spared_set]
+        result.redundant = [g for g in result.redundant if g.resolve() not in spared_set]
+        logger.info(
+            "Kept %d pre-existing genome(s) that curation would have dropped: %s.",
+            len(spared), ", ".join(sorted(spared_labels)),
+        )
     return result
 
 
@@ -372,10 +399,15 @@ def write_panel_tsv(
 _ROLE_LABEL = {
     "backbone": "backbone",
     "representative": "kept (parent)",
+    "sibling-kept": "kept (sibling, pre-existing)",
+    "redundant-kept": "kept (redundant, pre-existing)",
     "sibling-dropped": "dropped (sibling)",
     "redundant-dropped": "dropped (redundant)",
 }
-_ROLE_ORDER = {"backbone": 0, "representative": 1, "sibling-dropped": 2, "redundant-dropped": 3}
+_ROLE_ORDER = {
+    "backbone": 0, "representative": 1, "sibling-kept": 2, "redundant-kept": 3,
+    "sibling-dropped": 4, "redundant-dropped": 5,
+}
 
 
 def panel_table_html(table: list[dict], lineage_map: LineageMap | None = None) -> str:
````

`src/tessera/discover/run.py`:

````diff
diff --git a/src/tessera/discover/run.py b/src/tessera/discover/run.py
index 8019581..6b68fee 100644
--- a/src/tessera/discover/run.py
+++ b/src/tessera/discover/run.py
@@ -104,10 +104,16 @@ def find_references(params: FindRefParams, logger: logging.Logger) -> list[Candi
     _print_candidates(candidates, logger)
 
     if params.download is not None:
+        # What the download directory held before this run. It is often the user's own
+        # collection (`--download <collection>` is the documented way to grow one in
+        # place), so curation below may remove only what this run adds to it.
+        preexisting = (
+            collection_genomes(params.download) if params.download.is_dir() else []
+        )
         downloaded = _download(candidates, params.download, logger)
         _write_downloaded(params.output, downloaded, logger)
         if params.curate and downloaded:
-            _curate_download(params, query_label, query_row, logger)
+            _curate_download(params, query_label, query_row, preexisting, logger)
     elif candidates:
         logger.info(
             "Re-run with --download <collection_dir> to add the new references, "
@@ -287,13 +293,18 @@ def _download(
 
 
 def _curate_download(
-    params: FindRefParams, query_label: str, query_row: str, logger: logging.Logger
+    params: FindRefParams, query_label: str, query_row: str,
+    preexisting: list[Path], logger: logging.Logger,
 ) -> None:
-    """Drop the query's siblings and dereplicate the download directory in place.
+    """Drop the query's siblings and near-duplicates among this run's downloads.
 
     The backbone (the query's whole-genome anchor) is chosen from the existing
     ``--collection``, so a freshly-downloaded sibling cannot be mistaken for it. The
     query is reconstructed from its (de-gapped) MSA row, as skani needs a FASTA.
+
+    ``preexisting`` are the files the download directory held before this run; they are
+    compared against but never deleted, so pointing ``--download`` at the collection
+    itself cannot remove the user's own genomes.
     """
     from .panel import (
         curate_collection_dir,
@@ -311,18 +322,26 @@ def _curate_download(
         return
     qfasta = params.output / "query.degapped.fasta"
     qfasta.write_text(f">{query_label}\n{query_row.replace('-', '')}\n")
+    assert params.download is not None  # only called from the `download is not None` branch
+    # When --download is the collection itself, the collection now also holds this run's
+    # downloads; leave them out so a downloaded sibling cannot become the backbone.
+    held_before = {p.resolve() for p in preexisting}
+    new_files = {
+        p.resolve() for p in collection_genomes(params.download)
+        if p.resolve() not in held_before
+    }
     backbone = pick_backbone(
-        qfasta, collection_genomes(params.collection),
+        qfasta,
+        [g for g in collection_genomes(params.collection) if g.resolve() not in new_files],
         af_min=params.af_min, logger=logger,
     )
     if backbone is None:
         logger.warning("Could not determine a backbone from --collection; skipping curation.")
         return
-    assert params.download is not None  # only called from the `download is not None` branch
     curation = curate_collection_dir(
         qfasta, params.download, backbone,
         ani_margin=params.sibling_margin, af_min=params.af_min,
-        derep_ani=params.derep_ani, logger=logger,
+        derep_ani=params.derep_ani, protect=preexisting, logger=logger,
     )
     write_panel_tsv(params.output / "panel_lineages.tsv", curation.table, logger)
 
````

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_panel.py tests/unit/test_discover.py tests/unit/test_iterate.py -q`

Expected: all pass (the Step 2 command now gives `23 passed`).

- [ ] **Step 5: Run the fast gate**

Run: `ruff check src tests validation && pytest -m "not requires_binary" -q && mypy src`

Expected: `All checks passed!`, `597 passed, 3 deselected`, `Success: no issues found in 79 source files`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/cli/cmd_find_references.py src/tessera/discover/panel.py src/tessera/discover/run.py tests/unit/test_discover.py tests/unit/test_panel.py
git commit -m "$(cat <<'EOF'
Curate only downloaded genomes in find-references

With --download pointing at the collection itself, --curate deleted genomes that
were already there, and could pick a downloaded sibling as the backbone. Files
present before the run are now protected and excluded from backbone choice.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

---
### Task 6: A3-A5 -- fill loop: curate before building, align what was downloaded

**Spec item(s):** A3, A4, A5

**Files:**
- Modify: `src/tessera/discover/iterate.py` (module docstring; new `_curation_backbone`, `_curate_round`, `_labels`; `_grow_collection` rewritten)
- Modify: `src/tessera/cli/cmd_fill_references.py` (`--curate` help text)
- Test: `tests/unit/test_iterate.py` (append)

**Interfaces:**
- Consumes: `tessera.discover.panel.curate_collection_dir(query_fasta, collection, backbone, *, ani_margin, af_min, derep_ani, logger)` and `pick_backbone(query_fasta, genomes, *, af_min, logger)` as they are on `main` (this task does not need Task 5's `protect`).
- Produces (all in `tessera.discover.iterate`):
  - `_curation_backbone(params: FillParams, collection: Path, logger) -> Path | None`
  - `_curate_round(params, collection, backbone: Path | None, panel_rows: dict[str, dict], logger) -> None`
  - `_labels(collection: Path) -> set[str]`
  - `_grow_collection(...)` keeps its signature and return type `tuple[list[RoundResult], dict[str, dict], Path | None]`; the returned MSA may now be `<output>/final.msa.fasta`.

Three defects share one cause, the order of steps inside a round (`build -> scan -> search -> download -> curate`):

- **A3.** On a `max_rounds` exit the last round's downloads are in `collection/` but in no MSA.
- **A4.** Curation sits after two `break`s (no gaps; nothing downloaded), so it never runs on a supplied collection that already covers the query -- the sibling case.
- **A5.** The curation backbone is always auto-picked, so the genome named by `--reference` can be dropped as its twin.

New order: `[curate if owed] -> build -> scan -> stop checks -> search -> download`. Curation is "owed" before a build whenever the collection holds genomes that have not been through it: a collection the user supplied (before round 1) and each round's downloads (before the next build). After the loop, any owed curation runs, and if the collection is not what the last MSA was built from, one more MSA (`final.msa.fasta`) is built with no further search.

**Deviation from the spec, deliberately narrow.** The spec says curation runs "on round 1 before the first build". This plan does that only when the user supplied `--collection`. A collection seeded from scratch (`detect`, `build-panel`, `fill-references` without `-c`) is *not* curated before round 1, because seeding already filters siblings, and curating the seed would change what `detect` recruits: `filter_siblings` classes any genome covering >= 90 % of the query at an ANI within 1.5 points of the backbone as a twin, which on a panel under about 1.5 % divergent is every genome. `test_curate_does_not_run_on_a_freshly_seeded_collection` pins this; it passes before and after the change. For the same reason `_curate_round` warns when curation leaves fewer than two references.

Two of the new tests pass before the change and are guards, not regressions: `test_no_extra_build_when_the_last_round_adds_nothing` and `test_curate_does_not_run_on_a_freshly_seeded_collection`.

Test helpers: `_fake_skani` stubs skani at the `panel` layer so the real `pick_backbone` / `filter_siblings` / `curate_collection_dir` run; `caplog` tests use a logger outside the `tessera` namespace because the CLI tests set `propagate = False` on `tessera`, after which `caplog` no longer sees its records.

- [ ] **Step 1: Write the failing tests**

Append to the end of `tests/unit/test_iterate.py`:

```python
# --- the round structure: what is aligned, and when the panel is curated ----------

def _recording_build(monkeypatch, builds):
    """Replace build_msa with a stub that records which genomes each MSA was built from
    and writes them into the file, so the published panel can be inspected."""
    def fake_build(p, logger):
        members = sorted(f.name for f in p.collection.iterdir())
        builds.append((p.output.name, members))
        p.output.write_text("".join(f">{m}\nACGT\n" for m in ["q", *members]))
    monkeypatch.setattr(iterate, "build_msa", fake_build)


def _one_new_hit_per_round(monkeypatch):
    counter = iter(range(1, 50))

    def collect(*a, **k):
        return [Candidate(_gap(0.8), Hit(f"NEW{next(counter)}", "t", 90.0, 95.0, 1e-9), False)]

    def download(cands, dest, logger):
        for c in cands:
            (dest / f"{c.hit.accession}.fasta").write_text(f">{c.hit.accession}\nACGT\n")
        return cands

    monkeypatch.setattr(iterate, "collect_candidates", collect)
    monkeypatch.setattr(iterate, "_download", download)


def _fake_skani(monkeypatch, table):
    """Stub skani at the panel layer: ``table`` maps a genome label to (ANI %, AF %)."""
    from tessera.discover import panel

    def fake_ani(query, refs, logger):
        return {r: table[r.name.split(".")[0]] for r in refs}

    monkeypatch.setattr(panel, "skani_query_ani", fake_ani)
    monkeypatch.setattr(panel, "skani_available", lambda: True)
    monkeypatch.setattr(panel, "skder_available", lambda: False)
    monkeypatch.setattr(iterate, "skani_available", lambda: True)


def test_last_round_downloads_are_aligned(monkeypatch, tmp_path, logger):
    """On a max_rounds exit the final round's downloads used to sit in collection/ and be
    counted in the summary without ever reaching the published alignment."""
    query, coll, out = _setup(tmp_path)
    _common_mocks(monkeypatch, [([_gap(0.80)], 0.94), ([_gap(0.90)], 0.94)])
    builds: list[tuple[str, list[str]]] = []
    _recording_build(monkeypatch, builds)
    _one_new_hit_per_round(monkeypatch)

    trace = fill_references(
        FillParams(query=query, collection=coll, output=out, max_rounds=2), logger
    )

    assert [(r.round, r.added) for r in trace] == [(1, ["NEW1"]), (2, ["NEW2"])]
    final_members = ["NEW1.fasta", "NEW2.fasta", "refA.fasta"]
    assert sorted(p.name for p in (out / "collection").iterdir()) == final_members
    assert builds == [
        ("round1.msa.fasta", ["refA.fasta"]),
        ("round2.msa.fasta", ["NEW1.fasta", "refA.fasta"]),
        ("final.msa.fasta", final_members),
    ]
    published = (out / "panel.msa.fasta").read_text()
    assert ">NEW2.fasta" in published


def test_no_extra_build_when_the_last_round_adds_nothing(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    _common_mocks(monkeypatch, [([_gap(0.84)], 0.94), ([], 0.95)])
    builds: list[tuple[str, list[str]]] = []
    _recording_build(monkeypatch, builds)
    _one_new_hit_per_round(monkeypatch)

    fill_references(FillParams(query=query, collection=coll, output=out, max_rounds=2), logger)

    assert [name for name, _ in builds] == ["round1.msa.fasta", "round2.msa.fasta"]
    assert not (out / "final.msa.fasta").exists()


def test_curate_runs_on_a_supplied_collection_before_the_first_build(
    monkeypatch, tmp_path, logger
):
    """A whole-genome sibling in the starting collection leaves no coverage gap, so the
    loop converges in round 1. Curation used to sit after the download step and never ran
    -- in exactly the case it exists for."""
    query, coll, out = _setup(tmp_path)
    (coll / "SIBLING.fasta").write_text(">SIBLING\nACGT\n")
    (coll / "PARENT.fasta").write_text(">PARENT\nACGT\n")
    _common_mocks(monkeypatch, [([], 0.95)])
    builds: list[tuple[str, list[str]]] = []
    _recording_build(monkeypatch, builds)
    # refA is picked as backbone only if SIBLING does not outrank it; make the sibling a
    # whole-genome twin of refA (comparable ANI, full coverage) and PARENT a regional donor.
    _fake_skani(monkeypatch, {
        "refA": (97.0, 95.0), "SIBLING": (96.5, 95.0), "PARENT": (88.0, 30.0),
    })
    sections = {}
    monkeypatch.setattr(
        iterate, "run_recomb", lambda p, logger, **kw: sections.update(kw),
    )

    fill_references(FillParams(query=query, collection=coll, output=out, curate=True), logger)

    assert builds == [("round1.msa.fasta", ["PARENT.fasta", "refA.fasta"])]
    assert (coll / "SIBLING.fasta").exists()  # the user's own collection is untouched
    panel_tsv = (out / "panel_lineages.tsv").read_text()
    assert "SIBLING\tsibling-dropped" in panel_tsv
    assert "Reference panel" in [title for title, _ in sections["extra_sections"]]


def test_curate_does_not_run_on_a_freshly_seeded_collection(monkeypatch, tmp_path, logger):
    """Seeding already filters siblings (parents mode, pool selection). Curating the seed
    again before round 1 would change what `detect` recruits, so it is not done."""
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "ACGT" * 100 + "\n")
    out = tmp_path / "out"
    _common_mocks(monkeypatch, [([], 0.95)])
    monkeypatch.setattr(iterate, "skani_available", lambda: True)

    def seed(params, collection, query_records, exclude, logger):
        for name in ("S1", "S2"):
            (collection / f"{name}.fasta").write_text(f">{name}\nACGT\n")

    monkeypatch.setattr(iterate, "_seed_collection", seed)
    calls = []
    monkeypatch.setattr(iterate, "curate_collection_dir",
                        lambda *a, **k: calls.append("curate"))
    monkeypatch.setattr(iterate, "pick_backbone", lambda *a, **k: calls.append("backbone"))

    fill_references(FillParams(query=query, collection=None, output=out, curate=True), logger)

    assert calls == []


def test_curate_keeps_the_user_reference_as_backbone(monkeypatch, tmp_path, logger):
    """`--curate --reference X` used to auto-pick another backbone, drop X as its twin,
    and fail the next round with "Reference 'X' not found among the staged genomes"."""
    from tessera.core.io import select_reference

    query, coll, out = _setup(tmp_path)
    (coll / "MYREF.fasta").write_text(">MYREF\nACGT\n")
    _common_mocks(monkeypatch, [([_gap(0.80)], 0.94), ([_gap(0.90)], 0.94)])
    _one_new_hit_per_round(monkeypatch)
    # Auto-picking would choose refA (highest ANI); MYREF is then its whole-genome twin.
    _fake_skani(monkeypatch, {
        "refA": (97.0, 95.0), "MYREF": (96.5, 95.0),
        "NEW1": (88.0, 30.0), "NEW2": (87.0, 30.0),
    })
    builds: list[str] = []

    def fake_build(p, logger):  # resolve the reference the way the real build_msa does
        select_reference(sorted(p.collection.iterdir()), p.query, False, p.reference)
        builds.append(p.output.name)
        p.output.write_text(">q\nACGT\n")

    monkeypatch.setattr(iterate, "build_msa", fake_build)

    fill_references(
        FillParams(query=query, collection=coll, output=out, curate=True,
                   reference="MYREF", max_rounds=2),
        logger,
    )

    assert builds == ["round1.msa.fasta", "round2.msa.fasta", "final.msa.fasta"]
    assert (out / "collection" / "MYREF.fasta").exists()
    rows = dict(
        line.split("\t")[:2]
        for line in (out / "panel_lineages.tsv").read_text().splitlines()[1:]
    )
    assert rows["MYREF"] == "backbone"


def test_last_round_downloads_are_curated_before_the_final_build(monkeypatch, tmp_path, logger):
    """The final build sees a curated collection. Here the only download is a sibling, so
    after curation the collection is what round 1 was built from and no rebuild is needed."""
    query, coll, out = _setup(tmp_path)
    _common_mocks(monkeypatch, [([_gap(0.80)], 0.94)])
    builds: list[tuple[str, list[str]]] = []
    _recording_build(monkeypatch, builds)
    _one_new_hit_per_round(monkeypatch)
    # The one download (NEW1) is a sibling: closer than the backbone, whole-genome.
    _fake_skani(monkeypatch, {"refA": (92.0, 95.0), "NEW1": (99.0, 98.0)})

    trace = fill_references(
        FillParams(query=query, collection=coll, output=out, curate=True, max_rounds=1),
        logger,
    )

    assert builds == [("round1.msa.fasta", ["refA.fasta"])]
    assert not (out / "final.msa.fasta").exists()
    assert sorted(p.name for p in (out / "collection").iterdir()) == ["refA.fasta"]
    assert trace[0].added == []  # the sibling was downloaded, then curated away


def test_curate_accepts_a_reference_given_with_its_extension(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    (coll / "MYREF.fasta").write_text(">MYREF\nACGT\n")
    _common_mocks(monkeypatch, [([], 0.95)])
    _fake_skani(monkeypatch, {"refA": (97.0, 95.0), "MYREF": (96.5, 95.0)})

    fill_references(
        FillParams(query=query, collection=coll, output=out, curate=True,
                   reference="MYREF.fasta"),
        logger,
    )

    rows = dict(
        line.split("\t")[:2]
        for line in (out / "panel_lineages.tsv").read_text().splitlines()[1:]
    )
    assert rows == {"MYREF": "backbone", "refA": "sibling-dropped"}


def test_curate_warns_when_it_leaves_too_few_references(monkeypatch, tmp_path, caplog):
    """On a panel as close to the query as its own lineage, every genome but the backbone
    is classed as a sibling. Detection cannot run on one reference; say why."""
    import logging

    query, coll, out = _setup(tmp_path)
    for name in ("refB", "refC"):
        (coll / f"{name}.fasta").write_text(f">{name}\nACGT\n")
    _common_mocks(monkeypatch, [([], 0.95)])
    _fake_skani(monkeypatch, {
        "refA": (99.6, 99.0), "refB": (99.5, 99.0), "refC": (99.4, 99.0),
    })
    # Not under the "tessera" logger: the CLI tests switch its propagation off, and
    # caplog only sees records that reach the root logger.
    log = logging.getLogger("fill_loop_test")
    with caplog.at_level(logging.WARNING, logger="fill_loop_test"):
        fill_references(FillParams(query=query, collection=coll, output=out, curate=True), log)

    assert sorted(p.name for p in (out / "collection").iterdir()) == ["refA.fasta"]
    assert "Curation left 1 reference" in caplog.text
    assert "without --curate" in caplog.text
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_iterate.py -q`

Expected: `5 failed, 22 passed`. The failures are:

```
tests/unit/test_iterate.py::test_last_round_downloads_are_aligned
tests/unit/test_iterate.py::test_curate_runs_on_a_supplied_collection_before_the_first_build
tests/unit/test_iterate.py::test_curate_keeps_the_user_reference_as_backbone
tests/unit/test_iterate.py::test_curate_accepts_a_reference_given_with_its_extension
tests/unit/test_iterate.py::test_curate_warns_when_it_leaves_too_few_references
```

They fail because no `final.msa.fasta` is built, curation does not run on the supplied collection (so no `panel_lineages.tsv`), `MYREF` is deleted, and no warning is logged.

- [ ] **Step 3: Implement**

`src/tessera/cli/cmd_fill_references.py`:

````diff
diff --git a/src/tessera/cli/cmd_fill_references.py b/src/tessera/cli/cmd_fill_references.py
index bf7ece4..1849949 100644
--- a/src/tessera/cli/cmd_fill_references.py
+++ b/src/tessera/cli/cmd_fill_references.py
@@ -121,7 +121,9 @@ def fill_references(
     ),
     curate: bool = typer.Option(
         False, "--curate",
-        help="Drop the query's siblings and dereplicate each round (needs skani/skDER).",
+        help="Drop the query's siblings and dereplicate before each alignment build: the "
+        "supplied collection, then each round's downloads (needs skani/skDER). "
+        "--reference, when given, is the curation backbone and is never removed.",
     ),
     sibling_margin: float = typer.Option(
         3.0, "--sibling-margin",
````

`src/tessera/discover/iterate.py`:

````diff
diff --git a/src/tessera/discover/iterate.py b/src/tessera/discover/iterate.py
index 223179b..4895fa7 100644
--- a/src/tessera/discover/iterate.py
+++ b/src/tessera/discover/iterate.py
@@ -4,7 +4,9 @@ Each round rebuilds the MSA from the (growing) collection, scans it for coverage
 gaps, BLASTs the worst gaps against NCBI, and downloads the best new reference per
 gap into the collection. The loop stops when the gaps close, when no new reference
 can be found, when coverage stops improving (a stubborn residual is reported, not
-chased forever), or at ``max_rounds``.
+chased forever), or at ``max_rounds``. If the last round downloaded references, one
+further MSA is built from them, so the published panel always contains every
+reference the run reports.
 
 Because every round rebuilds the alignment, this needs an aligner binary and
 Entrez Direct, and it contacts NCBI over the network.
@@ -279,6 +281,58 @@ def fill_references(params: FillParams, logger: logging.Logger) -> list[RoundRes
     return trace
 
 
+def _curation_backbone(
+    params: FillParams, collection: Path, logger: logging.Logger
+) -> Path | None:
+    """The genome curation is anchored on: the user's ``--reference`` when given,
+    otherwise the query's closest whole-genome relative.
+
+    The backbone is never removed by curation, so honouring ``--reference`` here is what
+    keeps it in the collection for the next MSA build (which resolves the same name).
+    """
+    genomes = collection_genomes(collection)
+    if params.reference:
+        wanted = strip_sequence_extension(Path(params.reference).name)
+        for genome in genomes:
+            if strip_sequence_extension(genome.name) == wanted:
+                return genome
+        # Not in the collection: let build_msa report it with its own message.
+    return pick_backbone(params.query, genomes, af_min=params.af_min, logger=logger)
+
+
+def _curate_round(
+    params: FillParams,
+    collection: Path,
+    backbone: Path | None,
+    panel_rows: dict[str, dict],
+    logger: logging.Logger,
+) -> None:
+    """Drop siblings and near-duplicates from ``collection`` in place, recording roles."""
+    if backbone is None or not backbone.exists():
+        backbone = _curation_backbone(params, collection, logger)
+    if backbone is None:
+        return
+    curation = curate_collection_dir(
+        params.query, collection, backbone,
+        ani_margin=params.sibling_margin, af_min=params.af_min,
+        derep_ani=params.derep_ani, logger=logger,
+    )
+    for row in curation.table:
+        panel_rows[row["genome"]] = row
+    remaining = len(collection_genomes(collection))
+    if remaining < 2:
+        logger.warning(
+            "Curation left %d reference(s) in the panel; detection needs at least 2. "
+            "Every other genome was classed as a sibling of the query or a near-duplicate "
+            "of the backbone. On a panel this close to the query, run without --curate.",
+            remaining,
+        )
+
+
+def _labels(collection: Path) -> set[str]:
+    return {strip_sequence_extension(p.name) for p in collection_genomes(collection)}
+
+
 def _grow_collection(
     params: FillParams,
     collection: Path,
@@ -286,13 +340,22 @@ def _grow_collection(
     exclude: set[str],
     logger: logging.Logger,
 ) -> tuple[list[RoundResult], dict[str, dict], Path | None]:
-    """Run the build -> scan -> find -> download (-> curate) loop.
+    """Run the (curate ->) build -> scan -> find -> download loop.
 
     Each round rebuilds the MSA from the (growing) collection, scans for coverage
-    gaps, downloads the best new reference per gap, and (when ``--curate``)
-    dereplicates. Stops when the gaps close, when coverage stops improving, when
-    no new reference is found, or at ``max_rounds``. Mutates ``collection``;
-    returns the per-round trace, the curated panel-role table, and the last MSA.
+    gaps and downloads the best new reference per gap. Stops when the gaps close, when
+    coverage stops improving, when no new reference is found, or at ``max_rounds``.
+
+    With ``--curate`` the collection is dereplicated and cleared of the query's siblings
+    *before* a build whenever it holds genomes that have not been through that filter: a
+    collection the user supplied (before round 1), and anything a round downloaded
+    (before the next build). A freshly seeded collection is not curated before round 1 --
+    seeding applies its own sibling filter.
+
+    Whatever the exit, the alignment returned was built from the collection as it stands:
+    if the last round downloaded references, one more MSA (``final.msa.fasta``) is built
+    from them without a further search. Mutates ``collection``; returns the per-round
+    trace, the curated panel-role table, and the last MSA.
     """
     cov = CoverageParams.with_defaults(
         params.window_size, floor=params.coverage_floor, rel_drop=params.coverage_rel_drop,
@@ -303,8 +366,25 @@ def _grow_collection(
     # dropped genome keeps the role it had when removed, even after later rounds).
     panel_rows: dict[str, dict] = {}
     last_msa: Path | None = None
+    built_from: set[str] = set()  # labels the last MSA was built from
     prev_worst: float | None = None
+    # Curation owed before the next build, and the backbone to anchor it on.
+    curate_pending = params.curate and params.collection is not None
+    backbone: Path | None = None
+    pending_round: RoundResult | None = None  # the round whose downloads await curation
+
+    def curate_if_pending() -> None:
+        nonlocal curate_pending, pending_round
+        if not curate_pending:
+            return
+        _curate_round(params, collection, backbone, panel_rows, logger)
+        if pending_round is not None:
+            kept = _labels(collection)
+            pending_round.added = [a for a in pending_round.added if a in kept]
+        curate_pending, pending_round = False, None
+
     for rnd in range(1, params.max_rounds + 1):
+        curate_if_pending()
         msa = params.output / f"round{rnd}.msa.fasta"
         logger.info("=== Round %d: building MSA from %d reference(s) ===",
                     rnd, len(collection_genomes(collection)))
@@ -316,6 +396,7 @@ def _grow_collection(
             logger,
         )
         last_msa = msa
+        built_from = _labels(collection)
         result = compute_similarity(
             str(msa), query_label,
             window_size=params.window_size, window_step=params.window_step,
@@ -353,32 +434,34 @@ def _grow_collection(
         )
         # Pick the backbone from the pre-download (curated, sibling-free) collection
         # so a freshly-downloaded sibling cannot be mistaken for it.
-        backbone = None
         if params.curate:
-            backbone = pick_backbone(
-                params.query, collection_genomes(collection),
-                af_min=params.af_min, logger=logger,
-            )
+            backbone = _curation_backbone(params, collection, logger)
         downloaded = _download(candidates, collection, logger)
         rr.added = [c.hit.accession for c in downloaded]
         if not downloaded:
             logger.info("Stopping: no new references available to add.")
             break
-        if params.curate and backbone is not None:
-            curation = curate_collection_dir(
-                params.query, collection, backbone,
-                ani_margin=params.sibling_margin, af_min=params.af_min,
-                derep_ani=params.derep_ani, logger=logger,
-            )
-            for row in curation.table:
-                panel_rows[row["genome"]] = row
-            dropped = {c.hit.accession for c in downloaded} - {
-                strip_sequence_extension(p.name) for p in collection_genomes(collection)
-            }
-            rr.added = [a for a in rr.added if a not in dropped]
+        curate_pending, pending_round = params.curate, rr
     else:
         logger.info("Reached the maximum of %d round(s).", params.max_rounds)
 
+    # The last round may have downloaded references after its own MSA was built. Curate
+    # them like any other round's, then align whatever the collection now holds, so the
+    # published panel and the reported reference count describe the same set.
+    curate_if_pending()
+    if last_msa is not None and _labels(collection) != built_from:
+        final_msa = params.output / "final.msa.fasta"
+        logger.info("=== Final build: aligning %d reference(s) (no further search) ===",
+                    len(collection_genomes(collection)))
+        build_msa(
+            MsaParams(
+                query=params.query, collection=collection, output=final_msa,
+                aligner=params.aligner, reference=params.reference, threads=params.threads,
+            ),
+            logger,
+        )
+        last_msa = final_msa
+
     return trace, panel_rows, last_msa
 
 
````

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_iterate.py tests/unit/test_input_safety.py tests/unit/test_cli_keep_recombinant.py tests/unit/test_cli_nextclade.py -q`

Expected: all pass (the Step 2 command now gives `27 passed`).

- [ ] **Step 5: Run the fast gate**

Run: `ruff check src tests validation && pytest -m "not requires_binary" -q && mypy src`

Expected: `All checks passed!`, `605 passed, 3 deselected`, `Success: no issues found in 79 source files`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/cli/cmd_fill_references.py src/tessera/discover/iterate.py tests/unit/test_iterate.py
git commit -m "$(cat <<'EOF'
Curate before each build and align the last round's downloads in fill-references

On a max-rounds exit the final downloads were counted but never aligned; --curate
did nothing unless a round both found a gap and downloaded a reference; and
--curate could delete the genome named by --reference. Curation now runs before a
build whenever the collection holds unfiltered genomes, --reference anchors it, and
a final alignment is built when the last round added references. A freshly seeded
collection is still not curated before round 1.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

---
### Task 7: A8 -- MAF conversion keeps every backbone contig and every genome

**Spec item(s):** A8

**Files:**
- Modify: `src/tessera/converters/maf_to_fasta.py` (`maf_to_fasta` gains `ref_contigs`, `expected`, `logger`)
- Modify: `src/tessera/aligners/sibeliaz.py` (pass all three), `src/tessera/aligners/cactus.py` (pass `expected`, `logger`)
- Test: `tests/unit/test_converters.py` (append), `tests/unit/test_aligner_params.py` (append)

**Interfaces:**
- Consumes: `tessera.core.io.read_fasta(path) -> list[tuple[str, str]]`, `tessera.core.errors.OutputError`.
- Produces: `tessera.converters.maf_to_fasta.maf_to_fasta(maf_path, reference, out_path, name_map=None, exclude=None, ref_contigs: Sequence[tuple[str, int]] | None = None, expected: Sequence[str] | None = None, logger: logging.Logger | None = None) -> Path`.
  - `ref_contigs`: `(MAF source name, length)` per backbone contig, in FASTA file order. A backbone source in the MAF that is not listed raises `OutputError`.
  - `expected`: genome labels that must each have a row; one absent from the MAF becomes an all-gap row and is named in a warning on `logger`.
  - Row order: backbone first, then all other labels sorted.
  - This module now imports `OutputError`; Task 9 relies on that import.

The converter infers the backbone's layout and the set of genomes from the MAF alone. A MAF lists only contigs that fall in some block, and only genomes the aligner placed, so (i) contigs are laid out in sorted-name order rather than file order, (ii) a backbone contig no block covers is missing and later coordinates shift, (iii) a genome in no block has no row. Give the converter what the adapter knows: the backbone's contigs in file order and the full genome list.

For `sibeliaz` the MAF source names are the FASTA sequence IDs, so `read_fasta(reference)` supplies `ref_contigs` directly. For `cactus` only `expected` is passed: the sequence names `hal2maf` emits for a pangenome HAL could not be checked (the tool is not installed), so its contig layout is left as it is and documented in Task 10.

The `caplog` tests use a logger outside the `tessera` namespace (see Task 6's note).

- [ ] **Step 1: Write the failing tests**

Append to the end of `tests/unit/test_aligner_params.py`:

```python
def test_sibeliaz_lays_out_the_backbone_as_its_file_and_keeps_every_genome(
    monkeypatch, tmp_path: Path, caplog
) -> None:
    # Backbone: three contigs in the order contig_2, contig_10, contig_3. SibeliaZ gives
    # contig_10 no block (nothing is homologous to it) and gives the panel genome
    # `far` no block at all.
    ref = tmp_path / "ref.fasta"
    ref.write_text(">contig_2\nAAAA\n>contig_10\nGGGG\n>contig_3\nCCCC\n")
    qry = tmp_path / "qry.fasta"
    qry.write_text(">q1\nAAAACCCC\n")
    far = tmp_path / "far.fasta"
    far.write_text(">f1\nTTTTTTTT\n")

    def fake_run(caps, cmd, **kw):
        out_dir = Path(cmd[cmd.index("-o") + 1])
        (out_dir / "alignment.maf").write_text(
            "a\n"
            "s contig_3 0 4 + 4 CCCC\n"
            "s q1 4 4 + 8 CCCC\n"
            "\n"
            "a\n"
            "s contig_2 0 4 + 4 AAAA\n"
            "s q1 0 4 + 8 AAAA\n"
        )
        return ""

    monkeypatch.setattr(sz, "run_tool", fake_run)
    monkeypatch.setattr(sz, "_sibeliaz_invocation", lambda out_dir, logger: ["sibeliaz"])

    # Not under the "tessera" logger: the CLI tests switch its propagation off, and
    # caplog only sees records that reach the root logger.
    log = logging.getLogger("sibeliaz_adapter_test")
    with caplog.at_level(logging.WARNING, logger="sibeliaz_adapter_test"):
        result = sz.SibeliazAligner().align(
            [ref, qry, far], ref, tmp_path / "out", AlignParams(threads=1), log
        )

    rows: dict[str, str] = {}
    name = ""
    for line in result.msa_fasta.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            rows[name] = ""
        else:
            rows[name] += line
    assert rows == {
        "ref": "AAAA----CCCC",   # file order; the uncovered contig keeps its columns
        "far": "------------",   # placed in no block, but still a row
        "qry": "AAAA----CCCC",
    }
    assert "far" in caplog.text
```

Append to the end of `tests/unit/test_converters.py`:

```python
# --- multi-contig backbone layout and unplaced genomes ---------------------------

def test_maf_backbone_contigs_follow_file_order_when_given(tmp_path: Path) -> None:
    # The backbone FASTA lists contig_2 then contig_10. Sorted by name, contig_10 comes
    # first, which puts the backbone in a different order from its own file (and from
    # the minimap2 / mafft backends on the same input).
    maf = tmp_path / "order.maf"
    maf.write_text(
        "a\n"
        "s ref.contig_2 0 4 + 4 AAAA\n"
        "s qry.x 0 4 + 8 AAAA\n"
        "\n"
        "a\n"
        "s ref.contig_10 0 4 + 4 CCCC\n"
        "s qry.x 4 4 + 8 CCCC\n"
    )
    seqs = _read_fasta(maf_to_fasta(
        maf, "ref", tmp_path / "msa.fasta",
        ref_contigs=[("ref.contig_2", 4), ("ref.contig_10", 4)],
    ))
    assert seqs["ref"] == "AAAACCCC"
    assert seqs["qry"] == "AAAACCCC"


def test_maf_backbone_contig_without_a_block_keeps_its_columns(tmp_path: Path) -> None:
    # c2 has no homolog in any genome, so no MAF block mentions it. It is still part of
    # the backbone: dropping it would make the MSA 8 wide and shift c3 from 8-11 to 4-7.
    maf = tmp_path / "missing.maf"
    maf.write_text(
        "a\n"
        "s ref.c1 0 4 + 4 AAAA\n"
        "s qry.x 0 4 + 8 AAAA\n"
        "\n"
        "a\n"
        "s ref.c3 0 4 + 4 CCCC\n"
        "s qry.x 4 4 + 8 CCCC\n"
    )
    seqs = _read_fasta(maf_to_fasta(
        maf, "ref", tmp_path / "msa.fasta",
        ref_contigs=[("ref.c1", 4), ("ref.c2", 4), ("ref.c3", 4)],
    ))
    assert seqs["ref"] == "AAAA----CCCC"
    assert seqs["qry"] == "AAAA----CCCC"


def test_maf_rejects_a_backbone_contig_it_was_not_told_about(tmp_path: Path) -> None:
    import pytest

    from tessera.core.errors import OutputError

    maf = tmp_path / "extra.maf"
    maf.write_text("a\ns ref.c9 0 4 + 4 AAAA\ns qry.x 0 4 + 4 AAAA\n")
    with pytest.raises(OutputError, match="ref.c9"):
        maf_to_fasta(maf, "ref", tmp_path / "msa.fasta", ref_contigs=[("ref.c1", 4)])


def test_maf_genome_without_any_block_gets_an_all_gap_row(tmp_path: Path, caplog) -> None:
    # A panel member too divergent to be placed in any block used to have no row at all,
    # so the panel was silently one genome smaller than the collection.
    import logging

    maf = tmp_path / "vanish.maf"
    maf.write_text("a\ns r1 0 4 + 4 AAAA\ns q1 0 4 + 4 AAAT\n")
    name_map = {"r1": "ref", "q1": "qry", "p1": "divergent_panel"}
    # Not under the "tessera" logger: the CLI tests switch its propagation off, and
    # caplog only sees records that reach the root logger.
    log = logging.getLogger("maf_converter_test")
    with caplog.at_level(logging.WARNING, logger="maf_converter_test"):
        seqs = _read_fasta(maf_to_fasta(
            maf, "ref", tmp_path / "msa.fasta", name_map=name_map,
            expected=["ref", "qry", "divergent_panel"], logger=log,
        ))
    assert seqs == {"ref": "AAAA", "divergent_panel": "----", "qry": "AAAT"}
    assert "divergent_panel" in caplog.text
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_converters.py tests/unit/test_aligner_params.py -q`

Expected: `5 failed, 9 passed`. The failures are:

```
tests/unit/test_converters.py::test_maf_backbone_contigs_follow_file_order_when_given
tests/unit/test_converters.py::test_maf_backbone_contig_without_a_block_keeps_its_columns
tests/unit/test_converters.py::test_maf_rejects_a_backbone_contig_it_was_not_told_about
tests/unit/test_converters.py::test_maf_genome_without_any_block_gets_an_all_gap_row
tests/unit/test_aligner_params.py::test_sibeliaz_lays_out_the_backbone_as_its_file_and_keeps_every_genome
```

They fail because `maf_to_fasta` does not accept `ref_contigs` / `expected` / `logger` (TypeError), and the SibeliaZ adapter produces an 8-column MSA without the `far` row.

- [ ] **Step 3: Implement**

`src/tessera/aligners/cactus.py`:

````diff
diff --git a/src/tessera/aligners/cactus.py b/src/tessera/aligners/cactus.py
index 2b9aa1c..f34e22e 100644
--- a/src/tessera/aligners/cactus.py
+++ b/src/tessera/aligners/cactus.py
@@ -81,7 +81,14 @@ class CactusAligner(Aligner):
         # dotted query filename would otherwise be unfindable downstream).
         name_map = {_sample_name(g): g.stem for g in genomes}
         # Drop the Minigraph-Cactus backbone pseudo-genome so it is not a taxon.
-        maf_to_fasta(maf, ref_name, msa, name_map=name_map, exclude={"_MINIGRAPH_"})
+        # `expected` gives a genome that hal2maf placed in no block an all-gap row
+        # instead of leaving it out. The backbone contig order is not passed: the
+        # sequence names hal2maf emits for a pangenome HAL have not been verified
+        # against the input FASTA, so the layout stays as the MAF implies.
+        maf_to_fasta(
+            maf, ref_name, msa, name_map=name_map, exclude={"_MINIGRAPH_"},
+            expected=[g.stem for g in genomes], logger=logger,
+        )
         return AlignResult(msa_fasta=msa, native_format=hal)
 
 
````

`src/tessera/aligners/sibeliaz.py`:

````diff
diff --git a/src/tessera/aligners/sibeliaz.py b/src/tessera/aligners/sibeliaz.py
index 0ad80be..3d1568d 100644
--- a/src/tessera/aligners/sibeliaz.py
+++ b/src/tessera/aligners/sibeliaz.py
@@ -16,7 +16,7 @@ from pathlib import Path
 from ..converters.maf_to_fasta import maf_to_fasta
 from ..core.binaries import BinarySpec
 from ..core.errors import OutputError, UserInputError
-from ..core.io import normalize_reference
+from ..core.io import normalize_reference, read_fasta
 from ..core.plugins import ToolCapabilities
 from ..core.process import run_tool
 from .base import Aligner, AlignParams, AlignResult
@@ -112,7 +112,15 @@ class SibeliazAligner(Aligner):
         # genome filenames; build the seqid -> genome-stem map for the converter.
         name_map = _build_seqid_map(genomes)
         msa = out_dir / "msa.fasta"
-        maf_to_fasta(maf, reference.stem, msa, name_map=name_map)
+        # The MAF names only the backbone contigs that fall in a block, in no particular
+        # order, and only the genomes it placed. Hand the converter the backbone's own
+        # contig order and the full genome list so neither is inferred from the MAF.
+        maf_to_fasta(
+            maf, reference.stem, msa, name_map=name_map,
+            ref_contigs=[(seqid, len(seq)) for seqid, seq in read_fasta(reference)],
+            expected=[g.stem for g in genomes],
+            logger=logger,
+        )
         return AlignResult(msa_fasta=msa, native_format=maf)
 
 
````

`src/tessera/converters/maf_to_fasta.py`:

````diff
diff --git a/src/tessera/converters/maf_to_fasta.py b/src/tessera/converters/maf_to_fasta.py
index f5dd6d0..b81b3a9 100644
--- a/src/tessera/converters/maf_to_fasta.py
+++ b/src/tessera/converters/maf_to_fasta.py
@@ -17,9 +17,12 @@ reports positions that line up with the reference genome.
 
 from __future__ import annotations
 
+import logging
+from collections.abc import Sequence
 from dataclasses import dataclass
 from pathlib import Path
 
+from ..core.errors import OutputError
 from ..core.io import write_fasta_record
 
 _COMPLEMENTS = bytes.maketrans(
@@ -53,6 +56,9 @@ def maf_to_fasta(
     out_path: str | Path,
     name_map: dict[str, str] | None = None,
     exclude: set[str] | None = None,
+    ref_contigs: Sequence[tuple[str, int]] | None = None,
+    expected: Sequence[str] | None = None,
+    logger: logging.Logger | None = None,
 ) -> Path:
     """Project a MAF onto ``reference`` coordinates as an MSA-FASTA.
 
@@ -63,6 +69,18 @@ def maf_to_fasta(
 
     ``exclude`` drops genomes by label, e.g. ``{"_MINIGRAPH_"}`` to remove the
     Minigraph-Cactus backbone pseudo-genome so it is not emitted as a taxon.
+
+    ``ref_contigs`` gives the backbone's contigs as ``(MAF source name, length)`` in
+    the order of its FASTA file. The MAF alone cannot supply this: it lists only
+    contigs that fall in some block, in no particular order. With it, the backbone is
+    laid out as its file is, and a contig no block covers keeps its (all-gap) columns
+    instead of vanishing and shifting everything after it. Without it the contigs seen
+    in the MAF are laid out in sorted-name order.
+
+    ``expected`` lists every genome label that should have a row. A genome the aligner
+    placed in no block is absent from the MAF; it is written as an all-gap row and
+    named in a warning on ``logger``, so the panel is never silently smaller than the
+    collection.
     """
     maf_path = Path(maf_path)
     out_path = Path(out_path)
@@ -81,30 +99,49 @@ def maf_to_fasta(
     # each distinct reference source name is one contig (length = its src_size),
     # placed at a cumulative offset. Collapsing them into a single contig's
     # coordinate space would make later contigs overwrite earlier ones.
-    ref_contigs: dict[str, int] = {}  # source name -> contig length
+    seen_contigs: dict[str, int] = {}  # source name -> contig length
     species: set[str] = set()
     for block in blocks:
         for row in block:
             label = genome_of(row.name)
             species.add(label)
             if label == ref_key:
-                ref_contigs.setdefault(row.name, row.src_size)
+                seen_contigs.setdefault(row.name, row.src_size)
     species -= exclude
 
-    if not ref_contigs:
+    if not seen_contigs:
         raise ValueError(
             f"MAF projection onto reference '{ref_key}' found no reference rows. "
             f"Check that the reference label matches the MAF/name_map sequence names."
         )
 
+    if ref_contigs is not None:
+        unknown = sorted(set(seen_contigs) - {name for name, _ in ref_contigs})
+        if unknown:
+            raise OutputError(
+                f"{maf_path} aligns backbone sequence(s) {', '.join(unknown)} that are "
+                f"not in the backbone '{ref_key}' as staged. The aligner's sequence "
+                "names do not match the input FASTA."
+            )
+        layout = list(ref_contigs)
+    else:
+        layout = sorted(seen_contigs.items())
     ref_offsets: dict[str, int] = {}
     ref_length = 0
-    for name in sorted(ref_contigs):
+    for name, length in layout:
         ref_offsets[name] = ref_length
-        ref_length += ref_contigs[name]
-
+        ref_length += length
+
+    unplaced = sorted(set(expected or ()) - species - exclude - {ref_key})
+    if unplaced and logger is not None:
+        logger.warning(
+            "%d genome(s) share no alignment block with the backbone '%s' and are "
+            "written as all-gap rows: %s. They contribute nothing to the scan; they "
+            "may be too divergent for this aligner.",
+            len(unplaced), ref_key, ", ".join(unplaced),
+        )
     species.discard(ref_key)
-    ordered_species = [ref_key, *sorted(species)]
+    ordered_species = [ref_key, *sorted(species | set(unplaced))]
     out: dict[str, bytearray] = {s: bytearray(b"-" * ref_length) for s in ordered_species}
 
     for block in blocks:
````

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_converters.py tests/unit/test_aligner_params.py tests/unit/test_plugins.py -q`

Expected: all pass (the Step 2 command now gives `14 passed`).

- [ ] **Step 5: Run the fast gate**

Run: `ruff check src tests validation && pytest -m "not requires_binary" -q && mypy src`

Expected: `All checks passed!`, `610 passed, 3 deselected`, `Success: no issues found in 79 source files`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/aligners/cactus.py src/tessera/aligners/sibeliaz.py src/tessera/converters/maf_to_fasta.py tests/unit/test_aligner_params.py tests/unit/test_converters.py
git commit -m "$(cat <<'EOF'
Keep uncovered contigs and unplaced genomes in MAF-derived alignments

The MAF converter inferred the backbone layout and the genome set from the MAF,
which omits contigs and genomes that fall in no block. The sibeliaz adapter now
passes the backbone's contigs in file order, and both MAF backends pass the full
genome list, so a missing genome becomes an all-gap row with a warning.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

---
### Task 8: A9 -- progressiveMauve rows named from the staged label

**Spec item(s):** A9

**Files:**
- Modify: `src/tessera/aligners/progressivemauve.py` (`align_query` return value, `_concatenate`)
- Test: `tests/unit/test_aligner_backends.py` (append)

**Interfaces:**
- Consumes: `tessera.converters.xmfa_to_fasta.xmfa_to_fasta(xmfa_path, reference_name, flank, out_path, reference_length=None)`, which writes the reference row first and then one row per other sequence in the XMFA.
- Produces: `tessera.aligners.progressivemauve._concatenate(per_query: list[tuple[str, Path]], reference_label: str, out_path: Path) -> None`. Raises `OutputError` if a per-query FASTA does not hold exactly two records.

**This task cannot be verified against the real tool: progressiveMauve is not installed here.** The tests drive the adapter with a stubbed `run_tool` that writes the XMFA the way the existing test (`test_progressivemauve_concatenation_is_rectangular`) does, with the command-line paths as sequence names. The change does not depend on that detail any more, which is its point, but run one real alignment before relying on this backend (see Task 10).

`_concatenate` named each row `Path(name).stem`, where `name` is the path progressiveMauve echoed -- the *resolved* target of the staged symlink. Rows are now identified by position (reference first, query second) and labelled with the staged stem the adapter already holds.

- [ ] **Step 1: Write the failing tests**

Append to the end of `tests/unit/test_aligner_backends.py`:

```python
# --- progressiveMauve row names come from the staged labels ----------------------
def _fake_mauve(caps, cmd, **kw):
    """Write the XMFA the way progressiveMauve does: sequence names are the paths it
    was given on the command line (which the adapter passes resolved)."""
    cmd = [str(c) for c in cmd]
    xmfa = Path(cmd[cmd.index("--output") + 1])
    ref_path, query_path = cmd[-2], cmd[-1]
    ref_seq = "".join(
        line for line in Path(ref_path).read_text().splitlines() if not line.startswith(">")
    )
    qry_seq = "".join(
        line for line in Path(query_path).read_text().splitlines() if not line.startswith(">")
    )
    xmfa.write_text(
        f"#Sequence1File\t{ref_path}\n#Sequence2File\t{query_path}\n"
        f"> 1:1-8 + {ref_path}\n{ref_seq}\n> 2:1-8 + {query_path}\n{qry_seq}\n=\n"
    )
    return ""


def _mauve_rows(monkeypatch, staged: list[Path], out_dir: Path) -> list[tuple[str, str]]:
    monkeypatch.setattr(pm_mod, "run_tool", _fake_mauve)
    result = pm_mod.ProgressiveMauveAligner().align(
        staged, staged[0], out_dir, AlignParams(threads=1), _LOG
    )
    rows: list[tuple[str, str]] = []
    for line in result.msa_fasta.read_text().splitlines():
        if line.startswith(">"):
            rows.append((line[1:], ""))
        else:
            rows[-1] = (rows[-1][0], rows[-1][1] + line)
    return rows


def test_progressivemauve_names_rows_by_staged_label_not_symlink_target(
    monkeypatch, tmp_path: Path
) -> None:
    # Staged genomes are symlinks named <label>.fasta. Two of them point at files that
    # are both called genome.fna, and one target is named like the backbone's label.
    store = tmp_path / "store"
    for sub in ("x", "y", "z"):
        (store / sub).mkdir(parents=True)
    (store / "x" / "genome.fna").write_text(">r\nACGTACGT\n")
    (store / "y" / "genome.fna").write_text(">b\nACGTACGA\n")
    (store / "z" / "aRef.fna").write_text(">c\nTTTTTTTT\n")
    stage = tmp_path / "stage"
    stage.mkdir()
    staged = []
    for label, target in (("aRef", "x/genome.fna"), ("panelB", "y/genome.fna"),
                          ("panelC", "z/aRef.fna")):
        link = stage / f"{label}.fasta"
        link.symlink_to((store / target).resolve())
        staged.append(link)

    rows = _mauve_rows(monkeypatch, staged, tmp_path / "out")

    assert rows == [("aRef", "ACGTACGT"), ("panelB", "ACGTACGA"), ("panelC", "TTTTTTTT")]


def test_progressivemauve_tolerates_whitespace_in_the_path(monkeypatch, tmp_path: Path) -> None:
    stage = tmp_path / "dir with space"
    stage.mkdir()
    staged = []
    for label, seq in (("ref", "ACGTACGT"), ("qryA", "ACGTACGA"), ("qryB", "ACGTACGC")):
        path = stage / f"{label}.fasta"
        path.write_text(f">{label}\n{seq}\n")
        staged.append(path)

    rows = _mauve_rows(monkeypatch, staged, tmp_path / "out")

    assert [name for name, _ in rows] == ["ref", "qryA", "qryB"]


def test_progressivemauve_rejects_a_projection_that_is_not_pairwise(
    monkeypatch, tmp_path: Path
) -> None:
    from tessera.core.errors import OutputError

    def one_sequence_xmfa(caps, cmd, **kw):
        cmd = [str(c) for c in cmd]
        xmfa = Path(cmd[cmd.index("--output") + 1])
        ref_path = cmd[-2]
        xmfa.write_text(f"#Sequence1File\t{ref_path}\n> 1:1-8 + {ref_path}\nACGTACGT\n=\n")
        return ""

    monkeypatch.setattr(pm_mod, "run_tool", one_sequence_xmfa)
    genomes = []
    for label in ("ref", "qryA", "qryB"):
        path = tmp_path / f"{label}.fasta"
        path.write_text(f">{label}\nACGTACGT\n")
        genomes.append(path)
    with pytest.raises(OutputError, match=r"qryA\.fa.*found 1 record"):
        pm_mod.ProgressiveMauveAligner().align(
            genomes, genomes[0], tmp_path / "out", AlignParams(threads=1), _LOG
        )
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/unit/test_aligner_backends.py -q`

Expected: `3 failed, 11 passed`. The failures are:

```
tests/unit/test_aligner_backends.py::test_progressivemauve_names_rows_by_staged_label_not_symlink_target
tests/unit/test_aligner_backends.py::test_progressivemauve_tolerates_whitespace_in_the_path
tests/unit/test_aligner_backends.py::test_progressivemauve_rejects_a_projection_that_is_not_pairwise
```

They fail because rows come out named `genome`/`dir`, a genome whose target is named like the backbone is dropped, and a one-sequence projection is accepted.

- [ ] **Step 3: Implement**

`src/tessera/aligners/progressivemauve.py`:

````diff
diff --git a/src/tessera/aligners/progressivemauve.py b/src/tessera/aligners/progressivemauve.py
index 8d75881..c0f8b5c 100644
--- a/src/tessera/aligners/progressivemauve.py
+++ b/src/tessera/aligners/progressivemauve.py
@@ -17,6 +17,7 @@ from pathlib import Path
 
 from ..converters.xmfa_to_fasta import xmfa_to_fasta
 from ..core.binaries import BinarySpec
+from ..core.errors import OutputError
 from ..core.executors import parallel_map
 from ..core.io import normalize_reference, read_fasta, write_fasta_record
 from ..core.plugins import ToolCapabilities
@@ -68,7 +69,7 @@ class ProgressiveMauveAligner(Aligner):
         # resolved a yet-unexplained progressiveMauve error on some systems.
         workers = 1 if params.flag("single") else params.threads
 
-        def align_query(query: Path) -> Path:
+        def align_query(query: Path) -> tuple[str, Path]:
             stem = query.stem
             xmfa = xmfa_dir / f"{stem}.xmfa"
             fa = xmfa_dir / f"{stem}.fa"
@@ -79,25 +80,37 @@ class ProgressiveMauveAligner(Aligner):
                 log_prefix=f"progressivemauve:{stem}",
             )
             xmfa_to_fasta(xmfa, ref_arg, 0, fa, reference_length=ref_length)
-            return fa
+            return stem, fa
 
-        per_query_fastas = parallel_map(align_query, queries, workers, logger=logger)
+        per_query = parallel_map(align_query, queries, workers, logger=logger)
 
         msa = out_dir / "msa.fasta"
-        _concatenate(per_query_fastas, reference, msa)
+        _concatenate(per_query, reference.stem, msa)
         return AlignResult(msa_fasta=msa)
 
 
-def _concatenate(per_query_fastas: list[Path], reference: Path, out_path: Path) -> None:
-    """Write the reference row once, then each query row; leaf names are stems."""
-    ref_stem = reference.stem
-    written_ref = False
+def _concatenate(
+    per_query: list[tuple[str, Path]], reference_label: str, out_path: Path
+) -> None:
+    """Write the reference row once, then one row per query, named by staged label.
+
+    Each per-query FASTA holds exactly two records, in a fixed order: the reference
+    projection, then the query's. The names inside those files are the paths
+    progressiveMauve echoed from its command line -- resolved symlink targets, which can
+    carry an unrecognised extension, whitespace, or a basename shared with another
+    genome -- so rows are identified by position and labelled from the staged file,
+    never from those names.
+    """
     with open(out_path, "w") as out:
-        for fa in per_query_fastas:
-            for name, seq in read_fasta(fa):
-                leaf = Path(name).stem
-                if leaf == ref_stem:
-                    if written_ref:
-                        continue
-                    written_ref = True
-                write_fasta_record(out, leaf, seq)
+        for i, (label, fa) in enumerate(per_query):
+            records = read_fasta(fa)
+            if len(records) != 2:
+                raise OutputError(
+                    f"Expected a reference row and one query row in {fa} (the "
+                    f"projection of '{label}' onto '{reference_label}'), found "
+                    f"{len(records)} record(s). progressiveMauve's output is not a "
+                    "pairwise alignment."
+                )
+            if i == 0:
+                write_fasta_record(out, reference_label, records[0][1])
+            write_fasta_record(out, label, records[1][1])
````

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/unit/test_aligner_backends.py tests/unit/test_aligner_params.py -q`

Expected: all pass (the Step 2 command now gives `14 passed`).

- [ ] **Step 5: Run the fast gate**

Run: `ruff check src tests validation && pytest -m "not requires_binary" -q && mypy src`

Expected: `All checks passed!`, `613 passed, 3 deselected`, `Success: no issues found in 79 source files`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/aligners/progressivemauve.py tests/unit/test_aligner_backends.py
git commit -m "$(cat <<'EOF'
Name progressiveMauve rows from the staged genome label

Rows were named from the resolved path echoed in the XMFA, so symlinked
collections, unrecognised extensions and whitespace in a path gave wrong or
duplicate names. Rows are now taken by position and labelled from the staged file.
Checked with a simulated XMFA only; the tool was not available.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

---
### Task 9: A10 -- whitespace in sequence lines, and malformed input reported by file

**Spec item(s):** A10

**Files:**
- Modify: `src/tessera/core/io.py` (`read_fasta`; new `_has_sequence_whitespace`, `_write_clean`; `_stage_one`)
- Modify: `src/tessera/aligners/sibeliaz.py` (`_build_seqid_map`)
- Modify: `src/tessera/converters/maf_to_fasta.py` (truncated block), `xmfa_to_fasta.py` (missing reference), `mafft_merge.py` (empty output)
- Test: `tests/unit/test_input_safety.py` (append), `tests/unit/test_io.py` (append), `tests/integration/test_whitespace_genomes.py` (create; `requires_binary`)

**Interfaces:**
- Consumes: `OutputError` already imported in `maf_to_fasta.py` by Task 7. If this task is done before Task 7, add `from ..core.errors import OutputError` to that file.
- Produces:
  - `read_fasta` drops all whitespace from sequence lines and reads a nameless header (`>` or `> `) as `""`.
  - `tessera.core.io._has_sequence_whitespace(source: Path) -> bool`
  - `tessera.core.io._write_clean(src: Iterable[bytes], dst: BinaryIO) -> None`
  - `_stage_one` stages a cleaned copy (not a symlink) for a plain FASTA with whitespace in its sequence lines, and always cleans while decompressing a `.gz`.

**Deviation from the spec, found while verifying it.** The spec's design is "strip whitespace from sequence lines" in `read_fasta`. Done alone, that makes things worse for `minimap2`: Tessera would count a reference line as 60 bases while minimap2, reading the raw file, counts the trailing space as a 61st, and every column after the first line is shifted (measured: identity 0.25 over the second half of a 6 kb genome, with no error). The reader and the aligner must see the same sequence. So the cleaning also happens where every aligner's input passes -- staging: a genome with whitespace in its sequence lines is staged as a cleaned copy instead of a symlink, and a warning names it. `tests/integration/test_whitespace_genomes.py` checks both `mafft` and `minimap2` end to end; it fails on `main` for both (ragged rows).

A second, smaller deviation: the spec asks for a `UserInputError` on a header line `> `. `>` alone already reads as an unnamed record and `reassort` relies on that (it names such a record `segment_N`), so `> ` now reads the same way. The `sibeliaz` backend, which identifies genomes by sequence ID, rejects a nameless record with the file and line.

Error types: a truncated MAF block and an empty MAFFT result are tool output, so they raise `OutputError`; the XMFA check follows the existing guards in that file and raises `UserInputError`. All are `TesseraError`, so the CLI reports them cleanly.

- [ ] **Step 1: Write the failing tests**

Create `tests/integration/test_whitespace_genomes.py`:

```python
"""Genomes with whitespace in their sequence lines, through each pairwise aligner.

Needs the aligner binaries. A reference saved with a trailing space on every line used
to give a ragged MSA with mafft and, once the reader ignored the spaces, a silently
shifted one with minimap2 (which counts them as bases).
"""

from __future__ import annotations

import logging
import random
import shutil
from pathlib import Path

import pytest

from tessera.core.io import read_fasta
from tessera.msa.build import MsaParams, build_msa

_LOG = logging.getLogger("tessera.test")

pytestmark = pytest.mark.requires_binary


@pytest.mark.parametrize("aligner", ["mafft", "minimap2"])
def test_reference_with_trailing_spaces_aligns_in_register(aligner: str, tmp_path: Path) -> None:
    if shutil.which(aligner) is None:
        pytest.skip(f"{aligner} not installed")
    rng = random.Random(1)
    ref = "".join(rng.choice("ACGT") for _ in range(6000))
    qry = "".join(
        rng.choice([b for b in "ACGT" if b != c]) if rng.random() < 0.03 else c for c in ref
    )
    other = "".join(
        rng.choice([b for b in "ACGT" if b != c]) if rng.random() < 0.05 else c for c in ref
    )
    collection = tmp_path / "collection"
    collection.mkdir()
    spaced = "".join(ref[i:i + 60] + " \n" for i in range(0, len(ref), 60))
    (collection / "ref.fasta").write_text(f">r\n{spaced}")
    (collection / "other.fasta").write_text(f">o\n{other}\n")
    query = tmp_path / "query.fasta"
    query.write_text(f">q\n{qry}\n")

    msa = build_msa(
        MsaParams(query=query, collection=collection, output=tmp_path / "msa.fasta",
                  aligner=aligner, reference="ref", threads=1),
        _LOG,
    )
    rows = dict(read_fasta(msa))
    assert {len(row) for row in rows.values()} == {6000}
    tail_pairs = [
        (a, b) for a, b in zip(rows["ref"][3000:].upper(), rows["query"][3000:].upper(),
                               strict=True)
        if a in "ACGT" and b in "ACGT"
    ]
    assert len(tail_pairs) > 2500
    assert sum(a == b for a, b in tail_pairs) / len(tail_pairs) > 0.95
```

Append to the end of `tests/unit/test_input_safety.py`:

```python
# --- malformed records and tool output ---------------------------------------

def test_read_fasta_drops_whitespace_inside_sequence_lines(tmp_path: Path) -> None:
    """A trailing space or tab is not a base. Kept, it widened the backbone row and made
    the mafft MSA ragged (reference 6100 columns, query 6000)."""
    path = tmp_path / "ws.fasta"
    path.write_text(">ref desc\nACGT \nAC GT\t\n\nTTAA\n")
    assert read_fasta(path) == [("ref", "ACGTACGTTTAA")]


def test_read_fasta_accepts_a_header_with_no_name(tmp_path: Path) -> None:
    # ">" alone already read as an unnamed record; "> " raised IndexError instead.
    path = tmp_path / "blank.fasta"
    path.write_text("> \nACGT\n>\nTTAA\n")
    assert read_fasta(path) == [("", "ACGT"), ("", "TTAA")]


def test_sibeliaz_rejects_a_record_with_no_sequence_id(tmp_path: Path) -> None:
    path = tmp_path / "g1.fasta"
    path.write_text(">ok\nACGT\n> \nACGT\n")
    with pytest.raises(UserInputError, match=r"g1\.fasta.*line 3"):
        sibeliaz._build_seqid_map([path])


def test_truncated_maf_row_is_reported_against_the_file(tmp_path: Path) -> None:
    from tessera.converters.maf_to_fasta import maf_to_fasta

    maf = tmp_path / "cut.maf"
    maf.write_text("a\ns ref.c 0 8 + 8 ACGTACGT\ns qry.x 0 8 + 8 ACGT")
    with pytest.raises(OutputError, match=r"cut\.maf.*qry\.x"):
        maf_to_fasta(maf, "ref", tmp_path / "msa.fasta")


def test_xmfa_without_the_reference_is_reported_against_the_file(tmp_path: Path) -> None:
    from tessera.converters.xmfa_to_fasta import xmfa_to_fasta

    xmfa = tmp_path / "other.xmfa"
    xmfa.write_text(
        "#Sequence1File\t/data/a.fasta\n#Sequence2File\t/data/b.fasta\n"
        "> 1:1-4 + /data/a.fasta\nACGT\n> 2:1-4 + /data/b.fasta\nACGT\n=\n"
    )
    with pytest.raises(UserInputError, match=r"other\.xmfa.*/data/ref\.fasta"):
        xmfa_to_fasta(xmfa, "/data/ref.fasta", 0, tmp_path / "out.fasta")


def test_empty_mafft_output_is_a_tessera_error(tmp_path: Path) -> None:
    from tessera.converters.mafft_merge import merge_added_fragments

    empty = tmp_path / "empty.aln.fasta"
    empty.write_text("")
    with pytest.raises(OutputError, match=r"empty\.aln\.fasta"):
        merge_added_fragments(empty)
```

Append to the end of `tests/unit/test_io.py`:

```python
def test_stage_genomes_cleans_whitespace_out_of_sequence_lines(tmp_path: Path) -> None:
    """The aligner reads the staged file; Tessera reads the same genome through
    read_fasta. If only one of them ignores a trailing space the two disagree about the
    genome's length -- minimap2 counted the spaces and every later column was shifted."""
    coll = tmp_path / "coll"
    coll.mkdir()
    (coll / "clean.fasta").write_text(">c desc\nACGT\nTTAA\n")
    (coll / "spaces.fasta").write_text(">s desc\nACGT \nTT AA\t\n")
    (coll / "crlf.fasta").write_bytes(b">w desc\r\nACGT\r\nTTAA\r\n")
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGT\n")

    staged, _ = stage_genomes(query, coll, tmp_path / "stage", _LOG)
    by_stem = {p.stem: p for p in staged}

    assert by_stem["clean"].is_symlink()  # untouched input is still linked, not copied
    for stem, header in (("spaces", ">s desc"), ("crlf", ">w desc")):
        assert not by_stem[stem].is_symlink()
        assert by_stem[stem].read_text() == f"{header}\nACGT\nTTAA\n"
    # The user's own files are not modified.
    assert (coll / "spaces.fasta").read_text() == ">s desc\nACGT \nTT AA\t\n"


def test_stage_genomes_cleans_a_gzipped_genome(tmp_path: Path) -> None:
    coll = tmp_path / "coll"
    coll.mkdir()
    with gzip.open(coll / "z.fasta.gz", "wt") as fo:
        fo.write(">z\nACGT \n\nTTAA\n")
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGT\n")
    staged, _ = stage_genomes(query, coll, tmp_path / "stage", _LOG)
    assert {p.stem: p for p in staged}["z"].read_text() == ">z\nACGT\nTTAA\n"
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin" pytest tests/unit/test_input_safety.py tests/unit/test_io.py tests/integration/test_whitespace_genomes.py -q`

Expected: `10 failed, 34 passed`. The failures are:

```
tests/unit/test_input_safety.py::test_read_fasta_drops_whitespace_inside_sequence_lines
tests/unit/test_input_safety.py::test_read_fasta_accepts_a_header_with_no_name
tests/unit/test_input_safety.py::test_sibeliaz_rejects_a_record_with_no_sequence_id
tests/unit/test_input_safety.py::test_truncated_maf_row_is_reported_against_the_file
tests/unit/test_input_safety.py::test_xmfa_without_the_reference_is_reported_against_the_file
tests/unit/test_input_safety.py::test_empty_mafft_output_is_a_tessera_error
tests/unit/test_io.py::test_stage_genomes_cleans_whitespace_out_of_sequence_lines
tests/unit/test_io.py::test_stage_genomes_cleans_a_gzipped_genome
tests/integration/test_whitespace_genomes.py::test_reference_with_trailing_spaces_aligns_in_register[mafft]
tests/integration/test_whitespace_genomes.py::test_reference_with_trailing_spaces_aligns_in_register[minimap2]
```

They fail because whitespace is kept as sequence, `> ` raises IndexError, staging links the file as it is, and the converters raise bare IndexError / KeyError / ValueError.

- [ ] **Step 3: Implement**

`src/tessera/aligners/sibeliaz.py`:

````diff
diff --git a/src/tessera/aligners/sibeliaz.py b/src/tessera/aligners/sibeliaz.py
index 3d1568d..7884c4c 100644
--- a/src/tessera/aligners/sibeliaz.py
+++ b/src/tessera/aligners/sibeliaz.py
@@ -156,9 +156,16 @@ def _build_seqid_map(genomes) -> dict[str, str]:
     for genome in genomes:
         stem = genome.stem
         with open(genome) as fo:
-            for line in fo:
+            for lineno, line in enumerate(fo, start=1):
                 if line.startswith(">"):
-                    seqid = line[1:].split()[0]
+                    tokens = line[1:].split()
+                    if not tokens:
+                        raise UserInputError(
+                            f"{genome} has a record with no sequence ID (line {lineno}). "
+                            "The sibeliaz backend identifies genomes by sequence ID, so "
+                            "every record needs a name after '>'."
+                        )
+                    seqid = tokens[0]
                     owner = name_map.setdefault(seqid, stem)
                     if owner != stem:
                         # SibeliaZ names alignment rows by sequence ID alone, so two
````

`src/tessera/converters/maf_to_fasta.py`:

````diff
diff --git a/src/tessera/converters/maf_to_fasta.py b/src/tessera/converters/maf_to_fasta.py
index b81b3a9..39ac7cd 100644
--- a/src/tessera/converters/maf_to_fasta.py
+++ b/src/tessera/converters/maf_to_fasta.py
@@ -148,6 +148,17 @@ def maf_to_fasta(
         ref_row = next((r for r in block if genome_of(r.name) == ref_key), None)
         if ref_row is None:
             continue
+        # Every row of a block has the same aligned width. A shorter one means the file
+        # was cut off mid-write (a full disk, a killed aligner); indexing into it below
+        # would fail with a bare IndexError several frames from the file's name.
+        for row in block:
+            if len(row.text) != len(ref_row.text):
+                raise OutputError(
+                    f"Truncated alignment block in {maf_path}: row '{row.name}' has "
+                    f"{len(row.text)} column(s), the backbone row '{ref_row.name}' has "
+                    f"{len(ref_row.text)}. The aligner's output looks incomplete; check "
+                    "that it finished."
+                )
         contig_offset = ref_offsets[ref_row.name]
         if ref_row.strand == "-":
             # Reverse-complement the whole block into forward-reference orientation.
````

`src/tessera/converters/mafft_merge.py`:

````diff
diff --git a/src/tessera/converters/mafft_merge.py b/src/tessera/converters/mafft_merge.py
index 3abbcbb..17eecc8 100644
--- a/src/tessera/converters/mafft_merge.py
+++ b/src/tessera/converters/mafft_merge.py
@@ -13,6 +13,7 @@ from __future__ import annotations
 
 from pathlib import Path
 
+from ..core.errors import OutputError
 from ..core.io import read_fasta
 
 
@@ -24,7 +25,10 @@ def merge_added_fragments(aligned_path: str | Path) -> tuple[str, str]:
     """
     records = read_fasta(aligned_path)
     if not records:
-        raise ValueError(f"Empty MAFFT alignment: {aligned_path}")
+        raise OutputError(
+            f"MAFFT wrote an empty alignment at {aligned_path}. Check that mafft "
+            "completed and that the genome being added is a nucleotide FASTA."
+        )
     reference_row = records[0][1]
     width = len(reference_row)
     merged = bytearray(b"-" * width)
````

`src/tessera/converters/xmfa_to_fasta.py`:

````diff
diff --git a/src/tessera/converters/xmfa_to_fasta.py b/src/tessera/converters/xmfa_to_fasta.py
index 5821200..85bc44e 100644
--- a/src/tessera/converters/xmfa_to_fasta.py
+++ b/src/tessera/converters/xmfa_to_fasta.py
@@ -112,6 +112,13 @@ def xmfa_to_fasta(
                 seq[curr_pos : curr_pos + length_of_line] = stripped
                 curr_pos += length_of_line
 
+    if reference_name not in name2num:
+        listed = ", ".join(sorted(name2num)) or "none"
+        raise UserInputError(
+            f"{xmfa_path} does not list the reference {reference_name} among its "
+            f"sequence files (found: {listed}). The aligner's output does not belong "
+            "to this reference, or its header is incomplete."
+        )
     reference_num = name2num[reference_name]
 
     # Output width: the full reference length when known, otherwise the furthest
````

`src/tessera/core/io.py`:

````diff
diff --git a/src/tessera/core/io.py b/src/tessera/core/io.py
index c6dfa58..5950793 100644
--- a/src/tessera/core/io.py
+++ b/src/tessera/core/io.py
@@ -13,9 +13,9 @@ import gzip
 import logging
 import re
 import shutil
-from collections.abc import Sequence
+from collections.abc import Iterable, Sequence
 from pathlib import Path
-from typing import TextIO
+from typing import BinaryIO, TextIO
 
 from .errors import UserInputError
 
@@ -35,20 +35,25 @@ def read_fasta(path: str | Path) -> list[tuple[str, str]]:
 
     A ``.gz`` input is read through :mod:`gzip`: staging already accepts compressed
     genomes, so the readers that look inside a query must accept them too.
+
+    Whitespace inside a sequence line is dropped: a trailing space or tab is not a
+    base, and kept it lengthens the sequence (the backbone row then no longer matches
+    the rows aligned to it). A header with no name -- ``>`` or ``> `` -- reads as an
+    unnamed record (``""``).
     """
     records: list[tuple[str, str]] = []
     name: str | None = None
     seq: list[str] = []
     with _open_text(path) as fo:
         for line in fo:
-            line = line.rstrip("\r\n")
             if line.startswith(">"):
                 if name is not None:
                     records.append((name, "".join(seq)))
-                name = line[1:].split()[0] if len(line) > 1 else ""
+                tokens = line[1:].split()
+                name = tokens[0] if tokens else ""
                 seq = []
             else:
-                seq.append(line)
+                seq.append("".join(line.split()))
     if name is not None:
         records.append((name, "".join(seq)))
     return records
@@ -164,14 +169,58 @@ def _require_fasta(source: Path) -> None:
         )
 
 
+_WHITESPACE = re.compile(rb"\s")
+
+
+def _has_sequence_whitespace(source: Path) -> bool:
+    """True when a plain FASTA carries whitespace inside a sequence line.
+
+    A trailing space or tab, a space within the line, or a carriage return (CRLF line
+    endings). Blank lines and the header's own spaces do not count.
+    """
+    with open(source, "rb") as fo:
+        for line in fo:
+            if line.startswith(b">"):
+                continue
+            body = line[:-1] if line.endswith(b"\n") else line
+            if body and _WHITESPACE.search(body):
+                return True
+    return False
+
+
+def _write_clean(src: Iterable[bytes], dst: BinaryIO) -> None:
+    """Copy a FASTA, dropping whitespace from sequence lines (and blank lines)."""
+    for line in src:
+        if line.startswith(b">"):
+            dst.write(line.rstrip(b"\r\n") + b"\n")
+        else:
+            body = b"".join(line.split())
+            if body:
+                dst.write(body + b"\n")
+
+
 def _stage_one(source: Path, target_dir: Path, logger: logging.Logger) -> Path:
-    """Place one genome into ``target_dir`` as ``<label>.fasta``; return the path."""
+    """Place one genome into ``target_dir`` as ``<label>.fasta``; return the path.
+
+    The staged file is what the aligner reads, while Tessera reads the same genome
+    through :func:`read_fasta`, which ignores whitespace in sequence lines. The two
+    must agree on the genome's length, and not every aligner ignores a trailing space
+    (minimap2 counts it as a base), so a genome that carries such whitespace is staged
+    as a cleaned copy rather than a link to the original.
+    """
     label = strip_sequence_extension(source.name)
     target = target_dir / f"{label}.fasta"
     if source.name.endswith(".gz"):
         logger.debug("Decompressing %s -> %s", source, target)
         with gzip.open(source, "rb") as src, open(target, "wb") as dst:
-            shutil.copyfileobj(src, dst)
+            _write_clean(src, dst)
+    elif _has_sequence_whitespace(source):
+        logger.warning(
+            "%s has whitespace inside its sequence lines (trailing spaces, tabs or "
+            "Windows line endings); aligning a cleaned copy.", source,
+        )
+        with open(source, "rb") as src, open(target, "wb") as dst:
+            _write_clean(src, dst)
     else:
         logger.debug("Linking %s -> %s", source, target)
         target.symlink_to(source.resolve())
````

- [ ] **Step 4: Run the tests to verify they pass**

Run: `PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin" pytest tests/unit/test_input_safety.py tests/unit/test_io.py tests/integration/test_whitespace_genomes.py -q`

Expected: all pass (the Step 2 command now gives `44 passed`).

- [ ] **Step 5: Run the fast gate**

Run: `ruff check src tests validation && pytest -m "not requires_binary" -q && mypy src`

Expected: `All checks passed!`, `621 passed, 5 deselected`, `Success: no issues found in 79 source files`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/aligners/sibeliaz.py src/tessera/converters/maf_to_fasta.py src/tessera/converters/mafft_merge.py src/tessera/converters/xmfa_to_fasta.py src/tessera/core/io.py tests/integration/test_whitespace_genomes.py tests/unit/test_input_safety.py tests/unit/test_io.py
git commit -m "$(cat <<'EOF'
Clean whitespace from staged genomes and report malformed input by file

A FASTA with trailing spaces, tabs or CRLF line endings gave a ragged alignment.
Such a genome is now staged as a cleaned copy and the reader ignores whitespace in
sequence lines, so the aligner and Tessera agree on its length. A nameless header,
a truncated MAF block, an XMFA without the reference and an empty MAFFT result are
reported as errors that name the file.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

---
### Task 10: Docs, changelog and the full gate

**Spec item(s):** the documentation corrections named under A1, A3-A6, A8-A10.

**Files:**
- Modify: `docs/aligners.md`, `docs/reference-panels.md`, `CHANGELOG.md`

**Interfaces:**
- Consumes: the behaviour of Tasks 2-9.
- Produces: nothing code depends on.

The `--curate` help texts were already updated with the code in Tasks 5 and 6. What remains is the prose: `docs/aligners.md` says a fragmented query is "handled cleanly" by MAFFT (true only after Task 2) and that MAFFT "keeps" insertions (false: it runs with `--keeplength`); `docs/reference-panels.md` says curation runs "each round" and describes the reassort filter as "too little of their length".

- [ ] **Step 1: Apply the documentation and changelog changes**

`CHANGELOG.md`:

````diff
diff --git a/CHANGELOG.md b/CHANGELOG.md
index 39cbf1a..4f92d30 100644
--- a/CHANGELOG.md
+++ b/CHANGELOG.md
@@ -6,6 +6,53 @@ All notable changes to Tessera are recorded here. The format follows
 
 ## [Unreleased]
 
+### Fixed
+
+- **The MAFFT backend ignored strand.** A genome, or one contig of a draft assembly, on the
+  opposite strand to the backbone was aligned as given and came out at chance-level identity
+  (about 0.40 against 0.97 for the same genome in forward orientation), which the scan then
+  read as a divergent region. MAFFT now runs with `--adjustdirection`.
+- **`reassort` never applied its alignment-fraction filter.** The threshold was written as a
+  fraction (0.5) and compared with skani's percentage, so a tip aligning over a fifth of a
+  segment could outrank a full-length match and become the segment's nearest strain. The
+  threshold is now 50 %.
+- **`fill-references` did not align the last round's downloads.** On a `--max-rounds` exit the
+  final downloads were counted in `fill_summary.tsv`, the report and `lineages.tsv` but were
+  absent from `panel.msa.fasta`, so detection ran without them. One more alignment
+  (`final.msa.fasta`) is now built from them, without a further search.
+- **`fill-references --curate` did nothing unless a round both found a gap and downloaded a
+  reference.** A supplied collection that contains a whole-genome sibling of the query has
+  no coverage gap, so the loop converged before curation ran. Curation now runs before the
+  first alignment for a supplied collection and before each later build for that round's
+  downloads. A panel seeded from scratch is unchanged: seeding has its own sibling filter.
+- **`fill-references --curate --reference X` could delete X** and fail the next round with
+  "Reference 'X' not found". The given reference is now the curation backbone.
+- **`find-references --download <collection> --curate` deleted genomes that were already in
+  the collection.** Curation now removes only what the run downloaded; genomes already
+  present are reported as `sibling-kept` / `redundant-kept` in `panel_lineages.tsv`. A
+  downloaded sibling can no longer be picked as the curation backbone either.
+- **`reassort --scan-segments` could write outside the output directory.** A segment record
+  named `..` resolved the scan directory to the parent of the output and replaced its
+  `collection/`. Segment names are sanitised by one shared helper, and two names that
+  sanitise alike no longer share a scan directory.
+- **MAF-based backends (`sibeliaz`, `cactus`) could drop rows and columns.** A genome placed in
+  no alignment block had no row at all; it is now an all-gap row, named in a warning. With
+  `sibeliaz`, a multi-contig backbone was laid out in sorted-name order and lost any contig
+  that no block covered, shifting later coordinates; it is now laid out as its FASTA file is,
+  with uncovered contigs kept as gap columns.
+- **`progressivemauve` named rows after the resolved input path.** Symlinked collections,
+  unrecognised extensions and whitespace in a path gave wrong or duplicate row names, and a
+  genome whose link target shared the backbone's name was dropped. Rows are now named from
+  the staged genome label. (Checked against a simulated XMFA; the tool itself was not
+  available.)
+- **Whitespace in sequence lines was read as sequence.** A reference saved with trailing
+  spaces, tabs or Windows line endings produced a ragged alignment. Such a genome is now
+  aligned from a cleaned copy (a warning names it), and the FASTA reader ignores whitespace
+  in sequence lines.
+- A header line `> ` (no name), a truncated MAF block, an XMFA that does not list the
+  reference and an empty MAFFT result are reported as input/output errors that name the
+  file, instead of "Unexpected error".
+
 ## [1.2.0] - 2026-10-01
 
 ### Fixed
````

`docs/aligners.md`:

````diff
diff --git a/docs/aligners.md b/docs/aligners.md
index 55e9805..864d095 100644
--- a/docs/aligners.md
+++ b/docs/aligners.md
@@ -8,7 +8,7 @@ with `--aligner` and tune with repeatable `--aligner-arg key=value`.
 | Backend | Best for | Notes |
 |---|---|---|
 | `sibeliaz` (default) | Moderately divergent genomes, including rearrangements | Installs cleanly via conda; `kmer`, `abundance`, `bubble`, `filtermemory` |
-| `mafft` | Similar, largely collinear genomes | True base-level alignment, the canonical input for the window method; adds a fragmented query with `--addfragments`. `maxiterate`, `retree`, `op`, `ep`, `sixmerpair` |
+| `mafft` | Similar, largely collinear genomes | True base-level alignment, the canonical input for the window method; adds a fragmented query with `--addfragments` and reorients reverse-strand genomes or contigs (`--adjustdirection`). `maxiterate`, `retree`, `op`, `ep`, `sixmerpair` |
 | `minimap2` | Speed and assembly/contig queries | Fast assembly-to-reference projection; `preset` (default `asm20`, e.g. `asm10` for closer genomes) |
 | `progressivemauve` | Genomes with large rearrangements/inversions | Tolerant but slow, heavy, and not available as a conda build on all platforms; `seed_weight`, `single` |
 | `cactus` | Same-species pangenomes | Resource heavy (Toil/containers) |
@@ -17,8 +17,19 @@ with `--aligner` and tune with repeatable `--aligner-arg key=value`.
 data, reproduces `progressivemauve`'s recombination coordinates. For very similar,
 collinear genomes `mafft` gives the most faithful base-level signal and `minimap2`
 the fastest run (and the best fit for a fragmented query); `progressivemauve` remains
-an option for genomes with large rearrangements. Reference-anchored backends drop
-material inserted relative to the backbone; `mafft` keeps it as a true alignment.
+an option for genomes with large rearrangements. Every backend reports the alignment in
+backbone coordinates, so material inserted relative to the backbone is dropped -- `mafft`
+included, which runs with `--keeplength`.
+
+A multi-contig backbone is laid out as the concatenation of its contigs in the order of
+its FASTA file (`sibeliaz`, `mafft`, `minimap2`, `progressivemauve`); a contig that nothing
+aligns to keeps its columns, as gaps. The `cactus` backend lays contigs out in sorted-name
+order. A genome the aligner cannot place against the backbone at all is kept as an all-gap
+row and named in a warning, so the alignment always has one row per input genome.
+
+Genomes are staged before alignment. A FASTA with whitespace inside its sequence lines
+(trailing spaces, tabs, Windows line endings) is aligned from a cleaned copy and reported
+in a warning; the input file itself is not modified.
 
 Examples:
 
````

`docs/reference-panels.md`:

````diff
diff --git a/docs/reference-panels.md b/docs/reference-panels.md
index d1455fe..214ad96 100644
--- a/docs/reference-panels.md
+++ b/docs/reference-panels.md
@@ -105,7 +105,10 @@ It stops when the gaps close, when no new reference can be found, or when a roun
 longer improves the worst gap (so a genuinely hypervariable region is reported, not
 chased forever). Each round is recorded in `filled/fill_summary.tsv`, and the query's
 own record is auto-excluded from its FASTA header. This needs an aligner and Entrez
-Direct, and rebuilds the alignment every round.
+Direct, and rebuilds the alignment every round. If the last permitted round
+(`--max-rounds`) still downloaded references, one more alignment (`final.msa.fasta`) is
+built from them without a further search, so `panel.msa.fasta` holds every reference the
+run reports.
 
 Omit `--collection` to **start fresh** with no suggested references: the first round
 seeds the collection from an NCBI search, then the loop fills the remaining gaps as
@@ -277,9 +280,16 @@ SARS-CoV-2 sublineages (<1 % apart) alike. The curated `curated/collection/` and
 `panel_lineages.tsv` (each reference's role and ANI/coverage) are written; rebuild with
 `tessera msa` then `tessera recomb`.
 
-The same curation runs inside the fill loop with `fill-references --curate`, which
-keeps the growing panel diverse and sibling-free each round and adds a "Reference
-panel" section to the report. Both need skani (and skDER for dereplication):
+The same curation runs inside the fill loop with `fill-references --curate`, and adds a
+"Reference panel" section to the report. It runs before an alignment is built whenever the
+working collection holds genomes that have not been through it: the collection you
+supplied (before round 1), and each round's downloads (before the next build). A panel
+seeded from scratch is not curated before round 1, because seeding applies its own sibling
+filter. With `--reference`, that genome is the curation backbone and is never removed.
+`find-references --download ... --curate` curates only what that run downloaded: genomes
+already in the download directory are compared against but kept, and appear in
+`panel_lineages.tsv` as `sibling-kept` or `redundant-kept` if curation would otherwise
+have dropped them. Both need skani (and skDER for dereplication):
 `conda install -c bioconda skani skder`.
 
 ## Make a collection lineage-ready (`type-lineages`)
@@ -333,7 +343,7 @@ group); a `parent_group` column in `out/reassortment.tsv` records each segment's
 segments are linked only transitively (segment A shares a strain with B, B with C, but A and C share
 none) has no single spanning strain, so its `parent_strains` is empty and the text mosaic shows `?`
 for it; the verdict is still `clonal` because no pair disagrees. Segments
-below `--ani-floor` to every tip, or aligning over too little of their length, or with no resolvable
+below `--ani-floor` to every tip, or aligning over less than 50 % of their length, or with no resolvable
 dataset, are reported `unassigned` and excluded from the call. Because there is no influenza taxon
 alias, flu auto-typing needs the `nextclade` CLI (for `nextclade sort`) or explicit `--dataset
 SEGMENT=path` overrides; without either, each flu segment resolves to nothing and is left
````

- [ ] **Step 2: Run the full gate**

Run: `ruff check src tests validation && pytest -m "not requires_binary" -q && mypy src`
Expected: `All checks passed!`, `621 passed, 5 deselected`, `Success: no issues found in 79 source files`.

Run: `PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin" pytest -m requires_binary -q`
Expected: `5 passed, 621 deselected` (the existing MAFFT end-to-end test, two strand tests, two whitespace tests). A skip here means a binary is missing: report it as skipped, not passed.

Run: `COLUMNS=80 pytest tests/unit/test_cli_validation.py tests/unit/test_cli_commands.py -q`
Expected: `48 passed`. (Help text changed in Tasks 5 and 6; Rich wraps help to the terminal width, so check at CI's width.)

- [ ] **Step 3: Check the two shipped examples still give the same regions**

Run:

```bash
tessera recomb --msa example_data/divergent.msa.fasta --query query --output /tmp/a_div --window-size 300 --window-step 30 >/dev/null 2>&1
tessera recomb --msa example_data/cryptic_insert.msa.fasta --query query --output /tmp/a_cry >/dev/null 2>&1
cut -f1-6 /tmp/a_div/recombination_regions.tsv /tmp/a_cry/recombination_regions.tsv
```

Expected: one `parent_B` / `parent_A` region in each -- query 960-2010 in the divergent example and 4649-5349 in the cryptic one -- as on `main`. Nothing in this plan should change a region call; a difference means a task leaked into the scan.

- [ ] **Step 4: Optional checks that need tools or network not available when this plan was written**

Neither blocks the PR; record in the PR description which were run and what they showed.

- `python validation/run_reassort_benchmark.py` (needs network, `nextclade`, `skani`): Task 3 makes the alignment-fraction filter active for the first time, so the reassort precision/recall may move. Report the numbers before and after.
- One real `tessera msa --aligner progressivemauve` run on a collection of symlinked genomes, if the tool is available: Task 8 was verified against a simulated XMFA only. Check that the MSA has one row per genome, named by the collection's file names.

- [ ] **Step 5: Commit**

```bash
git add docs/aligners.md docs/reference-panels.md CHANGELOG.md
git commit -m "$(cat <<'EOF'
Document the data-safety and alignment fixes

Correct the statements that MAFFT keeps insertions and that curation runs each
round; describe backbone contig layout, unplaced genomes and cleaned staging; add
the changelog entries.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 6: Push and open the pull request**

```bash
git push -u origin fix-audit-data-safety
gh pr create --title "Fix data-loss and alignment defects from the post-1.2.0 audit (plan A)" --body "$(cat <<'EOF'
Implements plan A of docs/superpowers/specs/2026-10-01-post-1.2.0-audit-design.md (items A1-A10).

- MAFFT backend reorients reverse-strand genomes and contigs
- reassort alignment-fraction filter compared in percent
- fill-references aligns the last round's downloads; --curate runs before builds and honours --reference
- find-references --curate removes only what it downloaded
- reassort --scan-segments cannot write outside its output directory
- MAF backends keep uncovered contigs and unplaced genomes
- progressiveMauve rows named from the staged label (simulated XMFA only)
- whitespace in sequence lines cleaned at staging; malformed input reported by file

No region caller, default or output schema changes. Fast suite: 621 passed. requires_binary (mafft, minimap2): 5 passed.

Optional checks (reassort benchmark; a real progressiveMauve run): not run.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

Before running `gh pr create`, edit the two result lines in the body to what was actually observed: the pass counts from Step 2, and for the optional checks either "not run" or the numbers from Step 4.

Do not merge. Wait for CI, then hand the merge to the maintainer.
