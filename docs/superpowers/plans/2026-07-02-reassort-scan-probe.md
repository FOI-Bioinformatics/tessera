# reassort --scan-segments validation probe Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A durable, opt-in validation probe that builds a real-tip insert-shaped HA recombinant plus a clonal NA control, runs the shipped `reassort --scan-segments`, and reports whether the scan localizes the intragenic recombination.

**Architecture:** A new `validation/run_reassort_scan.py` with a pure `region_overlaps_span` scorer (CI-unit-tested) and a `main()` that imports tip-loading/splicing helpers from `run_hybrids`, drives `tessera.reassort.assign_segments(scan_segments=True)`, and prints an XPASS/XFAIL localization line plus reported attribution and specificity. Opt-in (needs the aligner env); not part of CI.

**Tech Stack:** Python 3.11+, skani + Nextclade + an aligner (external), pytest.

## Global Constraints

- No change to `src/` behavior; the probe drives the shipped `assign_segments`. The only CI-tested code is the pure scorer, unit-tested from `tests/unit`.
- XPASS is gated on localization only (a called region overlaps the true insert span); donor attribution and clonal-control specificity are reported, not gating.
- Opt-in: the end-to-end run needs skani + an aligner + network; it is not added to CI. Faithful reporting -- an HA miss or a flagged clonal control is stated, never retuned.
- Modest scientific language; no new runtime dependency.
- Ruff clean (`ruff check src tests validation`, line length 100, E/F/I/UP/B); CI coverage floor 80% holds.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

---

## File Structure

- `validation/run_reassort_scan.py` (create) -- the probe: `region_overlaps_span` (pure) + `main()`.
- `tests/unit/test_reassort_scan_probe.py` (create) -- CI unit test for `region_overlaps_span`.
- `validation/attribution-results.md` (modify) -- record the measured e2e.
- `validation/README.md` (modify) -- list the new opt-in harness.

---

## Task 1: The pure localization scorer

**Files:**
- Create: `validation/run_reassort_scan.py`
- Create: `tests/unit/test_reassort_scan_probe.py`

**Interfaces:**
- Produces: `region_overlaps_span(regions: list[dict], q_start: int, q_end: int) -> bool` -- True iff any present (non-`donor_absent`) region's `[query_start, query_end]` overlaps `[q_start, q_end]`.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_reassort_scan_probe.py`:

```python
"""Unit test for the reassort-scan probe's localization scorer (no binaries)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_reassort_scan.py"
_SPEC = importlib.util.spec_from_file_location("run_reassort_scan", _PATH)
rrs = importlib.util.module_from_spec(_SPEC)
sys.modules["run_reassort_scan"] = rrs
_SPEC.loader.exec_module(rrs)


def _rows(*triples):
    return [{"query_start": str(a), "query_end": str(b), "donor_absent": d}
            for a, b, d in triples]


def test_overlapping_region_is_localized():
    assert rrs.region_overlaps_span(_rows((600, 1200, "no")), 601, 1122) is True


def test_disjoint_region_is_not_localized():
    assert rrs.region_overlaps_span(_rows((0, 100, "no")), 601, 1122) is False


def test_touching_edge_counts_as_overlap():
    # region starts exactly at the span end -> touching interval counts.
    assert rrs.region_overlaps_span(_rows((1122, 1300, "no")), 601, 1122) is True


def test_donor_absent_region_is_excluded():
    assert rrs.region_overlaps_span(_rows((600, 1200, "yes")), 601, 1122) is False


def test_no_regions_is_not_localized():
    assert rrs.region_overlaps_span([], 601, 1122) is False
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/unit/test_reassort_scan_probe.py -q`
Expected: FAIL with `FileNotFoundError` / `spec_from_file_location` returning None (module does not exist yet).

- [ ] **Step 3: Create `run_reassort_scan.py` with the scorer**

Create `validation/run_reassort_scan.py`:

```python
"""Opt-in probe for ``reassort --scan-segments``: does the intragenic scan localize a real
recombinant?

Builds an HA query from two divergent H3N2 HA clade tips (a clade-B insert in a clade-A
backbone, real genomes with an exact known breakpoint) plus a clonal NA control, runs the
shipped ``assign_segments(scan_segments=True)``, and reports whether a called region overlaps
the true insert span (localization, the XPASS gate), whether the donor clade is named
(attribution, reported), and whether the clonal control is flagged (specificity, reported).

Needs skani, an aligner, and network access to the Nextclade datasets. Not part of CI. Run:

    export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
    python validation/run_reassort_scan.py
"""

from __future__ import annotations

import logging
import sys
import tempfile
from pathlib import Path

HA_DATASET = "nextstrain/flu/h3n2/ha/EPI1857216"
NA_DATASET = "nextstrain/flu/h3n2/na/EPI1857215"


def region_overlaps_span(regions: list[dict], q_start: int, q_end: int) -> bool:
    """True iff any present (non-``donor_absent``) region overlaps ``[q_start, q_end]``."""
    for r in regions:
        if r.get("donor_absent") == "yes":
            continue
        try:
            r_start, r_end = int(r["query_start"]), int(r["query_end"])
        except (KeyError, ValueError):
            continue
        if r_start <= q_end and r_end >= q_start:
            return True
    return False
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/unit/test_reassort_scan_probe.py -q`
Expected: PASS (5 passed).

- [ ] **Step 5: Lint**

Run: `ruff check validation/run_reassort_scan.py tests/unit/test_reassort_scan_probe.py`
Expected: `All checks passed!`

- [ ] **Step 6: Commit**

```bash
git add validation/run_reassort_scan.py tests/unit/test_reassort_scan_probe.py
git commit -m "Add the reassort-scan probe localization scorer

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Build the recombinant, run the scan, report

**Files:**
- Modify: `validation/run_reassort_scan.py` (add `main()` and helpers)
- Modify: `validation/attribution-results.md`
- Modify: `validation/README.md`

**Interfaces:**
- Consumes from `run_hybrids` (imported inside `main`): `_load_species(path, clade_key, logger) -> (genomes, reference, tips)`; `pick_parents(tips, reference, pinned, logger) -> (clade_a, clade_b, src_a, src_b)`; `make_hybrid(reference, src_a_muts, src_b_muts) -> (query, q_start, q_end)`; `_largest_clade(tips) -> (clade, central_acc)`; `reconstruct_gapped(reference, muts) -> str`; `parse_regions(path) -> list[dict]`; the exception `CaseSkipped`.
- Consumes from tessera: `assign_segments` (`tessera.reassort`), `write_fasta_record` (`tessera.core.io`).
- `tips[acc]` is `(clade, mutation_list)`, so a source accession's mutations are `tips[acc][1]`.

- [ ] **Step 1: Add the build/run/score/report logic**

Append to `validation/run_reassort_scan.py`:

```python
def _clonal_na(tips_na, reference_na, reconstruct_gapped, largest_clade):
    """A real single-clade NA genome (the largest clade's central tip) -> (clade, sequence)."""
    clade, src = largest_clade(tips_na)
    seq = reconstruct_gapped(reference_na, tips_na[src][1]).replace("-", "").upper()
    return clade, seq


def _run_probe(out_dir: Path, logger: logging.Logger) -> int:
    """Build the query, run the shipped scan, print the report. Always returns 0 (an opt-in
    probe never fails the caller); a SKIP is printed when no divergent HA pair is available."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import run_hybrids as rh
    from tessera.core.io import write_fasta_record
    from tessera.reassort import assign_segments

    _g, ref_ha, tips_ha = rh._load_species(HA_DATASET, None, logger)
    try:
        clade_a, clade_b, src_a, src_b = rh.pick_parents(tips_ha, ref_ha, [], logger)
    except rh.CaseSkipped as exc:
        print(f"[SKIP  ] HA localization  {exc}")
        return 0
    ha_query, q_start, q_end = rh.make_hybrid(ref_ha, tips_ha[src_a][1], tips_ha[src_b][1])

    _g2, ref_na, tips_na = rh._load_species(NA_DATASET, None, logger)
    na_clade, na_query = _clonal_na(tips_na, ref_na, rh.reconstruct_gapped, rh._largest_clade)

    query = out_dir / "query.fasta"
    with open(query, "w") as fo:
        write_fasta_record(fo, "HA", ha_query)
        write_fasta_record(fo, "NA", na_query)
    logger.info("HA insert: clade %s in clade %s backbone (true span %d-%d); NA clonal %s",
                clade_b, clade_a, q_start, q_end, na_clade)

    assign_segments(query, dataset_overrides={"HA": HA_DATASET, "NA": NA_DATASET},
                    scan_segments=True, output=out_dir, logger=logger)

    ha_regions = rh.parse_regions(out_dir / "HA" / "recombination_regions.tsv")
    na_regions = rh.parse_regions(out_dir / "NA" / "recombination_regions.tsv")
    localized = region_overlaps_span(ha_regions, q_start, q_end)
    attributed = any(clade_b in r.get("minor_parent", "") or clade_b in r.get("major_parent", "")
                     for r in ha_regions)
    na_present = [r for r in na_regions if r.get("donor_absent") != "yes"]

    print("\nTessera reassort --scan-segments probe (opt-in; needs the aligner env)")
    print("=" * 70)
    verdict = "XPASS " if localized else "XFAIL "
    print(f"[{verdict}] HA localization  insert clade {clade_b} (true span {q_start}-{q_end}); "
          f"{'a region overlaps it' if localized else 'no region overlaps it'}")
    print(f"[report] HA attribution   donor clade {clade_b} named in a region: "
          f"{'yes' if attributed else 'no'}")
    print(f"[report] NA specificity   clonal control ({na_clade}): "
          f"{len(na_present)} region(s) called")
    print("\n(opt-in probe; localization is the gate, attribution/specificity are reported)")
    return 0


def main(argv: list[str]) -> int:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    logger = logging.getLogger("reassort-scan-probe")
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(argv[0]) if argv else Path(tmp)
        out.mkdir(parents=True, exist_ok=True)
        return _run_probe(out, logger)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
```

- [ ] **Step 2: Lint and re-run the unit test**

Run: `ruff check validation/run_reassort_scan.py`
Expected: `All checks passed!`
Run: `python -m pytest tests/unit/test_reassort_scan_probe.py -q`
Expected: PASS (5 passed) -- the scorer is unchanged; the added `main`/`_run_probe` do not run under import.

- [ ] **Step 3: Opt-in end-to-end run and record**

Run (aligner env; not CI):
```bash
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
export NCBI_EMAIL="andreas.sjodin@gmail.com"
python validation/run_reassort_scan.py
```
Expected: an `HA localization` line (XPASS when a region overlaps the true insert span) plus the
reported attribution and NA-specificity lines. In `validation/attribution-results.md`, in the
`reassort --scan-segments` section, replace the note that a cleaner probe is "follow-up work" with
the measured result of this probe (the HA localization verdict, whether the donor clade was named,
and the NA clonal-control region count). State the measured result faithfully; if HA is not
localized, record what was measured and stop for review rather than retuning.

- [ ] **Step 4: Document the harness in `validation/README.md`**

Add a short entry describing `run_reassort_scan.py`: an opt-in probe that builds a real-tip
insert-shaped HA recombinant plus a clonal NA control and checks that `reassort --scan-segments`
localizes the intragenic recombination (localization-gated; attribution and specificity reported).
Note it needs the aligner env and is not part of CI. Match the surrounding entry style.

- [ ] **Step 5: Commit**

```bash
git add validation/run_reassort_scan.py validation/attribution-results.md validation/README.md
git commit -m "reassort-scan probe: build the recombinant, run the scan, record the e2e

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review

**Spec coverage:**
- Real-tip insert-shaped HA recombinant -> Task 2 `make_hybrid` on two `pick_parents` clades.
- Clonal NA control -> Task 2 `_clonal_na` (largest-clade central tip).
- Drives the shipped scan -> Task 2 `assign_segments(scan_segments=True, output=...)`.
- Localization-only XPASS gate -> Task 1 `region_overlaps_span`; Task 2 `verdict`.
- Attribution + NA specificity reported, not gating -> Task 2 report lines.
- CI unit test for the pure scorer -> Task 1.
- Opt-in e2e + faithful record -> Task 2 Steps 3-4.
- Dedicated file importing run_hybrids helpers -> Task 2 imports inside `_run_probe`.

**Placeholder scan:** every code step shows complete code; the only prose steps are the opt-in run
and the two doc edits (Task 2 Steps 3-4), which are a manual run and faithful documentation.

**Type consistency:** `region_overlaps_span(list[dict], int, int) -> bool` is identical in Task 1
(definition + test) and Task 2 (call site). `pick_parents` returns `(clade_a, clade_b, src_a, src_b)`
and `make_hybrid(reference, src_a_muts, src_b_muts) -> (query, q_start, q_end)`; the source
accessions index `tips`, and `tips[acc][1]` is the mutation list -- consistent with the Interfaces
block and the real `run_hybrids` signatures. `_load_species` returns `(genomes, reference, tips)`;
the probe uses `reference` and `tips` and ignores `genomes` (`_g`).

**Note for the implementer:** import `run_hybrids` inside `_run_probe` (not at module top) so the
Task 1 unit test can import `region_overlaps_span` without pulling in the heavy harness; the test
loads the module by file path via `importlib`.
