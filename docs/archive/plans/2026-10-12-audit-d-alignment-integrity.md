# Audit Plan D: Alignment Integrity and Adapter Resources Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the alignment that reaches the scan faithful to its inputs (a backbone row that is the backbone, contigs that cannot be confused, truncated or mismatched tool output that is reported rather than absorbed) and keep the adapters within their thread budget and failure-handling contract.

**Architecture:** Every task changes what alignment reaches the scan, or how an adapter fails, never how the scan calls regions. Checks are placed where the input is first seen: duplicate sequence IDs at staging (one place, every backend), backbone and contig-size checks in the adapter or converter that owns the layout, and process hygiene in the `run_tool` choke point. Two helpers are defined once in `aligners/base.py` and shared by the pairwise adapters: `split_thread_budget` (the MAFFT rule from commit `44e6b3f`, now also used by minimap2) and `warn_unplaced`.

**Tech Stack:** Python 3.11+, numpy, pytest, ruff, mypy. No new dependency.

**Spec:** `docs/superpowers/specs/2026-10-12-audit-swot-gap.md` (section 1.3, items C1-C7 and the "Suspected and minor" paragraph; section 3, G3 and G7, and the "Core side" test gaps).

## Global Constraints

- No new runtime dependency.
- Modest scientific language in code, docs and messages; reported numbers must be faithful (state what passes, what fails, what was skipped).
- "Could not test" must never be reported as "tested, found nothing".
- A new behaviour needs a test that fails without the change; if a test is a guard that passes before and after, the plan says so.
- Branch: create `fix-audit4-alignment` from `fix-audit-caller-behaviour` (or from `main` once that branch is merged). Commit messages explain the problem, the evidence and the trade-off, and end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Append the aligner env to PATH, never prepend (prepending swaps the Python interpreter): `export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"`.
- The dev Mac has 18 GB RAM: run one harness at a time.
- Each task ends with `ruff check src tests validation`, `mypy`, and `pytest -m "not requires_binary" -q` green before its commit.
- No Unicode in code or docs beyond what files already contain; ASCII punctuation.
- **Aligner binaries.** The spec says no aligner binary is installed. On the machine this plan was written on that is only partly true: the `recomfi-aln` env holds `mafft`, `minimap2` and `sibeliaz` (1.2.7); `progressiveMauve` and `cactus-pangenome` are not installed anywhere reachable. CI (`.github/workflows/ci.yml`) installs `mafft` only. Therefore: `requires_binary` tests run in CI with mafft alone (the minimap2, sibeliaz and progressiveMauve cases skip there); the minimap2 and sibeliaz cases run locally with the env appended to PATH; the progressiveMauve and cactus paths rest on fakes of the tool in every environment, and the PR says so. A real mafft run (the existing `requires_binary` end-to-end test) happens in CI; run it locally as well where the env exists.
- **No caller, default or region-calling change.** No task touches a caller or a default, so the `run_hybrids.py` gate of `CONTRIBUTING.md` is not required. The one change under `recomb/` (Task 6) splits a warning in two; `validation/run_specificity.py` is run in Task 12 as a check that its output is unchanged.
- **Shared files.** `src/tessera/recomb/run.py`, `src/tessera/core/io.py`, `docs/aligners.md` and `CHANGELOG.md` may also be edited by other fourth-audit plans. Rebase onto whichever lands first; every edit below is anchored on quoted text, not on line numbers alone.

## Review Focus

Input classes the spec implies but no example in it exercises, most likely to bite a user first. Each has a test in the task that owns the code.

1. **The default backend on macOS with an ordinary output path.** SibeliaZ's wrapper uses `paste --delimiters` and an `xargs -I` command longer than BSD's 255-byte limit, so on macOS every run wrote an empty MAF and the adapter blamed memory. *Task 3, `test_sibeliaz_bsd_wrapper_rewrites_paste_and_xargs` and `test_backbone_row_is_the_backbone_with_the_real_tool[sibeliaz]`.*
2. **An NCBI download, gzipped, whose records share an accession and differ only in their description.** The ID is the first token; the check must read gzip and split on whitespace. *Task 4, `test_a_repeated_id_is_found_in_a_gzipped_genome_by_its_first_token`.*
3. **A `--threads` budget that does not divide by the number of jobs** (five threads, two genomes). The total started must stay within the budget. *Task 5, `test_minimap2_never_starts_more_threads_than_the_budget`.*
4. **An all-N record beside an all-gap record.** Splitting the "unplaced" warning out of the alphabet warning must not swallow the alphabet warning for a placed record with no canonical base. *Task 6, `test_an_all_gap_row_is_unaligned_not_an_alphabet_problem` and `test_scan_names_an_all_gap_row_as_unaligned`.*
5. **A dangling symlink at the working-collection path.** `Path.exists()` is false for it, so an "is it a directory" check that starts from `exists()` misses it and `copytree` then fails with a bare `FileExistsError`. *Task 9, `test_a_non_directory_at_the_working_collection_path_is_a_user_error[dangling]`.*

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `src/tessera/aligners/base.py` | modify | `split_thread_budget`, `warn_unplaced` (shared by the pairwise adapters) |
| `src/tessera/aligners/progressivemauve.py` | modify | C1 backbone row from the reference; width check; unplaced warning; drop the `flank` argument |
| `src/tessera/aligners/sibeliaz.py` | modify | two more BSD wrapper patches; empty-MAF message (outside the spec, Task 3) |
| `src/tessera/aligners/minimap2.py` | modify | C3 shared thread budget; C5 unplaced warning |
| `src/tessera/aligners/mafft.py` | modify | use `split_thread_budget`; C5 unplaced warning |
| `src/tessera/aligners/cactus.py` | modify | refuse labels that collide once '.' becomes '_' |
| `src/tessera/core/io.py` | modify | C2 repeated sequence ID at staging; C7 non-directory working-collection path |
| `src/tessera/core/process.py`, `src/tessera/core/binaries.py` | modify | child processes get `stdin=DEVNULL` |
| `src/tessera/core/logging.py` | modify | close replaced handlers |
| `src/tessera/converters/maf_to_fasta.py` | modify | C6 short `s` row; backbone source size against the declared contig length |
| `src/tessera/converters/xmfa_to_fasta.py` | modify | remove the unreachable `flank` path |
| `src/tessera/discover/nextclade.py` | modify | C4 install the pool with `atomic_cache_dir`, with a manifest |
| `src/tessera/recomb/similarity.py`, `src/tessera/recomb/run.py` | modify | C5: name an all-gap row as unaligned, not as an alphabet problem |
| `tests/unit/test_backbone_invariant.py` | create | C1 invariant for all five backends against fakes; cactus sample names |
| `tests/integration/test_backbone_invariant_binaries.py` | create | the same invariant with the installed binaries (`requires_binary`) |
| `tests/unit/test_logging.py` | create | handler close |
| `tests/unit/test_io.py`, `test_aligner_backends.py`, `test_nextclade.py`, `test_similarity.py`, `test_run.py`, `test_converters.py`, `test_input_safety.py`, `test_process.py`, `test_binaries.py` | modify | new tests per task |
| `tests/unit/test_aligner_params.py`, `tests/unit/test_input_safety.py` | modify | two existing tests change (Task 6 stub return value; Task 11 call signature) |
| `docs/aligners.md`, `CHANGELOG.md` | modify | describe each change in the task that makes it |

---

### Task 1: Branch and baseline

Nothing is committed in this task.

**Files:** none.

**Interfaces:** Produces `validation/data/audit-d/baseline/spec.txt` (git-ignored) for Task 12 to diff against.

- [ ] **Step 1: Branch**

```bash
git switch fix-audit-caller-behaviour   # or main, once that branch is merged
git switch -c fix-audit4-alignment
pip install -e ".[dev]"
```

- [ ] **Step 2: Confirm the starting point**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected at `4b35e39`: `All checks passed!`, `Success: no issues found in 79 source files`, `757 passed, 5 deselected`. If the count differs, record the new baseline; each later task states how many tests it adds.

- [ ] **Step 3: Binary-gated baseline (where the env exists)**

```bash
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
pytest -m requires_binary -q -rs
```

Expected where the env exists: `5 passed, 757 deselected`. Without it, every case skips; write "requires_binary not run locally: no aligner env" in the PR and rely on CI for mafft.

- [ ] **Step 4: Specificity baseline**

```bash
mkdir -p validation/data/audit-d/baseline
python validation/run_specificity.py --reps 3 > validation/data/audit-d/baseline/spec.txt 2>&1
tail -15 validation/data/audit-d/baseline/spec.txt
```

No aligner, network or data is needed. Keep the summary lines for the PR.

---

### Task 2: The progressiveMauve backbone row is the backbone (audit items: C1)

`_concatenate` wrote the backbone row from the first query's projection, which is gapped wherever that query did not align. With two queries covering opposite halves of the backbone, half the backbone is lost, and which half depends on list order. The row stays full width, so nothing downstream notices. The fix writes the reference's own sequence, which the adapter already reads to pin the width. The invariant test covers all five backends against fakes of their tools; for mafft, minimap2, sibeliaz and cactus it is a guard that passes before and after.

**Files:**
- Modify: `src/tessera/aligners/progressivemauve.py` (lines 53-56 `ref_length`, line 88 the `_concatenate` call, lines 92-116 `_concatenate`, at `4b35e39`)
- Create: `tests/unit/test_backbone_invariant.py`
- Modify: `CHANGELOG.md` (new group under `### Fixed`)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `_concatenate(per_query: list[tuple[str, Path]], reference_label: str, reference_seq: str, out_path: Path) -> None`. Raises `OutputError` when a query row's width differs from `len(reference_seq)`.
- Produces (test module, reused by Tasks 3 and 11): `tests/unit/test_backbone_invariant.py` with `C1`, `C2`, `BACKBONE`, `_genomes(tmp_path) -> list[Path]`, `_read_msa(path) -> dict[str, str]`, `_LOG`.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_backbone_invariant.py`:

```python
"""The backbone row of every backend's alignment is the backbone's own sequence.

The backbone is both the coordinate anchor and a donor candidate, so a backbone row
that lost columns is a wrong alignment that still looks well-formed: full width,
rectangular, no ragged-row guard fires. Each backend is driven by a fake of its tool
that aligns two queries to opposite halves of a two-contig backbone, so no single
query covers the whole backbone. The fakes cover every backbone column, so for the
MAF-based backends (sibeliaz, cactus) equality is the right check here too.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from tessera.aligners import cactus as cactus_mod
from tessera.aligners import mafft as mafft_mod
from tessera.aligners import minimap2 as mm_mod
from tessera.aligners import progressivemauve as pm_mod
from tessera.aligners import sibeliaz as sz_mod
from tessera.aligners.base import AlignParams

_LOG = logging.getLogger("tessera.test")

# Two backbone contigs of 10 bp. Query A aligns only to c1, query B only to c2.
C1 = "ACGTACGTAA"
C2 = "TTGCATGCCG"
BACKBONE = C1 + C2


def _read_msa(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            out[name] = ""
        elif name is not None:
            out[name] += line
    return out


def _genomes(tmp_path: Path) -> list[Path]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    ref = tmp_path / "ref.fasta"
    ref.write_text(f">c1\n{C1}\n>c2\n{C2}\n")
    qa = tmp_path / "qryA.fasta"
    qa.write_text(f">a1\n{C1}\n")
    qb = tmp_path / "qryB.fasta"
    qb.write_text(f">b1\n{C2}\n")
    return [ref, qa, qb]


def _run_progressivemauve(monkeypatch, genomes: list[Path], out: Path) -> Path:
    def fake(caps, cmd, **kw):
        cmd = [str(c) for c in cmd]
        xmfa = Path(cmd[cmd.index("--output") + 1])
        ref_path, query_path = cmd[-2], cmd[-1]
        # progressiveMauve reads a multi-record FASTA as one concatenated sequence.
        lo, hi = (1, 10) if "qryA" in query_path else (11, 20)
        bases = BACKBONE[lo - 1:hi]
        xmfa.write_text(
            f"#Sequence1File\t{ref_path}\n#Sequence2File\t{query_path}\n"
            f"> 1:{lo}-{hi} + {ref_path}\n{bases}\n> 2:1-10 + {query_path}\n{bases}\n=\n"
        )
        return ""

    monkeypatch.setattr(pm_mod, "run_tool", fake)
    return pm_mod.ProgressiveMauveAligner().align(
        genomes, genomes[0], out, AlignParams(threads=1), _LOG
    ).msa_fasta


def _run_mafft(monkeypatch, genomes: list[Path], out: Path) -> Path:
    def fake(caps, cmd, **kw):
        added = str(cmd[cmd.index("--addfragments") + 1])
        row = C1 + "-" * 10 if "qryA" in added else "-" * 10 + C2
        kw["stdout_path"].write_text(f">ref\n{BACKBONE}\n>contig\n{row}\n")
        return ""

    monkeypatch.setattr(mafft_mod, "run_tool", fake)
    return mafft_mod.MafftAligner().align(
        genomes, genomes[0], out, AlignParams(threads=1), _LOG
    ).msa_fasta


def _run_minimap2(monkeypatch, genomes: list[Path], out: Path) -> Path:
    def fake(caps, cmd, **kw):
        query = str(cmd[-1])
        rname, seq = ("c1", C1) if "qryA" in query else ("c2", C2)
        kw["stdout_path"].write_text(
            "@SQ\tSN:c1\tLN:10\n@SQ\tSN:c2\tLN:10\n"
            f"q\t0\t{rname}\t1\t60\t10M\t*\t0\t0\t{seq}\t*\n"
        )
        return ""

    monkeypatch.setattr(mm_mod, "run_tool", fake)
    return mm_mod.Minimap2Aligner().align(
        genomes, genomes[0], out, AlignParams(threads=1), _LOG
    ).msa_fasta


def _run_sibeliaz(monkeypatch, genomes: list[Path], out: Path) -> Path:
    def fake(caps, cmd, **kw):
        out_dir = Path(cmd[cmd.index("-o") + 1])
        (out_dir / "alignment.maf").write_text(
            f"a\ns c1 0 10 + 10 {C1}\ns a1 0 10 + 10 {C1}\n\n"
            f"a\ns c2 0 10 + 10 {C2}\ns b1 0 10 + 10 {C2}\n"
        )
        return ""

    monkeypatch.setattr(sz_mod, "run_tool", fake)
    monkeypatch.setattr(sz_mod, "_sibeliaz_invocation", lambda out_dir, logger: ["sibeliaz"])
    return sz_mod.SibeliazAligner().align(
        genomes, genomes[0], out, AlignParams(threads=1), _LOG
    ).msa_fasta


def _run_cactus(monkeypatch, genomes: list[Path], out: Path) -> Path:
    def fake(caps, cmd, **kw):
        results = Path(cmd[cmd.index("--outDir") + 1])
        results.mkdir(parents=True, exist_ok=True)
        (results / "pangenome.full.hal").write_text("")
        return ""

    def fake_hal_to_maf(hal, ref_name, maf, logger, caps=None):
        # hal2maf names rows <sample>.<contig>.
        Path(maf).write_text(
            f"a\ns ref.c1 0 10 + 10 {C1}\ns qryA.a1 0 10 + 10 {C1}\n\n"
            f"a\ns ref.c2 0 10 + 10 {C2}\ns qryB.b1 0 10 + 10 {C2}\n"
        )
        return Path(maf)

    monkeypatch.setattr(cactus_mod, "run_tool", fake)
    monkeypatch.setattr(cactus_mod, "hal_to_maf", fake_hal_to_maf)
    return cactus_mod.CactusAligner().align(
        genomes, genomes[0], out, AlignParams(threads=1), _LOG
    ).msa_fasta


_BACKENDS = {
    "progressivemauve": _run_progressivemauve,
    "mafft": _run_mafft,
    "minimap2": _run_minimap2,
    "sibeliaz": _run_sibeliaz,
    "cactus": _run_cactus,
}


@pytest.mark.parametrize("backend", sorted(_BACKENDS))
@pytest.mark.parametrize("order", ["AB", "BA"])
def test_backbone_row_is_the_backbone_sequence(
    monkeypatch, tmp_path: Path, backend: str, order: str
) -> None:
    ref, qa, qb = _genomes(tmp_path / "in")
    genomes = [ref, qa, qb] if order == "AB" else [ref, qb, qa]
    msa = _read_msa(_BACKENDS[backend](monkeypatch, genomes, tmp_path / "out"))
    assert msa["ref"] == BACKBONE
    assert msa["qryA"] == C1 + "-" * 10
    assert msa["qryB"] == "-" * 10 + C2


def test_progressivemauve_rejects_a_projection_wider_than_the_backbone(
    monkeypatch, tmp_path: Path
) -> None:
    # The backbone row is now the reference itself, so a projection that reaches past
    # its end would make the rows ragged. It is reported against the file instead.
    from tessera.core.errors import OutputError

    def fake(caps, cmd, **kw):
        cmd = [str(c) for c in cmd]
        xmfa = Path(cmd[cmd.index("--output") + 1])
        ref_path, query_path = cmd[-2], cmd[-1]
        bases = BACKBONE + "ACGT"
        xmfa.write_text(
            f"#Sequence1File\t{ref_path}\n#Sequence2File\t{query_path}\n"
            f"> 1:1-24 + {ref_path}\n{bases}\n> 2:1-24 + {query_path}\n{bases}\n=\n"
        )
        return ""

    monkeypatch.setattr(pm_mod, "run_tool", fake)
    genomes = _genomes(tmp_path / "in")
    with pytest.raises(OutputError, match=r"qryA.*24 column.*20 bp"):
        pm_mod.ProgressiveMauveAligner().align(
            genomes, genomes[0], tmp_path / "out", AlignParams(threads=1), _LOG
        )
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_backbone_invariant.py -q`
Expected: `3 failed, 8 passed`. Observed while writing this plan:
- `test_backbone_row_is_the_backbone_sequence[AB-progressivemauve]`: `AssertionError: assert 'ACGTACGTAA----------' == 'ACGTACGTAATTGCATGCCG'`
- `[BA-progressivemauve]`: `assert '----------TTGCATGCCG' == 'ACGTACGTAATTGCATGCCG'` (the lost half follows the query order)
- `test_progressivemauve_rejects_a_projection_wider_than_the_backbone`: `Failed: DID NOT RAISE <class 'tessera.core.errors.OutputError'>`

The eight mafft, minimap2, sibeliaz and cactus cases pass before and after: they are guards.

- [ ] **Step 3: Implement**

In `src/tessera/aligners/progressivemauve.py`, replace

```python
        # Pin every per-query projection to the full reference length so the
        # rows concatenate into a rectangular MSA regardless of how far each
        # query aligns.
        ref_length = sum(len(seq) for _, seq in read_fasta(reference))
```

with

```python
        # Pin every per-query projection to the full reference length so the
        # rows concatenate into a rectangular MSA regardless of how far each
        # query aligns. A multi-record reference is read by progressiveMauve as one
        # concatenated sequence, so the backbone row is that concatenation.
        ref_seq = "".join(seq for _, seq in read_fasta(reference))
        ref_length = len(ref_seq)
```

replace

```python
        _concatenate(per_query, reference.stem, msa)
```

with

```python
        _concatenate(per_query, reference.stem, ref_seq, msa)
```

and replace the whole of `_concatenate` with

```python
def _concatenate(
    per_query: list[tuple[str, Path]], reference_label: str, reference_seq: str,
    out_path: Path,
) -> None:
    """Write the backbone row once, then one row per query, named by staged label.

    The backbone row is the reference's own sequence, not a projection. Each
    per-query projection carries a reference row gapped wherever that query did not
    align, so taking it from any one projection would lose the columns that query
    missed -- and which columns were lost would depend on the order of the queries.

    Each per-query FASTA holds exactly two records, in a fixed order: the reference
    projection, then the query's. The names inside those files are the paths
    progressiveMauve echoed from its command line -- resolved symlink targets, which can
    carry an unrecognised extension, whitespace, or a basename shared with another
    genome -- so rows are identified by position and labelled from the staged file,
    never from those names.
    """
    with open(out_path, "w") as out:
        write_fasta_record(out, reference_label, reference_seq)
        for label, fa in per_query:
            records = read_fasta(fa)
            if len(records) != 2:
                raise OutputError(
                    f"Expected a reference row and one query row in {fa} (the "
                    f"projection of '{label}' onto '{reference_label}'), found "
                    f"{len(records)} record(s). progressiveMauve's output is not a "
                    "pairwise alignment."
                )
            row = records[1][1]
            if len(row) != len(reference_seq):
                raise OutputError(
                    f"The projection of '{label}' in {fa} is {len(row)} column(s) wide; "
                    f"the backbone '{reference_label}' is {len(reference_seq)} bp. "
                    "progressiveMauve reported coordinates beyond the backbone."
                )
            write_fasta_record(out, label, row)
```

In `CHANGELOG.md`, under `## [Unreleased]` / `### Fixed`, insert a new group directly above the line `Data safety and alignment:`:

```markdown
Alignment integrity, from the fourth audit
(`docs/superpowers/specs/2026-10-12-audit-swot-gap.md`, section 1.3). These change what
alignment reaches the scan, not how the scan calls regions:

- **The progressiveMauve backbone row could lose columns.** It was taken from the first
  query's projection, which is gapped wherever that query did not align; with queries
  covering different parts of the backbone, the rest of the backbone came out as gaps,
  and which part depended on the order of the queries. The backbone row is now the
  backbone's own sequence. Checked against a fake of the tool only: progressiveMauve is
  not installed in CI or on the development machine.

```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_backbone_invariant.py tests/unit/test_aligner_backends.py tests/unit/test_aligner_params.py -q`
Expected: all pass.

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: clean, clean, `768 passed, 5 deselected` (757 + 11).

- [ ] **Step 6: Commit**

```bash
git add src/tessera/aligners/progressivemauve.py tests/unit/test_backbone_invariant.py CHANGELOG.md
git commit -m "$(cat <<'EOF'
Write the progressiveMauve backbone row from the backbone itself

Each per-query projection carries a reference row gapped wherever that query
did not align, and the adapter wrote the backbone row from the first
projection only. With two queries covering opposite halves of a 20 bp
backbone, half the backbone came back as gaps, and which half depended on
the order of the queries. The row stayed full width, so the ragged-row
guard never fired; the backbone is both the coordinate anchor and a donor
candidate, so this changes what the scan sees.

The adapter already read the reference to pin the width; it now writes that
sequence, and refuses a projection whose width differs from it rather than
emitting ragged rows. A new test asserts the backbone row equals the
backbone for all five backends against fakes of their tools, in both query
orders. Only progressiveMauve failed it; the other four cases are guards.
progressiveMauve itself is not installed here or in CI, so the fix is
checked against a fake XMFA only.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: SibeliaZ on macOS, and the backbone invariant with real binaries (outside the spec)

Found while writing this plan, by running the Task 2 invariant against the installed binaries. SibeliaZ 1.2.7's wrapper uses two more GNU-only constructs that the adapter's BSD patch list does not cover. `paste <(...) <(...) --delimiters ' '`: BSD `paste` reads `--delimiters` as a file name, so the MAF holds headers only. `xargs -I @ ... bash -c "block_align @ '$outfile'"`: BSD `xargs` limits an `-I` replacement to 255 bytes by default, and the command holds two absolute paths. Either way the adapter raised "SibeliaZ wrote an empty MAF ... most often spoa was OOM-killed". Reproduced on this Mac with two 3 kb genomes: the unpatched wrapper writes 0 alignment rows; with `paste -d ' '` it writes 2 for a short output path; for a 209-character output path it still writes 0 until `xargs -S 4096` is added, then 4. `sibeliaz` is the default backend, so on macOS the default `tessera msa` could not have produced an alignment with this SibeliaZ release. Whether this, rather than memory, explains the "sibeliaz OOMs on orthopox" note in `validation/README.md` was not tested (it needs the orthopox data); that note is not changed here.

A reviewer can reject this task independently of the others. If it is rejected, mark the `sibeliaz` case of the new binary test `xfail` on `sys.platform == "darwin"` with the reason above, and say so in the PR.

**Files:**
- Modify: `src/tessera/aligners/sibeliaz.py` (lines 24-43 `_BSD_PATCHES` and its comment; lines 104-111 the empty-MAF message, at `4b35e39`)
- Create: `tests/integration/test_backbone_invariant_binaries.py`
- Modify: `tests/unit/test_input_safety.py` (append)
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes: nothing from Task 2's code (the binary test reads rows itself).
- Produces: two more entries in `sibeliaz._BSD_PATCHES`; no signature change.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_input_safety.py`:

```python


def test_sibeliaz_bsd_wrapper_rewrites_paste_and_xargs(monkeypatch, tmp_path: Path) -> None:
    # BSD paste takes no long options and reads '--delimiters' as a file name, so the
    # wrapper's MAF step wrote headers only and the adapter reported an empty MAF
    # (blamed on memory) for every input on macOS.
    wrapper = tmp_path / "sibeliaz"
    wrapper.write_text(
        'find $outdir -name "*.tmp" -printf "%p\\n"\n'
        "paste <(echo \"$header\") <(echo \"$output\") --delimiters ' ' >> \"$alignment_output\"\n"
        "find $outdir -name \"*.tmp\" | xargs -I @ -P $threads bash -c \"block_align @\"\n"
    )
    monkeypatch.setattr(sibeliaz.sys, "platform", "darwin")
    monkeypatch.setattr(sibeliaz.shutil, "which", lambda name: str(wrapper))
    cmd = sibeliaz._sibeliaz_invocation(tmp_path, _LOG)
    patched = Path(cmd[-1]).read_text()
    assert "--delimiters" not in patched
    assert "paste -d ' ' <(echo \"$header\") <(echo \"$output\")" in patched
    # BSD xargs caps an -I replacement at 255 bytes; two absolute paths exceed it.
    assert "xargs -S 4096 -I @ -P $threads" in patched
```

Create `tests/integration/test_backbone_invariant_binaries.py`:

```python
"""The backbone row equals the backbone sequence, with the real aligner binaries.

The unit suite checks this invariant against fakes of each tool
(``tests/unit/test_backbone_invariant.py``). Here the installed backends align two
queries that each cover one contig of a two-contig backbone. Each case is skipped
when its binary is not on PATH.
"""

from __future__ import annotations

import logging
import random
import shutil
from pathlib import Path

import pytest

from tessera.aligners.base import AlignParams
from tessera.aligners.base import registry as aligner_registry

_LOG = logging.getLogger("tessera.test")

pytestmark = pytest.mark.requires_binary

_BINARIES = {
    "mafft": "mafft",
    "minimap2": "minimap2",
    "sibeliaz": "sibeliaz",
    "progressivemauve": "progressiveMauve",
}
# The MAF-based backends lay the backbone out from the blocks, so a backbone position no
# block covers is a gap in the backbone row (no other genome has a base there either).
# For them the invariant is that every base in the row is the backbone's own base at
# that position; the pairwise backends write the backbone itself.
_MAF_BASED = {"sibeliaz"}


def _mutate(seq: str, rng: random.Random, rate: float) -> str:
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def _read_rows(path: Path) -> dict[str, str]:
    rows: dict[str, str] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            rows[name] = ""
        elif name is not None:
            rows[name] += line
    return rows


@pytest.mark.parametrize("aligner", sorted(_BINARIES))
def test_backbone_row_is_the_backbone_with_the_real_tool(aligner: str, tmp_path: Path) -> None:
    if shutil.which(_BINARIES[aligner]) is None:
        pytest.skip(f"{_BINARIES[aligner]} not on PATH")
    rng = random.Random(7)
    c1 = "".join(rng.choice("ACGT") for _ in range(3000))
    c2 = "".join(rng.choice("ACGT") for _ in range(3000))
    ref = tmp_path / "ref.fasta"
    ref.write_text(f">c1\n{c1}\n>c2\n{c2}\n")
    qa = tmp_path / "qryA.fasta"
    qa.write_text(f">a1\n{_mutate(c1, rng, 0.01)}\n")
    qb = tmp_path / "qryB.fasta"
    qb.write_text(f">b1\n{_mutate(c2, rng, 0.01)}\n")

    result = aligner_registry.create(aligner).align(
        [ref, qa, qb], ref, tmp_path / "out", AlignParams(threads=2), _LOG
    )
    rows = _read_rows(result.msa_fasta)
    backbone = (c1 + c2).upper()
    assert {len(r) for r in rows.values()} == {len(backbone)}
    row = rows["ref"].upper()
    if aligner in _MAF_BASED:
        misplaced = [i for i, b in enumerate(row) if b != "-" and b != backbone[i]]
        assert misplaced == []
        assert row.count("-") < len(backbone) // 10  # most of the backbone is placed
    else:
        assert row == backbone
```

- [ ] **Step 2: Run them, expected failure**

```bash
pytest tests/unit/test_input_safety.py -k paste_and_xargs -q
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
pytest tests/integration/test_backbone_invariant_binaries.py -q -rs
```

Expected: the unit test fails with `AssertionError: assert '--delimiters' not in 'find $outdi...'`. The binary test, on macOS with the env: `1 failed, 2 passed, 1 skipped` -- `[sibeliaz]` fails with `OutputError: SibeliaZ wrote an empty MAF (no alignment rows) ...`, mafft and minimap2 pass, progressiveMauve skips (`progressiveMauve not on PATH`). On Linux the sibeliaz case is expected to pass before the change (GNU tools); in CI only mafft runs.

- [ ] **Step 3: Implement**

In `src/tessera/aligners/sibeliaz.py`, replace

```python
# SibeliaZ's bash wrapper uses GNU/Linux-only constructs in its alignment step
# (`free`, GNU `find -printf`, `stat -c`, `mktemp --suffix`). On macOS/BSD these
# fail silently, so blocks are found but no MAF is written. We run a patched copy
# of the wrapper on those platforms.
```

with

```python
# SibeliaZ's bash wrapper uses GNU/Linux-only constructs in its alignment step
# (`free`, GNU `find -printf`, `stat -c`, `mktemp --suffix`, `paste --delimiters`, and
# an `xargs -I` command longer than BSD's 255-byte default). On macOS/BSD these fail
# silently, so blocks are found but no MAF rows are written. We run a patched copy
# of the wrapper on those platforms.
```

and replace

```python
    ("ulimit $memory_min", ":"),
)
```

with

```python
    ("ulimit $memory_min", ":"),
    # BSD paste has no long options and reads '--delimiters' as a file name, so every
    # MAF row is lost and the MAF holds headers only (SibeliaZ 1.2.7).
    (
        "paste <(echo \"$header\") <(echo \"$output\") --delimiters ' '",
        "paste -d ' ' <(echo \"$header\") <(echo \"$output\")",
    ),
    # BSD xargs limits an -I replacement to 255 bytes; the command it builds holds two
    # absolute paths, so a long output directory stopped every block alignment.
    ("xargs -I @ -P $threads", "xargs -S 4096 -I @ -P $threads"),
)
```

In the same file, replace the empty-MAF message

```python
            raise OutputError(
                f"SibeliaZ wrote an empty MAF (no alignment rows) at {maf}. Its "
                "blocks were found but the per-block spoa alignment produced "
                "nothing -- most often spoa was OOM-killed on a large collinear "
                "block. Give the machine more memory or use the progressivemauve "
                "backend."
            )
```

with

```python
            raise OutputError(
                f"SibeliaZ wrote an empty MAF (no alignment rows) at {maf}. Its "
                "blocks were found but the per-block alignment step wrote nothing; "
                "the wrapper runs it with errors suppressed. Known causes are spoa "
                "being killed for lack of memory on a large collinear block, and a "
                "wrapper command this platform does not support. Give the machine "
                "more memory, or use the minimap2 or mafft backend."
            )
```

Append to the "Alignment integrity" group in `CHANGELOG.md` (immediately above the blank line before `Data safety and alignment:`):

```markdown
- **The `sibeliaz` backend wrote no alignment on macOS.** Its wrapper uses
  `paste --delimiters`, which BSD `paste` reads as a file name, and an `xargs -I` command
  longer than BSD's 255-byte default; either left the MAF without rows, and the adapter
  reported this as spoa running out of memory. The patched wrapper now rewrites both.
  Reproduced and checked with SibeliaZ 1.2.7 on macOS; Linux is unaffected.
```

- [ ] **Step 4: Run, expected pass**

Run the two commands of Step 2 again.
Expected: the unit test passes; the binary test gives `3 passed, 1 skipped` (progressiveMauve) on macOS with the env.

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: clean, clean, `769 passed, 9 deselected` (one fast test and four binary-gated cases added).

- [ ] **Step 6: Commit**

```bash
git add src/tessera/aligners/sibeliaz.py tests/unit/test_input_safety.py \
        tests/integration/test_backbone_invariant_binaries.py CHANGELOG.md
git commit -m "$(cat <<'EOF'
Patch two more GNU-only wrapper commands so SibeliaZ runs on macOS

Running the backbone invariant against the installed binaries showed that
the default backend wrote no alignment on macOS. SibeliaZ 1.2.7's wrapper
calls `paste ... --delimiters ' '`, which BSD paste reads as a file name,
and builds an `xargs -I` command from two absolute paths, which exceeds
BSD xargs's 255-byte replacement limit. Either left the MAF with headers
only, and the adapter reported it as spoa running out of memory. With two
3 kb genomes the unpatched wrapper wrote 0 alignment rows; with both
patches it wrote 4, also under a 209-character output path.

The patched copy of the wrapper now rewrites both commands, and the
empty-MAF message names a platform incompatibility as a cause next to
memory. The new binary-gated test checks the backbone row of each
installed backend: equal to the backbone for the pairwise backends, and
for sibeliaz no misplaced base (a MAF backend leaves positions no block
covers as gaps). Whether the recorded orthopox "OOM" had this cause was
not tested.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: Refuse a sequence ID repeated inside one genome (audit items: C2)

The sibeliaz, cactus and minimap2 layouts are dicts keyed on the sequence ID. Two records named `c1` in one file collapse to one entry: under sibeliaz the backbone row comes back half gapped and a donor block is written onto the wrong contig's columns; under minimap2 every query row shifts. The output is full width and well-formed. `msa/build.py:65` is the only path to `aligner.align`, and it stages first, so one check in `stage_genomes` covers every backend. Mafft and progressiveMauve do not use the ID, so for them this is a new refusal of input that used to align; a FASTA with a repeated ID is malformed for most tools, and the message says how to fix it.

**Files:**
- Modify: `src/tessera/core/io.py` (line 299 in `stage_genomes`; new function above `_WHITESPACE` at line 219, at `4b35e39`)
- Modify: `tests/unit/test_io.py` (append)
- Modify: `docs/aligners.md` (staging paragraph, lines 44-46), `CHANGELOG.md`

**Interfaces:**
- Consumes: nothing.
- Produces: `_reject_repeated_sequence_ids(source: Path) -> None` in `core/io.py`, raising `UserInputError` naming the file, the ID and both line numbers. Called for every collection file and the query before anything is staged.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_io.py`:

```python


# --- a sequence ID repeated inside one genome ------------------------------------
# The sibeliaz, cactus and minimap2 backends lay a genome's contigs out by sequence ID,
# in dicts keyed on that ID. Two records with one ID collapse to one entry: a donor
# block lands on the wrong contig's columns, or the backbone row comes back half gapped,
# and the alignment still looks well-formed.

def _collection_with(tmp_path: Path, refa_text: str) -> Path:
    collection = tmp_path / "collection"
    collection.mkdir()
    (collection / "refA.fasta").write_text(refa_text)
    (collection / "refB.fasta").write_text(">refB\nACGTACGT\n")
    return collection


def test_a_sequence_id_repeated_within_a_genome_is_rejected(tmp_path: Path) -> None:
    collection = _collection_with(tmp_path, ">c1\nACGT\n>c2\nGGGG\n>c1\nTTTT\n")
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGTACGT\n")
    with pytest.raises(UserInputError, match=r"refA\.fasta.*'c1'.*lines 1 and 5"):
        stage_genomes(query, collection, tmp_path / "genomes", _LOG)
    assert not (tmp_path / "genomes").exists()  # rejected before anything is staged


def test_a_sequence_id_repeated_within_the_query_is_rejected(tmp_path: Path) -> None:
    collection = _collection_with(tmp_path, ">refA\nACGTACGT\n")
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGT\n>q\nTTTT\n")
    with pytest.raises(UserInputError, match=r"query\.fasta.*'q'"):
        stage_genomes(query, collection, tmp_path / "genomes", _LOG)


def test_a_repeated_id_is_found_in_a_gzipped_genome_by_its_first_token(tmp_path: Path) -> None:
    # NCBI downloads arrive gzipped, and two records can share an accession while their
    # descriptions differ; the ID is the first token, which is what the aligners read.
    collection = _collection_with(tmp_path, ">refA\nACGTACGT\n")
    with gzip.open(collection / "refC.fasta.gz", "wt") as fo:
        fo.write(">NC_1.1 chromosome\nACGT\n>NC_1.1 plasmid pX\nTTTT\n")
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGTACGT\n")
    with pytest.raises(UserInputError, match=r"refC\.fasta\.gz.*'NC_1\.1'"):
        stage_genomes(query, collection, tmp_path / "genomes", _LOG)


def test_one_id_in_two_different_genomes_is_still_staged(tmp_path: Path) -> None:
    # Guard (passes before and after): only a repeat inside one file is rejected here.
    # An ID shared between genomes is the sibeliaz adapter's own check.
    collection = _collection_with(tmp_path, ">contig_1\nACGTACGT\n")
    (collection / "refB.fasta").write_text(">contig_1\nACGTACGA\n")
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGTACGT\n")
    staged, _ = stage_genomes(query, collection, tmp_path / "genomes", _LOG)
    assert len(staged) == 3
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_io.py -q`
Expected: `3 failed, 11 passed`; each of the three new rejection tests fails with `Failed: DID NOT RAISE <class 'tessera.core.errors.UserInputError'>`. The fourth new test is a guard and passes.

- [ ] **Step 3: Implement**

In `src/tessera/core/io.py`, in `stage_genomes`, replace

```python
    for source in [*collection_files, query]:
        _require_fasta(source)
```

with

```python
    for source in [*collection_files, query]:
        _require_fasta(source)
        _reject_repeated_sequence_ids(source)
```

and insert above `_WHITESPACE = re.compile(rb"\s")`:

```python
def _reject_repeated_sequence_ids(source: Path) -> None:
    """Refuse a genome file in which two records share a sequence ID.

    The ID is the header's first token, which is how the sibeliaz, cactus and minimap2
    backends identify a contig. They lay a genome's contigs out in tables keyed on it,
    so two records with one ID collapse to a single entry: a block aligned to one
    record is written onto the other's columns, or the backbone row loses one record
    entirely. The result is a full-width alignment with nothing to show it is wrong.
    """
    seen: dict[str, int] = {}
    opener = gzip.open if source.name.endswith(".gz") else open
    with opener(source, "rb") as fo:
        for lineno, line in enumerate(fo, start=1):
            if not line.startswith(b">"):
                continue
            tokens = line[1:].split()
            if not tokens:
                continue  # an unnamed record; the backends that need a name say so
            seqid = tokens[0].decode("utf-8", errors="replace")
            if seqid in seen:
                raise UserInputError(
                    f"{source} has two records with the sequence ID '{seqid}' (lines "
                    f"{seen[seqid]} and {lineno}). Contigs are placed in the alignment by "
                    "sequence ID, so the two cannot be told apart -- rename one."
                )
            seen[seqid] = lineno


```

In `docs/aligners.md`, replace

```markdown
Genomes are staged before alignment. A FASTA with whitespace inside its sequence lines
(trailing spaces, tabs, Windows line endings) is aligned from a cleaned copy and reported
in a warning; the input file itself is not modified.
```

with

```markdown
Genomes are staged before alignment. A FASTA with whitespace inside its sequence lines
(trailing spaces, tabs, Windows line endings) is aligned from a cleaned copy and reported
in a warning; the input file itself is not modified. A genome file in which two records
share a sequence ID (the first word of the header) is refused, naming the file and the
ID: several backends place contigs by that ID, and the two records would be written onto
each other's columns.
```

Append to the "Alignment integrity" group in `CHANGELOG.md`:

```markdown
- **Two records with one sequence ID in a genome were silently merged.** The `sibeliaz`,
  `cactus` and `minimap2` layouts are keyed on the ID, so the second record's blocks were
  written onto the first one's columns (or the backbone row lost a record) in an
  otherwise well-formed alignment. Such a file is now refused at staging, for every
  backend, naming the file, the ID and both lines.
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_io.py tests/unit/test_input_safety.py tests/integration/test_genome_io.py -q`
Expected: all pass.

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: clean, clean, `773 passed, 9 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/core/io.py tests/unit/test_io.py docs/aligners.md CHANGELOG.md
git commit -m "$(cat <<'EOF'
Refuse a genome file that repeats a sequence ID

The sibeliaz and cactus contig layouts and the minimap2 reference offsets
are dicts keyed on the sequence ID. A backbone with two records named c1
collapsed to one entry: the reference length came out short with a single
offset, the backbone row was half gapped, and a donor block aligned to the
first record was written onto the second one's columns. Under minimap2 the
reference row was right and every query row shifted. Each produced a
full-width, well-formed, wrong alignment with no warning.

Staging now reads each genome's headers (gzip included) and refuses a
repeated first token, naming the file, the ID and both lines. Staging is
the only path to an aligner, so one check covers every backend. The
trade-off: mafft and progressiveMauve do not use the ID, and a file that
used to align under them is now refused; the message says how to fix it.
An ID shared between two genomes is still allowed here and remains the
sibeliaz adapter's own check.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: Share the minimap2 thread budget (audit items: C3)

`minimap2.py` passed `-t <threads>` to each process and ran the pool with the same count, so `--threads 16` started 16 processes of 16 threads. Commit `44e6b3f` fixed the same pattern in MAFFT inline. The rule is moved into one helper in `aligners/base.py` and used by both adapters. The memory side named in the spec (each process indexes the backbone) was not measured; for viral backbones the index is small, so the thread count is the part this task fixes.

**Files:**
- Modify: `src/tessera/aligners/base.py` (new function above `class OutputKind` at line 27, at `4b35e39`)
- Modify: `src/tessera/aligners/minimap2.py` (line 23 import, line 62 `threads`, line 77 `parallel_map`)
- Modify: `src/tessera/aligners/mafft.py` (line 38 import, lines 126-138 the budget block)
- Modify: `tests/unit/test_aligner_backends.py` (append)
- Modify: `docs/aligners.md`, `CHANGELOG.md`

**Interfaces:**
- Consumes: nothing.
- Produces: `split_thread_budget(threads: int, tasks: int, job_limit: int | None = None) -> tuple[int, int]` in `tessera.aligners.base`, returning `(jobs, threads_per_job)` with `jobs * threads_per_job <= max(1, threads)`.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_aligner_backends.py`:

```python


# --- minimap2 shares its thread budget across concurrent jobs ----------------------
def _minimap2_runs(monkeypatch, tmp_path: Path, *, threads: int, n_queries: int
                   ) -> tuple[list[int], list[list[str]]]:
    """Run the minimap2 adapter against a fake `minimap2`; return (pool sizes, commands)."""
    commands: list[list[str]] = []
    workers: list[int] = []
    real_map = mm_mod.parallel_map

    def spy_map(fn, items, n, **kw):
        workers.append(n)
        return real_map(fn, items, n, **kw)

    def fake_run(caps, cmd, **kw):
        commands.append([str(c) for c in cmd])
        kw["stdout_path"].write_text(
            "@SQ\tSN:ref\tLN:12\nq\t0\tref\t1\t60\t12M\t*\t0\t0\tACGTACGTACGT\t*\n"
        )
        return ""

    monkeypatch.setattr(mm_mod, "parallel_map", spy_map)
    monkeypatch.setattr(mm_mod, "run_tool", fake_run)
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref\nACGTACGTACGT\n")
    genomes = [ref]
    for i in range(n_queries):
        qry = tmp_path / f"qry{i}.fasta"
        qry.write_text(f">qry{i}\nACGTACGTACGT\n")
        genomes.append(qry)
    mm_mod.Minimap2Aligner().align(genomes, ref, tmp_path / "out",
                                   AlignParams(threads=threads), _LOG)
    return workers, commands


def test_minimap2_splits_threads_across_concurrent_jobs(monkeypatch, tmp_path: Path) -> None:
    """--threads 4 used to start four minimap2 processes of four threads each."""
    workers, commands = _minimap2_runs(monkeypatch, tmp_path, threads=4, n_queries=5)
    assert workers == [4]
    assert len(commands) == 5
    assert {cmd[cmd.index("-t") + 1] for cmd in commands} == {"1"}


def test_minimap2_gives_a_single_job_all_threads(monkeypatch, tmp_path: Path) -> None:
    workers, commands = _minimap2_runs(monkeypatch, tmp_path, threads=4, n_queries=1)
    assert workers == [1]
    assert commands[0][commands[0].index("-t") + 1] == "4"


def test_minimap2_never_starts_more_threads_than_the_budget(monkeypatch, tmp_path: Path) -> None:
    # A budget that does not divide evenly: five threads over two queries. Two jobs of
    # two threads stay within it; two jobs of five, or of three, would not.
    workers, commands = _minimap2_runs(monkeypatch, tmp_path, threads=5, n_queries=2)
    per_job = {int(cmd[cmd.index("-t") + 1]) for cmd in commands}
    assert workers == [2]
    assert per_job == {2}
    assert workers[0] * max(per_job) <= 5
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_aligner_backends.py -k minimap2 -q`
Expected: `3 failed, 2 passed`. Observed: `assert {'4'} == {'1'}` (five jobs of four threads); `assert [4] == [1]` (a pool of four for one genome); `assert [5] == [2]`. The existing `test_minimap2_forwards_preset_and_threads` (three threads, one query, `-t 3`) passes before and after.

- [ ] **Step 3: Implement**

In `src/tessera/aligners/base.py`, insert above `class OutputKind(Enum):`:

```python
def split_thread_budget(
    threads: int, tasks: int, job_limit: int | None = None
) -> tuple[int, int]:
    """Divide a ``--threads`` budget between concurrent tool runs.

    Returns ``(jobs, threads_per_job)``: at most one job per task, never more jobs than
    threads, and never fewer than one of either. ``job_limit`` caps the number of jobs
    further (a memory bound, say). ``jobs * threads_per_job`` never exceeds the budget,
    so ``--threads 16`` starts at most sixteen threads in all rather than sixteen
    processes of sixteen threads each.
    """
    budget = max(1, threads)
    jobs = max(1, min(budget, tasks))
    if job_limit is not None:
        jobs = max(1, min(jobs, job_limit))
    return jobs, budget // jobs


```

In `src/tessera/aligners/minimap2.py`, replace

```python
from .base import Aligner, AlignParams, AlignResult
```

with

```python
from .base import Aligner, AlignParams, AlignResult, split_thread_budget
```

replace

```python
        preset = str(params.extra.get("preset", "asm20"))
        threads = str(max(1, params.threads))
        queries = [g for g in genomes if g != reference]
```

with

```python
        preset = str(params.extra.get("preset", "asm20"))
        queries = [g for g in genomes if g != reference]
        # The thread budget is shared between the concurrent mapping jobs; giving each
        # job all of it started --threads processes of --threads threads each.
        jobs, per_job = split_thread_budget(params.threads, len(queries))
        threads = str(per_job)
```

and replace

```python
        rows = parallel_map(map_genome, queries, params.threads, logger=logger)
```

with

```python
        rows = parallel_map(map_genome, queries, jobs, logger=logger)
```

In `src/tessera/aligners/mafft.py`, replace

```python
from .base import Aligner, AlignParams, AlignResult
```

with

```python
from .base import Aligner, AlignParams, AlignResult, split_thread_budget
```

and replace

```python
        budget = max(1, params.threads)
        jobs = max(1, min(budget, len(queries)))
        if not sixmer:
            limit = _multipair_job_limit(len(ref_seq), max(longest_added, 1))
            if limit < jobs:
                logger.info(
                    "Running %d MAFFT job(s) at a time instead of %d to bound memory.",
                    limit, jobs,
                )
                jobs = limit
        threads = str(budget // jobs)
```

with

```python
        jobs, per_job = split_thread_budget(params.threads, len(queries))
        if not sixmer:
            limit = _multipair_job_limit(len(ref_seq), max(longest_added, 1))
            if limit < jobs:
                logger.info(
                    "Running %d MAFFT job(s) at a time instead of %d to bound memory.",
                    limit, jobs,
                )
                jobs, per_job = split_thread_budget(params.threads, len(queries), limit)
        threads = str(per_job)
```

In `docs/aligners.md`, insert after the MAFFT memory paragraph (the one ending ``pass `sixmerpair=true`.``) a new paragraph:

```markdown
Threads with `minimap2`: genomes are mapped concurrently, and `--threads` is shared
between the concurrent `minimap2` processes rather than given to each, as with `mafft`.
```

Append to the "Alignment integrity" group in `CHANGELOG.md`:

```markdown
- **`minimap2` multiplied its thread budget.** `--threads N` started N concurrent
  `minimap2` processes with N threads each. The budget is now shared between them, by the
  same rule the `mafft` backend uses.
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_aligner_backends.py tests/unit/test_backbone_invariant.py -q`
Expected: all pass, including the eight existing MAFFT thread and memory tests.

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: clean, clean, `776 passed, 9 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/aligners/base.py src/tessera/aligners/minimap2.py \
        src/tessera/aligners/mafft.py tests/unit/test_aligner_backends.py \
        docs/aligners.md CHANGELOG.md
git commit -m "$(cat <<'EOF'
Share the minimap2 thread budget between concurrent jobs

The minimap2 adapter passed -t <threads> to every process and also ran
the mapping pool with <threads> workers, so --threads 16 started sixteen
processes of sixteen threads. MAFFT had the same pattern and was fixed in
44e6b3f; minimap2, the backend used for the 200 kb orthopox and mpox
cases, was not.

The MAFFT rule (at most one job per genome, threads divided between the
jobs, a further cap where memory requires it) moves into
split_thread_budget in aligners/base.py and both adapters call it, so the
rule is written once. MAFFT's behaviour is unchanged; its existing thread
and memory tests pass as before. The per-process index memory the audit
also named was not measured.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Name a genome the aligner could not place (audit items: C5)

`docs/aligners.md` says a genome the aligner cannot place is named in a warning. Only `maf_to_fasta` (sibeliaz, cactus) does so. Under minimap2, mafft and progressiveMauve the genome becomes an all-gap row in silence and surfaces later in `run.py` as "only 0% A/C/G/T outside gaps ... Check its alphabet", which points at the wrong cause. This task adds the warning to the three pairwise adapters through one helper, and in the scan reports an all-gap row as unaligned, separately from the alphabet warning. The scan change is a warning only: no row is added, removed or treated differently.

**Files:**
- Modify: `src/tessera/aligners/base.py` (new function next to `split_thread_budget`)
- Modify: `src/tessera/aligners/minimap2.py`, `src/tessera/aligners/mafft.py` (imports; before `msa = out_dir / "msa.fasta"`)
- Modify: `src/tessera/aligners/progressivemauve.py` (`_concatenate` returns the query rows; the call site)
- Modify: `src/tessera/recomb/similarity.py` (lines 121-141 `low_canonical_records`; new function above `_read_alignment` at line 144, at `4b35e39`)
- Modify: `src/tessera/recomb/run.py` (line 40 import; lines 287-293 the warning loop; new helper above `_discover_lineage_tsv` at line 195)
- Modify: `tests/unit/test_aligner_backends.py`, `tests/unit/test_similarity.py`, `tests/unit/test_run.py` (append)
- Modify: `tests/unit/test_aligner_params.py` line 50 (existing test changes, see Step 3)
- Modify: `docs/aligners.md` (lines 41-42), `CHANGELOG.md`

**Interfaces:**
- Consumes: `split_thread_budget` imports from Task 5 (the import lines below include it).
- Produces: `warn_unplaced(rows: Sequence[tuple[str, str]], backbone: str, logger: logging.Logger) -> None` in `tessera.aligners.base`; `_concatenate(...) -> list[tuple[str, str]]` (the query rows); `unaligned_records(rows: dict[str, np.ndarray]) -> list[str]` in `tessera.recomb.similarity`; `low_canonical_records` no longer lists all-gap rows; `_warn_unusable_rows(rows, logger: logging.Logger) -> None` in `tessera.recomb.run`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_aligner_backends.py`:

```python


# --- a genome the aligner did not place is named -----------------------------------
# Not under the "tessera" logger: the CLI tests switch its propagation off, and caplog
# only sees records that reach the root logger.
_WARN_LOG = logging.getLogger("aligner_unplaced_test")


def _ref_and_queries(tmp_path: Path) -> list[Path]:
    paths = []
    for name, seq in (("ref", "ACGTACGTACGT"), ("near", "ACGTACGTACGT"), ("far", "TTTTTTTTTTTT")):
        path = tmp_path / f"{name}.fasta"
        path.write_text(f">{name}\n{seq}\n")
        paths.append(path)
    return paths


def test_minimap2_names_a_genome_it_could_not_map(monkeypatch, tmp_path: Path, caplog) -> None:
    def fake_run(caps, cmd, **kw):
        if "far" in str(cmd[-1]):
            body = "far\t4\t*\t0\t0\t*\t*\t0\t0\tTTTTTTTTTTTT\t*\n"  # unmapped
        else:
            body = "near\t0\tref\t1\t60\t12M\t*\t0\t0\tACGTACGTACGT\t*\n"
        kw["stdout_path"].write_text("@SQ\tSN:ref\tLN:12\n" + body)
        return ""

    monkeypatch.setattr(mm_mod, "run_tool", fake_run)
    genomes = _ref_and_queries(tmp_path)
    with caplog.at_level(logging.WARNING, logger=_WARN_LOG.name):
        result = mm_mod.Minimap2Aligner().align(genomes, genomes[0], tmp_path / "out",
                                                AlignParams(threads=1), _WARN_LOG)
    assert _read_msa(result.msa_fasta)["far"] == "-" * 12
    assert "1 genome(s)" in caplog.text
    assert "far" in caplog.text and "near" not in caplog.text


def test_mafft_names_a_genome_it_could_not_place(monkeypatch, tmp_path: Path, caplog) -> None:
    def fake_run(caps, cmd, **kw):
        row = "-" * 12 if "far" in str(cmd[cmd.index("--addfragments") + 1]) else "ACGTACGTACGT"
        kw["stdout_path"].write_text(f">ref\nACGTACGTACGT\n>added\n{row}\n")
        return ""

    monkeypatch.setattr(mafft_mod, "run_tool", fake_run)
    genomes = _ref_and_queries(tmp_path)
    with caplog.at_level(logging.WARNING, logger=_WARN_LOG.name):
        mafft_mod.MafftAligner().align(genomes, genomes[0], tmp_path / "out",
                                       AlignParams(threads=1), _WARN_LOG)
    assert "far" in caplog.text and "near" not in caplog.text


def test_progressivemauve_names_a_genome_in_no_block(monkeypatch, tmp_path: Path, caplog) -> None:
    def fake_run(caps, cmd, **kw):
        cmd = [str(c) for c in cmd]
        xmfa = Path(cmd[cmd.index("--output") + 1])
        ref_path, query_path = cmd[-2], cmd[-1]
        header = f"#Sequence1File\t{ref_path}\n#Sequence2File\t{query_path}\n"
        if "far" in query_path:  # the reference aligned to nothing
            xmfa.write_text(header + f"> 1:1-12 + {ref_path}\nACGTACGTACGT\n=\n")
        else:
            xmfa.write_text(header + f"> 1:1-12 + {ref_path}\nACGTACGTACGT\n"
                            f"> 2:1-12 + {query_path}\nACGTACGTACGT\n=\n")
        return ""

    monkeypatch.setattr(pm_mod, "run_tool", fake_run)
    genomes = _ref_and_queries(tmp_path)
    with caplog.at_level(logging.WARNING, logger=_WARN_LOG.name):
        pm_mod.ProgressiveMauveAligner().align(genomes, genomes[0], tmp_path / "out",
                                               AlignParams(threads=1), _WARN_LOG)
    assert "far" in caplog.text and "near" not in caplog.text
```

Append to `tests/unit/test_similarity.py`:

```python


def test_an_all_gap_row_is_unaligned_not_an_alphabet_problem(tmp_path: Path) -> None:
    """An all-gap row means the aligner placed nothing of that genome; reporting it
    as an unreadable alphabet sent the user looking in the wrong place."""
    from tessera.recomb.similarity import unaligned_records

    msa = write_fasta(tmp_path / "m.fasta", {
        "good": "ACGTACGTAC",
        "unplaced": "----------",
        "unknown": "NNNNNNNNNN",   # placed, but no canonical base: still an alphabet flag
    })
    rows = _read_alignment(str(msa))
    assert unaligned_records(rows) == ["unplaced"]
    flagged = dict(low_canonical_records(rows, min_fraction=0.5))
    assert "unplaced" not in flagged
    assert flagged["unknown"] == 0.0
```

Append to `tests/unit/test_run.py`:

```python


def test_scan_names_an_all_gap_row_as_unaligned(tmp_path: Path, caplog) -> None:
    import numpy as np

    from tessera.recomb.run import _warn_unusable_rows

    rows = {
        "good": np.frombuffer(b"ACGTACGTAC", dtype=np.uint8),
        "unplaced": np.frombuffer(b"----------", dtype=np.uint8),
        "junk": np.frombuffer(b"NNNNNNNNNN", dtype=np.uint8),
    }
    log = logging.getLogger("scan_rows_test")
    with caplog.at_level(logging.WARNING, logger="scan_rows_test"):
        _warn_unusable_rows(rows, log)
    lines = caplog.text.splitlines()
    unplaced = [line for line in lines if "'unplaced'" in line]
    junk = [line for line in lines if "'junk'" in line]
    assert len(unplaced) == 1 and "no aligned base" in unplaced[0]
    assert "alphabet" not in unplaced[0]
    assert len(junk) == 1 and "alphabet" in junk[0]
```

- [ ] **Step 2: Run them, expected failure**

Run: `pytest tests/unit/test_aligner_backends.py tests/unit/test_similarity.py tests/unit/test_run.py -q`
Expected: `5 failed`. Observed: `assert '1 genome(s)' in ''` (minimap2); `assert ('far' in '')` (mafft, progressiveMauve); `ImportError: cannot import name 'unaligned_records' from 'tessera.recomb.similarity'`; `ImportError: cannot import name '_warn_unusable_rows' from 'tessera.recomb.run'`.

- [ ] **Step 3: Implement**

In `src/tessera/aligners/base.py`, insert after `split_thread_budget` (above `class OutputKind(Enum):`):

```python
def warn_unplaced(
    rows: Sequence[tuple[str, str]], backbone: str, logger: logging.Logger
) -> None:
    """Name every genome whose row holds no aligned base.

    The pairwise backends write a genome they could not place against the backbone as
    an all-gap row, so the alignment keeps one row per input genome. Such a row
    contributes nothing to the scan; without a warning the panel is silently smaller
    than the collection. The MAF-based backends say the same in their converter.
    """
    unplaced = sorted(label for label, row in rows if not row.strip("-"))
    if unplaced:
        logger.warning(
            "%d genome(s) were not aligned to any part of the backbone '%s' and are "
            "written as all-gap rows: %s. They contribute nothing to the scan; they "
            "may be too divergent for this aligner.",
            len(unplaced), backbone, ", ".join(unplaced),
        )


```

(`logging` and `Sequence` are already imported in `base.py`.)

In both `src/tessera/aligners/minimap2.py` and `src/tessera/aligners/mafft.py`, replace

```python
from .base import Aligner, AlignParams, AlignResult, split_thread_budget
```

with

```python
from .base import (
    Aligner,
    AlignParams,
    AlignResult,
    split_thread_budget,
    warn_unplaced,
)
```

and in each, replace

```python
        msa = out_dir / "msa.fasta"
        with open(msa, "w") as out:
            write_fasta_record(out, reference.stem, ref_seq)
```

with

```python
        warn_unplaced(rows, reference.stem, logger)
        msa = out_dir / "msa.fasta"
        with open(msa, "w") as out:
            write_fasta_record(out, reference.stem, ref_seq)
```

In `src/tessera/aligners/progressivemauve.py`, replace

```python
from .base import Aligner, AlignParams, AlignResult
```

with

```python
from .base import Aligner, AlignParams, AlignResult, warn_unplaced
```

replace

```python
        _concatenate(per_query, reference.stem, ref_seq, msa)
```

with

```python
        rows = _concatenate(per_query, reference.stem, ref_seq, msa)
        warn_unplaced(rows, reference.stem, logger)
```

and in `_concatenate` (as written in Task 2) make these three replacements:

```diff
 def _concatenate(
     per_query: list[tuple[str, Path]], reference_label: str, reference_seq: str,
     out_path: Path,
-) -> None:
+) -> list[tuple[str, str]]:
@@
-    never from those names.
+    never from those names. Returns the query rows as ``(label, row)``.
     """
+    rows: list[tuple[str, str]] = []
     with open(out_path, "w") as out:
@@
             write_fasta_record(out, label, row)
+            rows.append((label, row))
+    return rows
```

**Existing test that changes.** `tests/unit/test_aligner_params.py:50` stubs `_concatenate` with `lambda *a, **k: None`; the adapter now passes the return value to `warn_unplaced`, so the stub must return an empty list of rows. Replace

```python
    monkeypatch.setattr(pm, "_concatenate", lambda *a, **k: None)
```

with

```python
    monkeypatch.setattr(pm, "_concatenate", lambda *a, **k: [])  # no query rows
```

Without this change that test fails with `TypeError: 'NoneType' object is not iterable`.

In `src/tessera/recomb/similarity.py`, in `low_canonical_records`, replace

```python
    otherwise fail silently. Gaps are excluded from the denominator -- they are
    expected in an alignment and say nothing about the alphabet.
    """
    flagged: list[tuple[str, float]] = []
    for label, seq in rows.items():
        ungapped = int(np.count_nonzero(seq != GAP))
        if ungapped == 0:
            flagged.append((label, 0.0))
            continue
```

with

```python
    otherwise fail silently. Gaps are excluded from the denominator -- they are
    expected in an alignment and say nothing about the alphabet. An all-gap record
    has no alphabet to judge; :func:`unaligned_records` reports it instead.
    """
    flagged: list[tuple[str, float]] = []
    for label, seq in rows.items():
        ungapped = int(np.count_nonzero(seq != GAP))
        if ungapped == 0:
            continue
```

and insert above `def _read_alignment(msa_path: str)`:

```python
def unaligned_records(rows: dict[str, np.ndarray]) -> list[str]:
    """Records with no aligned base at all (all-gap rows), sorted by label.

    An aligner writes a genome it could not place against the backbone this way; the
    record then contributes nothing to the scan.
    """
    return sorted(label for label, seq in rows.items() if not np.any(seq != GAP))


```

In `src/tessera/recomb/run.py`, replace

```python
    low_canonical_records,
)
```

with

```python
    low_canonical_records,
    unaligned_records,
)
```

replace

```python
    # Only canonical bases count toward a window, so a record that is mostly something
    # else contributes almost nothing. Say so rather than let it drop out in silence.
    for label, fraction in low_canonical_records(bp_result.rows):
        logger.warning(
            "Sequence '%s' is only %.0f%% A/C/G/T outside gaps; it will contribute "
            "little to the scan. Check its alphabet.", label, 100.0 * fraction,
        )
```

with

```python
    _warn_unusable_rows(bp_result.rows, logger)
```

and insert above `def _discover_lineage_tsv(`:

```python
def _warn_unusable_rows(rows, logger: logging.Logger) -> None:
    """Name the records that will contribute little or nothing to the scan.

    Only canonical bases count toward a window, so such a record drops out in silence
    unless it is named. An all-gap row and an unreadable alphabet have different
    causes -- the aligner placed nothing of the genome, or the reader did not
    understand its letters -- so they are reported separately.
    """
    for label in unaligned_records(rows):
        logger.warning(
            "Sequence '%s' has no aligned base (an all-gap row): the aligner did not "
            "place it against the backbone, so it contributes nothing to the scan.",
            label,
        )
    for label, fraction in low_canonical_records(rows):
        logger.warning(
            "Sequence '%s' is only %.0f%% A/C/G/T outside gaps; it will contribute "
            "little to the scan. Check its alphabet.", label, 100.0 * fraction,
        )


```

In `docs/aligners.md`, replace

```markdown
order. A genome the aligner cannot place against the backbone at all is kept as an all-gap
row and named in a warning, so the alignment always has one row per input genome.
```

with

```markdown
order. A genome the aligner cannot place against the backbone at all is kept as an all-gap
row and named in a warning by every backend, so the alignment always has one row per
input genome. The scan names such a row again, as unaligned, when it reads the alignment.
```

Append to the "Alignment integrity" group in `CHANGELOG.md`:

```markdown
- **A genome the `minimap2`, `mafft` or `progressivemauve` backend could not place was
  dropped in silence.** It became an all-gap row with no warning (only the MAF-based
  backends named it), and the scan later reported it as an alphabet problem. These
  backends now name such genomes, and the scan reports an all-gap row as unaligned,
  separately from a record whose letters it could not read.
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_aligner_backends.py tests/unit/test_similarity.py tests/unit/test_run.py tests/unit/test_aligner_params.py tests/unit/test_backbone_invariant.py -q`
Expected: all pass.

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: clean, clean, `781 passed, 9 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/aligners/base.py src/tessera/aligners/minimap2.py \
        src/tessera/aligners/mafft.py src/tessera/aligners/progressivemauve.py \
        src/tessera/recomb/similarity.py src/tessera/recomb/run.py \
        tests/unit/test_aligner_backends.py tests/unit/test_similarity.py \
        tests/unit/test_run.py tests/unit/test_aligner_params.py \
        docs/aligners.md CHANGELOG.md
git commit -m "$(cat <<'EOF'
Name a genome the pairwise backends could not place

docs/aligners.md promised that a genome the aligner cannot place is named
in a warning. Only the MAF converter (sibeliaz, cactus) did so. Under
minimap2, mafft and progressiveMauve the genome became an all-gap row in
silence and surfaced in the scan as "only 0% A/C/G/T outside gaps ...
Check its alphabet", which points at the wrong cause.

The three pairwise adapters now call one helper, warn_unplaced, that names
every all-gap row. In the scan, an all-gap row is reported as unaligned and
no longer counted as a low-canonical record; a placed record with no
canonical base still gets the alphabet warning. Only warnings change: no
row is added, removed or treated differently by the callers.

One existing test stubbed progressiveMauve's _concatenate to return None;
it now returns the query rows, so the stub returns an empty list.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 7: Truncated and mismatched MAF input is reported (audit items: C6, Core-side test gap "a source size disagreeing with the declared contig length")

Two silent paths in `maf_to_fasta.py`. A row cut before its seventh field was dropped by `_iter_blocks` (`if len(parts) >= 7`), so the width guard written for truncation never saw it, and the genome became an all-gap row reported as divergent. And when the adapter supplies the backbone layout (`ref_contigs`), a block placed by the MAF's own `srcSize` was never compared with the declared length. Observed before the fix: with `c1` declared 4 bp and a MAF that says 6 bp, the backbone row came out as `----GGTT`: `c1`'s columns were empty and the block at `c1:4-6` was written over the first two columns of `c2`, whose bases are `TTTT`.

**Files:**
- Modify: `src/tessera/converters/maf_to_fasta.py` (line 129 `layout = list(ref_contigs)`; lines 204-219 `_iter_blocks`, at `4b35e39`)
- Modify: `tests/unit/test_converters.py` (append)
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes: nothing.
- Produces: `maf_to_fasta` raises `OutputError` for an `s` line with fewer than seven fields (naming the file and the line) and, when `ref_contigs` is given, for a backbone contig whose MAF `srcSize` differs from its declared length.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_converters.py`:

```python


def test_maf_row_cut_before_its_text_is_reported_against_the_file(tmp_path: Path) -> None:
    # A row cut off before its seventh field used to be dropped while parsing, so the
    # width guard never saw it and the genome came out as an all-gap row -- read as
    # divergence rather than as a truncated file.
    import pytest

    from tessera.core.errors import OutputError

    maf = tmp_path / "cut.maf"
    maf.write_text("a\ns ref.c 0 8 + 8 ACGTACGT\ns qry.x 0 8 + 8\n")
    with pytest.raises(OutputError, match=r"cut\.maf.*line 3.*qry\.x"):
        maf_to_fasta(maf, "ref", tmp_path / "msa.fasta")


def test_maf_backbone_source_size_must_match_the_declared_contig(tmp_path: Path) -> None:
    # The MAF says contig c1 is 6 bp; the backbone FASTA says 4. Laid out by the FASTA,
    # a block at c1:4-6 would land on c2's columns.
    import pytest

    from tessera.core.errors import OutputError

    maf = tmp_path / "size.maf"
    maf.write_text(
        "a\ns c1 4 2 + 6 GG\ns q1 0 2 + 2 GG\n\n"
        "a\ns c2 0 4 + 4 TTTT\ns q1 0 4 + 4 TTTT\n"
    )
    name_map = {"c1": "ref", "c2": "ref", "q1": "qry"}
    with pytest.raises(OutputError, match=r"size\.maf.*c1 \(6 bp in the MAF, 4 bp"):
        maf_to_fasta(maf, "ref", tmp_path / "msa.fasta", name_map=name_map,
                     ref_contigs=[("c1", 4), ("c2", 4)])
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_converters.py -q`
Expected: `2 failed, 9 passed`; both fail with `Failed: DID NOT RAISE <class 'tessera.core.errors.OutputError'>`. The existing `test_truncated_maf_row_is_reported_against_the_file` in `test_input_safety.py` (a seven-field row with short text) passes before and after.

- [ ] **Step 3: Implement**

In `src/tessera/converters/maf_to_fasta.py`, replace

```python
        layout = list(ref_contigs)
    else:
```

with

```python
        declared = dict(ref_contigs)
        resized = sorted(
            f"{name} ({size} bp in the MAF, {declared[name]} bp in the FASTA)"
            for name, size in seen_contigs.items() if size != declared[name]
        )
        if resized:
            # Laid out by the FASTA, a block placed by the MAF's own length would land
            # on the next contig's columns.
            raise OutputError(
                f"{maf_path} gives backbone sequence lengths that differ from the "
                f"backbone '{ref_key}' as staged: {', '.join(resized)}. The aligner's "
                "output does not belong to this input."
            )
        layout = list(ref_contigs)
    else:
```

and replace the start of `_iter_blocks`

```python
def _iter_blocks(maf_path: Path):
    block: list[_Row] = []
    with open(maf_path) as fo:
        for line in fo:
            if line.startswith("a"):
                if block:
                    yield block
                block = []
            elif line.startswith("s"):
                parts = line.split()
                # s src start size strand srcSize text
                if len(parts) >= 7:
                    block.append(_Row(
                        name=parts[1], start=int(parts[2]), size=int(parts[3]),
                        strand=parts[4], src_size=int(parts[5]), text=parts[6],
                    ))
```

with

```python
def _iter_blocks(maf_path: Path):
    block: list[_Row] = []
    with open(maf_path) as fo:
        for lineno, line in enumerate(fo, start=1):
            if line.startswith("a"):
                if block:
                    yield block
                block = []
            elif line.startswith("s"):
                parts = line.split()
                # s src start size strand srcSize text
                if len(parts) < 7:
                    # Dropping the row would leave its genome as an all-gap row, which
                    # reads as divergence rather than as a file cut off mid-write.
                    raise OutputError(
                        f"Truncated alignment row in {maf_path} at line {lineno}: "
                        f"{line.strip()!r} has {len(parts)} field(s), a MAF 's' line "
                        "has 7. The aligner's output looks incomplete; check that it "
                        "finished."
                    )
                block.append(_Row(
                    name=parts[1], start=int(parts[2]), size=int(parts[3]),
                    strand=parts[4], src_size=int(parts[5]), text=parts[6],
                ))
```

Append to the "Alignment integrity" group in `CHANGELOG.md`:

```markdown
- **Two MAF inputs were absorbed instead of reported.** A MAF row cut off before its
  sequence text was dropped, so its genome became an all-gap row that read as divergence;
  and a backbone contig whose length in the MAF differed from the staged FASTA had its
  blocks written onto the next contig's columns. Both now stop the run with a message
  naming the file.
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_converters.py tests/unit/test_input_safety.py tests/unit/test_backbone_invariant.py -q`
Expected: all pass.

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: clean, clean, `783 passed, 9 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/converters/maf_to_fasta.py tests/unit/test_converters.py CHANGELOG.md
git commit -m "$(cat <<'EOF'
Report a truncated MAF row and a backbone length mismatch

The MAF reader kept an 's' line only when it had seven fields, so a row
cut off before its text was dropped while parsing and the width guard
written for truncated output never saw it. Its genome became an all-gap
row, reported as sharing no block with the backbone -- divergence rather
than a truncated file. Such a row now stops the run, naming the file and
the line.

When the adapter supplies the backbone's contig layout, a block was placed
using the MAF's own srcSize without comparing it with the declared length.
With c1 declared 4 bp and a MAF saying 6 bp, a block at c1:4-6 was written
onto c2's columns. The two lengths are now compared and a mismatch is
reported. Both cases were driven with constructed MAFs; neither has been
seen from a real aligner.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 8: Install the Nextclade pool with `atomic_cache_dir` (audit items: C4)

`build_pool` staged into a temporary directory and installed it with a bare `os.replace`, which fails on a non-empty destination (`OSError: [Errno 66] Directory not empty` on macOS, `ENOTEMPTY` on Linux), surfacing as "Unexpected error". Two triggers: a concurrent run on the same dataset, and a cache directory holding only files the reader does not recognise (read as absent, so every run rebuilt and then failed). `core/cache.py:143-182` (`atomic_cache_dir`) already handles both: a pool another run installed first is kept when it carries a manifest, and a directory without one is replaced. The pool now writes that manifest. The reader keeps `manifest_required=False`, so pools cached by earlier releases (which were already installed atomically) are still used.

**Files:**
- Modify: `src/tessera/discover/nextclade.py` (line 29 import; lines 256-258 local imports; lines 276-355 the build and install, at `4b35e39`)
- Modify: `tests/unit/test_nextclade.py` (append)
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes: `atomic_cache_dir`, `write_cache_manifest` from `tessera.core.cache` (existing).
- Produces: `build_pool` signature unchanged. The installed pool directory holds `tessera_cache.json` with keys `dataset`, `tag`, `per_clade_consensus`, `genomes`, `tessera_version`, `written (UTC)`.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_nextclade.py`:

```python


# --- installing the pool into the cache ------------------------------------------
def test_build_pool_installs_over_a_directory_of_stray_files(monkeypatch, tmp_path, logger):
    # A cache directory holding only files the reader does not recognise reads as
    # absent, so the pool is rebuilt -- and a bare rename onto that non-empty directory
    # failed with an OSError on every run.
    _patch_files(monkeypatch)
    cache = tmp_path / "c"
    cache.mkdir()
    (cache / "notes.txt").write_text("not a genome\n")
    genomes = nc.build_pool(_dataset(), cache_dir=cache, logger=logger)
    assert {g.stem for g in genomes} >= {"ACC1", "ACC2", "EX1"}
    assert all(g.parent == cache for g in genomes)


def test_build_pool_keeps_a_pool_a_concurrent_run_installed_first(
    monkeypatch, tmp_path, logger
):
    # Two runs on one dataset: the other run installs its pool while this one is still
    # building. Its copy answers the same request, so this run uses it rather than fail.
    _patch_files(monkeypatch)
    cache = tmp_path / "c"
    real_write = nc._write_genome
    raced = {"done": False}

    def write_then_race(out_dir, label, note, seq, floor):
        if not raced["done"]:
            raced["done"] = True
            nc.build_pool(_dataset(), cache_dir=cache, logger=logger)  # the other run
        return real_write(out_dir, label, note, seq, floor)

    monkeypatch.setattr(nc, "_write_genome", write_then_race)
    genomes = nc.build_pool(_dataset(), cache_dir=cache, logger=logger)
    assert {g.stem for g in genomes} >= {"ACC1", "ACC2", "EX1"}
    assert not [p for p in tmp_path.iterdir() if p.name != "c"]  # no staging dir left


def test_build_pool_marks_the_installed_pool_complete(monkeypatch, tmp_path, logger):
    from tessera.core.cache import read_cache_manifest

    _patch_files(monkeypatch)
    cache = tmp_path / "c"
    nc.build_pool(_dataset(), cache_dir=cache, logger=logger)
    manifest = read_cache_manifest(cache)
    assert manifest is not None
    assert manifest["dataset"] == "x/y" and manifest["tag"] == "t"
    assert manifest["genomes"] == 3
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_nextclade.py -q`
Expected: `3 failed, 15 passed`. Observed: the first two fail with `OSError: [Errno 66] Directory not empty: '.../tmpXXXX' -> '.../c'`; the third with `assert None is not None`.

- [ ] **Step 3: Implement**

In `src/tessera/discover/nextclade.py`, replace

```python
from ..core.cache import cached_genomes
```

with

```python
from ..core.cache import atomic_cache_dir, cached_genomes, write_cache_manifest
```

In `build_pool`, replace

```python
    import os
    import shutil as _shutil
    import tempfile

    from ..recomb.consensus import consensus_sequence
```

with

```python
    from .. import __version__
    from ..recomb.consensus import consensus_sequence
```

replace

```python
    cache_dir.parent.mkdir(parents=True, exist_ok=True)
    build_dir = Path(tempfile.mkdtemp(dir=cache_dir.parent))
    try:
        written: list[Path] = []
```

with

```python
    # Staged and installed by atomic_cache_dir, which also settles a collision: a pool
    # that a concurrent run installed first (it carries a manifest) is kept, and a
    # directory without one -- stray files the reader does not recognise -- is replaced.
    with atomic_cache_dir(cache_dir) as build_dir:
        written: list[Path] = []
```

(the body of the former `try:` block keeps its indentation unchanged) and replace

```python
        if not written:
            raise UserInputError(
                f"Nextclade dataset '{dataset.path}' yielded no usable pool genomes."
            )
        os.replace(build_dir, cache_dir)  # atomic on the same filesystem
    except BaseException:
        _shutil.rmtree(build_dir, ignore_errors=True)
        raise
```

with

```python
        if not written:
            raise UserInputError(
                f"Nextclade dataset '{dataset.path}' yielded no usable pool genomes."
            )
        write_cache_manifest(
            build_dir, dataset=dataset.path, tag=dataset.tag,
            per_clade_consensus=per_clade_consensus, genomes=len(written),
            tessera_version=__version__,
        )
```

`atomic_cache_dir` removes the staging tree on any exception, as the removed `except BaseException` did, so the `UserInputError` path is unchanged.

Append to the "Alignment integrity" group in `CHANGELOG.md`:

```markdown
- **The Nextclade pool cache could not be installed over an existing directory.** Two runs
  building the same pool at once, or a cache directory holding only unrecognised files,
  made the install fail with "Unexpected error". The pool is now installed the way the
  NCBI Virus panel is: a pool another run finished first is kept, a directory without a
  completeness record is replaced, and the pool carries such a record. Pools cached by
  earlier releases are still used.
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_nextclade.py tests/unit/test_input_safety.py tests/unit/test_lineage_assign.py tests/unit/test_reassort_scan.py tests/unit/test_reassort_assign.py -q`
Expected: all pass.

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: clean, clean, `786 passed, 9 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/discover/nextclade.py tests/unit/test_nextclade.py CHANGELOG.md
git commit -m "$(cat <<'EOF'
Install the Nextclade pool cache with atomic_cache_dir

build_pool staged the pool in a temporary directory and installed it with
a bare os.replace, which fails on a non-empty destination. Two runs
building one dataset's pool at once, or a cache directory holding only
files the reader does not recognise (read as absent, so every run rebuilt
and then failed), surfaced as "Unexpected error: Directory not empty".

core/cache.py already has the helper for this, used by the NCBI Virus
panel: a copy another run finished first is kept if it carries a
manifest, and a directory without one is replaced. The pool now goes
through it and writes that manifest (dataset, tag, consensus mode, genome
count, version). The reader still accepts a pool without a manifest, so
caches written by earlier releases, which were installed atomically, are
reused rather than fetched again.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 9: A non-directory at the working-collection path (audit items: C7)

`_clear_working_copy` already refuses a non-directory that Tessera did not create. When a previous run's marker is present, a plain file at `<output>/collection` reached `shutil.rmtree` and raised `NotADirectoryError`, and a symlink raised `OSError: Cannot call rmtree on a symbolic link`. A dangling symlink was missed by `dest.exists()` and failed later in `copytree` with `FileExistsError`. All three now raise `UserInputError`. A symlink to a directory is refused rather than followed: clearing it would empty a directory elsewhere. The spec's description ("calls rmtree without checking that the destination is a directory") is slightly broader than the code: the unmarked case was already a user error.

**Files:**
- Modify: `src/tessera/core/io.py` (lines 148-154 in `_clear_working_copy`, at `4b35e39`)
- Modify: `tests/unit/test_input_safety.py` (append)
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes: nothing.
- Produces: `copy_collection` and `new_working_collection` raise `UserInputError` ("exists and is not a directory") for a file or any symlink at the destination.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_input_safety.py`:

```python


@pytest.mark.parametrize("kind", ["file", "symlink", "dangling"])
def test_a_non_directory_at_the_working_collection_path_is_a_user_error(
    tmp_path: Path, kind: str
) -> None:
    # A previous run's marker is present, so the path counts as Tessera's own; what
    # sits there is not a directory. Clearing it raised a bare NotADirectoryError (an
    # OSError for a symlink, a FileExistsError for a dangling one), reported as
    # "Unexpected error".
    source = _collection(tmp_path / "src")
    out = tmp_path / "out"
    out.mkdir()
    (out / ".tessera-working-collection").write_text("")
    dest = out / "collection"
    if kind == "file":
        dest.write_text("not a directory\n")
    elif kind == "dangling":
        dest.symlink_to(tmp_path / "gone", target_is_directory=True)
    else:
        elsewhere = _collection(tmp_path / "elsewhere", names=("keep",))
        dest.symlink_to(elsewhere, target_is_directory=True)
    with pytest.raises(UserInputError, match=r"collection.*not a directory"):
        copy_collection(source, dest)
    assert dest.exists() or dest.is_symlink()
    if kind == "symlink":
        assert (tmp_path / "elsewhere" / "keep.fasta").exists()
```

- [ ] **Step 2: Run it, expected failure**

Run: `pytest tests/unit/test_input_safety.py -k non_directory -q`
Expected: `3 failed`. Observed: `NotADirectoryError: [Errno 20] Not a directory: PosixPath('.../out/collection')`; `OSError: Cannot call rmtree on a symbolic link`; `FileExistsError: [Errno 17] File exists: '.../out/collection'`.

- [ ] **Step 3: Implement**

In `src/tessera/core/io.py`, in `_clear_working_copy`, replace

```python
    if not dest.exists():
        return
    parent = dest.parent
    ours = (parent / _WORKING_COPY_MARKER).exists() or any(
        (parent / name).exists() for name in _RUN_ARTEFACTS
    )
    if not ours and (not dest.is_dir() or any(dest.iterdir())):
```

with

```python
    if dest.is_symlink() or (dest.exists() and not dest.is_dir()):
        # Only a directory can be cleared; a symlink is refused even when it points at
        # one, as clearing it would empty a directory somewhere else.
        raise UserInputError(
            f"{dest} exists and is not a directory. A run keeps its working copy of the "
            "references there; move it elsewhere or choose a different output "
            "directory."
        )
    if not dest.exists():
        return
    parent = dest.parent
    ours = (parent / _WORKING_COPY_MARKER).exists() or any(
        (parent / name).exists() for name in _RUN_ARTEFACTS
    )
    if not ours and any(dest.iterdir()):
```

(`not dest.is_dir()` in the last condition can no longer be true, so it is dropped.)

Append to the "Alignment integrity" group in `CHANGELOG.md`:

```markdown
- **A file or symlink at `<output>/collection` gave a raw error.** When a previous run's
  marker was present, clearing the working collection failed with `NotADirectoryError`,
  an `OSError` for a symlink, or `FileExistsError` for a dangling one. Each is now a user
  error that names the path; a symlink is refused rather than followed.
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_input_safety.py tests/integration/test_genome_io.py tests/unit/test_iterate.py tests/unit/test_panel.py -q`
Expected: all pass, including the three existing `copy_collection` tests.

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: clean, clean, `789 passed, 9 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/core/io.py tests/unit/test_input_safety.py CHANGELOG.md
git commit -m "$(cat <<'EOF'
Refuse a non-directory at the working-collection path

The working collection is cleared at the start of every run. When a
previous run's marker was present, a plain file at <output>/collection
reached shutil.rmtree and raised NotADirectoryError, a symlink raised
"Cannot call rmtree on a symbolic link", and a dangling symlink slipped
past exists() and failed in copytree with FileExistsError. Each reached
the user as "Unexpected error".

All three now raise a user error naming the path. A symlink to a directory
is refused rather than followed: clearing through it would empty a
directory elsewhere. The case without a marker was already refused.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 10: Child processes get no standard input; replaced log handlers are closed (audit items: "Suspected and minor": stdin and the logging handler)

Child processes inherited the parent's standard input (`process.py:96,116`, `binaries.py:47`). A tool that reads stdin, or prompts, waits on the terminal; alignment has no default timeout, so the run would hang rather than fail. No current call site is known to trigger it, so the test makes the condition explicit: standard input is replaced with a pipe that never reaches end of file, as a terminal does. Separately, `configure_logging` cleared handlers without closing them, leaving the previous run log open in any process that configures logging more than once. The two are folded because both are process hygiene with no effect on results.

**Files:**
- Modify: `src/tessera/core/process.py` (lines 3-6 docstring; lines 96-99 and 116-118 `subprocess.run`, at `4b35e39`)
- Modify: `src/tessera/core/binaries.py` (lines 47-52 `subprocess.run`)
- Modify: `src/tessera/core/logging.py` (lines 19-24)
- Modify: `tests/unit/test_process.py`, `tests/unit/test_binaries.py` (append)
- Create: `tests/unit/test_logging.py`
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes: nothing.
- Produces: every `subprocess.run` in `core/` passes `stdin=subprocess.DEVNULL`; `configure_logging` closes the handlers it removes. No signature change.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_process.py`:

```python


# --- standard input -----------------------------------------------------------------

@pytest.fixture
def stdin_left_open():
    """Make this process's standard input a pipe that never reaches end of file, as an
    interactive terminal does. A child that inherits it and reads it waits forever."""
    import os

    read_end, write_end = os.pipe()
    saved = os.dup(0)
    os.dup2(read_end, 0)
    try:
        yield
    finally:
        os.dup2(saved, 0)
        for fd in (saved, read_end, write_end):
            os.close(fd)


@pytest.mark.parametrize("to_file", [False, True])
def test_a_tool_that_reads_stdin_gets_end_of_file(
    stdin_left_open, tmp_path: Path, to_file: bool
) -> None:
    # Without stdin redirected the child inherits the terminal and waits on it; with no
    # timeout on alignment, the run would hang rather than fail.
    cmd = _python("import sys", "print(len(sys.stdin.read()))")
    out = tmp_path / "out.txt"
    output = run_tool(CAPS, cmd, logger=LOG, timeout=5, stdout_path=out if to_file else None)
    assert (out.read_text() if to_file else output).strip() == "0"
```

Append to `tests/unit/test_binaries.py`:

```python


def test_version_probe_does_not_inherit_stdin(monkeypatch, on_path) -> None:
    # A tool that prompts on its version option would otherwise wait on the terminal
    # for the full 30 s probe timeout.
    import subprocess

    from tessera.core import binaries

    fake_binary(on_path, "prompty", "prompty 1.2.3")
    seen: dict = {}
    real_run = subprocess.run

    def spy(*args, **kwargs):
        seen.update(kwargs)
        return real_run(*args, **kwargs)

    monkeypatch.setattr(binaries.subprocess, "run", spy)
    assert check_binaries((BinarySpec("prompty"),)) == {"prompty": "1.2.3"}
    assert seen.get("stdin") is subprocess.DEVNULL
```

Create `tests/unit/test_logging.py`:

```python
"""Logging configuration: a reconfigured run log releases the file it held."""

from __future__ import annotations

import logging
from pathlib import Path

from tessera.core.logging import configure_logging


def _run_log_handler(logger: logging.Logger) -> logging.FileHandler:
    (handler,) = [h for h in logger.handlers if isinstance(h, logging.FileHandler)]
    return handler


def test_reconfiguring_closes_the_previous_run_log(tmp_path: Path) -> None:
    # A process that configures logging more than once (the test suite, a notebook)
    # replaced the file handler without closing it, leaving one open file per call.
    first = _run_log_handler(configure_logging(tmp_path / "a"))
    try:
        configure_logging(tmp_path / "b")
        assert first.stream is None  # FileHandler.close() releases the stream
    finally:
        configure_logging(None)
    assert (tmp_path / "a" / "tessera.log").exists()
```

- [ ] **Step 2: Run them, expected failure**

Run: `pytest tests/unit/test_process.py tests/unit/test_binaries.py tests/unit/test_logging.py -q`
Expected: `4 failed`, after about 10 s (two 5 s timeouts). Observed: both `test_a_tool_that_reads_stdin_gets_end_of_file` cases fail with `ToolExecutionError: command failed (exit -1)` (the timeout); `test_version_probe_does_not_inherit_stdin` with `assert None is -3` (`-3` is `subprocess.DEVNULL`); `test_reconfiguring_closes_the_previous_run_log` with `assert <_io.TextIOWrapper name='.../a/tessera.log' ...> is None`.

- [ ] **Step 3: Implement**

In `src/tessera/core/process.py`, replace the module docstring's first paragraph

```python
Every adapter runs its binary through :func:`run_tool` rather than calling
``subprocess`` directly. Commands are passed as argument lists and executed
without a shell, so user-supplied paths and filenames cannot inject shell
metacharacters. A non-zero exit raises :class:`ToolExecutionError`.
```

with

```python
Every adapter runs its binary through :func:`run_tool` rather than calling
``subprocess`` directly. Commands are passed as argument lists and executed
without a shell, so user-supplied paths and filenames cannot inject shell
metacharacters. A non-zero exit raises :class:`ToolExecutionError`. Standard input
is ``/dev/null``: no tool here is fed on stdin, and one that waited on an inherited
terminal would hang the run, as alignment has no default timeout.
```

replace

```python
                proc = subprocess.run(
                    cmd, stdout=out, stderr=subprocess.PIPE, text=True,
                    timeout=timeout, env=env,
                )
```

with

```python
                proc = subprocess.run(
                    cmd, stdin=subprocess.DEVNULL, stdout=out, stderr=subprocess.PIPE,
                    text=True, timeout=timeout, env=env,
                )
```

and replace

```python
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=timeout, env=env
            )
```

with

```python
            proc = subprocess.run(
                cmd, stdin=subprocess.DEVNULL, capture_output=True, text=True,
                timeout=timeout, env=env,
            )
```

In `src/tessera/core/binaries.py`, replace

```python
            [name, *version_args],
            capture_output=True,
```

with

```python
            [name, *version_args],
            stdin=subprocess.DEVNULL,
            capture_output=True,
```

In `src/tessera/core/logging.py`, replace

```python
    Idempotent: re-configuring replaces existing handlers so repeated calls (e.g.
    in tests) do not duplicate log lines.
    """
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
```

with

```python
    Idempotent: re-configuring replaces existing handlers so repeated calls (e.g.
    in tests) do not duplicate log lines. The replaced handlers are closed, so a
    previous run log does not stay open.
    """
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(logging.DEBUG)
    for handler in logger.handlers:
        handler.close()
    logger.handlers.clear()
```

(`StreamHandler.close()` does not close `sys.stderr`.)

Append to the "Alignment integrity" group in `CHANGELOG.md`:

```markdown
- **External tools inherited standard input.** A tool that read stdin or prompted would
  have waited on the terminal, and alignment has no default timeout. Every external tool,
  and every version probe, now runs with standard input closed. No current tool is known
  to have done this. Reconfiguring logging also now closes the previous run log.
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_process.py tests/unit/test_binaries.py tests/unit/test_logging.py tests/unit/test_cli_commands.py -q`
Expected: all pass.

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: clean, clean, `793 passed, 9 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/core/process.py src/tessera/core/binaries.py src/tessera/core/logging.py \
        tests/unit/test_process.py tests/unit/test_binaries.py tests/unit/test_logging.py \
        CHANGELOG.md
git commit -m "$(cat <<'EOF'
Close stdin for external tools and close replaced log handlers

run_tool and the version probe started child processes without
redirecting standard input, so a tool that read stdin or prompted would
wait on the inherited terminal; alignment deliberately has no default
timeout, so the run would hang rather than fail. No current call site is
known to trigger it. The test makes the condition explicit by giving the
test process a stdin pipe that never reaches end of file: before the
change the child waited until the 5 s test timeout, after it the child
reads end of file at once.

configure_logging cleared its handlers without closing them, so a process
that configured logging more than once kept every earlier run log open.
The handlers are now closed before they are replaced.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 11: Cactus sample-name collisions; remove the unreachable XMFA flank path (audit items: "Suspected and minor": cactus sample names, `xmfa_to_fasta` flank)

Two small changes to code that cannot be run here (cactus and progressiveMauve are not installed), folded because both are adapter-converter tidying a reviewer would take or leave together.

1. `cactus._sample_name` maps '.' to '_', so labels `MN908947.3` and `MN908947_3` (distinct, so staging accepts them) become one Cactus sample name, one seqFile line and one key in the converter's name map. The spec marks this SUSPECTED; what cactus does with a repeated sample name was not tested. The guard refuses the input before cactus is started, which needs no binary to test. It sits in the cactus adapter rather than in staging: the collision exists only for this backend, and core does not know about adapters.
2. `xmfa_to_fasta(..., flank, ...)` has a deletion-flank branch that no caller reaches: the only caller (`progressivemauve.py:82`) and the only test (`test_input_safety.py:368`) pass `0`. The parameter and the branch are removed. This is dead-code removal: no test can fail before and pass after; the existing XMFA tests and the Task 2 invariant are the guard.

**Files:**
- Modify: `src/tessera/aligners/cactus.py` (line 19 import; line 46 in `align`; new function above `_find_hal` at line 100, at `4b35e39`)
- Modify: `src/tessera/converters/xmfa_to_fasta.py` (lines 1-9 docstring; lines 34-45 signature and docstring; lines 144-186 and 193-201 the flank branches)
- Modify: `src/tessera/aligners/progressivemauve.py` (the `xmfa_to_fasta` call)
- Modify: `tests/unit/test_backbone_invariant.py` (append)
- Modify: `tests/unit/test_input_safety.py` line 368 (existing test changes: drops the `0` argument)
- Modify: `docs/aligners.md` (backend table, `cactus` row), `CHANGELOG.md`

**Interfaces:**
- Consumes: `tests/unit/test_backbone_invariant.py` helpers from Task 2 (`_LOG`, `cactus_mod`, `AlignParams`, `pytest`).
- Produces: `cactus._reject_colliding_sample_names(genomes: Sequence[Path]) -> None`; `xmfa_to_fasta(xmfa_path: str | Path, reference_name: str, out_path: str | Path, reference_length: int | None = None) -> Path` (the `flank` positional argument is gone).

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_backbone_invariant.py`:

```python


def test_cactus_refuses_labels_that_collide_once_dots_are_replaced(
    monkeypatch, tmp_path: Path
) -> None:
    # Cactus sample names may not contain '.', so the adapter maps it to '_'. Two
    # accessions differing only there would then share one sample name and one row.
    from tessera.core.errors import UserInputError

    def must_not_run(caps, cmd, **kw):
        raise AssertionError("cactus-pangenome must not be started")

    monkeypatch.setattr(cactus_mod, "run_tool", must_not_run)
    genomes = []
    for label in ("ref", "MN908947.3", "MN908947_3"):
        path = tmp_path / f"{label}.fasta"
        path.write_text(f">{label}\nACGT\n")
        genomes.append(path)
    with pytest.raises(UserInputError, match=r"MN908947\.3.*MN908947_3.*MN908947_3"):
        cactus_mod.CactusAligner().align(
            genomes, genomes[0], tmp_path / "out", AlignParams(threads=1), _LOG
        )
```

In `tests/unit/test_input_safety.py`, in `test_xmfa_without_the_reference_is_reported_against_the_file`, replace

```python
        xmfa_to_fasta(xmfa, "/data/ref.fasta", 0, tmp_path / "out.fasta")
```

with

```python
        xmfa_to_fasta(xmfa, "/data/ref.fasta", tmp_path / "out.fasta")
```

- [ ] **Step 2: Run them, expected failure**

Run: `pytest tests/unit/test_backbone_invariant.py tests/unit/test_input_safety.py -q`
Expected: `2 failed`. Observed: the cactus test fails with `AssertionError: cactus-pangenome must not be started`; `test_xmfa_without_the_reference_is_reported_against_the_file` fails with `TypeError: xmfa_to_fasta() missing 1 required positional argument: 'out_path'` (a signature-only failure, expected until Step 3).

- [ ] **Step 3: Implement**

In `src/tessera/aligners/cactus.py`, replace

```python
from ..core.errors import OutputError
```

with

```python
from ..core.errors import OutputError, UserInputError
```

replace

```python
        genomes, reference = normalize_reference(genomes, reference, tool="Cactus")
        out_dir.mkdir(parents=True, exist_ok=True)
```

with

```python
        genomes, reference = normalize_reference(genomes, reference, tool="Cactus")
        _reject_colliding_sample_names(genomes)
        out_dir.mkdir(parents=True, exist_ok=True)
```

and insert above `def _find_hal(results: Path) -> Path | None:`:

```python
def _reject_colliding_sample_names(genomes: Sequence[Path]) -> None:
    """Refuse two genomes whose labels become one Cactus sample name.

    Staging already refuses two inputs with one label, but ``_sample_name`` maps '.' to
    '_', so ``MN908947.3`` and ``MN908947_3`` -- distinct labels -- would share a
    sample name, a seqFile line and an MSA row.
    """
    seen: dict[str, str] = {}
    for genome in genomes:
        sample = _sample_name(genome)
        if sample in seen and seen[sample] != genome.stem:
            raise UserInputError(
                f"Genomes '{seen[sample]}' and '{genome.stem}' both become the Cactus "
                f"sample name '{sample}' (Cactus does not allow '.' in sample names). "
                "Rename one, or use another --aligner."
            )
        seen[sample] = genome.stem


```

(`Sequence` is already imported in `cactus.py`.)

In `src/tessera/converters/xmfa_to_fasta.py`, make these replacements:

```diff
 subprocess by the legacy ``phylo.py``. Behaviour is preserved: the reference
-sequence defines the coordinate system, reference gaps are removed, optional
-deletion flanks are screened, and every sequence is projected onto the reference
-coordinates. Strings are handled as ``bytearray`` / ``bytes`` throughout (the
-original relied on Python 2 ``str``/``bytearray`` equivalence); regex scanning
-uses a latin-1 view so byte positions are 1:1 with characters.
+sequence defines the coordinate system, reference gaps are removed, and every
+sequence is projected onto the reference coordinates. The original's optional
+deletion-flank screen is not carried over: no caller used it. Strings are handled
+as ``bytearray`` / ``bytes`` throughout (the original relied on Python 2
+``str``/``bytearray`` equivalence); regex scanning uses a latin-1 view so byte
+positions are 1:1 with characters.
```

```diff
 def xmfa_to_fasta(
     xmfa_path: str | Path,
     reference_name: str,
-    flank: int,
     out_path: str | Path,
     reference_length: int | None = None,
 ) -> Path:
     """Convert a progressiveMauve XMFA file to a reference-anchored FASTA.
 
     ``reference_name`` must match a ``#SequenceNFile`` entry in the XMFA.
-    ``flank`` (>=0) screens deletion flanks by that many bases.
```

```diff
         # Find reference gaps and remove those columns from every sequence.
         search_pos = 0
-        list_of_gaps: list[list[int]] = []
         ref_view = a_gen[alignment][reference_num]["seq"].decode("latin-1")
```

```diff
                 for sequence in a_gen[alignment]:
                     a_gen[alignment][sequence]["seq"][start:end] = bytearray()
-                if flank > 0:
-                    for g in list_of_gaps:
-                        g[0] = g[0] - end + start
-                        g[1] = g[1] - end + start
-                    list_of_gaps.append([start, start])
-
-        if flank > 0:
-            search_pos = 0
-            for sequence in a_gen[alignment]:
-                if sequence == reference_num:
-                    continue
-                seq_view = a_gen[alignment][sequence]["seq"].decode("latin-1")
-                while search_pos < length_of_reference:
-                    gap_hit = _PATTERN_GAP.search(seq_view, search_pos)
-                    if gap_hit:
-                        list_of_gaps.append([gap_hit.start(), gap_hit.end()])
-                        search_pos = gap_hit.end()
-                    else:
-                        break
-            for non_ref_gap in list_of_gaps:
-                for sequence in a_gen[alignment]:
-                    seq = a_gen[alignment][sequence]["seq"]
-                    new_start = max(0, non_ref_gap[0] - flank)
-                    new_end = min(non_ref_gap[1] + flank, len(seq))
-                    seq[new_start:new_end] = bytearray(b"-" * (new_end - new_start))
 
     # Project each alignment block onto reference coordinates.
```

```diff
         end = int(ref["p2"])
-        if flank > 0:
-            start += flank
-            end -= flank
         forward = ref["sign"] == "+"
@@
-            fragment = bytes(seq[flank:-flank]) if flank > 0 else bytes(seq)
+            fragment = bytes(seq)
```

In `src/tessera/aligners/progressivemauve.py`, replace

```python
            xmfa_to_fasta(xmfa, ref_arg, 0, fa, reference_length=ref_length)
```

with

```python
            xmfa_to_fasta(xmfa, ref_arg, fa, reference_length=ref_length)
```

In `docs/aligners.md`, replace the table row

```markdown
| `cactus` | Same-species pangenomes | Resource heavy (Toil/containers) |
```

with

```markdown
| `cactus` | Same-species pangenomes | Resource heavy (Toil/containers); two genome labels that differ only in `.` and `_` are refused, as Cactus sample names cannot contain `.` |
```

Append to the "Alignment integrity" group in `CHANGELOG.md`:

```markdown
- **`cactus`: two labels differing only in `.` and `_` would share one row.** Cactus
  sample names cannot contain `.`, and the adapter maps it to `_`; such inputs are now
  refused before Cactus starts. Checked without the binary, which is not installed here.
```

- [ ] **Step 4: Run, expected pass**

Run: `pytest tests/unit/test_backbone_invariant.py tests/unit/test_input_safety.py tests/unit/test_aligner_backends.py -q && grep -n "flank" src/tessera/converters/xmfa_to_fasta.py`
Expected: all pass; the grep finds only the docstring sentence.

- [ ] **Step 5: Full checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: clean, clean, `794 passed, 9 deselected`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/aligners/cactus.py src/tessera/converters/xmfa_to_fasta.py \
        src/tessera/aligners/progressivemauve.py tests/unit/test_backbone_invariant.py \
        tests/unit/test_input_safety.py docs/aligners.md CHANGELOG.md
git commit -m "$(cat <<'EOF'
Refuse colliding Cactus sample names; drop the unused XMFA flank path

The cactus adapter maps '.' to '_' in sample names, as Cactus reserves '.'
for haplotype suffixes. Staging accepts MN908947.3 and MN908947_3 as two
labels, but they would become one sample name, one seqFile line and one
key in the converter's name map. What cactus does with that was not
tested (it is not installed here); the adapter now refuses such input
before starting it, which needs no binary to check.

xmfa_to_fasta carried a deletion-flank branch from the Python 2 original
that no caller reached: the adapter and the one test both passed 0. The
parameter and the branch are removed. No behaviour changes; the existing
XMFA tests and the backbone invariant are the guard, and the one test
that passed the 0 drops it.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 12: Final checks and pull request

**Files:** none changed, unless a check fails.

**Interfaces:** none.

- [ ] **Step 1: Full suite and static checks**

Run: `ruff check src tests validation && mypy && pytest -m "not requires_binary" -q`
Expected: `All checks passed!`, `Success: no issues found in 79 source files`, `794 passed, 9 deselected` (from 757 and 5).

- [ ] **Step 2: Binary-gated tests**

```bash
export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
pytest -m requires_binary -q -rs
```

Expected where the env holds mafft, minimap2 and sibeliaz: `8 passed, 1 skipped` (progressiveMauve not on PATH). Observed on the candidate while writing this plan: exactly that, on macOS. In CI only mafft is installed, so the minimap2 and sibeliaz cases skip there.

- [ ] **Step 3: Specificity harness, unchanged output**

```bash
mkdir -p validation/data/audit-d/final
python validation/run_specificity.py --reps 3 > validation/data/audit-d/final/spec.txt 2>&1
diff validation/data/audit-d/baseline/spec.txt validation/data/audit-d/final/spec.txt
```

Expected: no difference. Observed while writing this plan, base and candidate alike: `TOTAL 0/12 (0%, CI 0%-24%)` false regions at `--min-methods 2`, and positive control `detected 3/3 | correct donor 3/3 | median breakpoint error 13 bp`. This harness feeds already-aligned data, so of this plan's changes only the Task 6 warning split runs in it. Any difference is a regression: find the task with `git bisect` and fix or revert it.

- [ ] **Step 4: Pull request**

Push and open the PR only when the maintainer asks (`fix-audit-caller-behaviour` is itself unmerged). The description lists, per task, the audit item and the test that fails without it; states that progressiveMauve and cactus are checked only against fakes; states which binary-gated cases ran where (Step 2); and gives the Step 3 result. It ends with the Claude Code attribution line the session's instructions give for pull requests (the robot emoji followed by `Generated with [Claude Code](https://claude.com/claude-code)`).

## Self-review notes

**Coverage of the spec.** C1 (Task 2), C2 (Task 4), C3 (Task 5), C4 (Task 8), C5 (Task 6), C6 (Task 7), C7 (Task 9); minors: stdin (Task 10), logging handler (Task 10), cactus sample names (Task 11), XMFA flank (Task 11). Core-side test gaps from section 3: backbone row per backend (Tasks 2 and 3), two records sharing an ID (Task 4), the minimap2 worker count (Task 5), the short-MAF-row branch (Task 7), a source size disagreeing with the declared contig length (Task 7), the pool install under concurrency and with stray files (Task 8), a non-directory at the working collection path (Task 9).

**Where the spec is not accurate.**
- "No external aligner binary is installed": on the development Mac the `recomfi-aln` env holds mafft, minimap2 and sibeliaz 1.2.7; progressiveMauve and cactus are absent. CI installs mafft.
- G3 asks for "the concatenated backbone row equals the reference sequence" for every backend. That holds for the pairwise backends. A MAF-based backend (sibeliaz, cactus) lays the backbone out from blocks, so a position no block covers is a gap in the backbone row by design; with real SibeliaZ the last 17 bp of a 3 kb contig were uncovered. No other genome has a base in such a column, so the scan cannot compare anything there. The binary test therefore checks "no misplaced base" for MAF backends; the unit fakes cover every column, so equality holds there.
- C7: the unmarked case (a file at the path with no marker) was already a `UserInputError`; the raw errors occur only with the marker present, and for a dangling symlink.
- C3: the per-process memory of the minimap2 index was not measured; for viral backbones it is small, and the change is about threads.

**Found while writing the plan, outside the spec.** The default `sibeliaz` backend wrote an empty MAF on macOS for every input tried, because of two GNU-only wrapper commands (Task 3). `validation/README.md` and the project notes record "sibeliaz OOMs on the orthopox example on a laptop"; whether that was this incompatibility rather than memory was not tested and is left as an open question. Suggested check after Task 3, by the maintainer: run the orthopox example with `--aligner sibeliaz` on macOS, one run at a time.

**Deferred, with reasons.**
- `collection_genomes` skips genomes in subdirectories. Fixing it needs a decision (warn, recurse, or refuse), and the function has no logger; the spec marks it minor and optional.
- Two unnamed records (`>` with no ID) in one genome are not refused by Task 4: the sibeliaz adapter already refuses an unnamed record, and the message for an empty ID would be unclear. Under minimap2 two unnamed backbone records would share one offset key; not reproduced.
- Alignment still has no default timeout (`process.py` documents why). With stdin closed (Task 10), the hang the spec describes no longer depends on it.
- The MAF converter lays cactus contigs out in sorted-name order; recorded in the previous audit and in `docs/aligners.md`, and not changed here.

**Existing tests that change.** `tests/unit/test_aligner_params.py:50` (Task 6): the `_concatenate` stub returns `[]` instead of `None`, because the return value is now used. `tests/unit/test_input_safety.py:368` (Task 11): the `xmfa_to_fasta` call drops the `0` flank argument.

**Verification done while writing.** Every new test was run in a detached worktree of `4b35e39` against the unmodified source (28 of the new fast tests failed for the reasons recorded in each Step 2; the eight non-progressiveMauve invariant cases and the one-ID-in-two-genomes test passed, as guards should) and against the cumulative candidate (`794 passed, 9 deselected`; ruff and mypy clean; `requires_binary` 8 passed, 1 skipped on macOS). The tasks were applied cumulatively, not one commit at a time; the per-task counts above are derived from the number of tests each task adds.
