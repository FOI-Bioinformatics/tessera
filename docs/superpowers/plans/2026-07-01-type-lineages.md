# `tessera type-lineages` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A `tessera type-lineages` command that assigns a lineage to every genome in an NCBI/local collection (metadata -> nearest-reference -> de-novo clustering) and writes the shared `lineages.tsv`, making any collection lineage-ready.

**Architecture:** A new `discover/lineage_assign.py` holds the three assignment methods and the orchestration; a shared `typing.is_recombinant_lineage` flags recombinant names; a thin `cli/cmd_type_lineages.py` wraps it. Heavy reuse of existing typing, Nextclade, and skani helpers; the only external tool is the already-optional `skani`.

**Tech Stack:** Python 3.11+, Typer CLI, `skani` (optional external binary, already used), pytest.

## Global Constraints

- Dependency-free by design: no new runtime dependency. `skani` is already an optional external tool; native reference placement (no `nextclade` binary required).
- Reuse over reimplement: `build_lineage_map`/`titles_from_collection`/`write_lineage_map` (`recomb/typing.py`), `resolve_dataset`/`build_pool`/`first_header` (`discover/nextclade.py` + `recomb/typing.py`), `skani_query_ani`/`dereplicate` (`discover/panel.py`).
- Command-only scope: assign + write `lineages.tsv`. No change to panel selection, detection, or the report this cycle.
- Modest scientific language; faithful reporting. Commit messages end with the co-author trailer; branch before committing on the default branch.
- CI unit tests must not need external binaries (stub `skani_query_ani`); the end-to-end command run is opt-in.

---

### Task 1: shared `is_recombinant_lineage`

**Files:**
- Modify: `src/tessera/recomb/typing.py`
- Modify: `validation/run_hybrids.py` (re-point its local helper)
- Test: `tests/unit/test_typing_recombinant.py` (create)

**Interfaces:**
- Produces: `is_recombinant_lineage(name: str) -> bool` in `recomb/typing.py`.

- [ ] **Step 1: Write the failing test**

```python
# tests/unit/test_typing_recombinant.py
"""The shared recombinant-lineage name test."""

from __future__ import annotations

from tessera.recomb.typing import is_recombinant_lineage


def test_recombinant_names_detected():
    assert is_recombinant_lineage("CRF01_AE")
    assert is_recombinant_lineage("URF1")
    assert is_recombinant_lineage("XBB")          # Pango recombinant (X + letter)
    assert is_recombinant_lineage("recombinant")  # literal, case-insensitive
    assert is_recombinant_lineage("RECOMBINANT")


def test_non_recombinant_names_pass():
    assert not is_recombinant_lineage("A.1")
    assert not is_recombinant_lineage("B.1.1.7")
    assert not is_recombinant_lineage("A.D.1.8")
    assert not is_recombinant_lineage("clade 9")
    assert not is_recombinant_lineage("")
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/unit/test_typing_recombinant.py -q`
Expected: FAIL -- `ImportError: cannot import name 'is_recombinant_lineage'`.

- [ ] **Step 3: Implement in `recomb/typing.py`**

Add near the top-level helpers (after the imports / `_LINEAGE_TOKEN` block):

```python
# Recombinant lineage names (HIV CRF/URF, Pango X-lineages): a clean parental panel
# excludes them, since they carry both parents' segments and mask the true parents.
_RECOMBINANT_LINEAGE = re.compile(r"^(CRF|URF)\d|^X[A-Z]|recombinant", re.IGNORECASE)


def is_recombinant_lineage(name: str) -> bool:
    """True when a lineage name denotes a recombinant (CRF/URF, a Pango X-lineage,
    or the literal 'recombinant')."""
    return bool(_RECOMBINANT_LINEAGE.search(name))
```

(`import re` is already present in `typing.py`; confirm and do not duplicate it.)

- [ ] **Step 4: Re-point the harness to the shared helper**

In `validation/run_hybrids.py`, delete the local `_RECOMBINANT_CLADE` regex and `is_recombinant_clade` definition (`:187-192`) and replace with a re-export so existing call sites keep working:

```python
from tessera.recomb.typing import is_recombinant_lineage as is_recombinant_clade
```

- [ ] **Step 5: Run tests + lint**

Run: `python -m pytest tests/unit/test_typing_recombinant.py -q` (expect 2 passed) and
`python -m pytest -m "not requires_binary" -q` (expect all pass -- the harness re-point is import-safe) and
`ruff check src/tessera/recomb/typing.py validation/run_hybrids.py tests/unit/test_typing_recombinant.py`.

- [ ] **Step 6: Commit**

```bash
git add src/tessera/recomb/typing.py validation/run_hybrids.py tests/unit/test_typing_recombinant.py
git commit -m "Add shared is_recombinant_lineage; re-point the harness to it

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: de-novo clustering (`cluster_by_ani` + `_ani_edges`)

**Files:**
- Create: `src/tessera/discover/lineage_assign.py`
- Test: `tests/unit/test_lineage_assign.py` (create)

**Interfaces:**
- Consumes: `skani_query_ani(query_fasta, refs, logger) -> dict[Path, tuple[float, float]]` (`discover/panel.py`).
- Produces: `cluster_by_ani(genomes: list[Path], edges: list[tuple[Path, Path]]) -> dict[str, list[Path]]`; `_ani_edges(genomes: list[Path], *, threshold: float, logger) -> list[tuple[Path, Path]]`; module constants `DEFAULT_REF_ANI_FLOOR = 90.0`, `DEFAULT_CLUSTER_ANI = 95.0`.

- [ ] **Step 1: Write the failing test (pure clustering)**

```python
# tests/unit/test_lineage_assign.py
"""Lineage assignment: de-novo clustering (pure) and nearest-reference (stubbed skani)."""

from __future__ import annotations

from pathlib import Path

from tessera.discover import lineage_assign as la


def _g(name: str) -> Path:
    return Path(f"/pool/{name}.fasta")


def test_cluster_by_ani_connected_components():
    a, b, c, d = _g("a"), _g("b"), _g("c"), _g("d")
    # a-b-c linked (transitively one lineage); d isolated (its own)
    clusters = la.cluster_by_ani([a, b, c, d], [(a, b), (b, c)])
    groups = sorted(sorted(p.name for p in members) for members in clusters.values())
    assert groups == [["a.fasta", "b.fasta", "c.fasta"], ["d.fasta"]]
    assert all(k.startswith("denovo_") for k in clusters)


def test_cluster_by_ani_all_singletons_when_no_edges():
    a, b = _g("a"), _g("b")
    clusters = la.cluster_by_ani([a, b], [])
    assert sorted(len(v) for v in clusters.values()) == [1, 1]
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/unit/test_lineage_assign.py -q`
Expected: FAIL -- `ModuleNotFoundError: No module named 'tessera.discover.lineage_assign'`.

- [ ] **Step 3: Implement the module (clustering + edges)**

```python
# src/tessera/discover/lineage_assign.py
"""Assign lineage labels to a genome collection (the ``type-lineages`` command).

Three methods, most-authoritative first: (1) metadata/header mining, (2) nearest
neighbour to a Nextclade dataset's clade-labelled tips, (3) de-novo ANI clustering.
Produces the shared ``lineages.tsv`` so any collection becomes lineage-ready. The only
external tool is the already-optional ``skani``.
"""

from __future__ import annotations

import logging
from pathlib import Path

from ..core.io import strip_sequence_extension
from .panel import skani_query_ani

DEFAULT_REF_ANI_FLOOR = 90.0   # a genome joins a reference clade only above this ANI
DEFAULT_CLUSTER_ANI = 95.0     # de-novo: genomes at/above this ANI share a lineage


def cluster_by_ani(
    genomes: list[Path], edges: list[tuple[Path, Path]]
) -> dict[str, list[Path]]:
    """Connected components of the ANI graph -> de-novo lineages.

    ``edges`` are the genome pairs whose ANI cleared the clustering threshold. Every
    genome is a node; each connected component becomes one ``denovo_N`` lineage,
    numbered in order of first appearance. Returns ``{lineage_label: [genome, ...]}``.
    """
    parent: dict[Path, Path] = {g: g for g in genomes}

    def find(x: Path) -> Path:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        if a in parent and b in parent:
            parent[find(a)] = find(b)

    components: dict[Path, list[Path]] = {}
    for g in genomes:
        components.setdefault(find(g), []).append(g)
    return {f"denovo_{i}": members for i, members in enumerate(components.values(), start=1)}


def _ani_edges(
    genomes: list[Path], *, threshold: float, logger: logging.Logger
) -> list[tuple[Path, Path]]:
    """All-vs-all genome pairs with ANI >= ``threshold`` (pairwise ``skani_query_ani``;
    one skani call per genome vs the rest -- fine for a one-off typing pass)."""
    edges: list[tuple[Path, Path]] = []
    for i, g in enumerate(genomes):
        others = genomes[:i] + genomes[i + 1:]
        if not others:
            continue
        ani = skani_query_ani(g, others, logger)
        edges += [(g, other) for other, (a, _af) in ani.items() if a >= threshold]
    return edges
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `python -m pytest tests/unit/test_lineage_assign.py -q`
Expected: PASS (2 passed).

- [ ] **Step 5: Lint and commit**

```bash
ruff check src/tessera/discover/lineage_assign.py tests/unit/test_lineage_assign.py
git add src/tessera/discover/lineage_assign.py tests/unit/test_lineage_assign.py
git commit -m "Add de-novo ANI clustering for lineage assignment

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: nearest-reference assignment

**Files:**
- Modify: `src/tessera/discover/lineage_assign.py`
- Test: `tests/unit/test_lineage_assign.py` (extend)

**Interfaces:**
- Consumes: `skani_query_ani` (as in Task 2); `strip_sequence_extension`.
- Produces: `assign_by_nearest_reference(genomes: list[Path], labeled_tips: dict[Path, str], *, ani_floor: float, logger) -> dict[str, str]` (returns `{genome_label: clade}`, omitting genomes below the floor).

- [ ] **Step 1: Write the failing test (stubbed skani)**

Add to `tests/unit/test_lineage_assign.py`:

```python
import logging


def test_assign_by_nearest_reference(monkeypatch):
    g1, g2 = _g("q1"), _g("q2")
    tipA, tipB = _g("refA"), _g("refB")
    labeled = {tipA: "A.1", tipB: "B.2"}

    # q1 is closest to tipA (95%) and clears the floor; q2's best (85%) is below it.
    fake = {
        g1: {tipA: (95.0, 90.0), tipB: (80.0, 88.0)},
        g2: {tipA: (85.0, 90.0), tipB: (82.0, 80.0)},
    }
    monkeypatch.setattr(la, "skani_query_ani", lambda q, refs, logger: fake[q])

    out = la.assign_by_nearest_reference([g1, g2], labeled, ani_floor=90.0,
                                         logger=logging.getLogger("t"))
    assert out == {"q1": "A.1"}          # q1 assigned; q2 below floor -> omitted
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/unit/test_lineage_assign.py::test_assign_by_nearest_reference -q`
Expected: FAIL -- `AttributeError: ... has no attribute 'assign_by_nearest_reference'`.

- [ ] **Step 3: Implement**

Append to `src/tessera/discover/lineage_assign.py`:

```python
def assign_by_nearest_reference(
    genomes: list[Path], labeled_tips: dict[Path, str], *, ani_floor: float,
    logger: logging.Logger,
) -> dict[str, str]:
    """Assign each genome the clade of its nearest labelled reference tip.

    ``labeled_tips`` maps a reference-tip FASTA path to its clade. For each genome the
    tip with the highest whole-genome ANI wins, but only when it clears ``ani_floor``;
    genomes below the floor are left unassigned (absent from the result). Returns
    ``{genome_label: clade}``.
    """
    tips = list(labeled_tips)
    assigned: dict[str, str] = {}
    if not tips:
        return assigned
    for g in genomes:
        ani = skani_query_ani(g, tips, logger)
        best = max(ani, key=lambda t: ani[t][0], default=None)
        if best is not None and ani[best][0] >= ani_floor:
            assigned[strip_sequence_extension(g.name)] = labeled_tips[best]
    return assigned
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `python -m pytest tests/unit/test_lineage_assign.py -q`
Expected: PASS (3 passed).

- [ ] **Step 5: Lint and commit**

```bash
ruff check src/tessera/discover/lineage_assign.py tests/unit/test_lineage_assign.py
git add src/tessera/discover/lineage_assign.py tests/unit/test_lineage_assign.py
git commit -m "Add nearest-reference clade assignment

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: reference tips + `assign_lineages` orchestration

**Files:**
- Modify: `src/tessera/discover/lineage_assign.py`
- Test: `tests/unit/test_lineage_assign.py` (extend)

**Interfaces:**
- Consumes: `build_lineage_map` / `titles_from_collection` / `first_header` (`recomb/typing.py`); `resolve_dataset` / `build_pool` (`discover/nextclade.py`); `dereplicate` (`discover/panel.py`); `assign_by_nearest_reference`, `_ani_edges`, `cluster_by_ani` (Tasks 2-3).
- Produces: `_reference_tips(*, query, nextclade_dataset, email, cache_dir, logger) -> dict[Path, str]`; `assign_lineages(genomes, *, user_lineage_map=None, taxon=None, nextclade_dataset=None, ref_ani_floor=DEFAULT_REF_ANI_FLOOR, cluster_ani=DEFAULT_CLUSTER_ANI, email=None, cache_dir=None, logger) -> list[tuple[str, str, str]]`.

- [ ] **Step 1: Write the failing test (ladder precedence, sub-steps stubbed)**

Add to `tests/unit/test_lineage_assign.py`:

```python
def test_assign_lineages_ladder_precedence(monkeypatch):
    m1, ref1, dn1, dn2 = _g("m1"), _g("ref1"), _g("dn1"), _g("dn2")
    genomes = [m1, ref1, dn1, dn2]

    # (1) metadata types only m1
    monkeypatch.setattr(la, "build_lineage_map",
                        lambda **kw: [("m1", "B.1", "title")])
    monkeypatch.setattr(la, "titles_from_collection", lambda files: {})
    # (2) a reference exists; nearest-neighbour types ref1
    monkeypatch.setattr(la, "_reference_tips", lambda **kw: {_g("tip"): "A.1"})
    monkeypatch.setattr(la, "assign_by_nearest_reference",
                        lambda genomes, tips, *, ani_floor, logger: {"ref1": "A.1"})
    # (3) de-novo clusters the remaining dn1, dn2 together
    monkeypatch.setattr(la, "_ani_edges", lambda genomes, *, threshold, logger: [(dn1, dn2)])

    rows = la.assign_lineages(genomes, taxon="Test virus", logger=logging.getLogger("t"))
    by_label = {label: (lin, src) for label, lin, src in rows}
    assert by_label["m1"] == ("B.1", "title")
    assert by_label["ref1"] == ("A.1", "nextclade-nn")
    assert by_label["dn1"][1] == "denovo" and by_label["dn2"][1] == "denovo"
    assert by_label["dn1"][0] == by_label["dn2"][0]      # same de-novo lineage
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/unit/test_lineage_assign.py::test_assign_lineages_ladder_precedence -q`
Expected: FAIL -- `AttributeError: ... has no attribute 'assign_lineages'`.

- [ ] **Step 3: Implement**

Append to `src/tessera/discover/lineage_assign.py` (add the imports at the top of the file, not inline):

```python
# add to the module's imports:
from ..recomb.typing import (
    build_lineage_map,
    first_header,
    titles_from_collection,
)
from .nextclade import build_pool, resolve_dataset
from .panel import dereplicate


def _reference_tips(
    *, query: Path, nextclade_dataset: str | None, email: str | None,
    cache_dir: Path | None, logger: logging.Logger,
) -> dict[Path, str]:
    """Resolve a Nextclade dataset and return its clade-labelled tips, dereplicated per
    clade. Empty when no dataset resolves. ``resolve_dataset`` auto-detects the dataset
    from ``query`` (a representative collection genome) unless ``nextclade_dataset`` is
    given. Each tip FASTA header is ``>{accession} {clade}``, so the clade is the
    header's remainder."""
    from ..core.cache import nextclade_cache

    try:
        # resolve_dataset(query, override, *, email, logger): the override (an explicit
        # dataset path) wins; otherwise the dataset is auto-detected from the query genome.
        dataset = resolve_dataset(query, nextclade_dataset, email=email, logger=logger)
    except Exception as exc:  # noqa: BLE001 - no dataset is a normal, non-fatal outcome
        logger.info("No Nextclade dataset for reference typing (%s); skipping.", exc)
        return {}
    cache = cache_dir if cache_dir is not None else nextclade_cache(dataset.path, dataset.tag)
    tips = build_pool(dataset, cache_dir=cache, logger=logger)
    clade_of: dict[Path, str] = {}
    for tip in tips:
        header = first_header(tip)
        parts = header.split(None, 1)
        clade = parts[1].strip() if len(parts) > 1 else ""
        if clade and clade not in ("NA", "example"):
            clade_of[tip] = clade
    # one representative per clade so nearest-neighbour is not dominated by a big clade
    label_by_name = {strip_sequence_extension(p.name): c for p, c in clade_of.items()}
    reps, _ = dereplicate(list(clade_of), clade_of=label_by_name, logger=logger)
    return {p: clade_of[p] for p in reps if p in clade_of}


def assign_lineages(
    genomes: list[Path], *, user_lineage_map: Path | None = None,
    taxon: str | None = None, nextclade_dataset: str | None = None,
    ref_ani_floor: float = DEFAULT_REF_ANI_FLOOR, cluster_ani: float = DEFAULT_CLUSTER_ANI,
    email: str | None = None, cache_dir: Path | None = None, logger: logging.Logger,
) -> list[tuple[str, str, str]]:
    """Assign a lineage to every genome via the metadata -> nearest-reference ->
    de-novo ladder. Returns ``(label, lineage, source)`` rows for ``write_lineage_map``."""
    rows = build_lineage_map(
        user_tsv=user_lineage_map,
        title_by_label=titles_from_collection(genomes),
        organism=taxon,
    )
    typed = {label for label, _lin, _src in rows}
    untyped = [g for g in genomes if strip_sequence_extension(g.name) not in typed]

    if untyped:
        labeled = _reference_tips(query=untyped[0], nextclade_dataset=nextclade_dataset,
                                  email=email, cache_dir=cache_dir, logger=logger)
        if labeled:
            nn = assign_by_nearest_reference(untyped, labeled, ani_floor=ref_ani_floor,
                                             logger=logger)
            rows += [(label, clade, "nextclade-nn") for label, clade in nn.items()]
            typed |= set(nn)
            untyped = [g for g in untyped if strip_sequence_extension(g.name) not in typed]

    if len(untyped) >= 2:
        edges = _ani_edges(untyped, threshold=cluster_ani, logger=logger)
        for lineage, members in cluster_by_ani(untyped, edges).items():
            rows += [(strip_sequence_extension(g.name), lineage, "denovo") for g in members]
    elif untyped:
        rows.append((strip_sequence_extension(untyped[0].name), "denovo_1", "denovo"))

    rows.sort()
    logger.info("Typed %d/%d genome(s) (%d untyped remain).",
                len({r[0] for r in rows}), len(genomes),
                len(genomes) - len({r[0] for r in rows}))
    return rows
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `python -m pytest tests/unit/test_lineage_assign.py -q`
Expected: PASS (4 passed).

- [ ] **Step 5: Lint + full fast suite + commit**

Run: `ruff check src/tessera/discover/lineage_assign.py` and `python -m pytest -m "not requires_binary" -q` (expect all pass).

```bash
git add src/tessera/discover/lineage_assign.py tests/unit/test_lineage_assign.py
git commit -m "Add assign_lineages orchestration (metadata -> nearest-ref -> de-novo)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: the `type-lineages` CLI command + docs

**Files:**
- Create: `src/tessera/cli/cmd_type_lineages.py`
- Modify: `src/tessera/cli/main.py` (register the subcommand)
- Modify: `docs/reference-panels.md` (document the command)
- Test: none new (thin Typer wiring; the logic is covered by Tasks 1-4). Verify via `--help`.

**Interfaces:**
- Consumes: `assign_lineages` (Task 4); `write_lineage_map`, `LINEAGES_TSV` (`recomb/typing.py`); the `app`/`get_logger`/`stage_errors` CLI helpers (`cli/main.py`).

- [ ] **Step 1: Create the command**

```python
# src/tessera/cli/cmd_type_lineages.py
"""The ``tessera type-lineages`` command: assign lineages to a genome collection.

Types an NCBI/local set of genomes (header/metadata, then nearest-neighbour to a
Nextclade dataset, then de-novo ANI clustering) and writes ``lineages.tsv``, so the
collection is lineage-ready for detection and panel building. Needs skani; the
nearest-neighbour step contacts the Nextclade dataset server on first use.
"""

from __future__ import annotations

import os
from pathlib import Path

import typer

from .main import app, get_logger, stage_errors


@app.command(name="type-lineages")
def type_lineages(
    collection: Path = typer.Option(
        ..., "-c", "--collection", help="Directory of genome FASTAs to type."
    ),
    output: Path = typer.Option(
        ..., "-o", "--output", help="Output directory (writes lineages.tsv)."
    ),
    taxon: str | None = typer.Option(
        None, "--taxon", help="Taxon for the Nextclade reference step (auto-detected if omitted)."
    ),
    nextclade_dataset: str | None = typer.Option(
        None, "--nextclade-dataset", help="Nextclade dataset path (else auto-detected)."
    ),
    lineage_map: Path | None = typer.Option(
        None, "--lineage-map", help="User TSV (accession<TAB>genotype) taking precedence."
    ),
    ref_ani_floor: float = typer.Option(
        90.0, "--ref-ani-floor", help="Min ANI to inherit a nearest reference tip's clade."
    ),
    cluster_ani: float = typer.Option(
        95.0, "--cluster-ani", help="De-novo: genomes at/above this ANI share a lineage."
    ),
    email: str | None = typer.Option(
        None, "--email", help="Contact email for NCBI (or set NCBI_EMAIL)."
    ),
) -> None:
    """Assign a lineage to every genome in a collection and write lineages.tsv."""
    from ..discover.lineage_assign import assign_lineages
    from ..recomb.typing import LINEAGES_TSV, write_lineage_map

    logger = get_logger(output)
    with stage_errors(logger):
        genomes = sorted(
            p for p in collection.iterdir()
            if p.is_file() and p.suffix.lower() in (".fasta", ".fa", ".fna")
        )
        if not genomes:
            raise typer.BadParameter(f"No FASTA genomes found in {collection}")
        rows = assign_lineages(
            genomes, user_lineage_map=lineage_map, taxon=taxon,
            nextclade_dataset=nextclade_dataset, ref_ani_floor=ref_ani_floor,
            cluster_ani=cluster_ani, email=email or os.environ.get("NCBI_EMAIL"),
            logger=logger,
        )
        output.mkdir(parents=True, exist_ok=True)
        write_lineage_map(output / LINEAGES_TSV, rows)
        logger.info("Wrote %d lineage row(s): %s", len(rows), output / LINEAGES_TSV)
```

- [ ] **Step 2: Register the subcommand**

In `src/tessera/cli/main.py`, add `cmd_type_lineages` to the side-effect import block (`:107-115`), keeping the list alphabetical:

```python
from . import (  # noqa: E402,F401
    cmd_build_panel,
    cmd_curate_panel,
    cmd_detect,
    cmd_fill_references,
    cmd_find_references,
    cmd_msa,
    cmd_recomb,
    cmd_type_lineages,
)
```

- [ ] **Step 3: Verify the command registers**

Run: `python -m pytest -m "not requires_binary" -q` (expect all pass) and
`ruff check src tests validation` (clean) and, via the installed console entry point,
`tessera type-lineages --help` (expect the options to render; do NOT use `python -m tessera.cli.main`).

- [ ] **Step 4: Document it**

In `docs/reference-panels.md`, add a short subsection under the typing/lineage discussion:

```
## Make a collection lineage-ready (`type-lineages`)

`tessera type-lineages -c genomes/ -o out/` assigns a lineage to every genome and writes
`out/lineages.tsv`, so an NCBI or local set is typed like a Nextclade pool. Each genome takes
the best label available: a header/metadata designation, else the clade of its nearest
Nextclade-dataset reference tip (skani ANI, above `--ref-ani-floor`), else a de-novo lineage
from ANI clustering (`--cluster-ani`). Downstream `recomb`/`detect`/`fill-references` read the
`lineages.tsv` beside the MSA or output, so the typed names carry through. Needs skani; the
nearest-neighbour step fetches the Nextclade dataset on first use.
```

- [ ] **Step 5: Commit**

```bash
git add src/tessera/cli/cmd_type_lineages.py src/tessera/cli/main.py docs/reference-panels.md
git commit -m "Add the tessera type-lineages command

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 6: Opt-in end-to-end check (needs skani; not CI)**

```bash
export PATH="$HOME/miniforge3/envs/recomfi-aln/bin:$PATH"
tessera type-lineages -c <a small genome dir> -o /tmp/tl_out --taxon "<taxon>"
column -t /tmp/tl_out/lineages.tsv | head
```

Expected: a `lineages.tsv` where genomes with known lineages are labelled (source `title`/
`nextclade-nn`) and the rest carry `denovo_*`. Spot-check a few labels against expectations.

---

## Notes for the implementer
- Tasks 1-5 Steps up to commit are CI-verifiable (the skani calls are stubbed in tests). Only Task 5 Step 6 (end-to-end) needs `skani` + network.
- Do not reimplement typing/Nextclade/skani helpers -- import and reuse them.
- `assign_lineages` monkeypatch targets in tests reference the module-level names
  (`la.build_lineage_map`, `la._reference_tips`, `la.assign_by_nearest_reference`,
  `la._ani_edges`), so those must be module-level imports/defs (not inline), matching the code above.
