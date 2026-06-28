# Nextclade Seed-Source Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `--seed-source nextclade` to recomfi that auto-detects a Nextclade dataset from the query, reconstructs its reference-tree tips into a candidate genome pool, and feeds that pool through the existing regional-selection path.

**Architecture:** A new self-contained module `discover/nextclade.py` resolves a dataset (override / `nextclade sort` / BLAST+alias) and reconstructs its tree tips + examples into clade-labeled genomes. `iterate.py` gains a thin `nextclade` branch that reuses the existing `select_regional` -> collection -> MSA -> HMM flow unchanged. A per-`path@tag` on-disk cache avoids re-fetching.

**Tech Stack:** Python 3.11+, stdlib `urllib`/`json`/`re`, recomfi's existing `core`/`discover`/`recomb` packages, Typer CLI, pytest with mocked `urlopen`/`run_tool`.

## Global Constraints

- Python >= 3.11; line length <= 100 (ruff `E,F,I,UP,B`, B008 ignored).
- No new Python dependencies — stdlib only (`urllib.request`, `json`, `re`, `hashlib`).
- The `nextclade` CLI is OPTIONAL: probed with `shutil.which`, never required.
- No Unicode in any file (user rule); modest scientific language in code/docs.
- Nextclade v3 file URL pattern: `https://data.clades.nextstrain.org/v3/{path}/{tag}/{filename}`; index at `https://data.clades.nextstrain.org/v3/index.json`.
- Errors raise `recomfi.core.errors.UserInputError` (never `sys.exit`).
- Tests are binary-free: mock `urlopen`/`run_tool`/`detect_taxon`; never require `nextclade`.
- Detection precedence: explicit `--candidate-pool` (local) > nextclade > blast.
- Dataset detection precedence: `--nextclade-dataset` override > `nextclade sort` > BLAST+alias.

---

### Task 1: `nextclade_cache` cache key

**Files:**
- Modify: `src/recomfi/core/cache.py` (add function after `ncbi_virus_cache`)
- Test: `tests/unit/test_cache.py` (new file)

**Interfaces:**
- Consumes: existing `cache_root`, `_slug` in `core/cache.py`.
- Produces: `nextclade_cache(path: str, tag: str, *, override: str | Path | None = None) -> Path`

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_cache.py`:

```python
"""On-disk cache key helpers."""

from __future__ import annotations

from pathlib import Path

from recomfi.core.cache import nextclade_cache


def test_nextclade_cache_path_is_under_override_and_stable(tmp_path: Path):
    a = nextclade_cache("community/neherlab/hiv-1/hxb2", "2026-04-14--11-55-23Z", override=tmp_path)
    b = nextclade_cache("community/neherlab/hiv-1/hxb2", "2026-04-14--11-55-23Z", override=tmp_path)
    assert a == b  # deterministic
    assert tmp_path in a.parents
    assert a.parent.name == "nextclade"


def test_nextclade_cache_tag_changes_key(tmp_path: Path):
    old = nextclade_cache("nextstrain/sars-cov-2/XBB", "2024-01-01--00-00-00Z", override=tmp_path)
    new = nextclade_cache("nextstrain/sars-cov-2/XBB", "2026-01-01--00-00-00Z", override=tmp_path)
    assert old != new  # a newer tag is a distinct cache dir
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_cache.py -v`
Expected: FAIL with `ImportError: cannot import name 'nextclade_cache'`

- [ ] **Step 3: Write minimal implementation**

In `src/recomfi/core/cache.py`, after the `ncbi_virus_cache` function, add:

```python
def nextclade_cache(path: str, tag: str, *, override: str | Path | None = None) -> Path:
    """Cache directory for a reconstructed Nextclade pool, keyed by ``path@tag``."""
    key = hashlib.sha1(f"{path}@{tag}".encode()).hexdigest()[:12]  # noqa: S324 - non-cryptographic
    return cache_root(override) / "nextclade" / f"{_slug(path)}_{key}"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_cache.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add src/recomfi/core/cache.py tests/unit/test_cache.py
git commit -m "feat: add nextclade_cache key helper"
```

---

### Task 2: Reconstruction and labeling helpers

**Files:**
- Create: `src/recomfi/discover/nextclade.py`
- Test: `tests/unit/test_nextclade.py` (new file)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces (module-private, used by Task 4):
  - `_reconstruct_sequence(reference: str, nuc_muts: list[str]) -> str` (ungapped, upper-cased)
  - `_clade_of(node_attrs: dict) -> str`
  - `_accession_of(node: dict) -> str`

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_nextclade.py`:

```python
"""Nextclade pool source: reconstruction, labeling, resolution, caching."""

from __future__ import annotations

from recomfi.discover import nextclade as nc


def test_reconstruct_applies_subs_and_deletions_and_strips_gaps():
    # reference ACGT; A1G sets pos1 -> G; T4- deletes pos4 -> stripped.
    assert nc._reconstruct_sequence("ACGT", ["A1G", "T4-"]) == "GCG"


def test_reconstruct_ignores_out_of_range_and_malformed():
    assert nc._reconstruct_sequence("ACGT", ["A9G", "garbage", "C2T"]) == "ATGT"


def test_clade_priority_prefers_pango_then_universal_membership():
    assert nc._clade_of({"Nextclade_pango": {"value": "BA.2"},
                         "clade_membership": {"value": "21L"}}) == "BA.2"
    assert nc._clade_of({"clade_membership": "H2"}) == "H2"  # bare scalar
    assert nc._clade_of({"author": {"value": "x"}}) == "NA"  # no clade key


def test_accession_prefers_attr_then_name_and_is_sanitised():
    assert nc._accession_of({"node_attrs": {"accession": {"value": "AY/03 7026"}},
                             "name": "ignored"}) == "AY_03_7026"
    assert nc._accession_of({"node_attrs": {}, "name": "JN417241"}) == "JN417241"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_nextclade.py -v`
Expected: FAIL with `ModuleNotFoundError` or `AttributeError` (module/functions absent)

- [ ] **Step 3: Write minimal implementation**

Create `src/recomfi/discover/nextclade.py`:

```python
"""Build a candidate genome pool from a Nextclade dataset.

Resolves a dataset from the query (an explicit override, else ``nextclade sort``
when the binary is present, else recomfi's BLAST taxon detection mapped to a
dataset by keyword), then reconstructs every reference-tree tip by applying its
root-to-tip nucleotide mutations to the dataset reference. The reconstructed tips
(plus the dataset's example sequences) form a diverse, clade-labeled pool that the
normal regional-selection path turns into a focused recombination panel.

A Nextclade dataset is a clade-typing reference tree, not a set of clean
recombinant-parental endpoints, so this pool suits broad screens better than the
recovery of a specific recombinant's parental mosaic (see the HIV-1 CRF01_AE
finding in the project memory). Everything here is stdlib HTTP + text, with the
``nextclade`` binary used only, and optionally, to improve dataset detection.
"""

from __future__ import annotations

import json
import logging
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen

from ..core.binaries import BinarySpec
from ..core.errors import UserInputError
from ..core.io import write_fasta_record
from ..core.plugins import ToolCapabilities
from ..core.process import run_tool

_BASE_URL = "https://data.clades.nextstrain.org/v3"
_INDEX_URL = f"{_BASE_URL}/index.json"

NEXTCLADE = ToolCapabilities(
    name="nextclade", conda=("bioconda::nextclade",),
    required_binaries=(BinarySpec("nextclade", version_args=("--version",)),),
)

# A nuc mutation: ref-base, 1-based position, alt-base ('-' = deletion).
_MUT = re.compile(r"^([A-Za-z])(\d+)([A-Za-z-])$")
# Clade-label node attributes, most specific first; clade_membership is universal.
_CLADE_KEYS = ("Nextclade_pango", "clade_nextstrain", "clade_membership", "subclade",
               "LANL_subtype")


def _attr(value):
    """A node attribute is either ``{"value": x}`` or a bare scalar."""
    return value.get("value") if isinstance(value, dict) else value


def _clade_of(node_attrs: dict) -> str:
    for key in _CLADE_KEYS:
        if key in node_attrs:
            val = _attr(node_attrs[key])
            if val:
                return str(val)
    return "NA"


def _accession_of(node: dict) -> str:
    acc = _attr((node.get("node_attrs") or {}).get("accession")) or node.get("name", "")
    return re.sub(r"[^\w.]", "_", str(acc))


def _reconstruct_sequence(reference: str, nuc_muts: list[str]) -> str:
    """Apply substitutions and deletions to ``reference``; return the ungapped genome."""
    chars = list(reference)
    for mut in nuc_muts:
        m = _MUT.match(mut)
        if not m:
            continue
        pos, alt = int(m.group(2)), m.group(3)
        i = pos - 1
        if 0 <= i < len(chars):
            chars[i] = alt  # '-' for a deletion
    return "".join(chars).replace("-", "").upper()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_nextclade.py -v`
Expected: PASS (4 passed)

- [ ] **Step 5: Commit**

```bash
git add src/recomfi/discover/nextclade.py tests/unit/test_nextclade.py
git commit -m "feat: nextclade tip reconstruction and clade labeling helpers"
```

---

### Task 3: Index fetch and dataset resolution

**Files:**
- Modify: `src/recomfi/discover/nextclade.py`
- Test: `tests/unit/test_nextclade.py`

**Interfaces:**
- Consumes: `_BASE_URL`, `_INDEX_URL`, `NEXTCLADE` from Task 2; `pool.detect_taxon`.
- Produces:
  - `@dataclass NextcladeDataset(path: str, tag: str, files: dict)`
  - `fetch_index() -> dict[str, dict]` (path -> index entry)
  - `_dataset_file_url(dataset: NextcladeDataset, role: str) -> str`
  - `resolve_dataset(query: Path, override: str | None, *, email: str | None, logger) -> NextcladeDataset`
  - `_ALIAS: dict[str, str]`

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_nextclade.py`:

```python
import json as _json
from pathlib import Path

import pytest

from recomfi.core.errors import UserInputError

_FAKE_INDEX = {
    "collections": [{
        "datasets": [
            {"path": "community/neherlab/hiv-1/hxb2",
             "versions": [{"tag": "2026-04-14--11-55-23Z"}, {"tag": "2024-01-01--00-00-00Z"}],
             "files": {"reference": "reference.fasta", "treeJson": "tree.json",
                       "examples": "sequences.fasta"}},
            {"path": "nextstrain/sars-cov-2/XBB",
             "versions": [{"tag": "2025-09-09--12-13-13Z"}],
             "files": {"reference": "reference.fasta", "treeJson": "tree.json",
                       "examples": "sequences.fasta"}},
        ],
    }],
}


class _Resp:
    def __init__(self, data: bytes):
        self._data = data
    def read(self):
        return self._data
    def __enter__(self):
        return self
    def __exit__(self, *a):
        return False


def _patch_index(monkeypatch):
    monkeypatch.setattr(nc, "urlopen",
                        lambda url, timeout=0: _Resp(_json.dumps(_FAKE_INDEX).encode()))


def test_fetch_index_flattens_and_keeps_newest_tag(monkeypatch):
    _patch_index(monkeypatch)
    idx = nc.fetch_index()
    assert set(idx) == {"community/neherlab/hiv-1/hxb2", "nextstrain/sars-cov-2/XBB"}
    assert idx["community/neherlab/hiv-1/hxb2"]["versions"][0]["tag"] == "2026-04-14--11-55-23Z"


def test_resolve_override_validates_against_index(monkeypatch, tmp_path, logger):
    _patch_index(monkeypatch)
    q = tmp_path / "q.fasta"; q.write_text(">q\nACGT\n")
    ds = nc.resolve_dataset(q, "nextstrain/sars-cov-2/XBB", email=None, logger=logger)
    assert ds.path == "nextstrain/sars-cov-2/XBB"
    assert ds.tag == "2025-09-09--12-13-13Z"
    assert ds.files["treeJson"] == "tree.json"


def test_resolve_override_unknown_path_raises(monkeypatch, tmp_path, logger):
    _patch_index(monkeypatch)
    q = tmp_path / "q.fasta"; q.write_text(">q\nACGT\n")
    with pytest.raises(UserInputError, match="not a Nextclade dataset"):
        nc.resolve_dataset(q, "nextstrain/nope", email=None, logger=logger)


def test_resolve_blast_fallback_maps_organism_to_dataset(monkeypatch, tmp_path, logger):
    _patch_index(monkeypatch)
    monkeypatch.setattr(nc.shutil, "which", lambda name: None)  # no nextclade binary
    monkeypatch.setattr(nc, "detect_taxon",
                        lambda q, *, email, logger: "Human immunodeficiency virus 1")
    q = tmp_path / "q.fasta"; q.write_text(">q\nACGT\n")
    ds = nc.resolve_dataset(q, None, email=None, logger=logger)
    assert ds.path == "community/neherlab/hiv-1/hxb2"


def test_resolve_blast_fallback_unmapped_raises(monkeypatch, tmp_path, logger):
    _patch_index(monkeypatch)
    monkeypatch.setattr(nc.shutil, "which", lambda name: None)
    monkeypatch.setattr(nc, "detect_taxon", lambda q, *, email, logger: "Unknown mystery virus")
    q = tmp_path / "q.fasta"; q.write_text(">q\nACGT\n")
    with pytest.raises(UserInputError, match="--nextclade-dataset"):
        nc.resolve_dataset(q, None, email=None, logger=logger)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_nextclade.py -k "index or resolve" -v`
Expected: FAIL (`fetch_index` / `resolve_dataset` undefined)

- [ ] **Step 3: Write minimal implementation**

Append to `src/recomfi/discover/nextclade.py`:

```python
from .pool import detect_taxon  # placed with the other imports at the top of the file

# Keyword (lower-cased organism substring) -> Nextclade dataset path, for the
# BLAST-detection fallback. Covers the core pathogens; extend as needed.
_ALIAS: dict[str, str] = {
    "human immunodeficiency virus": "community/neherlab/hiv-1/hxb2",
    "severe acute respiratory syndrome coronavirus 2": "nextstrain/sars-cov-2/XBB",
    "monkeypox": "nextstrain/mpox/all-clades",
    "mpox": "nextstrain/mpox/all-clades",
    "measles": "nextstrain/measles/genome/WHO-2012",
    "respiratory syncytial virus": "nextstrain/rsv/a/EPI_ISL_412866",
    "dengue": "nextstrain/dengue/all",
    "mumps": "nextstrain/mumps/genome",
    "rubella": "nextstrain/rubella/genome",
    "west nile": "nextstrain/wnv/all-lineages",
    "yellow fever": "nextstrain/yellow-fever/prM-E",
    "ebola": "nextstrain/orthoebolavirus/ebov",
    "zaire ebolavirus": "nextstrain/orthoebolavirus/ebov",
    "metapneumovirus": "nextstrain/hmpv/all-clades/NC_039199",
}


@dataclass
class NextcladeDataset:
    path: str
    tag: str
    files: dict


def fetch_index() -> dict[str, dict]:
    """Download and flatten the Nextclade v3 index to ``{path: entry}``."""
    try:
        with urlopen(_INDEX_URL, timeout=60) as resp:  # noqa: S310 - fixed https URL
            data = json.loads(resp.read())
    except (URLError, OSError, json.JSONDecodeError) as exc:
        raise UserInputError(f"Could not fetch the Nextclade dataset index: {exc}") from exc
    out: dict[str, dict] = {}

    def walk(coll: dict) -> None:
        for child in coll.get("collections", []):
            walk(child)
        for ds in coll.get("datasets", []):
            if ds.get("path") and ds.get("versions"):
                out[ds["path"]] = ds

    walk(data)
    return out


def _dataset_file_url(dataset: NextcladeDataset, role: str) -> str:
    filename = dataset.files[role]
    return f"{_BASE_URL}/{dataset.path}/{dataset.tag}/{filename}"


def _entry_to_dataset(path: str, entry: dict) -> NextcladeDataset:
    return NextcladeDataset(path=path, tag=entry["versions"][0]["tag"],
                            files=entry.get("files", {}))


def _sort_dataset(query: Path, index: dict[str, dict], logger) -> str | None:
    """Best dataset path from ``nextclade sort``; None on any failure (caller falls back)."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "sort.tsv"
        try:
            run_tool(NEXTCLADE,
                     ["nextclade", "sort", str(query), "--output-results-tsv", str(out)],
                     logger=logger, log_prefix="nextclade-sort")
        except Exception as exc:  # noqa: BLE001 - any failure becomes a clean fallback
            logger.warning("nextclade sort failed (%s); falling back to BLAST detection.", exc)
            return None
        if not out.exists():
            return None
        rows = out.read_text().splitlines()
    if len(rows) < 2:
        return None
    header = rows[0].split("\t")
    if "dataset" not in header:
        return None
    col = header.index("dataset")
    for line in rows[1:]:
        fields = line.split("\t")
        if len(fields) > col and fields[col] in index:
            return fields[col]
    return None


def resolve_dataset(
    query: Path, override: str | None, *, email: str | None, logger: logging.Logger
) -> NextcladeDataset:
    """Resolve the Nextclade dataset for ``query`` (override / sort / BLAST+alias)."""
    index = fetch_index()
    if override:
        if override not in index:
            raise UserInputError(
                f"'{override}' is not a Nextclade dataset. Browse the index at "
                f"{_INDEX_URL} and pass a valid --nextclade-dataset path."
            )
        logger.info("Using the Nextclade dataset '%s' (from --nextclade-dataset).", override)
        return _entry_to_dataset(override, index[override])

    if shutil.which("nextclade") is not None:
        path = _sort_dataset(query, index, logger)
        if path:
            logger.info("nextclade sort selected the dataset '%s'.", path)
            return _entry_to_dataset(path, index[path])

    taxon = detect_taxon(query, email=email, logger=logger).lower()
    for keyword, path in _ALIAS.items():
        if keyword in taxon and path in index:
            logger.info("Mapped '%s' to the Nextclade dataset '%s'.", taxon, path)
            return _entry_to_dataset(path, index[path])
    raise UserInputError(
        f"Could not map the query's taxon ('{taxon}') to a Nextclade dataset. "
        f"Pass --nextclade-dataset <path> explicitly (see {_INDEX_URL})."
    )
```

Move the `from .pool import detect_taxon` line up to the import block at the top of the file (keep imports together; the inline comment above is only to show where it belongs).

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_nextclade.py -v`
Expected: PASS (all reconstruction + index + resolve tests pass)

- [ ] **Step 5: Commit**

```bash
git add src/recomfi/discover/nextclade.py tests/unit/test_nextclade.py
git commit -m "feat: resolve a Nextclade dataset from the query (override/sort/blast)"
```

---

### Task 4: `build_pool` (reconstruct + cache)

**Files:**
- Modify: `src/recomfi/discover/nextclade.py`
- Test: `tests/unit/test_nextclade.py`

**Interfaces:**
- Consumes: `NextcladeDataset`, `_dataset_file_url`, `_reconstruct_sequence`, `_clade_of`, `_accession_of`, `urlopen`.
- Produces: `build_pool(dataset: NextcladeDataset, *, cache_dir: Path, logger) -> list[Path]`

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_nextclade.py`:

```python
_REF = ">NC ref\n" + "ACGTACGTAC" * 1 + "\n"  # 10 bp reference

_TREE = {
    "tree": {
        "name": "root",
        "branch_attrs": {"mutations": {"nuc": ["A1G"]}},
        "children": [
            {"name": "TIP1",
             "branch_attrs": {"mutations": {"nuc": ["C2T"]}},
             "node_attrs": {"accession": "ACC1", "clade_membership": {"value": "X"}},
             "children": []},
            {"name": "TIP2",
             "branch_attrs": {"mutations": {"nuc": ["T4-"]}},
             "node_attrs": {"accession": "ACC2", "clade_membership": "Y"},
             "children": []},
        ],
    },
}

_EXAMPLES = ">EX1 example\n" + "ACGTACGTAC" + "\n"


def _patch_files(monkeypatch):
    payloads = {
        "reference.fasta": _REF.encode(),
        "tree.json": _json.dumps(_TREE).encode(),
        "sequences.fasta": _EXAMPLES.encode(),
    }

    def fake_urlopen(url, timeout=0):
        for name, data in payloads.items():
            if url.endswith(name):
                return _Resp(data)
        raise AssertionError(f"unexpected url {url}")

    monkeypatch.setattr(nc, "urlopen", fake_urlopen)


def _dataset():
    return nc.NextcladeDataset(
        path="x/y", tag="t",
        files={"reference": "reference.fasta", "treeJson": "tree.json",
               "examples": "sequences.fasta"},
    )


def test_build_pool_reconstructs_tips_examples_and_labels(monkeypatch, tmp_path, logger):
    _patch_files(monkeypatch)
    genomes = nc.build_pool(_dataset(), cache_dir=tmp_path / "c", logger=logger)
    by_name = {g.stem: g.read_text() for g in genomes}
    assert {"ACC1", "ACC2", "EX1"} <= set(by_name)
    # TIP1: A1G then C2T on ACGTACGTAC -> GTGTACGTAC, header carries clade X
    assert ">ACC1 X" in by_name["ACC1"]
    assert "GTGTACGTAC" in by_name["ACC1"].replace("\n", "")
    # TIP2: A1G then T4 deleted -> GCGACGTAC (gap stripped), clade Y
    assert ">ACC2 Y" in by_name["ACC2"]


def test_build_pool_uses_cache_on_second_call(monkeypatch, tmp_path, logger):
    _patch_files(monkeypatch)
    cache = tmp_path / "c"
    nc.build_pool(_dataset(), cache_dir=cache, logger=logger)

    def boom(url, timeout=0):
        raise AssertionError("network must not be touched on a cache hit")

    monkeypatch.setattr(nc, "urlopen", boom)
    again = nc.build_pool(_dataset(), cache_dir=cache, logger=logger)
    assert {g.stem for g in again} >= {"ACC1", "ACC2", "EX1"}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_nextclade.py -k build_pool -v`
Expected: FAIL (`build_pool` undefined)

- [ ] **Step 3: Write minimal implementation**

Append to `src/recomfi/discover/nextclade.py`:

```python
from ..core.cache import cached_genomes  # add to the top import block


def _download_text(dataset: NextcladeDataset, role: str, logger: logging.Logger) -> str:
    url = _dataset_file_url(dataset, role)
    try:
        with urlopen(url, timeout=120) as resp:  # noqa: S310 - fixed https URL
            return resp.read().decode("utf-8", "replace")
    except (URLError, OSError) as exc:
        raise UserInputError(f"Could not download {url}: {exc}") from exc


def _read_reference(text: str) -> str:
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))


def _write_genome(out_dir: Path, label: str, note: str, seq: str, floor: int) -> Path | None:
    if len(seq) < floor:  # a truncated reconstruction
        return None
    path = out_dir / f"{label}.fasta"
    with open(path, "w") as fo:
        write_fasta_record(fo, f"{label} {note}", seq)
    return path


def build_pool(
    dataset: NextcladeDataset, *, cache_dir: Path, logger: logging.Logger
) -> list[Path]:
    """Return the Nextclade pool genomes for ``dataset`` (cached per ``path@tag``).

    Genomes are reconstructed into a temporary directory and atomically moved into
    ``cache_dir`` on success, so a download interrupted midway never leaves a
    partial directory that a later run would mistake for a complete cache.
    """
    import os
    import shutil as _shutil
    import tempfile

    existing = cached_genomes(cache_dir)
    if existing:
        logger.info("Using the cached Nextclade pool for '%s' (%d genome(s)): %s",
                    dataset.path, len(existing), cache_dir)
        return existing

    reference = _read_reference(_download_text(dataset, "reference", logger))
    if not reference:
        raise UserInputError(f"Nextclade dataset '{dataset.path}' has an empty reference.")
    floor = max(1, len(reference) // 2)
    tree = json.loads(_download_text(dataset, "treeJson", logger))["tree"]

    cache_dir.parent.mkdir(parents=True, exist_ok=True)
    build_dir = Path(tempfile.mkdtemp(dir=cache_dir.parent))
    try:
        written: list[Path] = []
        seen: set[str] = set()

        def walk(node: dict, inherited: list[str]) -> None:
            nuc = (node.get("branch_attrs", {}).get("mutations", {}) or {}).get("nuc", [])
            path_muts = inherited + nuc
            children = node.get("children", [])
            if children:
                for child in children:
                    walk(child, path_muts)
                return
            acc = _accession_of(node)
            if not acc or acc in seen:
                return
            seen.add(acc)
            clade = _clade_of(node.get("node_attrs", {}) or {})
            out = _write_genome(build_dir, acc, clade,
                                _reconstruct_sequence(reference, path_muts), floor)
            if out is not None:
                written.append(out)

        walk(tree, [])
        logger.info("Reconstructed %d Nextclade tree-tip genome(s).", len(written))

        # Add example sequences; the header's first token (often subtype.country...)
        # gives a label and clade-ish note, deduped against the reconstructed tips.
        if "examples" in dataset.files:
            examples = 0
            for header, seq in read_fasta_text(_download_text(dataset, "examples", logger)):
                acc = re.sub(r"[^\w.]", "_", header.split()[0]) if header else ""
                if not acc or acc in seen:
                    continue
                seen.add(acc)
                note = header.split(".")[0] if "." in header else "example"
                out = _write_genome(build_dir, acc, note, seq.upper(), floor)
                if out is not None:
                    written.append(out)
                    examples += 1
            logger.info("Added %d Nextclade example genome(s).", examples)

        if not written:
            raise UserInputError(
                f"Nextclade dataset '{dataset.path}' yielded no usable pool genomes."
            )
        os.replace(build_dir, cache_dir)  # atomic on the same filesystem
    except BaseException:
        _shutil.rmtree(build_dir, ignore_errors=True)
        raise

    genomes = cached_genomes(cache_dir)  # paths now under the final cache dir
    logger.info("Nextclade pool: %d genome(s) -> %s", len(genomes), cache_dir)
    return genomes


def read_fasta_text(text: str) -> list[tuple[str, str]]:
    """Parse FASTA from an in-memory string into ``(header, sequence)`` records."""
    records: list[tuple[str, str]] = []
    name: str | None = None
    seq: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq)))
            name = line[1:].strip()
            seq = []
        elif name is not None:
            seq.append(line.strip())
    if name is not None:
        records.append((name, "".join(seq)))
    return records
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_nextclade.py -v`
Expected: PASS (all nextclade tests)

- [ ] **Step 5: Verify lint and full module**

Run: `ruff check src/recomfi/discover/nextclade.py && pytest tests/unit/test_nextclade.py -q`
Expected: no lint errors; all tests pass

- [ ] **Step 6: Commit**

```bash
git add src/recomfi/discover/nextclade.py tests/unit/test_nextclade.py
git commit -m "feat: build a cached Nextclade genome pool from tree tips and examples"
```

---

### Task 5: Wire the nextclade source into the fill loop

**Files:**
- Modify: `src/recomfi/discover/iterate.py`
- Test: `tests/unit/test_iterate.py`

**Interfaces:**
- Consumes: `nextclade.resolve_dataset`, `nextclade.build_pool`, `core.cache.nextclade_cache`, existing `_select_from`, `_copy_into`, `_seed_from_pool`.
- Produces: `FillParams.nextclade_dataset: str | None`; a `nextclade` branch in `_seed_from_pool`; `_fetch_nextclade(params, logger) -> list[Path]`.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_iterate.py`:

```python
def test_seed_source_nextclade_routes_through_pool_selection(monkeypatch, tmp_path, logger):
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "ACGT" * 100 + "\n")
    out = tmp_path / "out"
    _common_mocks(monkeypatch, [([], 0.95)])  # converge immediately, no NCBI rounds
    monkeypatch.setattr(iterate, "read_fasta", lambda p: [("q", "ACGT" * 100)])

    pool_dir = tmp_path / "pool"
    pool_dir.mkdir()
    (pool_dir / "REF1.fasta").write_text(">REF1 A1\nACGT\n")

    captured = {}

    def fake_fetch(params, logger):
        captured["dataset"] = params.nextclade_dataset
        return [pool_dir / "REF1.fasta"]

    def fake_select(params, genomes, logger):
        from recomfi.discover.pool import PoolSelection
        return PoolSelection(selected=list(genomes))

    monkeypatch.setattr(iterate, "_fetch_nextclade", fake_fetch)
    monkeypatch.setattr(iterate, "_select_from", fake_select)

    fill_references(
        FillParams(query=query, collection=None, output=out, seed_source="nextclade",
                   nextclade_dataset="nextstrain/sars-cov-2/XBB"),
        logger,
    )
    assert captured["dataset"] == "nextstrain/sars-cov-2/XBB"
    assert (out / "collection" / "REF1.fasta").exists()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_iterate.py -k nextclade -v`
Expected: FAIL (`FillParams` has no `nextclade_dataset`, or `_fetch_nextclade` undefined)

- [ ] **Step 3: Write minimal implementation**

In `src/recomfi/discover/iterate.py`:

(a) Add the field to `FillParams` (next to `candidate_pool` / `taxon`):

```python
    nextclade_dataset: str | None = None  # Nextclade dataset path; None = auto-detect
```

(b) In `_seed_collection`, broaden the pool-source membership test:

```python
    if params.seed_source in ("local", "ncbi-virus", "nextclade"):
        _seed_from_pool(params, collection, logger)
        return
```

(c) In `_seed_from_pool`, add the nextclade branch before the local/ncbi branch:

```python
    if params.seed_source == "nextclade" and not force_ncbi:
        genomes = _fetch_nextclade(params, logger)
    elif params.seed_source == "local" and not force_ncbi:
        if params.candidate_pool is None:
            raise UserInputError("--seed-source local needs --candidate-pool <dir>.")
        genomes = iter_pool_genomes(params.candidate_pool)
    else:  # ncbi-virus (or auto-switched into it)
        genomes = _fetch_diverse(params, logger)
    _copy_into(_select_from(params, genomes, logger).selected, collection, logger)
```

(d) Add the fetch helper near `_fetch_diverse`:

```python
def _fetch_nextclade(params: FillParams, logger: logging.Logger) -> list[Path]:
    """Reconstruct a candidate pool from a Nextclade dataset (cached per path@tag)."""
    from ..core.cache import nextclade_cache
    from .nextclade import build_pool, resolve_dataset

    dataset = resolve_dataset(
        params.query, params.nextclade_dataset, email=params.email, logger=logger
    )
    cache = nextclade_cache(dataset.path, dataset.tag, override=params.cache_dir)
    return build_pool(dataset, cache_dir=cache, logger=logger)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_iterate.py -k nextclade -v`
Expected: PASS

- [ ] **Step 5: Run the full discover suite for regressions**

Run: `pytest tests/unit/test_iterate.py tests/unit/test_pool.py -q`
Expected: all pass

- [ ] **Step 6: Commit**

```bash
git add src/recomfi/discover/iterate.py tests/unit/test_iterate.py
git commit -m "feat: route seed-source nextclade through regional pool selection"
```

---

### Task 6: CLI options

**Files:**
- Modify: `src/recomfi/cli/cmd_fill_references.py`
- Modify: `src/recomfi/cli/cmd_detect.py`
- Modify: `src/recomfi/cli/cmd_build_panel.py`
- Test: `tests/unit/test_cli_nextclade.py` (new file)

**Interfaces:**
- Consumes: `FillParams.nextclade_dataset`, `seed_source="nextclade"` from Task 5; the Typer `app`.
- Produces: `--nextclade-dataset` on all three commands; `--nextclade/--no-nextclade` on `detect` and `build-panel`; `nextclade` added to the `fill-references` `--seed-source` choices.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_cli_nextclade.py`:

```python
"""CLI: the nextclade seed-source options reach FillParams."""

from __future__ import annotations

from typer.testing import CliRunner

from recomfi.cli import cmd_build_panel, cmd_detect, cmd_fill_references  # noqa: F401
from recomfi.cli.main import app

runner = CliRunner()


def _capture(monkeypatch):
    captured = {}
    import recomfi.discover.iterate as iterate

    def fake_fill(params, logger):
        captured["params"] = params
        return []

    # all three commands import fill_references from discover.iterate at call time
    monkeypatch.setattr(iterate, "fill_references", fake_fill)
    return captured


def test_detect_nextclade_flag_sets_seed_source(monkeypatch, tmp_path):
    captured = _capture(monkeypatch)
    q = tmp_path / "q.fasta"; q.write_text(">q\nACGT\n")
    result = runner.invoke(app, ["detect", "-q", str(q), "-o", str(tmp_path / "o"),
                                 "--nextclade"])
    assert result.exit_code == 0, result.output
    assert captured["params"].seed_source == "nextclade"


def test_detect_nextclade_dataset_override(monkeypatch, tmp_path):
    captured = _capture(monkeypatch)
    q = tmp_path / "q.fasta"; q.write_text(">q\nACGT\n")
    result = runner.invoke(app, ["detect", "-q", str(q), "-o", str(tmp_path / "o"),
                                 "--nextclade-dataset", "nextstrain/sars-cov-2/XBB"])
    assert result.exit_code == 0, result.output
    assert captured["params"].seed_source == "nextclade"
    assert captured["params"].nextclade_dataset == "nextstrain/sars-cov-2/XBB"


def test_fill_references_accepts_nextclade_seed_source(monkeypatch, tmp_path):
    captured = _capture(monkeypatch)
    q = tmp_path / "q.fasta"; q.write_text(">q\nACGT\n")
    result = runner.invoke(app, ["fill-references", "-q", str(q), "-o", str(tmp_path / "o"),
                                 "--seed-source", "nextclade"])
    assert result.exit_code == 0, result.output
    assert captured["params"].seed_source == "nextclade"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_cli_nextclade.py -v`
Expected: FAIL (unknown options / seed-source choice rejected)

- [ ] **Step 3: Implement the CLI changes**

In `src/recomfi/cli/cmd_fill_references.py`:
- In the `_require_choice(seed_source, {...}, "--seed-source")` call, add `"nextclade"`:

```python
        _require_choice(seed_source, {"blast", "local", "ncbi-virus", "nextclade"}, "--seed-source")
```

- Add a new option in the signature (next to `taxon`):

```python
    nextclade_dataset: str | None = typer.Option(
        None, "--nextclade-dataset",
        help="Nextclade dataset path for --seed-source nextclade "
        "(e.g. nextstrain/sars-cov-2/XBB). Auto-detected from the query if omitted.",
    ),
```

- Pass it into `FillParams(...)`:

```python
            taxon=taxon, nextclade_dataset=nextclade_dataset,
```

In `src/recomfi/cli/cmd_detect.py` and `src/recomfi/cli/cmd_build_panel.py`, add two options (next to `candidate_pool`):

```python
    nextclade: bool = typer.Option(
        False, "--nextclade/--no-nextclade",
        help="Recruit the panel from a Nextclade dataset auto-detected from the query.",
    ),
    nextclade_dataset: str | None = typer.Option(
        None, "--nextclade-dataset",
        help="Nextclade dataset path (implies --nextclade; e.g. nextstrain/sars-cov-2/XBB).",
    ),
```

In both `detect` and `build_panel`, replace the `seed_source=` line in the `FillParams(...)` construction. Current line:

```python
            seed_source="local" if candidate_pool else "blast",
```

becomes:

```python
            seed_source=_seed_source(candidate_pool, nextclade, nextclade_dataset),
            nextclade_dataset=nextclade_dataset,
```

and add this helper at the top of each module (after the imports, before the command function):

```python
def _seed_source(candidate_pool, nextclade: bool, nextclade_dataset: str | None) -> str:
    """Local pool wins, then Nextclade, then the default BLAST recruitment."""
    if candidate_pool:
        return "local"
    if nextclade or nextclade_dataset:
        return "nextclade"
    return "blast"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_cli_nextclade.py -v`
Expected: PASS (3 passed)

- [ ] **Step 5: Lint and full CLI regression**

Run: `ruff check src/recomfi/cli && pytest tests/unit/test_cli_nextclade.py -q`
Expected: no lint errors; tests pass

- [ ] **Step 6: Commit**

```bash
git add src/recomfi/cli/cmd_fill_references.py src/recomfi/cli/cmd_detect.py \
        src/recomfi/cli/cmd_build_panel.py tests/unit/test_cli_nextclade.py
git commit -m "feat: expose --nextclade / --nextclade-dataset on detect, build-panel, fill-references"
```

---

### Task 7: Documentation

**Files:**
- Modify: `README.md` (the seed-source section, after the `ncbi-virus` bullet around the "finite pool" paragraph)

**Interfaces:** none (docs only).

- [ ] **Step 1: Add the Nextclade subsection**

In `README.md`, after the `- `ncbi-virus` —` bullet in the "finite pool" `--seed-source` list, add a third bullet:

```markdown
- `nextclade` — a pool reconstructed from a Nextclade dataset
  (https://docs.nextstrain.org/projects/nextclade/en/stable/user/datasets.html).
  recomfi auto-detects the dataset from the query (`nextclade sort` when the CLI is
  installed, otherwise BLAST taxon detection mapped to a dataset), or you pass it
  explicitly with `--nextclade-dataset <path>` (e.g. `nextstrain/sars-cov-2/XBB`,
  `community/neherlab/hiv-1/hxb2`). Every reference-tree tip is reconstructed from
  the dataset reference plus its mutations and labelled by clade, so the report
  names parents by clade. Fetched pools are cached per dataset version.
```

Then add a short usage block after the existing `--seed-source local` example:

```markdown
```
recomfi detect --query CRF01_AE.fasta --output out/ --nextclade
recomfi fill-references --query q.fasta --output out/ --seed-source nextclade \
    --nextclade-dataset nextstrain/sars-cov-2/XBB
```
A Nextclade dataset is a clade-typing reference tree, so the pool spans clade
representatives rather than clean recombinant-parental endpoints; it suits a broad
screen. Single-segment datasets (most flu segments) yield segment-length pool
genomes, which is logged. The `nextclade` CLI is optional (it only sharpens dataset
auto-detection).
```

- [ ] **Step 2: Verify the docs render and reference real flags**

Run: `grep -n "seed-source nextclade\|--nextclade" README.md`
Expected: the new bullet and usage block are present and reference `--nextclade` / `--nextclade-dataset` / `--seed-source nextclade`.

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "docs: document the nextclade seed-source"
```

---

### Task 8: Full-suite verification

**Files:** none (verification only).

- [ ] **Step 1: Run the whole binary-free suite**

Run: `pytest -m "not requires_binary" -q`
Expected: all pass (existing ~110 tests + the new cache/nextclade/iterate/cli tests)

- [ ] **Step 2: Run the linter over everything touched**

Run: `ruff check src tests`
Expected: no errors

- [ ] **Step 3: Smoke-check the CLI help surfaces the new options**

Run: `python -m recomfi.cli.main detect --help | grep -i nextclade`
Expected: `--nextclade` and `--nextclade-dataset` appear

- [ ] **Step 4: Commit any lint fixes**

```bash
git add -A
git commit -m "chore: lint and suite green for the nextclade seed-source" || echo "nothing to commit"
```

---

## Notes for the implementer

- The `nextclade` binary is never required by tests; detection failure must fall through to BLAST, and `build_pool` must work entirely over mocked `urlopen`.
- Keep `discover/nextclade.py` self-contained — `iterate.py` only calls `resolve_dataset` and `build_pool` (via `_fetch_nextclade`).
- Reconstruction is reference-coordinate: substitutions and deletions only; insertions relative to the reference are intentionally dropped.
- Module imports are added in Task 2 but several are first *used* in Tasks 3-4; the first explicit `ruff check` is in Task 4 Step 5, by which point every import is referenced.
- The existing `_fetch_diverse`/`_seed_from_pool` structure is the template; the nextclade branch should read as a sibling of the ncbi-virus branch.
