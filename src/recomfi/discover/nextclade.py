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
from ..core.cache import cached_genomes
from ..core.errors import UserInputError
from ..core.io import write_fasta_record
from ..core.plugins import ToolCapabilities
from ..core.process import run_tool
from .pool import detect_taxon

_BASE_URL = "https://data.clades.nextstrain.org/v3"
_INDEX_URL = f"{_BASE_URL}/index.json"

NEXTCLADE = ToolCapabilities(
    name="nextclade", conda=("bioconda::nextclade",),
    required_binaries=(BinarySpec("nextclade", version_args=("--version",)),),
)

# A nuc mutation: ref-base, 1-based position, alt-base ('-' = deletion).
_MUT = re.compile(r"^([A-Za-z])(\d+)([A-Za-z-])$")
# Clade-label node attributes in priority order. The dataset-specific informative
# labels come first -- Nextclade_pango / clade_nextstrain (SARS-CoV-2), LANL_subtype
# (HIV) -- then clade_membership, the universal key present in most datasets but often
# the coarser label; subclade (flu) is the final fallback.
_CLADE_KEYS = ("Nextclade_pango", "clade_nextstrain", "LANL_subtype", "clade_membership",
               "subclade")


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

        # Iterative DFS: avoids hitting Python's recursion limit on very deep trees.
        stack: list[tuple[dict, list[str]]] = [(tree, [])]
        while stack:
            node, inherited = stack.pop()
            nuc = (node.get("branch_attrs", {}).get("mutations", {}) or {}).get("nuc", [])
            path_muts = inherited + nuc
            children = node.get("children", [])
            if children:
                for child in children:
                    stack.append((child, path_muts))
                continue
            acc = _accession_of(node)
            if not acc or acc in seen:
                continue
            seen.add(acc)
            clade = _clade_of(node.get("node_attrs", {}) or {})
            out = _write_genome(build_dir, acc, clade,
                                _reconstruct_sequence(reference, path_muts), floor)
            if out is not None:
                written.append(out)

        logger.info("Reconstructed %d Nextclade tree-tip genome(s).", len(written))

        # Add example sequences; the header's first token (often subtype.country...)
        # gives a label and clade-ish note, deduped against the reconstructed tips.
        # A download failure here is non-fatal: the tips-only pool is still valid.
        if "examples" in dataset.files:
            try:
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
            except (UserInputError, URLError, OSError) as exc:
                logger.warning(
                    "Could not fetch Nextclade example sequences (%s); "
                    "using the reconstructed tips only.",
                    exc,
                )

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
