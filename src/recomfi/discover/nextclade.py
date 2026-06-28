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
from .pool import detect_taxon

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
