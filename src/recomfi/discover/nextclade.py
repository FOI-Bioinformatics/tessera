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
