"""Download candidate reference genomes from NCBI via Entrez Direct (efetch).

The in-package counterpart of ``validation/fetch.py``: routes through
``core.process.run_tool`` and declares ``efetch`` as a binary so a missing tool
is reported the same way as a missing aligner.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.errors import OutputError
from ..core.io import strip_sequence_extension
from ..core.plugins import ToolCapabilities
from ..core.process import run_tool

EFETCH = ToolCapabilities(
    name="efetch",
    conda=("bioconda::entrez-direct",),
    required_binaries=(BinarySpec("efetch", version_args=()),),
)


def efetch_available() -> bool:
    return shutil.which("efetch") is not None


def efetch_fasta(accession: str, collection_dir: Path, logger: logging.Logger) -> Path:
    """Download ``accession`` as FASTA into ``collection_dir`` as ``<label>.fasta``.

    The label is the accession with sequence extensions stripped, so the file
    stages cleanly as a collection genome later.
    """
    collection_dir.mkdir(parents=True, exist_ok=True)
    label = strip_sequence_extension(accession)
    dest = collection_dir / f"{label}.fasta"
    run_tool(
        EFETCH,
        ["efetch", "-db", "nuccore", "-id", accession, "-format", "fasta"],
        logger=logger, log_prefix=f"efetch:{accession}", stdout_path=dest,
    )
    if not dest.exists() or dest.stat().st_size == 0 or not dest.read_text().startswith(">"):
        raise OutputError(f"efetch returned no FASTA for accession '{accession}'.")
    return dest
