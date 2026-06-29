"""HAL -> MAF wrapper (Cactus output normalization).

Cactus emits a HAL alignment; ``hal2maf`` projects it to MAF against a reference
genome, which :mod:`tessera.converters.maf_to_fasta` then turns into an
MSA-FASTA.
"""

from __future__ import annotations

import logging
from pathlib import Path

from ..core.plugins import ToolCapabilities
from ..core.process import run_tool


def hal_to_maf(
    hal_path: str | Path,
    reference: str,
    out_path: str | Path,
    logger: logging.Logger,
    caps: ToolCapabilities | None = None,
) -> Path:
    """Project a HAL to MAF with ``hal2maf``.

    ``--noAncestors`` drops ancestral nodes; the Minigraph-Cactus
    ``_MINIGRAPH_`` backbone pseudo-genome is filtered downstream in
    :func:`maf_to_fasta`.
    """
    hal_path = Path(hal_path)
    out_path = Path(out_path)
    caps = caps or ToolCapabilities(name="hal2maf")
    run_tool(
        caps,
        ["hal2maf", "--refGenome", reference, "--noAncestors", hal_path, out_path],
        logger=logger,
        log_prefix="hal2maf",
    )
    return out_path
