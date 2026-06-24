"""Build a reference-anchored MSA from a query and a reference collection.

Stages the genomes, selects a backbone reference, runs the chosen aligner
backend, and writes the resulting MSA-FASTA to the requested output path. The
leaf names in the MSA are the genome labels (filenames without extension); the
query's label is what the recombination step is later asked to analyse.
"""

from __future__ import annotations

import logging
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from ..aligners.base import AlignParams
from ..aligners.base import registry as aligner_registry
from ..core.errors import OutputError
from ..core.io import select_reference, stage_genomes

DEFAULT_ALIGNER = "progressivemauve"


@dataclass
class MsaParams:
    query: Path
    collection: Path
    output: Path
    aligner: str = DEFAULT_ALIGNER
    reference: str | None = None
    query_as_backbone: bool = False
    threads: int = 1
    extra: dict = field(default_factory=dict)


def build_msa(params: MsaParams, logger: logging.Logger) -> Path:
    """Run the MSA build described by ``params`` and return the output path."""
    aligner = aligner_registry.create(params.aligner)
    versions = aligner.preflight()
    logger.info("Using aligner '%s' (%s)", params.aligner, _format_versions(versions))

    output = Path(params.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="recomfi_msa_", dir=output.parent) as tmp:
        tmpdir = Path(tmp)
        genomes_dir = tmpdir / "genomes"
        genomes, query_staged = stage_genomes(
            params.query, params.collection, genomes_dir, logger
        )
        reference = select_reference(
            genomes, query_staged, params.query_as_backbone, params.reference
        )
        logger.info("Reference (backbone): %s", reference.stem)

        align_params = AlignParams(
            threads=params.threads,
            reference=reference,
            extra=dict(params.extra),
        )
        result = aligner.align(genomes, reference, tmpdir / "align", align_params, logger)

        msa = result.msa_fasta
        if not (msa.exists() and msa.stat().st_size > 0):
            raise OutputError(
                f"Aligner '{params.aligner}' did not produce a non-empty MSA. "
                "Verify the backend and its dependencies are installed."
            )
        # Move the MSA out before the temporary directory is removed.
        Path(msa).replace(output)

    logger.info("MSA written to %s", output)
    return output


def _format_versions(versions: dict[str, str]) -> str:
    if not versions:
        return "version unknown"
    return ", ".join(f"{name} {ver}" for name, ver in versions.items())
