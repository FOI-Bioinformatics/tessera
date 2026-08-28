"""Build a reference-anchored MSA from a query and a reference collection.

Stages the genomes, selects a backbone reference, runs the chosen aligner
backend, and writes the resulting MSA-FASTA to the requested output path. The
leaf names in the MSA are the genome labels (filenames without extension); the
query's label is what the recombination step is later asked to analyse.

Alongside the MSA a provenance sidecar records which aligner built it, at which
version, with which arguments and backbone. An alignment is an interpretation of
the input, not a fact about it, so a result that cannot name the aligner that
produced it cannot be reproduced by a reader. :func:`provenance_path` gives the
sidecar's location for a given MSA; ``recomb`` picks it up automatically.
"""

from __future__ import annotations

import json
import logging
import tempfile
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from .. import __version__
from ..aligners.base import AlignParams
from ..aligners.base import registry as aligner_registry
from ..core.errors import OutputError
from ..core.io import select_reference, stage_genomes

DEFAULT_ALIGNER = "sibeliaz"


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


def provenance_path(msa: Path) -> Path:
    """Where the provenance sidecar for ``msa`` lives (``<msa stem>.provenance.json``)."""
    return Path(msa).with_suffix(".provenance.json")


def build_msa(params: MsaParams, logger: logging.Logger) -> Path:
    """Run the MSA build described by ``params`` and return the output path.

    Also writes the provenance sidecar beside the MSA; see :func:`provenance_path`.
    """
    aligner = aligner_registry.create(params.aligner)
    versions = aligner.preflight()
    logger.info("Using aligner '%s' (%s)", params.aligner, _format_versions(versions))

    output = Path(params.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="tessera_msa_", dir=output.parent) as tmp:
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
        _write_provenance(
            output, params=params, versions=versions,
            backbone=reference.stem, n_genomes=len(genomes),
        )

    logger.info("MSA written to %s", output)
    return output


def _write_provenance(
    output: Path,
    *,
    params: MsaParams,
    versions: dict[str, str],
    backbone: str,
    n_genomes: int,
) -> Path:
    """Record how this alignment was built, beside the alignment itself."""
    record = {
        "tessera version": __version__,
        "date (UTC)": datetime.now(UTC).strftime("%Y-%m-%d %H:%M:%S"),
        "aligner": params.aligner,
        "aligner versions": dict(versions),
        "aligner args": dict(params.extra),
        "backbone": backbone,
        "genomes": n_genomes,
        "query": str(params.query),
        "collection": str(params.collection),
        "threads": params.threads,
    }
    path = provenance_path(output)
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return path


def read_provenance(msa: Path) -> dict | None:
    """The provenance recorded for ``msa``, or ``None`` when there is no usable sidecar.

    Absent for an alignment a user built themselves, which is expected -- the caller
    reports what it has rather than requiring the sidecar.
    """
    try:
        loaded = json.loads(provenance_path(msa).read_text())
    except (OSError, json.JSONDecodeError):
        return None
    return loaded if isinstance(loaded, dict) else None


def _format_versions(versions: dict[str, str]) -> str:
    if not versions:
        return "version unknown"
    return ", ".join(f"{name} {ver}" for name, ver in versions.items())
