"""Aligner interface.

An aligner takes genome FASTAs (and usually a reference) and produces a
canonical reference-anchored multiple-sequence-alignment FASTA. Raw tool formats
(XMFA, MAF, HAL) are normalized to MSA-FASTA by the converters package so the
downstream recombination scan is identical regardless of aligner.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

from ..core.plugins import Registry, ToolCapabilities, preflight

registry: Registry[Aligner] = Registry("recomfi.aligners")


class OutputKind(Enum):
    MSA_FASTA = "msa_fasta"


@dataclass
class AlignParams:
    threads: int = 1
    reference: Path | None = None
    extra: dict = field(default_factory=dict)


@dataclass
class AlignResult:
    msa_fasta: Path
    native_format: Path | None = None  # raw HAL/MAF/XMFA kept for provenance


class Aligner(ABC):
    capabilities: ToolCapabilities
    output_kind: OutputKind = OutputKind.MSA_FASTA

    def preflight(self) -> dict[str, str]:
        return preflight(self.capabilities)

    @abstractmethod
    def align(
        self,
        genomes: Sequence[Path],
        reference: Path | None,
        out_dir: Path,
        params: AlignParams,
        logger: logging.Logger,
    ) -> AlignResult:
        raise NotImplementedError
