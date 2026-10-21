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

from ..core.errors import UserInputError
from ..core.plugins import Registry, ToolCapabilities, preflight

registry: Registry[Aligner] = Registry("tessera.aligners")

_TRUTHY = frozenset({"1", "true", "yes", "on"})
_FALSY = frozenset({"0", "false", "no", "off", ""})


def split_thread_budget(
    threads: int, tasks: int, job_limit: int | None = None
) -> tuple[int, int]:
    """Divide a ``--threads`` budget between concurrent tool runs.

    Returns ``(jobs, threads_per_job)``: at most one job per task, never more jobs than
    threads, and never fewer than one of either. ``job_limit`` caps the number of jobs
    further (a memory bound, say). ``jobs * threads_per_job`` never exceeds the budget,
    so ``--threads 16`` starts at most sixteen threads in all rather than sixteen
    processes of sixteen threads each.
    """
    budget = max(1, threads)
    jobs = max(1, min(budget, tasks))
    if job_limit is not None:
        jobs = max(1, min(jobs, job_limit))
    return jobs, budget // jobs


def warn_unplaced(
    rows: Sequence[tuple[str, str]], backbone: str, logger: logging.Logger
) -> None:
    """Name every genome whose row holds no aligned base.

    The pairwise backends write a genome they could not place against the backbone as
    an all-gap row, so the alignment keeps one row per input genome. Such a row
    contributes nothing to the scan; without a warning the panel is silently smaller
    than the collection. The MAF-based backends say the same in their converter.
    """
    unplaced = sorted(label for label, row in rows if not row.strip("-"))
    if unplaced:
        logger.warning(
            "%d genome(s) were not aligned to any part of the backbone '%s' and are "
            "written as all-gap rows: %s. They contribute nothing to the scan; they "
            "may be too divergent for this aligner.",
            len(unplaced), backbone, ", ".join(unplaced),
        )


class OutputKind(Enum):
    MSA_FASTA = "msa_fasta"


@dataclass
class AlignParams:
    threads: int = 1
    reference: Path | None = None
    extra: dict = field(default_factory=dict)

    def flag(self, key: str, default: bool = False) -> bool:
        """Interpret a boolean-style ``--aligner-arg`` extra (e.g. ``single=true``).

        Extras arrive as strings, so a plain truthiness test would read
        ``single=false`` as ``True``. Values are matched case-insensitively
        against a small true/false vocabulary; an unrecognised value raises
        :class:`UserInputError` rather than silently defaulting.
        """
        if key not in self.extra:
            return default
        value = str(self.extra[key]).strip().lower()
        if value in _TRUTHY:
            return True
        if value in _FALSY:
            return False
        raise UserInputError(
            f"--aligner-arg {key}={self.extra[key]!r} is not a boolean; "
            f"use one of true/false (or yes/no, 1/0)."
        )


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
