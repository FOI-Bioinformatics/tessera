"""Aligner adapters, discovered via the ``recomfi.aligners`` entry-point group."""

from __future__ import annotations

from .base import Aligner, AlignParams, AlignResult, OutputKind, registry

__all__ = ["Aligner", "AlignParams", "AlignResult", "OutputKind", "registry"]
