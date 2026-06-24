"""Shared pytest fixtures."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

EXAMPLE_DATA = Path(__file__).resolve().parents[1] / "example_data"


@pytest.fixture
def logger() -> logging.Logger:
    log = logging.getLogger("recomfi.test")
    log.addHandler(logging.NullHandler())
    return log


@pytest.fixture
def example_data() -> Path:
    return EXAMPLE_DATA


def write_fasta(path: Path, records: dict[str, str]) -> Path:
    """Write a simple (aligned) FASTA from a label -> sequence mapping."""
    with open(path, "w") as fo:
        for name, seq in records.items():
            fo.write(f">{name}\n{seq}\n")
    return path
