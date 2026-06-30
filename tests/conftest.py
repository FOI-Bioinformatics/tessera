"""Shared pytest fixtures."""

from __future__ import annotations

import logging
import random
from pathlib import Path

import pytest

EXAMPLE_DATA = Path(__file__).resolve().parents[1] / "example_data"


@pytest.fixture
def logger() -> logging.Logger:
    log = logging.getLogger("tessera.test")
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


def recombinant_msa(tmp_path: Path, *, recombinant: bool) -> Path:
    """Synthetic four-record MSA: an A backbone with an optional B insert.

    Shared by the triplet callers' tests (3seq, maxchi, bootscan). When
    ``recombinant`` is set, sites [2000, 4000) of the query are copied from
    parent B so the caller should recover an A/B recombinant.
    """
    rng = random.Random(7)
    base = "".join(rng.choice("ACGT") for _ in range(6000))

    def mut(seq: str, frac: float) -> str:
        chars = list(seq)
        for i in range(len(chars)):
            if rng.random() < frac:
                chars[i] = rng.choice("ACGT")
        return "".join(chars)

    a, b, other = mut(base, 0.03), mut(base, 0.03), mut(base, 0.10)
    query = list(a)
    if recombinant:
        query[2000:4000] = list(b[2000:4000])  # A backbone, B insert
    return write_fasta(tmp_path / "m.fasta",
                       {"query": "".join(query), "A": a, "B": b, "other": other})
