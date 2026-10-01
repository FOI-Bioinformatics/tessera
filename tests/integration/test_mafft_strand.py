"""MAFFT backend on reverse-strand input (needs the mafft binary).

An assembly carries no promise about strand: a whole genome, or one contig of a draft,
may be the reverse complement of the backbone. The row must still be a real alignment.
"""

from __future__ import annotations

import logging
import random
import shutil
from pathlib import Path

import pytest

from tessera.aligners.base import AlignParams
from tessera.aligners.mafft import MafftAligner
from tessera.core.io import read_fasta

_LOG = logging.getLogger("tessera.test")

pytestmark = pytest.mark.requires_binary

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _genomes() -> tuple[str, str]:
    """A 6 kb random reference and a query 3 % diverged from it."""
    rng = random.Random(1)
    ref = "".join(rng.choice("ACGT") for _ in range(6000))
    qry = "".join(
        rng.choice([b for b in "ACGT" if b != c]) if rng.random() < 0.03 else c for c in ref
    )
    return ref, qry


def _identity_by_third(ref_row: str, qry_row: str) -> list[float]:
    out = []
    for lo in range(0, len(ref_row), 2000):
        pairs = [
            (a, b) for a, b in zip(ref_row[lo:lo + 2000].upper(), qry_row[lo:lo + 2000].upper(),
                                   strict=True)
            if a in "ACGT" and b in "ACGT"
        ]
        out.append(sum(a == b for a, b in pairs) / len(pairs) if pairs else 0.0)
    return out


def _align(tmp_path: Path, query_fasta: str) -> tuple[str, str]:
    ref, _ = _genomes()
    (tmp_path / "ref.fasta").write_text(f">ref\n{ref}\n")
    (tmp_path / "qry.fasta").write_text(query_fasta)
    result = MafftAligner().align(
        [tmp_path / "ref.fasta", tmp_path / "qry.fasta"], tmp_path / "ref.fasta",
        tmp_path / "out", AlignParams(threads=1), _LOG,
    )
    msa = dict(read_fasta(result.msa_fasta))
    return msa["ref"], msa["qry"]


@pytest.fixture(autouse=True)
def _needs_mafft() -> None:
    if shutil.which("mafft") is None:
        pytest.skip("mafft not installed")


def test_reverse_complemented_query_aligns(tmp_path: Path) -> None:
    _, qry = _genomes()
    ref_row, qry_row = _align(tmp_path, f">q\n{_revcomp(qry)}\n")
    assert len(ref_row) == len(qry_row) == 6000
    assert all(identity > 0.95 for identity in _identity_by_third(ref_row, qry_row))


def test_draft_query_with_one_reversed_contig_aligns(tmp_path: Path) -> None:
    _, qry = _genomes()
    contigs = f">c1\n{qry[:2000]}\n>c2\n{_revcomp(qry[2000:4000])}\n>c3\n{qry[4000:]}\n"
    ref_row, qry_row = _align(tmp_path, contigs)
    assert len(ref_row) == len(qry_row) == 6000
    assert all(identity > 0.95 for identity in _identity_by_third(ref_row, qry_row))
