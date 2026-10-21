"""Genomes with whitespace in their sequence lines, through each pairwise aligner.

Needs the aligner binaries. A reference saved with a trailing space on every line used
to give a ragged MSA with mafft and, once the reader ignored the spaces, a silently
shifted one with minimap2 (which counts them as bases).
"""

from __future__ import annotations

import logging
import random
import shutil
from pathlib import Path

import pytest

from tessera.core.io import read_fasta
from tessera.msa.build import MsaParams, build_msa

_LOG = logging.getLogger("tessera.test")

pytestmark = pytest.mark.requires_binary


@pytest.mark.parametrize("aligner", ["mafft", "minimap2"])
def test_reference_with_trailing_spaces_aligns_in_register(aligner: str, tmp_path: Path) -> None:
    if shutil.which(aligner) is None:
        pytest.skip(f"{aligner} not installed")
    rng = random.Random(1)
    ref = "".join(rng.choice("ACGT") for _ in range(6000))
    qry = "".join(
        rng.choice([b for b in "ACGT" if b != c]) if rng.random() < 0.03 else c for c in ref
    )
    other = "".join(
        rng.choice([b for b in "ACGT" if b != c]) if rng.random() < 0.05 else c for c in ref
    )
    collection = tmp_path / "collection"
    collection.mkdir()
    spaced = "".join(ref[i:i + 60] + " \n" for i in range(0, len(ref), 60))
    (collection / "ref.fasta").write_text(f">r\n{spaced}")
    (collection / "other.fasta").write_text(f">o\n{other}\n")
    query = tmp_path / "query.fasta"
    query.write_text(f">q\n{qry}\n")

    msa = build_msa(
        MsaParams(query=query, collection=collection, output=tmp_path / "msa.fasta",
                  aligner=aligner, reference="ref", threads=1),
        _LOG,
    )
    rows = dict(read_fasta(msa))
    assert {len(row) for row in rows.values()} == {6000}
    tail_pairs = [
        (a, b) for a, b in zip(rows["ref"][3000:].upper(), rows["query"][3000:].upper(),
                               strict=True)
        if a in "ACGT" and b in "ACGT"
    ]
    assert len(tail_pairs) > 2500
    assert sum(a == b for a, b in tail_pairs) / len(tail_pairs) > 0.95
