"""The backbone row equals the backbone sequence, with the real aligner binaries.

The unit suite checks this invariant against fakes of each tool
(``tests/unit/test_backbone_invariant.py``). Here the installed backends align two
queries that each cover one contig of a two-contig backbone. Each case is skipped
when its binary is not on PATH.
"""

from __future__ import annotations

import logging
import random
import shutil
from pathlib import Path

import pytest

from tessera.aligners.base import AlignParams
from tessera.aligners.base import registry as aligner_registry

_LOG = logging.getLogger("tessera.test")

pytestmark = pytest.mark.requires_binary

_BINARIES = {
    "mafft": "mafft",
    "minimap2": "minimap2",
    "sibeliaz": "sibeliaz",
    "progressivemauve": "progressiveMauve",
}
# The MAF-based backends lay the backbone out from the blocks, so a backbone position no
# block covers is a gap in the backbone row (no other genome has a base there either).
# For them the invariant is that every base in the row is the backbone's own base at
# that position; the pairwise backends write the backbone itself.
_MAF_BASED = {"sibeliaz"}


def _mutate(seq: str, rng: random.Random, rate: float) -> str:
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def _read_rows(path: Path) -> dict[str, str]:
    rows: dict[str, str] = {}
    name = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:]
            rows[name] = ""
        elif name is not None:
            rows[name] += line
    return rows


@pytest.mark.parametrize("aligner", sorted(_BINARIES))
def test_backbone_row_is_the_backbone_with_the_real_tool(aligner: str, tmp_path: Path) -> None:
    if shutil.which(_BINARIES[aligner]) is None:
        pytest.skip(f"{_BINARIES[aligner]} not on PATH")
    rng = random.Random(7)
    c1 = "".join(rng.choice("ACGT") for _ in range(3000))
    c2 = "".join(rng.choice("ACGT") for _ in range(3000))
    ref = tmp_path / "ref.fasta"
    ref.write_text(f">c1\n{c1}\n>c2\n{c2}\n")
    qa = tmp_path / "qryA.fasta"
    qa.write_text(f">a1\n{_mutate(c1, rng, 0.01)}\n")
    qb = tmp_path / "qryB.fasta"
    qb.write_text(f">b1\n{_mutate(c2, rng, 0.01)}\n")

    result = aligner_registry.create(aligner).align(
        [ref, qa, qb], ref, tmp_path / "out", AlignParams(threads=2), _LOG
    )
    rows = _read_rows(result.msa_fasta)
    backbone = (c1 + c2).upper()
    assert {len(r) for r in rows.values()} == {len(backbone)}
    row = rows["ref"].upper()
    if aligner in _MAF_BASED:
        misplaced = [i for i, b in enumerate(row) if b != "-" and b != backbone[i]]
        assert misplaced == []
        assert row.count("-") < len(backbone) // 10  # most of the backbone is placed
    else:
        assert row == backbone
