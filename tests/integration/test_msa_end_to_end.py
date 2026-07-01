"""End-to-end MSA build + recomb on the bundled example data.

Raw genomes are derived on the fly from the example alignment (each record ungapped),
then re-aligned by each supported aligner and scanned for recombination. Parametrised
over the aligners; each case is skipped when its binary is not installed.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import pytest

from tessera.msa.build import MsaParams, build_msa
from tessera.recomb.run import RecombParams, run_recomb

_LOG = logging.getLogger("tessera.test")

pytestmark = pytest.mark.requires_binary

EXAMPLE = Path(__file__).resolve().parents[2] / "example_data"

# Aligner key -> the executable that must be on PATH for that aligner to run. The
# bundled example data is small synthetic genomes, so mafft (the default per-genome
# aligner) is exercised here; the whole-bacterial-genome aligners (progressiveMauve,
# sibeliaz) and the long-genome minimap2 path suit larger inputs and are covered by
# the opt-in validation harnesses instead.
ALIGNER_BINARIES = {
    "mafft": "mafft",
}


def _derive_genomes(msa: Path, out_dir: Path) -> Path:
    """Split an alignment into raw (ungapped) genomes: a ``query.fasta`` plus a
    ``collection/`` of the reference genomes. Returns the query path."""
    records: dict[str, list[str]] = {}
    label: str | None = None
    for line in msa.read_text().splitlines():
        if line.startswith(">"):
            label = line[1:].strip().split()[0]
            records[label] = []
        elif label is not None:
            records[label].append(line.strip())
    seqs = {k: "".join(v).replace("-", "").upper() for k, v in records.items()}

    collection = out_dir / "collection"
    collection.mkdir(parents=True, exist_ok=True)
    query = out_dir / "query.fasta"
    query.write_text(f">query\n{seqs['query']}\n")
    for name, seq in seqs.items():
        if name != "query":
            (collection / f"{name}.fasta").write_text(f">{name}\n{seq}\n")
    return query


@pytest.mark.parametrize("aligner", sorted(ALIGNER_BINARIES))
def test_msa_then_recomb_on_example_data(aligner: str, tmp_path: Path) -> None:
    if shutil.which(ALIGNER_BINARIES[aligner]) is None:
        pytest.skip(f"{ALIGNER_BINARIES[aligner]} not installed")

    # The divergent example: a parent_A backbone carrying a parent_B insert, parents
    # ~11 % apart -- a robust mosaic that any aligner should let the scan recover.
    query = _derive_genomes(EXAMPLE / "divergent.msa.fasta", tmp_path)
    collection = tmp_path / "collection"
    msa_out = tmp_path / "msa.fasta"

    build_msa(
        MsaParams(
            query=query, collection=collection, output=msa_out,
            aligner=aligner, threads=1,
        ),
        _LOG,
    )
    assert msa_out.exists() and msa_out.stat().st_size > 0

    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=msa_out, output=out, query="query",
                     window_size=300, window_step=30),
        _LOG,
    )
    winners = (out / "window_winners.tsv").read_text().splitlines()[1:]
    ranked = [line.split("\t")[0] for line in winners]
    # The query's two parents rank highest: parent_A (backbone) and parent_B (donor).
    assert "parent_A" in ranked[:3]
    assert "parent_B" in ranked[:3]
    # A recombinant region should be called, with parent_B the donor over the insert.
    regions = (out / "recombination_regions.tsv").read_text().splitlines()[1:]
    assert any("parent_B" in row for row in regions)
