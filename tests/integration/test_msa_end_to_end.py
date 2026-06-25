"""End-to-end MSA build + recomb on the bundled example data.

Parametrised over the supported MSA aligners; each case is skipped when its
binary is not installed.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import pytest

from recomfi.msa.build import MsaParams, build_msa
from recomfi.recomb.run import RecombParams, run_recomb

_LOG = logging.getLogger("recomfi.test")

pytestmark = pytest.mark.requires_binary

EXAMPLE = Path(__file__).resolve().parents[2] / "example_data"

# Aligner key -> the executable that must be on PATH for that aligner to run.
ALIGNER_BINARIES = {
    "progressivemauve": "progressiveMauve",
    "sibeliaz": "sibeliaz",
}


@pytest.mark.parametrize("aligner", sorted(ALIGNER_BINARIES))
def test_msa_then_recomb_on_example_data(aligner: str, tmp_path: Path) -> None:
    if shutil.which(ALIGNER_BINARIES[aligner]) is None:
        pytest.skip(f"{ALIGNER_BINARIES[aligner]} not installed")

    query = EXAMPLE / "cowpox_with_variolaInsert.fasta.gz"
    collection = EXAMPLE / "collection"
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
        RecombParams(msa=msa_out, output=out, query="cowpox_with_variolaInsert"),
        _LOG,
    )
    winners = (out / "window_winners.tsv").read_text().splitlines()[1:]
    ranked = [line.split("\t")[0] for line in winners]
    # The synthetic cowpox+variola query should rank these two backbones highest.
    assert "cowpox_KC813504" in ranked[:3]
    assert "variola" in ranked[:3]
    # a recombinant region should be called over the variola insert
    assert (out / "recombination_regions.tsv").exists()
