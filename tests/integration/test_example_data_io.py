"""Real-data IO on the bundled example_data/ (no external binaries).

Exercises the staging and header-parsing paths against the actual orthopoxvirus
files: a gzipped, lowercase, mixed-header collection plus an 8-contig query.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from tessera.aligners.sibeliaz import _build_seqid_map
from tessera.core.io import select_reference, stage_genomes

_LOG = logging.getLogger("tessera.test")

EXAMPLE = Path(__file__).resolve().parents[2] / "example_data"
COLLECTION = EXAMPLE / "collection"
QUERY = EXAMPLE / "cowpox_with_variolaInsert.fasta.gz"

EXPECTED_LABELS = {
    "camelpox", "cowpox", "cowpox_KC813504", "monkeypox",
    "taterapox", "vaccinia", "variola", "cowpox_with_variolaInsert",
}


@pytest.fixture
def staged(tmp_path: Path):
    genomes, query = stage_genomes(QUERY, COLLECTION, tmp_path / "genomes", _LOG)
    return genomes, query


def test_stage_genomes_decompresses_and_names(staged) -> None:
    genomes, query = staged
    # seven references plus the query, all normalised to <label>.fasta
    assert len(genomes) == 8
    assert {g.stem for g in genomes} == EXPECTED_LABELS
    assert all(g.suffix == ".fasta" for g in genomes)
    # gzip content was decompressed to plain FASTA text
    assert query.read_text(errors="ignore").startswith(">")


def test_query_contig_count_preserved(staged) -> None:
    _, query = staged
    headers = [ln for ln in query.read_text().splitlines() if ln.startswith(">")]
    assert len(headers) == 8  # the eight assembly scaffolds


def test_reference_selection(staged) -> None:
    genomes, query = staged
    # default picks a collection genome, not the query
    default_ref = select_reference(genomes, query, query_as_backbone=False)
    assert default_ref != query
    assert default_ref.stem in EXPECTED_LABELS - {"cowpox_with_variolaInsert"}
    # --query-as-backbone picks the query itself
    assert select_reference(genomes, query, query_as_backbone=True) == query
    # an explicit reference resolves by label
    assert select_reference(genomes, query, False, "cowpox_KC813504").stem == "cowpox_KC813504"


def test_sibeliaz_seqid_map_handles_mixed_headers(staged) -> None:
    genomes, _ = staged
    name_map = _build_seqid_map(genomes)
    # BV-BRC pipe header, NCBI accession header, and the eight scaffold headers
    assert name_map["accn|NC_003391"] == "camelpox"
    assert name_map["KC813504.1"] == "cowpox_KC813504"
    assert name_map["accn|NC_001611"] == "variola"
    for i in range(8):
        assert name_map[f"scaffold_{i}"] == "cowpox_with_variolaInsert"
