"""Real-data IO: genome staging and header parsing (no external binaries).

Exercises the staging and header-parsing paths against a self-contained fixture that
mirrors the awkward inputs seen in practice: a gzipped collection whose records carry
mixed header styles (BV-BRC ``accn|...`` pipes and NCBI ``ACCESSION.version``), plus a
multi-contig assembly query. The fixture is built in a temp dir so the test does not
depend on any shipped data.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import pytest

from tessera.aligners.sibeliaz import _build_seqid_map
from tessera.core.io import select_reference, stage_genomes

_LOG = logging.getLogger("tessera.test")

# filename stem -> (record header, sequence). The filename is the label; the header's
# first token is the seqid the sibeliaz map must resolve back to that label.
_COLLECTION = {
    "virus_a": (">accn|NC_111111 Virus A, complete genome", "acgtacgtacgtacgtaccg"),
    "virus_b": (">AB000001.1 Virus B strain X", "acgtacgtacgtacgtaccg"),
    "virus_c": (">accn|NC_333333 Virus C, complete genome", "acgtacgtacgtacgtacca"),
}
_QUERY_CONTIGS = 4  # an assembly query split across scaffolds
EXPECTED_LABELS = {"virus_a", "virus_b", "virus_c", "assembly"}


def _write_gz(path: Path, text: str) -> None:
    with gzip.open(path, "wt") as fo:
        fo.write(text)


@pytest.fixture
def staged(tmp_path: Path):
    collection = tmp_path / "collection"
    collection.mkdir()
    for label, (header, seq) in _COLLECTION.items():
        _write_gz(collection / f"{label}.fasta.gz", f"{header}\n{seq}\n")
    query = tmp_path / "assembly.fasta.gz"
    _write_gz(query, "".join(f">scaffold_{i}\nacgtacgtac\n" for i in range(_QUERY_CONTIGS)))
    genomes, staged_query = stage_genomes(query, collection, tmp_path / "genomes", _LOG)
    return genomes, staged_query


def test_stage_genomes_decompresses_and_names(staged) -> None:
    genomes, query = staged
    # three references plus the query, all normalised to <label>.fasta
    assert len(genomes) == 4
    assert {g.stem for g in genomes} == EXPECTED_LABELS
    assert all(g.suffix == ".fasta" for g in genomes)
    # gzip content was decompressed to plain FASTA text
    assert query.read_text(errors="ignore").startswith(">")


def test_query_contig_count_preserved(staged) -> None:
    _, query = staged
    headers = [ln for ln in query.read_text().splitlines() if ln.startswith(">")]
    assert len(headers) == _QUERY_CONTIGS


def test_reference_selection(staged) -> None:
    genomes, query = staged
    # default picks a collection genome, not the query
    default_ref = select_reference(genomes, query, query_as_backbone=False)
    assert default_ref != query
    assert default_ref.stem in EXPECTED_LABELS - {"assembly"}
    # --query-as-backbone picks the query itself
    assert select_reference(genomes, query, query_as_backbone=True) == query
    # an explicit reference resolves by label
    assert select_reference(genomes, query, False, "virus_b").stem == "virus_b"


def test_sibeliaz_seqid_map_handles_mixed_headers(staged) -> None:
    genomes, _ = staged
    name_map = _build_seqid_map(genomes)
    # BV-BRC pipe header, NCBI accession.version header, and the scaffold headers
    assert name_map["accn|NC_111111"] == "virus_a"
    assert name_map["AB000001.1"] == "virus_b"
    assert name_map["accn|NC_333333"] == "virus_c"
    for i in range(_QUERY_CONTIGS):
        assert name_map[f"scaffold_{i}"] == "assembly"
