"""MAF -> reference-anchored FASTA projection."""

from __future__ import annotations

from pathlib import Path

from recomfi.converters.maf_to_fasta import maf_to_fasta


def _read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            out[name] = ""
        elif name is not None:
            out[name] += line
    return out


def test_maf_projects_onto_reference_and_drops_ref_gaps(tmp_path: Path) -> None:
    # Two blocks. In block 2 the reference has a gap column that must be dropped,
    # so each sequence contributes only its reference-aligned positions.
    maf = tmp_path / "a.maf"
    maf.write_text(
        "a\n"
        "s ref.chr 0 4 + 8 ACGT\n"
        "s qry.chr 0 4 + 8 ACGA\n"
        "\n"
        "a\n"
        "s ref.chr 4 4 + 8 AC-G\n"
        "s qry.chr 4 4 + 8 ACTG\n"
    )
    out = maf_to_fasta(maf, "ref", tmp_path / "msa.fasta")
    seqs = _read_fasta(out)
    # ref keeps its non-gap columns: ACGT (block1) + ACG (block2) = ACGTACG
    assert seqs["ref"] == "ACGTACG"
    # qry projected onto ref columns: ACGA (block1) + ACG (block2, gap col dropped)
    assert seqs["qry"] == "ACGAACG"
    assert len(seqs["ref"]) == len(seqs["qry"])
