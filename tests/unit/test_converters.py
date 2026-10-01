"""MAF -> reference-anchored FASTA projection."""

from __future__ import annotations

from pathlib import Path

from tessera.converters.maf_to_fasta import maf_to_fasta


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
    # Full reference length is 8 (srcSize); blocks cover positions 0-6, so the
    # uncovered reference position 7 is a trailing gap for every sequence.
    # ref: ACGT (block1, pos 0-3) + ACG (block2, pos 4-6) + gap (pos 7)
    assert seqs["ref"] == "ACGTACG-"
    # qry projected onto ref columns: ACGA (block1) + ACG (block2, gap col dropped)
    assert seqs["qry"] == "ACGAACG-"
    assert len(seqs["ref"]) == len(seqs["qry"]) == 8


def test_maf_reverse_strand_reference_is_placed_forward(tmp_path: Path) -> None:
    # A block whose reference row is on the '-' strand must be reverse-complemented
    # into forward-reference coordinates, not kept reversed.
    maf = tmp_path / "rc.maf"
    # reference srcSize 8; this block covers forward positions 4-7.
    #   '-' strand, start=0 size=4 -> fstart = 8 - 0 - 4 = 4
    #   ref text on '-' strand "ACGT" -> forward revcomp "ACGT"
    #   qry text "AAGT" -> forward revcomp "ACTT"
    maf.write_text(
        "a\n"
        "s ref.chr 0 4 - 8 ACGT\n"
        "s qry.chr 0 4 + 4 AAGT\n"
    )
    out = maf_to_fasta(maf, "ref", tmp_path / "msa.fasta")
    seqs = _read_fasta(out)
    assert len(seqs["ref"]) == 8
    # forward positions 4-7 carry the reverse-complemented block; 0-3 are gaps
    assert seqs["ref"] == "----ACGT"
    assert seqs["qry"] == "----ACTT"


def test_maf_multicontig_reference_concatenates_contigs(tmp_path: Path) -> None:
    # A reference split across two contigs (ref.c1, ref.c2) must be laid out as a
    # concatenation, not collapsed into one contig's coordinate space. Each contig
    # is 4 bp, so the reference coordinate system is 8 bp wide.
    maf = tmp_path / "mc.maf"
    maf.write_text(
        "a\n"
        "s ref.c1 0 4 + 4 ACGT\n"
        "s qry.x 0 4 + 8 ACGT\n"
        "\n"
        "a\n"
        "s ref.c2 0 4 + 4 TTTT\n"
        "s qry.x 4 4 + 8 TTTT\n"
    )
    seqs = _read_fasta(maf_to_fasta(maf, "ref", tmp_path / "msa.fasta"))
    assert len(seqs["ref"]) == len(seqs["qry"]) == 8
    assert seqs["ref"] == "ACGTTTTT"  # c1 at 0-3, c2 at 4-7
    assert seqs["qry"] == "ACGTTTTT"


def test_maf_name_map_relabels_to_real_stems(tmp_path: Path) -> None:
    # Mirrors the cactus fix: dot-sanitised MAF sample names map back to the real
    # genome stems (which may themselves contain dots) for the recomb lookup.
    maf = tmp_path / "a.maf"
    maf.write_text(
        "a\n"
        "s cowpox_KC813504.chr 0 4 + 4 ACGT\n"
        "s sample_1.chr 0 4 + 4 ACGA\n"
    )
    name_map = {"cowpox_KC813504": "cowpox_KC813504", "sample_1": "sample.1"}
    out = maf_to_fasta(maf, "cowpox_KC813504", tmp_path / "msa.fasta", name_map=name_map)
    seqs = _read_fasta(out)
    assert set(seqs) == {"cowpox_KC813504", "sample.1"}  # dotted stem restored
    assert seqs["sample.1"] == "ACGA"


# --- multi-contig backbone layout and unplaced genomes ---------------------------

def test_maf_backbone_contigs_follow_file_order_when_given(tmp_path: Path) -> None:
    # The backbone FASTA lists contig_2 then contig_10. Sorted by name, contig_10 comes
    # first, which puts the backbone in a different order from its own file (and from
    # the minimap2 / mafft backends on the same input).
    maf = tmp_path / "order.maf"
    maf.write_text(
        "a\n"
        "s ref.contig_2 0 4 + 4 AAAA\n"
        "s qry.x 0 4 + 8 AAAA\n"
        "\n"
        "a\n"
        "s ref.contig_10 0 4 + 4 CCCC\n"
        "s qry.x 4 4 + 8 CCCC\n"
    )
    seqs = _read_fasta(maf_to_fasta(
        maf, "ref", tmp_path / "msa.fasta",
        ref_contigs=[("ref.contig_2", 4), ("ref.contig_10", 4)],
    ))
    assert seqs["ref"] == "AAAACCCC"
    assert seqs["qry"] == "AAAACCCC"


def test_maf_backbone_contig_without_a_block_keeps_its_columns(tmp_path: Path) -> None:
    # c2 has no homolog in any genome, so no MAF block mentions it. It is still part of
    # the backbone: dropping it would make the MSA 8 wide and shift c3 from 8-11 to 4-7.
    maf = tmp_path / "missing.maf"
    maf.write_text(
        "a\n"
        "s ref.c1 0 4 + 4 AAAA\n"
        "s qry.x 0 4 + 8 AAAA\n"
        "\n"
        "a\n"
        "s ref.c3 0 4 + 4 CCCC\n"
        "s qry.x 4 4 + 8 CCCC\n"
    )
    seqs = _read_fasta(maf_to_fasta(
        maf, "ref", tmp_path / "msa.fasta",
        ref_contigs=[("ref.c1", 4), ("ref.c2", 4), ("ref.c3", 4)],
    ))
    assert seqs["ref"] == "AAAA----CCCC"
    assert seqs["qry"] == "AAAA----CCCC"


def test_maf_rejects_a_backbone_contig_it_was_not_told_about(tmp_path: Path) -> None:
    import pytest

    from tessera.core.errors import OutputError

    maf = tmp_path / "extra.maf"
    maf.write_text("a\ns ref.c9 0 4 + 4 AAAA\ns qry.x 0 4 + 4 AAAA\n")
    with pytest.raises(OutputError, match="ref.c9"):
        maf_to_fasta(maf, "ref", tmp_path / "msa.fasta", ref_contigs=[("ref.c1", 4)])


def test_maf_genome_without_any_block_gets_an_all_gap_row(tmp_path: Path, caplog) -> None:
    # A panel member too divergent to be placed in any block used to have no row at all,
    # so the panel was silently one genome smaller than the collection.
    import logging

    maf = tmp_path / "vanish.maf"
    maf.write_text("a\ns r1 0 4 + 4 AAAA\ns q1 0 4 + 4 AAAT\n")
    name_map = {"r1": "ref", "q1": "qry", "p1": "divergent_panel"}
    # Not under the "tessera" logger: the CLI tests switch its propagation off, and
    # caplog only sees records that reach the root logger.
    log = logging.getLogger("maf_converter_test")
    with caplog.at_level(logging.WARNING, logger="maf_converter_test"):
        seqs = _read_fasta(maf_to_fasta(
            maf, "ref", tmp_path / "msa.fasta", name_map=name_map,
            expected=["ref", "qry", "divergent_panel"], logger=log,
        ))
    assert seqs == {"ref": "AAAA", "divergent_panel": "----", "qry": "AAAT"}
    assert "divergent_panel" in caplog.text


def test_maf_genome_aligned_only_away_from_the_backbone_is_named(tmp_path: Path, caplog) -> None:
    # SibeliaZ also emits blocks between non-backbone genomes. A genome seen only in such
    # blocks projects to an all-gap row just like one seen in no block, and must be named
    # in the same warning.
    import logging

    maf = tmp_path / "offref.maf"
    maf.write_text(
        "a\ns r1 0 4 + 4 AAAA\ns a1 0 4 + 4 AAAT\n\n"
        "a\ns a1 0 4 + 4 AAAT\ns b1 0 4 + 4 AAAT\n"
    )
    name_map = {"r1": "ref", "a1": "A", "b1": "B"}
    log = logging.getLogger("maf_converter_test")
    with caplog.at_level(logging.WARNING, logger="maf_converter_test"):
        seqs = _read_fasta(maf_to_fasta(
            maf, "ref", tmp_path / "msa.fasta", name_map=name_map,
            expected=["ref", "A", "B", "C"], logger=log,
        ))
    assert seqs == {"ref": "AAAA", "A": "AAAT", "B": "----", "C": "----"}
    assert "2 genome(s)" in caplog.text
    assert "B, C" in caplog.text
