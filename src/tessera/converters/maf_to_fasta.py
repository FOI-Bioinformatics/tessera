"""MAF -> reference-anchored MSA-FASTA.

Projects a MAF (SibeliaZ, MULTIZ, or Cactus via hal2maf) onto the coordinate
system of a chosen reference, placing every block at its **true reference
coordinate** in a full-reference-length alignment -- the same model as
:mod:`tessera.converters.xmfa_to_fasta`. Each reference base maps to one column,
positions no block covers are left as gaps, and blocks whose reference row is on
the ``-`` strand are reverse-complemented into forward-reference orientation.

This matters: SibeliaZ routinely reports a large fraction of blocks with the
reference on the ``-`` strand, so a converter that merely concatenated
reference-covered columns in start order (ignoring strand) scrambled coordinates
by tens of kilobases. Placing blocks at their forward-reference coordinate keeps
the output coordinates equal to the reference's own, so the recombination scan
reports positions that line up with the reference genome.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from ..core.errors import OutputError
from ..core.io import write_fasta_record

_COMPLEMENTS = bytes.maketrans(
    b"acgtrymkbdhvACGTRYMKBDHV", b"tgcayrkmvhdbTGCAYRKMVHDB"
)


def _revcomp(text: str) -> str:
    """Reverse-complement an aligned string; gaps stay gaps."""
    return text.encode("latin-1").translate(_COMPLEMENTS).decode("latin-1")[::-1]


@dataclass
class _Row:
    name: str
    start: int
    size: int
    strand: str
    src_size: int
    text: str


def _species(name: str) -> str:
    """MAF source names are often ``genome.contig``; key on the genome part."""
    return name.split(".")[0]


def maf_to_fasta(
    maf_path: str | Path,
    reference: str,
    out_path: str | Path,
    name_map: dict[str, str] | None = None,
    exclude: set[str] | None = None,
    ref_contigs: Sequence[tuple[str, int]] | None = None,
    expected: Sequence[str] | None = None,
    logger: logging.Logger | None = None,
) -> Path:
    """Project a MAF onto ``reference`` coordinates as an MSA-FASTA.

    ``name_map`` maps a MAF source name (often a contig/sequence ID) to a genome
    label. When given, sequences are grouped by genome (needed for tools like
    SibeliaZ whose MAF uses raw sequence IDs); otherwise the ``genome.contig``
    convention is assumed (Cactus via hal2maf).

    ``exclude`` drops genomes by label, e.g. ``{"_MINIGRAPH_"}`` to remove the
    Minigraph-Cactus backbone pseudo-genome so it is not emitted as a taxon.

    ``ref_contigs`` gives the backbone's contigs as ``(MAF source name, length)`` in
    the order of its FASTA file. The MAF alone cannot supply this: it lists only
    contigs that fall in some block, in no particular order. With it, the backbone is
    laid out as its file is, and a contig no block covers keeps its (all-gap) columns
    instead of vanishing and shifting everything after it. Without it the contigs seen
    in the MAF are laid out in sorted-name order.

    ``expected`` lists every genome label that should have a row. A genome the aligner
    placed in no block is absent from the MAF; it is written as an all-gap row. Such a
    genome, and one aligned only in blocks that do not include the backbone, is named in
    a warning on ``logger``, so the panel is never silently smaller than the collection.
    """
    maf_path = Path(maf_path)
    out_path = Path(out_path)
    exclude = exclude or set()

    def genome_of(src: str) -> str:
        if name_map is not None:
            return name_map.get(src, name_map.get(src.split(".")[0], _species(src)))
        return _species(src)

    ref_key = name_map.get(reference, reference) if name_map else _species(reference)

    blocks = list(_iter_blocks(maf_path))

    # A multi-contig reference is laid out as a concatenation of its contigs:
    # each distinct reference source name is one contig (length = its src_size),
    # placed at a cumulative offset. Collapsing them into a single contig's
    # coordinate space would make later contigs overwrite earlier ones.
    seen_contigs: dict[str, int] = {}  # source name -> contig length
    species: set[str] = set()
    placed: set[str] = set()  # genomes sharing at least one block with the backbone
    for block in blocks:
        labels = [genome_of(row.name) for row in block]
        species.update(labels)
        if ref_key in labels:
            placed.update(labels)
        for row, label in zip(block, labels, strict=True):
            if label == ref_key:
                seen_contigs.setdefault(row.name, row.src_size)
    species -= exclude

    if not seen_contigs:
        raise ValueError(
            f"MAF projection onto reference '{ref_key}' found no reference rows. "
            f"Check that the reference label matches the MAF/name_map sequence names."
        )

    if ref_contigs is not None:
        unknown = sorted(set(seen_contigs) - {name for name, _ in ref_contigs})
        if unknown:
            raise OutputError(
                f"{maf_path} aligns backbone sequence(s) {', '.join(unknown)} that are "
                f"not in the backbone '{ref_key}' as staged. The aligner's sequence "
                "names do not match the input FASTA."
            )
        layout = list(ref_contigs)
    else:
        layout = sorted(seen_contigs.items())
    ref_offsets: dict[str, int] = {}
    ref_length = 0
    for name, length in layout:
        ref_offsets[name] = ref_length
        ref_length += length

    # A genome in no block at all, and one aligned only in blocks that lack the backbone,
    # both project to an all-gap row; name both.
    unplaced = sorted((set(expected or ()) | species) - placed - exclude - {ref_key})
    if unplaced and logger is not None:
        logger.warning(
            "%d genome(s) share no alignment block with the backbone '%s' and are "
            "written as all-gap rows: %s. They contribute nothing to the scan; they "
            "may be too divergent for this aligner.",
            len(unplaced), ref_key, ", ".join(unplaced),
        )
    species.discard(ref_key)
    ordered_species = [ref_key, *sorted(species | set(unplaced))]
    out: dict[str, bytearray] = {s: bytearray(b"-" * ref_length) for s in ordered_species}

    for block in blocks:
        ref_row = next((r for r in block if genome_of(r.name) == ref_key), None)
        if ref_row is None:
            continue
        # Every row of a block has the same aligned width. A shorter one means the file
        # was cut off mid-write (a full disk, a killed aligner); indexing into it below
        # would fail with a bare IndexError several frames from the file's name.
        for row in block:
            if len(row.text) != len(ref_row.text):
                raise OutputError(
                    f"Truncated alignment block in {maf_path}: row '{row.name}' has "
                    f"{len(row.text)} column(s), the backbone row '{ref_row.name}' has "
                    f"{len(ref_row.text)}. The aligner's output looks incomplete; check "
                    "that it finished."
                )
        contig_offset = ref_offsets[ref_row.name]
        if ref_row.strand == "-":
            # Reverse-complement the whole block into forward-reference orientation.
            fstart = ref_row.src_size - ref_row.start - ref_row.size
            rows = [(genome_of(r.name), _revcomp(r.text)) for r in block]
            ref_text = _revcomp(ref_row.text)
        else:
            fstart = ref_row.start
            rows = [(genome_of(r.name), r.text) for r in block]
            ref_text = ref_row.text
        fstart += contig_offset

        pos = fstart
        for col, ref_char in enumerate(ref_text):
            if ref_char == "-":
                continue  # insertion relative to the reference: dropped
            if 0 <= pos < ref_length:
                for label, text in rows:
                    if label in out and out[label][pos] == ord("-"):
                        ch = text[col]
                        if ch != "-":
                            out[label][pos] = ord(ch)
            pos += 1

    if not any(c != ord("-") for c in out[ref_key]):
        raise ValueError(
            f"MAF projection onto reference '{ref_key}' produced an empty alignment "
            f"(no usable blocks). Check that the reference label matches the "
            f"MAF/name_map sequence names, or that the genomes share alignable regions."
        )

    with open(out_path, "w") as fo:
        for s in ordered_species:
            write_fasta_record(fo, s, out[s].decode("latin-1"))
    return out_path


def _iter_blocks(maf_path: Path):
    block: list[_Row] = []
    with open(maf_path) as fo:
        for line in fo:
            if line.startswith("a"):
                if block:
                    yield block
                block = []
            elif line.startswith("s"):
                parts = line.split()
                # s src start size strand srcSize text
                if len(parts) >= 7:
                    block.append(_Row(
                        name=parts[1], start=int(parts[2]), size=int(parts[3]),
                        strand=parts[4], src_size=int(parts[5]), text=parts[6],
                    ))
            elif line.strip() == "" and block:
                yield block
                block = []
        if block:
            yield block
