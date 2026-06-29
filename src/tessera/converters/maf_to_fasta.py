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

from dataclasses import dataclass
from pathlib import Path

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
) -> Path:
    """Project a MAF onto ``reference`` coordinates as an MSA-FASTA.

    ``name_map`` maps a MAF source name (often a contig/sequence ID) to a genome
    label. When given, sequences are grouped by genome (needed for tools like
    SibeliaZ whose MAF uses raw sequence IDs); otherwise the ``genome.contig``
    convention is assumed (Cactus via hal2maf).

    ``exclude`` drops genomes by label, e.g. ``{"_MINIGRAPH_"}`` to remove the
    Minigraph-Cactus backbone pseudo-genome so it is not emitted as a taxon.
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
    ref_contigs: dict[str, int] = {}  # source name -> contig length
    species: set[str] = set()
    for block in blocks:
        for row in block:
            label = genome_of(row.name)
            species.add(label)
            if label == ref_key:
                ref_contigs.setdefault(row.name, row.src_size)
    species -= exclude

    if not ref_contigs:
        raise ValueError(
            f"MAF projection onto reference '{ref_key}' found no reference rows. "
            f"Check that the reference label matches the MAF/name_map sequence names."
        )

    ref_offsets: dict[str, int] = {}
    ref_length = 0
    for name in sorted(ref_contigs):
        ref_offsets[name] = ref_length
        ref_length += ref_contigs[name]

    species.discard(ref_key)
    ordered_species = [ref_key, *sorted(species)]
    out: dict[str, bytearray] = {s: bytearray(b"-" * ref_length) for s in ordered_species}

    for block in blocks:
        ref_row = next((r for r in block if genome_of(r.name) == ref_key), None)
        if ref_row is None:
            continue
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
