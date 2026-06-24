"""MAF -> reference-anchored MSA-FASTA.

Projects a MAF (SibeliaZ, MULTIZ, or Cactus via hal2maf) onto the coordinate
system of a chosen reference: within each alignment block, columns where the
reference has a gap are dropped, every other sequence contributes its
gap-trimmed row, and sequences absent from a block are padded with gaps. Blocks
are concatenated in reference order. The result has no insertions relative to
the reference -- a fixed-width MSA.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class _Row:
    name: str
    ref_start: int
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

    blocks = list(_iter_blocks(maf_path))
    # The reference is passed as a genome label (e.g. a filename stem), which may
    # itself contain dots. With a name_map the row names are contig IDs that map
    # to those same labels, so resolve the reference to a label WITHOUT
    # version-stripping: map it only if it is itself a contig key, otherwise use
    # it verbatim.
    ref_key = name_map.get(reference, reference) if name_map else _species(reference)

    species: set[str] = set()
    for block in blocks:
        for row in block:
            species.add(genome_of(row.name))
    species -= exclude
    species.discard(ref_key)
    ordered_species = [ref_key, *sorted(species)]

    # collect, per species, its concatenated aligned text in reference order
    pieces: dict[str, list[str]] = {s: [] for s in ordered_species}

    placed_blocks: list[tuple[int, dict[str, str]]] = []
    for block in blocks:
        ref_row = next((r for r in block if genome_of(r.name) == ref_key), None)
        if ref_row is None:
            continue
        keep = [i for i, ch in enumerate(ref_row.text) if ch != "-"]
        if not keep:
            continue
        width = len(keep)
        block_cols: dict[str, str] = {}
        for s in ordered_species:
            row = next((r for r in block if genome_of(r.name) == s), None)
            if row is None:
                block_cols[s] = "-" * width
            else:
                block_cols[s] = "".join(row.text[i] for i in keep)
        placed_blocks.append((ref_row.ref_start, block_cols))

    for _, block_cols in sorted(placed_blocks, key=lambda x: x[0]):
        for s in ordered_species:
            pieces[s].append(block_cols[s])

    width = sum(len(p) for p in pieces[ref_key]) if ref_key in pieces else 0
    if width == 0:
        raise ValueError(
            f"MAF projection onto reference '{ref_key}' produced a zero-length "
            f"alignment (no usable blocks). Check that the reference label matches "
            f"the MAF/name_map sequence names, or that the genomes share alignable "
            f"regions (whole-genome aligners need collinearity not present across "
            f"highly divergent inputs)."
        )

    with open(out_path, "w") as fo:
        for s in ordered_species:
            seq = "".join(pieces[s])
            fo.write(f">{s}\n")
            for pos in range(0, len(seq), 80):
                fo.write(seq[pos : pos + 80] + "\n")
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
                    block.append(_Row(name=parts[1], ref_start=int(parts[2]), text=parts[6]))
            elif line.strip() == "" and block:
                yield block
                block = []
        if block:
            yield block
