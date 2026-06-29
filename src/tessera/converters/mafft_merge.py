"""Merge a MAFFT ``--addfragments --keeplength`` alignment into one genome row.

``mafft --addfragments <genome> --keeplength <reference>`` returns an alignment
whose length equals the reference's: the first record is the (unchanged)
reference, followed by one record per added contig, each already projected onto
reference coordinates (insertions relative to the reference are dropped by
``--keeplength``). A draft genome contributes several contig rows; this merges
them into a single reference-anchored row by taking, per column, the first
non-gap base across the contigs.
"""

from __future__ import annotations

from pathlib import Path

from ..core.io import read_fasta


def merge_added_fragments(aligned_path: str | Path) -> tuple[str, str]:
    """Return ``(reference_row, merged_genome_row)`` from a keeplength alignment.

    The reference is the first record; all remaining records are the added
    contigs of one genome and are merged column-wise (first non-gap wins).
    """
    records = read_fasta(aligned_path)
    if not records:
        raise ValueError(f"Empty MAFFT alignment: {aligned_path}")
    reference_row = records[0][1]
    width = len(reference_row)
    merged = bytearray(b"-" * width)
    for _, seq in records[1:]:
        for i, base in enumerate(seq):
            if i >= width:
                break
            if base != "-" and merged[i] == ord("-"):
                merged[i] = ord(base)
    return reference_row, merged.decode("ascii")
