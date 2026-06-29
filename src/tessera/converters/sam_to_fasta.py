"""Project a SAM alignment onto reference coordinates (minimap2 output).

minimap2 ``-a`` aligns each query contig to the reference and emits SAM with the
query SEQ already on the forward reference strand. Walking the CIGAR places each
query base at its reference position, dropping insertions relative to the
reference (consistent with the other reference-anchored converters) and leaving
gaps where the reference is deleted/unaligned. All of a genome's contigs project
onto one row in reference coordinates.
"""

from __future__ import annotations

import re
from pathlib import Path

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

_FLAG_UNMAPPED = 0x4
_FLAG_SECONDARY = 0x100


def project_sam_to_row(
    sam_path: str | Path, rname_to_offset: dict[str, int], total_len: int
) -> str:
    """Return the genome's row (length ``total_len``) in reference coordinates.

    ``rname_to_offset`` maps each reference contig name to its start offset in the
    concatenated reference row (so a multi-contig reference is supported).
    Primary and supplementary alignments are placed; secondary and unmapped
    records are skipped.
    """
    row = bytearray(b"-" * total_len)
    with open(sam_path) as fo:
        for line in fo:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            if flag & (_FLAG_UNMAPPED | _FLAG_SECONDARY):
                continue
            rname, pos, cigar, seq = fields[2], int(fields[3]), fields[5], fields[9]
            if cigar == "*" or seq == "*" or rname not in rname_to_offset:
                continue

            ref_i = rname_to_offset[rname] + pos - 1
            q_i = 0
            for length_str, op in _CIGAR.findall(cigar):
                length = int(length_str)
                if op in "M=X":
                    for k in range(length):
                        rk = ref_i + k
                        if 0 <= rk < total_len and q_i + k < len(seq) and row[rk] == ord("-"):
                            row[rk] = ord(seq[q_i + k])
                    ref_i += length
                    q_i += length
                elif op in "IS":  # insertion / soft-clip: consume query only
                    q_i += length
                elif op in "DN":  # deletion / skip: consume reference only
                    ref_i += length
                # H (hard-clip) and P (padding) consume neither
    return row.decode("ascii")
