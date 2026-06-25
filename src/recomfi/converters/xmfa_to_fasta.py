"""XMFA -> FASTA converter (Python 3 port of x2fa.py).

Port of ``x2fa.py`` v9 by Adrian Larkeryd, originally invoked as a Python 2
subprocess by the legacy ``phylo.py``. Behaviour is preserved: the reference
sequence defines the coordinate system, reference gaps are removed, optional
deletion flanks are screened, and every sequence is projected onto the reference
coordinates. Strings are handled as ``bytearray`` / ``bytes`` throughout (the
original relied on Python 2 ``str``/``bytearray`` equivalence); regex scanning
uses a latin-1 view so byte positions are 1:1 with characters.
"""

from __future__ import annotations

import re
from pathlib import Path

# Complement table over the IUPAC subset the original handled.
_COMPLEMENTS = bytes.maketrans(b"acgtrymkbdhvACGTRYMKBDHV", b"tgcayrkmvhdbTGCAYRKMVHDB")

_PATTERN_START = re.compile(r"^>\s*(\d+):(\d+)-(\d+) ([+-])")
_PATTERN_SEQ_NAME = re.compile(r"#Sequence(\d+)File")
_PATTERN_COMMENT = re.compile(r"#")
_PATTERN_GAP = re.compile("-+")


def reverse_complement(dna: bytes) -> bytes:
    """Complement and reverse a DNA byte string."""
    return dna.translate(_COMPLEMENTS)[::-1]


def xmfa_to_fasta(
    xmfa_path: str | Path,
    reference_name: str,
    flank: int,
    out_path: str | Path,
    reference_length: int | None = None,
) -> Path:
    """Convert a progressiveMauve XMFA file to a reference-anchored FASTA.

    ``reference_name`` must match a ``#SequenceNFile`` entry in the XMFA.
    ``flank`` (>=0) screens deletion flanks by that many bases.

    ``reference_length`` pins the output width to the full reference genome
    length. Without it the width is the furthest aligned reference position in
    this XMFA, which varies per query; pinning it keeps every per-query
    projection on the same coordinate system so they concatenate into a
    rectangular MSA.
    """
    xmfa_path = Path(xmfa_path)
    out_path = Path(out_path)

    a_gen: dict[int, dict[int, dict]] = {}
    rm_h: dict[int, dict[int, int]] = {}
    name2num: dict[str, int] = {}
    num2name: dict[int, str] = {}

    alignment_number = 0
    curr_seq: int | None = None
    curr_pos = 0
    a_gen[alignment_number] = {}
    rm_h[alignment_number] = {}

    with open(xmfa_path) as xmfa:
        for line in xmfa:
            if _PATTERN_START.search(line):
                curr_seq = int(line.split(":")[0].split(" ")[1])
                curr_pos = 0
                a_gen[alignment_number][curr_seq] = {}

                startend = line.split(" ")[1].split(":")[1].split("-")
                p1, p2 = int(startend[0]), int(startend[1])
                lo, hi = (p1, p2) if p1 < p2 else (p2, p1)
                a_gen[alignment_number][curr_seq]["p1"] = lo
                a_gen[alignment_number][curr_seq]["p2"] = hi
                a_gen[alignment_number][curr_seq]["sign"] = line.split(" ")[2]
                a_gen[alignment_number][curr_seq]["seq"] = bytearray(hi - lo)
            elif line.strip() == "=":
                alignment_number += 1
                a_gen[alignment_number] = {}
                rm_h[alignment_number] = {}
            elif _PATTERN_SEQ_NAME.search(line):
                num = int(line.split("Sequence")[1].split("File")[0])
                name = line.split("\t")[1].strip()
                name2num[name] = num
                num2name[num] = name
            elif _PATTERN_COMMENT.search(line):
                pass
            else:
                if curr_seq is None:
                    continue
                stripped = line.strip().encode("utf-8")
                length_of_line = len(line) - 1
                seq = a_gen[alignment_number][curr_seq]["seq"]
                seq[curr_pos : curr_pos + length_of_line] = stripped
                curr_pos += length_of_line

    reference_num = name2num[reference_name]

    # Output width: the full reference length when known, otherwise the furthest
    # aligned reference position in this XMFA. The max() guards against a p2 that
    # somehow exceeds the supplied length, which would otherwise extend one row.
    length_of_reference = 0
    for alignment in a_gen:
        ref = a_gen[alignment].get(reference_num)
        if ref is not None and ref["p2"] > length_of_reference:
            length_of_reference = ref["p2"]
    if reference_length is not None:
        length_of_reference = max(length_of_reference, reference_length)

    outseqs: dict[int, bytearray] = {
        num: bytearray(b"-" * length_of_reference) for num in num2name
    }

    for alignment in list(a_gen.keys()):
        if reference_num not in a_gen[alignment]:
            del a_gen[alignment]
            continue

        # Find reference gaps and remove those columns from every sequence.
        search_pos = 0
        list_of_gaps: list[list[int]] = []
        ref_view = a_gen[alignment][reference_num]["seq"].decode("latin-1")
        while search_pos < len(a_gen[alignment][reference_num]["seq"]):
            gap_hit = _PATTERN_GAP.search(ref_view, search_pos)
            if gap_hit:
                rm_h[alignment][gap_hit.start()] = gap_hit.end() - gap_hit.start()
                search_pos = gap_hit.end()
            else:
                break

        for pos in reversed(sorted(rm_h[alignment].keys())):
            if rm_h[alignment][pos]:
                start = pos
                end = pos + rm_h[alignment][pos]
                for sequence in a_gen[alignment]:
                    a_gen[alignment][sequence]["seq"][start:end] = bytearray()
                if flank > 0:
                    for g in list_of_gaps:
                        g[0] = g[0] - end + start
                        g[1] = g[1] - end + start
                    list_of_gaps.append([start, start])

        if flank > 0:
            search_pos = 0
            for sequence in a_gen[alignment]:
                if sequence == reference_num:
                    continue
                seq_view = a_gen[alignment][sequence]["seq"].decode("latin-1")
                while search_pos < length_of_reference:
                    gap_hit = _PATTERN_GAP.search(seq_view, search_pos)
                    if gap_hit:
                        list_of_gaps.append([gap_hit.start(), gap_hit.end()])
                        search_pos = gap_hit.end()
                    else:
                        break
            for non_ref_gap in list_of_gaps:
                for sequence in a_gen[alignment]:
                    seq = a_gen[alignment][sequence]["seq"]
                    new_start = max(0, non_ref_gap[0] - flank)
                    new_end = min(non_ref_gap[1] + flank, len(seq))
                    seq[new_start:new_end] = bytearray(b"-" * (new_end - new_start))

    # Project each alignment block onto reference coordinates.
    for alignment in a_gen:
        ref = a_gen[alignment][reference_num]
        start = int(ref["p1"]) - 1
        end = int(ref["p2"])
        if flank > 0:
            start += flank
            end -= flank
        forward = ref["sign"] == "+"
        for sequence in a_gen[alignment]:
            seq = a_gen[alignment][sequence]["seq"]
            if not (start >= 0 and end > 0 and seq):
                continue
            fragment = bytes(seq[flank:-flank]) if flank > 0 else bytes(seq)
            if not forward:
                fragment = reverse_complement(fragment)
            outseqs[sequence][start:end] = fragment

    with open(out_path, "w") as fo:
        _write_sequence(fo, num2name[reference_num], outseqs[reference_num])
        for sequence in outseqs:
            if sequence != reference_num:
                _write_sequence(fo, num2name[sequence], outseqs[sequence])
    return out_path


def _write_sequence(fo, name: str, seq: bytearray, width: int = 80) -> None:
    fo.write(f">{name}\n")
    text = bytes(seq).decode("latin-1")
    for pos in range(0, len(text), width):
        fo.write(text[pos : pos + width] + "\n")
