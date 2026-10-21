"""Realistic reference-anchored MSA factories for tests.

Two deterministic (seeded) generators, both shaped after real runs rather than
uniform-random noise.

:func:`make_recombination_msa` is the positive case, shaped like the orthopoxvirus
example dataset; :func:`make_gap_degraded_msa` is the negative case behind issues
#81 and #82 -- a clonal query on a near-identical, clade-structured panel whose
nearest reference is misaligned over one locus, where any region reported is a
false positive.

The orthopoxvirus-shaped MSA carries:

- a backbone reference the query matches closely (``cowpox_KC813504`` in the real
  data, median similarity ~1.0);
- seven references at *graded* divergence matching the README medians
  (``cowpox_KC813504`` ~1.0, ``variola``/``camelpox``/``taterapox`` ~0.98,
  ``vaccinia``/``cowpox``/``monkeypox`` ~0.97);
- a recombinant donor block where the query equals ``variola`` (so the major
  parent is the backbone and one large region is called for variola);
- optional indels (gap columns), an N-run, and a fragmented multi-contig query
  whose inter-contig regions are gaps in the query row (mirroring the eight
  scaffolds), which makes query coordinates diverge from MSA columns.

Scale is reduced (~20 kbp by default) so tests run fast while preserving the
qualitative structure of the real data.
"""

from __future__ import annotations

import random
from dataclasses import dataclass

_BASES = "ACGT"

# Per-reference substitutions/site relative to the backbone, graded to mirror the
# README median similarities (1 - divergence).
GRADED_DIVERGENCE: dict[str, float] = {
    "cowpox_KC813504": 0.001,  # backbone-near: the query is essentially this
    "variola": 0.02,
    "camelpox": 0.02,
    "taterapox": 0.02,
    "vaccinia": 0.03,
    "cowpox": 0.03,
    "monkeypox": 0.03,
}


@dataclass
class RealisticMSA:
    records: dict[str, str]
    backbone: str
    donor: str
    donor_block: tuple[int, int] | None
    length: int


def _mutate(seq: str, rate: float, rng: random.Random) -> str:
    if rate <= 0:
        return seq
    chars = list(seq)
    for i, base in enumerate(chars):
        if rng.random() < rate:
            chars[i] = rng.choice([b for b in _BASES if b != base])
    return "".join(chars)


def make_recombination_msa(
    *,
    seed: int = 0,
    length: int = 20000,
    backbone: str = "cowpox_KC813504",
    donor: str = "variola",
    donor_block: tuple[int, int] | None = None,
    recombinant: bool = True,
    n_contigs: int = 1,
    contig_gap: int = 500,
    n_run: tuple[int, int] | None = None,
    uppercase: bool = True,
) -> RealisticMSA:
    """Build a realistic reference-anchored MSA (label -> aligned sequence)."""
    rng = random.Random(seed)
    base = "".join(rng.choice(_BASES) for _ in range(length))

    # References: each is the backbone mutated at its graded rate.
    references = {label: _mutate(base, rate, rng) for label, rate in GRADED_DIVERGENCE.items()}

    # Query: equals the backbone outside the donor block; equals the donor inside.
    query = list(base)
    if recombinant:
        if donor_block is None:
            donor_block = (int(length * 0.35), int(length * 0.60))
        d0, d1 = donor_block
        query[d0:d1] = list(references[donor][d0:d1])
    else:
        donor_block = None

    # A couple of indels: a short deletion in one reference and one in the query
    # (away from the donor block), represented as gap columns.
    references["camelpox"] = _with_gap(references["camelpox"], int(length * 0.10), 30)
    query = list(_with_gap("".join(query), int(length * 0.80), 30))

    # Optional N-run in the query (robustness: must not count as identity).
    if n_run is not None:
        s, ln = n_run
        query[s:s + ln] = list("N" * ln)

    # Fragment the query into contigs: blank out inter-contig regions as gaps in
    # the query row only (references stay intact).
    if n_contigs > 1:
        seg = length // n_contigs
        for k in range(1, n_contigs):
            g0 = k * seg - contig_gap // 2
            query[g0:g0 + contig_gap] = list("-" * contig_gap)

    records = {backbone: references[backbone]}
    for label, seq in references.items():
        if label != backbone:
            records[label] = seq
    records["query"] = "".join(query)

    if not uppercase:
        records = {k: v.lower() for k, v in records.items()}

    return RealisticMSA(
        records=records, backbone=backbone, donor=donor,
        donor_block=donor_block, length=length,
    )


def _with_gap(seq: str, start: int, gap_len: int) -> str:
    return seq[:start] + "-" * gap_len + seq[start + gap_len:]


# Default gap rates for the gap-degraded fixture below, taken from the mpox locus
# that prompted issues #81 and #82: 8.65% gaps in the major parent at the locus
# against 1.64% genome-wide.
LOCUS_GAP_RATE = 0.087
BACKGROUND_GAP_RATE = 0.016


@dataclass
class GapDegradedMSA:
    records: dict[str, str]
    nearest: str  # the query's sister tip: nearest reference, degraded at the locus
    locus: tuple[int, int]  # MSA columns over which that reference is misaligned
    locus_gap_rate: float  # its realized gap rate inside the locus
    background_gap_rate: float  # its realized gap rate genome-wide
    length: int
    donor_block: None = None  # the query is clonal: nothing was copied into it


def _sprinkle_gaps(
    seq: str, rate: float, rng: random.Random, span: tuple[int, int] | None = None
) -> str:
    """Gap out individual columns at ``rate``, modelling an ordinary indel load."""
    chars = list(seq)
    lo, hi = span if span is not None else (0, len(chars))
    for i in range(lo, hi):
        if rng.random() < rate:
            chars[i] = "-"
    return "".join(chars)


def _misalign(
    seq: str, locus: tuple[int, int], gap_rate: float, shift: int, rng: random.Random
) -> str:
    """Rewrite ``locus`` as a frame-shifted, gap-riddled copy of itself.

    A locally poor alignment is not a clean deletion: the residues placed in those
    columns come from elsewhere in the same sequence, so they mismatch the query at
    roughly the background base composition while the gap rate rises. Shifting
    rather than randomizing keeps the row's own composition intact, so the record
    still looks like the organism it came from.
    """
    lo, hi = locus
    donor = seq[lo + shift:hi + shift] or seq[lo:hi]
    chars = list(seq)
    for offset, col in enumerate(range(lo, hi)):
        chars[col] = "-" if rng.random() < gap_rate else donor[offset % len(donor)]
    return "".join(chars)


def _gap_rate(seq: str, span: tuple[int, int] | None = None) -> float:
    lo, hi = span if span is not None else (0, len(seq))
    window = seq[lo:hi]
    return window.count("-") / len(window) if window else 0.0


def make_gap_degraded_msa(
    *,
    seed: int = 1,
    length: int = 20000,
    locus: tuple[int, int] = (8000, 9500),
    locus_gap_rate: float = LOCUS_GAP_RATE,
    background_gap_rate: float = BACKGROUND_GAP_RATE,
    shift: int = 37,
    tip_divergence: float = 0.0015,
    clade_divergence: float = 0.004,
    deep_divergence: float = 0.012,
) -> GapDegradedMSA:
    """A clonal query on a low-divergence panel whose nearest reference is locally
    misaligned -- the geometry behind issues #81 and #82.

    Unlike :func:`make_recombination_msa` this is clade-structured rather than a
    star: ``query`` and ``sister`` are tips inside one clade, so ``sister`` wins the
    scan genome-wide on the derived sites they share. ``cousin`` sits in the same
    clade one step further out, three ``outclade`` tips sit in a second clade and
    ``distant`` is an outgroup, putting the whole panel within ~1% of the query --
    the near-identical regime where single columns decide a window.

    Nothing is copied into the query from anywhere, so **every region a caller
    reports on this MSA is a false positive**. The one defect planted is in
    ``sister``'s row over ``locus``, which is frame-shifted and extra-gapped as a
    locally poor alignment would be. The intended reading is that the locus is badly
    aligned; the reading a caller currently reaches is that the query is a mosaic of
    ``sister`` and whichever reference stayed intact there.
    """
    rng = random.Random(seed)
    ancestor = "".join(rng.choice(_BASES) for _ in range(length))

    # Clade A: the query and the two references that share its derived sites.
    stem_a = _mutate(ancestor, clade_divergence, rng)
    query = _mutate(stem_a, tip_divergence, rng)
    records = {
        "sister": _mutate(stem_a, tip_divergence, rng),
        "cousin": _mutate(stem_a, tip_divergence * 2, rng),
    }
    # Clade B plus an outgroup, so the panel has somewhere else to place a window.
    stem_b = _mutate(ancestor, clade_divergence, rng)
    for i in range(3):
        records[f"outclade{i}"] = _mutate(stem_b, tip_divergence * (i + 1), rng)
    records["distant"] = _mutate(ancestor, deep_divergence, rng)

    # The ordinary indel load every row in a real MSA carries, so the locus has to
    # stand out against a background rather than against a gap-free alignment.
    query = _sprinkle_gaps(query, background_gap_rate, rng)
    records = {
        label: _sprinkle_gaps(seq, background_gap_rate, rng)
        for label, seq in records.items()
    }
    background = _gap_rate(records["sister"])
    records["sister"] = _misalign(records["sister"], locus, locus_gap_rate, shift, rng)

    return GapDegradedMSA(
        records={"query": query, **records},
        nearest="sister",
        locus=locus,
        locus_gap_rate=_gap_rate(records["sister"], locus),
        background_gap_rate=background,
        length=length,
    )
