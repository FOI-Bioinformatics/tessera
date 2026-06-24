"""Realistic reference-anchored MSA factory for tests.

Builds a deterministic (seeded) MSA shaped like the real orthopoxvirus example
dataset rather than uniform-random noise:

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
