"""A near-identical sibling of the query must not be reported as the major parent.

Issue #112. The HMM sets the query's whole-genome siblings aside before it segments, and
the log says they are excluded from the competition. The base-pair vote that
``reconcile_major`` consults has no such step: a sibling that is a near-copy of the query
wins every window, so the HMM's major (the true backbone) "wins no window" and was
replaced by the vote winner -- the sibling just excluded. The one region reported then
read ``minor=refB major=sib`` with the minor parent less similar than the major.

The geometry is the one from the issue: a refA backbone carrying a refB-derived tract,
and a sibling that is the query plus 0.05 % noise (so it carries the tract too).
"""

from __future__ import annotations

import csv
import logging
import random
from pathlib import Path

from tessera.recomb.run import RecombParams, run_recomb

from ..conftest import write_fasta

LENGTH = 6000
TRACT = (2000, 3500)


def _mutate(seq: str, rate: float, rng: random.Random) -> str:
    out = []
    for c in seq:
        if rng.random() < rate:
            out.append(rng.choice([b for b in "ACGT" if b != c]))
        else:
            out.append(c)
    return "".join(out)


def _records(*, with_clone: bool, seed: int = 21) -> dict[str, str]:
    rng = random.Random(seed)
    root = "".join(rng.choice("ACGT") for _ in range(LENGTH))
    a, b = _mutate(root, 0.03, rng), _mutate(root, 0.03, rng)
    query = a[: TRACT[0]] + b[TRACT[0]: TRACT[1]] + a[TRACT[1]:]
    query = _mutate(query, 0.005, rng)
    records = [
        ("refA", _mutate(a, 0.02, rng)),
        ("refB", _mutate(b, 0.02, rng)),
        ("out1", _mutate(root, 0.06, rng)),
        ("out2", _mutate(root, 0.06, rng)),
        ("query", query),
    ]
    if with_clone:
        records.append(("sib", _mutate(query, 0.0005, rng)))
    return dict(records)


def _run(tmp_path: Path, records) -> tuple[list[dict[str, str]], str]:
    msa = write_fasta(tmp_path / "msa.fasta", records)
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=msa, output=out, query="query", window_size=500,
                     window_step=100, plot_format="png"),
        logging.getLogger("tessera.test"),
    )
    with open(out / "recombination_regions.tsv", encoding="utf-8", newline="") as fo:
        rows = list(csv.DictReader(fo, delimiter="\t"))
    return [r for r in rows if r["donor_absent"] == "no"], str(out)


def test_a_clone_sibling_is_not_named_as_the_major_parent(tmp_path: Path) -> None:
    rows, _ = _run(tmp_path, _records(with_clone=True))
    assert rows, "the refB tract should still be found"
    assert all(r["major_parent"] != "sib" for r in rows)
    assert all(r["minor_parent"] != "sib" for r in rows)
    assert {r["major_parent"] for r in rows} == {"refA"}
    assert {r["minor_parent"] for r in rows} == {"refB"}


def test_a_reported_region_never_has_a_less_similar_donor_than_backbone(tmp_path: Path) -> None:
    rows, _ = _run(tmp_path, _records(with_clone=True))
    for r in rows:
        assert float(r["mean_sim_minor"]) >= float(r["mean_sim_major"])


def test_the_panel_without_a_clone_is_unchanged(tmp_path: Path) -> None:
    """The control: the same tract, no sibling, is called against refA as before."""
    rows, _ = _run(tmp_path, _records(with_clone=False))
    assert rows
    assert {r["major_parent"] for r in rows} == {"refA"}
    assert {r["minor_parent"] for r in rows} == {"refB"}
