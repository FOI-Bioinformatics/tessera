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


def _run(tmp_path: Path, records, **overrides) -> tuple[list[dict[str, str]], str]:
    msa = write_fasta(tmp_path / "msa.fasta", records)
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=msa, output=out, query="query", window_size=500,
                     window_step=100, plot_format="png", **overrides),
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


# --- the site callers compete without the clone as well (issue #120) -----------------

def _spans(rows: list[dict[str, str]]) -> set[tuple[str, str, str, str]]:
    return {(r["minor_parent"], r["major_parent"], r["query_start"], r["query_end"])
            for r in rows}


def test_the_tract_is_called_by_more_than_one_caller(tmp_path: Path) -> None:
    """Before: only the HMM set the clone aside, so 3SEQ, MaxChi and Bootscan competed it
    and found nothing; the one region listed ``hmm`` alone."""
    rows, _ = _run(tmp_path, _records(with_clone=True))
    assert rows
    assert all(len(r["methods"].split(",")) >= 2 for r in rows), [r["methods"] for r in rows]


def test_the_tract_survives_a_gate_of_two(tmp_path: Path) -> None:
    rows, _ = _run(tmp_path, _records(with_clone=True), min_methods=2)
    assert {r["minor_parent"] for r in rows} == {"refB"}
    assert {r["major_parent"] for r in rows} == {"refA"}


def test_a_site_caller_alone_also_sets_the_clone_aside(tmp_path: Path) -> None:
    """Sibling detection does not depend on the HMM being among the selected callers."""
    rows, _ = _run(tmp_path, _records(with_clone=True), methods=("3seq",))
    assert {r["minor_parent"] for r in rows} == {"refB"}
    assert {r["major_parent"] for r in rows} == {"refA"}


def test_the_order_of_the_selected_callers_does_not_change_the_result(tmp_path: Path) -> None:
    (tmp_path / "f").mkdir()
    (tmp_path / "b").mkdir()
    forward, _ = _run(tmp_path / "f", _records(with_clone=True), methods=("hmm", "3seq"))
    backward, _ = _run(tmp_path / "b", _records(with_clone=True), methods=("3seq", "hmm"))
    assert _spans(forward) == _spans(backward)


def test_without_sibling_exclusion_the_clone_still_masks_the_event(tmp_path: Path) -> None:
    """Guard: ``--no-exclude-siblings`` keeps its meaning for every caller."""
    rows, _ = _run(tmp_path, _records(with_clone=True), exclude_siblings=False)
    assert rows == []


def test_a_clonal_query_with_a_copy_of_itself_reports_nothing(tmp_path: Path) -> None:
    """Guard: a sibling is set aside only when removing it exposes a mosaic, so a query
    with no recombinant tract is untouched however many copies of it the panel holds."""
    records = _records(with_clone=False)
    clonal = records["refA"]
    records["query"] = _mutate(clonal, 0.005, random.Random(5))
    records["sib"] = _mutate(records["query"], 0.0005, random.Random(6))
    rows, _ = _run(tmp_path, records)
    assert rows == []


def _twins(tmp_path: Path, records, *, regions_found: bool = True):
    from tessera.recomb.analyze import analyze
    from tessera.recomb.regions import RegionParams, call_regions, masking_twins
    from tessera.recomb.similarity import compute_similarity

    msa = write_fasta(tmp_path / "msa.fasta", records)
    result = compute_similarity(str(msa), "query", window_size=500, window_step=100)
    params = RegionParams.with_defaults(500, method="hmm")
    regions, _, dropped = call_regions(result, analyze(result), 500, params)
    return masking_twins(result, 500, params, regions if regions_found else [], dropped), dropped


def test_the_clone_is_a_masking_twin(tmp_path: Path) -> None:
    twins, dropped = _twins(tmp_path, _records(with_clone=True))
    assert [e.label for e in dropped] == ["sib"]
    assert [e.label for e in twins] == ["sib"]


def test_a_sibling_decision_the_hmm_did_not_back_with_a_region_is_not_inherited(
    tmp_path: Path,
) -> None:
    """On a clonal query the HMM's guard can fire by chance and drop a reference without
    reporting anything; the other callers must not change because of it."""
    twins, dropped = _twins(tmp_path, _records(with_clone=True), regions_found=False)
    assert dropped  # the HMM did set it aside
    assert twins == []
