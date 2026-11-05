"""Every reported region is at least 50 bp unless --min-region says otherwise (issue #121).

On aligned genomes with indels the site callers report runs of a few bp to a few tens of bp
that sit on local alignment structure, some with q-values far below 0.001. In the indel
simulations a 50 bp floor localised 37 of 72 tracts against 38 without one, and cut the
regions away from any true tract from 73 to 12 (validation/indel-specificity-results.md).
"""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

from tessera.recomb import run as run_module
from tessera.recomb.ensemble import filter_by_length
from tessera.recomb.run import DEFAULT_MIN_REGION_BP, RecombParams, run_recomb

from ..conftest import write_fasta


def _region(length: int) -> SimpleNamespace:
    """``filter_by_length`` reads only ``length_bp``."""
    return SimpleNamespace(length_bp=length)


def test_the_default_floor_is_fifty_bp() -> None:
    assert DEFAULT_MIN_REGION_BP == 50


def test_a_region_just_under_the_floor_is_removed_and_one_at_it_is_kept() -> None:
    regions = [_region(49), _region(50)]
    kept, _, removed = filter_by_length(regions, [{}, {}], DEFAULT_MIN_REGION_BP)
    assert [r.length_bp for r in kept] == [50]
    assert removed == 1


@pytest.fixture
def floors_seen(monkeypatch: pytest.MonkeyPatch) -> list[int]:
    seen: list[int] = []
    real = run_module.filter_by_length

    def spy(regions, breakdown, min_length):
        seen.append(min_length)
        return real(regions, breakdown, min_length)

    monkeypatch.setattr(run_module, "filter_by_length", spy)
    return seen


def _msa(tmp_path: Path) -> Path:
    import random
    rng = random.Random(3)
    root = "".join(rng.choice("ACGT") for _ in range(3000))

    def mutate(seq: str, rate: float) -> str:
        return "".join(rng.choice("ACGT") if rng.random() < rate else c for c in seq)

    records = {"refA": mutate(root, 0.03), "refB": mutate(root, 0.03),
               "out1": mutate(root, 0.06), "out2": mutate(root, 0.06)}
    records["query"] = mutate(records["refA"], 0.005)
    return write_fasta(tmp_path / "msa.fasta", records)


def _run(tmp_path: Path, **overrides) -> None:
    run_recomb(RecombParams(msa=_msa(tmp_path), output=tmp_path / "out", query="query",
                            window_size=500, window_step=100, plot_format="png", **overrides),
               logging.getLogger("tessera.test"))


def test_without_min_region_the_default_floor_is_applied(
    tmp_path: Path, floors_seen: list[int]
) -> None:
    _run(tmp_path)
    assert floors_seen == [DEFAULT_MIN_REGION_BP]


def test_an_explicit_min_region_replaces_the_default_floor(
    tmp_path: Path, floors_seen: list[int]
) -> None:
    _run(tmp_path, min_region=1)
    assert floors_seen == [1]
