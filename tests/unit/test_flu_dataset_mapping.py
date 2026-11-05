"""Influenza segment taxa map to the per-subtype, per-segment Nextclade datasets (#128).

``tessera reassort`` resolves a dataset for every segment from the taxon of a BLAST probe.
The alias table is keyed by organism, so a flu segment ("influenza a virus (a/.../2022(h3n2))
segment 6 neuraminidase (na) gene") mapped to nothing and the run ended UNDETERMINED.
"""

from __future__ import annotations

import logging

import pytest

from tessera.core.errors import UserInputError
from tessera.discover import nextclade
from tessera.discover.nextclade import _flu_dataset, _map_taxon

INDEX = {p: {"path": p, "versions": [{"tag": "t1"}], "files": {}} for p in (
    "nextstrain/flu/h3n2/ha/EPI1857216", "nextstrain/flu/h3n2/ha/CY163680",
    "nextstrain/flu/h3n2/na/EPI1857215", "nextstrain/flu/h3n2/pb2",
    "nextstrain/flu/h3n2/mp", "nextstrain/flu/h1n1pdm/ha/MW626062",
    "nextstrain/flu/h1n1pdm/na/MW626056", "nextstrain/flu/h1n1/ha",
    "nextstrain/flu/vic/ha/KX058884", "nextstrain/flu/yam/ha/JN993010",
)}


@pytest.mark.parametrize(("taxon", "expected"), [
    ("influenza a virus (a/hong kong sar/1/1968(h3n2)) segment 4 hemagglutinin (ha) gene",
     "nextstrain/flu/h3n2/ha/EPI1857216"),
    ("influenza a virus (a/ezhou/1812/2022(h3n2)) segment 6 neuraminidase (na) gene",
     "nextstrain/flu/h3n2/na/EPI1857215"),
    ("influenza a virus (a/x/1/2022(h3n2)) segment 1 polymerase pb2 (pb2) gene",
     "nextstrain/flu/h3n2/pb2"),
    ("influenza a virus (a/x/1/2022(h3n2)) segment 7 matrix protein 2 (m2) and matrix "
     "protein 1 (m1) genes", "nextstrain/flu/h3n2/mp"),
    ("influenza a virus (a/california/07/2009(h1n1)) segment 4 hemagglutinin (ha) gene",
     "nextstrain/flu/h1n1pdm/ha/MW626062"),
    ("influenza a virus (a/brisbane/59/2007(h1n1)) segment 4 hemagglutinin (ha) gene",
     "nextstrain/flu/h1n1/ha"),
    ("influenza b virus (b/victoria/2/87) segment 4 hemagglutinin (ha) gene",
     "nextstrain/flu/vic/ha/KX058884"),
    ("influenza b virus (b/yamagata/16/88) segment 4 hemagglutinin (ha) gene",
     "nextstrain/flu/yam/ha/JN993010"),
])
def test_a_flu_segment_taxon_maps_to_its_dataset(taxon: str, expected: str) -> None:
    assert _flu_dataset(taxon, INDEX) == expected


@pytest.mark.parametrize("taxon", [
    "influenza a virus (a/x/1/2022(h3n2)) segment 4 hemagglutinin (ha) gene "
    .replace("(h3n2)", ""),                                     # no subtype
    "influenza a virus (a/x/1/2022(h3n2)) hypothetical protein",    # no segment
    "influenza b virus (b/x/1/2022) segment 4 hemagglutinin (ha) gene",  # lineage unknown
    "influenza a virus (a/x/1/2022(h7n9)) segment 4 hemagglutinin (ha) gene",  # no dataset
    "dengue virus 2",
])
def test_a_name_missing_the_subtype_or_segment_is_refused_not_guessed(taxon: str) -> None:
    assert _flu_dataset(taxon, INDEX) is None


def test_map_taxon_uses_the_flu_mapping_and_logs_it(caplog: pytest.LogCaptureFixture) -> None:
    logger = logging.getLogger("tessera.test.flu")
    with caplog.at_level(logging.INFO, logger=logger.name):
        ds = _map_taxon("influenza a virus (a/x/1/2022(h3n2)) segment 6 neuraminidase (na) gene",
                        INDEX, logger)
    assert ds.path == "nextstrain/flu/h3n2/na/EPI1857215"
    assert "Mapped" in caplog.text


def test_the_refusal_names_the_option_the_command_takes() -> None:
    logger = logging.getLogger("tessera.test.flu")
    with pytest.raises(UserInputError, match=r"--dataset SEGMENT=<path>"):
        _map_taxon("sars-cov-2 is not flu", {}, logger, option="--dataset SEGMENT=<path>")
    with pytest.raises(UserInputError, match="--nextclade-dataset"):
        _map_taxon("nothing", {}, logger)


def test_resolve_dataset_passes_the_option_into_the_refusal(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(nextclade, "fetch_index", lambda: {})
    monkeypatch.setattr(nextclade.shutil, "which", lambda name: None)
    monkeypatch.setattr(nextclade, "detect_taxon", lambda *a, **k: "Unknown virus")
    with pytest.raises(UserInputError, match=r"--dataset SEGMENT=<path>"):
        nextclade.resolve_dataset(
            tmp_path / "q.fa", None, email=None, logger=logging.getLogger("t"),
            option="--dataset SEGMENT=<path>",
        )


def test_a_long_parent_list_is_shortened_for_the_log() -> None:
    from tessera.cli.cmd_reassort import _few

    assert _few(["a", "b", "c"]) == "a,b,c"
    assert _few([f"s{i}" for i in range(100)]) == "s0,s1,s2 and 97 more"
    assert _few([]) == ""
