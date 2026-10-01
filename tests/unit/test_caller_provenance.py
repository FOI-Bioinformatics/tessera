"""The run provenance names the callers that ran and the settings that shape the result.

`run_provenance.json` and the report's "Run parameters" table are what a reader uses to
reproduce a run. A caller described under another caller's name, or a setting that
changes which regions are reported but is not recorded, makes the record misleading.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tessera.recomb.regions import CALLERS
from tessera.recomb.run import RecombParams, caller_description, run_recomb


def _params(example_data: Path, output: Path, **kwargs) -> RecombParams:
    return RecombParams(
        msa=example_data / "divergent.msa.fasta", output=output, query="query",
        window_size=300, window_step=30, plot_format="png", **kwargs,
    )


@pytest.mark.parametrize("method", CALLERS)
def test_each_caller_is_described_under_its_own_name(example_data, tmp_path, method) -> None:
    description = caller_description(method, _params(example_data, tmp_path))
    assert description.startswith(method)


@pytest.mark.parametrize("method", [m for m in CALLERS if m != "heuristic"])
def test_only_the_heuristic_caller_is_called_heuristic(example_data, tmp_path, method) -> None:
    # maxchi, bootscan, geneconv and barcode used to fall through to the heuristic text.
    assert "heuristic" not in caller_description(method, _params(example_data, tmp_path))


def test_default_run_names_its_four_callers(example_data, tmp_path, logger) -> None:
    out = tmp_path / "out"
    run_recomb(_params(example_data, out), logger)
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert "heuristic" not in run["caller"]
    for name in ("hmm", "3seq", "maxchi", "bootscan"):
        assert name in run["caller"]


def test_provenance_records_the_settings_that_change_what_is_reported(
    example_data, tmp_path, logger
) -> None:
    out = tmp_path / "default"
    run_recomb(_params(example_data, out), logger)
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["min methods (agreement gate)"] == "1"
    assert run["sibling exclusion"] == "on"
    assert run["lineage clustering"] == "on"
    assert run["donor re-attribution"] == "off"

    out = tmp_path / "changed"
    run_recomb(
        _params(example_data, out, min_methods=2, exclude_siblings=False,
                cluster_lineages=False, reattribute_donors=True, reattribute_margin=0.05),
        logger,
    )
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["min methods (agreement gate)"] == "2"
    assert run["sibling exclusion"] == "off"
    assert run["lineage clustering"] == "off"
    assert run["donor re-attribution"] == "on (margin 0.05)"
