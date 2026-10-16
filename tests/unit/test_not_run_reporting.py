"""A test or caller that gave no result is reported as "not run", with its reason.

Audit 4, items B1 and B7: ``--no-phi`` made the report say the alignment had too few
informative sites, every region read ``parent_free_support = no``, and the provenance said
nothing; a single-state HMM read as a clean negative.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from tessera.recomb.run import RecombParams, run_recomb


def _rows(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(), delimiter="\t"))


def _divergent(example_data: Path, tmp_path: Path, logger, **kwargs) -> Path:
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=example_data / "divergent.msa.fasta", output=out, query="query",
                     window_size=300, window_step=30, plot_format="png", **kwargs),
        logger,
    )
    return out


def test_no_phi_is_reported_as_not_run(example_data, tmp_path, logger) -> None:
    out = _divergent(example_data, tmp_path, logger, phi=False)
    report = (out / "report.html").read_text()
    # divergent.msa.fasta has 159 informative columns; the old text blamed the alignment.
    assert "Too few informative sites" not in report
    assert "disabled with --no-phi" in report
    assert "not evidence against recombination" in report

    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["recombination signal (PHI)"] == "not run (disabled with --no-phi)"

    regions = _rows(out / "recombination_regions.tsv")
    assert regions, "the divergent example has a called region"
    assert {r["parent_free_support"] for r in regions} == {"not run"}  # was "no"
    methods = _rows(out / "recombination_methods.tsv")
    assert {r["parent_free_support"] for r in methods} == {"not run"}


def test_no_phi_on_a_single_caller_run_still_says_not_run(
    example_data, tmp_path, logger
) -> None:
    """Review focus 1: with one caller there is no method table, so the regions TSV is
    the only machine-readable place that says the test did not run."""
    out = _divergent(example_data, tmp_path, logger, phi=False, methods=("hmm",))
    assert not (out / "recombination_methods.tsv").exists()
    regions = _rows(out / "recombination_regions.tsv")
    assert regions and {r["parent_free_support"] for r in regions} == {"not run"}


def test_phi_that_ran_still_reports_yes_or_no(example_data, tmp_path, logger) -> None:
    """Guard (passes before and after): with the test on, the columns keep yes/no."""
    out = _divergent(example_data, tmp_path, logger, phi_window=20)
    regions = _rows(out / "recombination_regions.tsv")
    assert {r["parent_free_support"] for r in regions} <= {"yes", "no"}
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["recombination signal (PHI)"].startswith("p=")


def test_too_few_informative_sites_is_recorded_in_the_provenance(
    tmp_path: Path, logger
) -> None:
    """Below ten informative columns the test cannot run; the provenance said nothing."""
    import random

    rng = random.Random(1)
    base = [rng.choice("ACGT") for _ in range(3000)]
    seqs = {name: list(base) for name in ("query", "A", "B", "C")}
    for col in (500, 1500, 2500):  # three informative columns, split 2/2
        alt = "A" if base[col] != "A" else "C"
        seqs["query"][col] = seqs["A"][col] = alt
    msa = tmp_path / "few.fasta"
    msa.write_text("".join(f">{k}\n{''.join(v)}\n" for k, v in seqs.items()))
    out = tmp_path / "out"
    run_recomb(RecombParams(msa=msa, output=out, query="query", plot_format="png"), logger)
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["recombination signal (PHI)"].startswith("not run (fewer than 10 informative")
    assert "not run" in (out / "report.html").read_text()
