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
    return list(csv.DictReader(path.open(encoding="utf-8"), delimiter="\t"))


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


def _one_lineage_panel(tmp_path: Path) -> Path:
    """Five references 0.25 % apart and a query carrying a real 6 kb tract from one.

    Under base-pair windowing every pair of references stays above the 98.5 % pooling
    floor in every window, so lineage clustering pools them into a single HMM state.
    """
    import random

    rng = random.Random(3)
    anc = [rng.choice("ACGT") for _ in range(20000)]

    def mut(seq: list[str], rate: float) -> list[str]:
        return [rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c
                for c in seq]

    refs = {f"ref_{i}": mut(anc, 0.0025) for i in range(5)}
    query = list(refs["ref_0"])
    query[8000:14000] = refs["ref_3"][8000:14000]
    msa = tmp_path / "pool.fasta"
    msa.write_text(">query\n" + "".join(query) + "\n"
                   + "".join(f">{k}\n{''.join(v)}\n" for k, v in refs.items()))
    return msa


def test_single_state_hmm_is_reported_as_not_run(tmp_path: Path, logger) -> None:
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=_one_lineage_panel(tmp_path), output=out, query="query",
                     plot_format="png", methods=("hmm",), informative_sites=False),
        logger,
    )
    assert _rows(out / "recombination_regions.tsv") == []  # unchanged: nothing called
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["callers not run"].startswith("hmm (lineage clustering pooled every reference")
    report = (out / "report.html").read_text()
    assert "No recombination among the present references" not in report
    assert "not a negative result" in report
    assert "--no-cluster-lineages" in report


def test_single_state_hmm_is_not_run_in_the_method_table(tmp_path: Path, logger) -> None:
    """In an ensemble the other callers still call the tract; the HMM column must not
    read as a caller that looked and found nothing."""
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=_one_lineage_panel(tmp_path), output=out, query="query",
                     plot_format="png", informative_sites=False),
        logger,
    )
    rows = _rows(out / "recombination_methods.tsv")
    assert rows, "the site callers find the 6 kb tract"
    assert {r["hmm"] for r in rows} == {"not run"}  # was "no"
    regions = _rows(out / "recombination_regions.tsv")
    assert all("hmm" not in r["methods"] for r in regions)


def test_barcode_and_one_state_hmm_are_both_reported_not_run(tmp_path: Path, logger) -> None:
    """Review focus 2: two callers not run in one run, each with its own reason.
    Neither counts toward the agreement gate."""
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=_one_lineage_panel(tmp_path), output=out, query="query",
                     plot_format="png", informative_sites=False,
                     methods=("hmm", "3seq", "barcode"), min_methods=2),
        logger,
    )
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["callers not run"].startswith("hmm (lineage clustering")
    assert "; barcode (it needs typed references" in run["callers not run"]
    assert run["min methods (agreement gate)"] == "1 (requested 2; 1 caller could test)"
