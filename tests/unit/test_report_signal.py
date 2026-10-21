"""The parent-free signal section states the alpha the run used, and says "not
testable" when the PHI test could not have rejected."""

from __future__ import annotations

import json
import logging

from tessera.recomb.diagnostics import RecombinationSignal
from tessera.recomb.report_html import _signal_html
from tessera.recomb.report_text import write_profile_tsv
from tessera.recomb.run import RecombParams, run_recomb


def _signal(phi_p: float | None, n_informative: int = 500) -> RecombinationSignal:
    return RecombinationSignal(
        n_informative=n_informative, phi_p=phi_p, phi_observed=0.1, phi_window=100,
        rmin=2, rmin_intervals=[(100, 200), (900, 1000)], profile=[(150, 150, 0.2)],
    )


def test_signal_section_says_not_testable() -> None:
    html = _signal_html(_signal(None, n_informative=54), alpha=0.05)
    assert "not testable" in html
    assert "p =" not in html
    assert "no significant recombination signal" not in html
    assert "54 informative sites" in html
    # 54 sites and a window of 53 ranks is also untestable (every pair is inside it), so
    # "do not exceed the window" would be false there; describe the actual condition.
    assert "do not exceed" not in html
    assert "every pair" in html
    assert "--phi-window" in html  # tells the reader what to change


def test_profile_tsv_header_marks_an_untestable_phi(tmp_path) -> None:
    write_profile_tsv(_signal(None, n_informative=54), tmp_path, logging.getLogger("tessera"))
    header = (tmp_path / "recombination_profile.tsv").read_text().splitlines()[0]
    assert header.split("\t")[1] == "NA"
    assert "not testable" in header


def test_run_reports_an_untestable_phi(example_data, tmp_path, logger) -> None:
    # cryptic_insert has 54 informative sites: fewer than the default window of 100.
    out = tmp_path / "cryptic"
    run_recomb(
        RecombParams(msa=example_data / "cryptic_insert.msa.fasta", output=out,
                     query="query", plot_format="png"),
        logger,
    )
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["recombination signal (PHI)"].startswith("not testable")
    assert "not testable" in (out / "report.html").read_text()


def test_report_states_the_alpha_the_run_used(example_data, tmp_path, logger) -> None:
    """The section always printed "alpha 0.05" and judged significance at 0.05, while
    the per-region PHI flag used --alpha."""
    out = tmp_path / "divergent"
    run_recomb(
        RecombParams(msa=example_data / "divergent.msa.fasta", output=out, query="query",
                     window_size=300, window_step=30, plot_format="png", alpha=0.2,
                     phi_window=20),
        logger,
    )
    report = (out / "report.html").read_text()
    assert "(alpha 0.2;" in report
    assert "(alpha 0.05;" not in report
