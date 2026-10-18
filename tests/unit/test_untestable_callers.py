"""Callers that cannot test do not count toward the agreement gate; --merge-gap reaches the
coverage diagnostic.

Audit 4, items B7 (the gate half, handed over by plan F) and B6 (the coverage half).
"""

from __future__ import annotations

import json
import logging
import random
from pathlib import Path

from tessera.recomb import run as run_module
from tessera.recomb.run import RecombParams, run_recomb


def _rows(path: Path) -> list[dict[str, str]]:
    lines = [x for x in path.read_text().splitlines() if x.strip()]
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, x.split("\t"), strict=False)) for x in lines[1:]]


def _one_lineage_panel(tmp_path: Path) -> Path:
    """Five references 0.25 % apart and a query carrying a real 6 kb tract from one.

    Under base-pair windowing every pair stays above the 98.5 % pooling floor in every
    window, so lineage clustering pools them into a single HMM state. (The same panel
    as plan F's ``test_not_run_reporting._one_lineage_panel``.)
    """
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


def test_a_one_state_hmm_does_not_count_toward_the_agreement_gate(tmp_path, logger) -> None:
    # hmm + 3seq at --min-methods 2. The HMM is left one state and can call nothing, so
    # a gate of two could only be met by a caller that tested nothing: every 3SEQ region
    # was suppressed and the run read as a negative.
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=_one_lineage_panel(tmp_path), output=out, query="query",
                     plot_format="png", methods=("hmm", "3seq"), min_methods=2,
                     informative_sites=False),
        logger,
    )
    regions = _rows(out / "recombination_regions.tsv")
    called = [r for r in regions if r.get("donor_absent") != "yes"]
    assert called, "3SEQ's region is reported once the gate is lowered to the callers that ran"
    assert all(r["methods"] == "3seq" for r in called)
    gate = json.loads((out / "run_provenance.json").read_text())["run"][
        "min methods (agreement gate)"]
    assert gate == "1 (requested 2; 1 caller could test)"


def test_a_testable_hmm_still_counts(tmp_path, logger) -> None:
    # Guard: with clustering off the HMM has five states and counts as before.
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=_one_lineage_panel(tmp_path), output=out, query="query",
                     plot_format="png", methods=("hmm", "3seq"), min_methods=2,
                     informative_sites=False, cluster_lineages=False),
        logger,
    )
    gate = json.loads((out / "run_provenance.json").read_text())["run"][
        "min methods (agreement gate)"]
    assert gate == "2"


def test_merge_gap_reaches_the_coverage_diagnostic(tmp_path, logger, monkeypatch) -> None:
    seen: list[int] = []
    real = run_module.call_coverage_gaps

    def spy(result, window_size, params):
        seen.append(params.merge_gap)
        return real(result, window_size, params)

    monkeypatch.setattr(run_module, "call_coverage_gaps", spy)
    msa = _one_lineage_panel(tmp_path)
    for label, merge_gap in (("default", None), ("set", 250)):
        run_recomb(RecombParams(msa=msa, output=tmp_path / label, query="query",
                                plot_format="png", methods=("3seq",),
                                merge_gap=merge_gap), logger)
    assert seen == [1000, 250]  # the window size by default, the option when given


def _warnings(logger: logging.Logger) -> tuple[list[str], logging.Handler]:
    """Collect WARNING-level messages from ``logger`` (the fixture logger may not
    propagate to pytest's caplog, so attach a handler of our own)."""
    seen: list[str] = []

    class Collect(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            if record.levelno >= logging.WARNING:
                seen.append(record.getMessage())

    handler = Collect()
    logger.addHandler(handler)
    return seen, handler


def test_the_lowered_gate_warning_says_could_test_not_ran(tmp_path, logger) -> None:
    # The HMM executes and decodes one state, so it ran but could not test. The warning
    # and the provenance must not call it a caller that ran: only the provenance string
    # was pinned, so this pins the warning.
    seen, handler = _warnings(logger)
    try:
        run_recomb(
            RecombParams(msa=_one_lineage_panel(tmp_path), output=tmp_path / "out",
                         query="query", plot_format="png", methods=("hmm", "3seq"),
                         min_methods=2, informative_sites=False),
            logger,
        )
    finally:
        logger.removeHandler(handler)
    lowered = [m for m in seen if "--min-methods 2 cannot be met" in m]
    assert len(lowered) == 1, seen
    assert "1 caller(s) could test" in lowered[0]
    assert "ran" not in lowered[0]


def test_a_one_state_hmm_changes_nothing_with_the_default_callers(tmp_path, logger) -> None:
    # The confinement claim: the clamp moves only when the callers that could test number
    # fewer than --min-methods. With the default four callers and --min-methods 2 the HMM
    # drops to three, the min() still saturates, and nothing about the gate changes.
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=_one_lineage_panel(tmp_path), output=out, query="query",
                     plot_format="png", min_methods=2, informative_sites=False),
        logger,
    )
    run = json.loads((out / "run_provenance.json").read_text())["run"]
    assert run["min methods (agreement gate)"] == "2"  # not lowered, and no "requested"
    assert run["callers not run"].startswith("hmm (")  # the HMM still could not test
    called = [r for r in _rows(out / "recombination_regions.tsv")
              if r.get("donor_absent") != "yes"]
    assert called and all("hmm" not in r["methods"] for r in called)


def _two_gap_panel(tmp_path: Path) -> Path:
    """A query with two stretches no reference covers, 4000 columns apart."""
    rng = random.Random(11)
    anc = [rng.choice("ACGT") for _ in range(15000)]
    refs = {f"ref_{i}": [rng.choice([b for b in "ACGT" if b != c]) if rng.random() < 0.01
                         else c for c in anc] for i in range(3)}
    query = list(refs["ref_0"])
    for lo, hi in ((2500, 5500), (9500, 12500)):
        query[lo:hi] = [rng.choice("ACGT") for _ in range(hi - lo)]
    msa = tmp_path / "gaps.fasta"
    msa.write_text(">query\n" + "".join(query) + "\n"
                   + "".join(f">{k}\n{''.join(v)}\n" for k, v in refs.items()))
    return msa


def test_merge_gap_merges_the_coverage_gaps_it_is_given(tmp_path, logger) -> None:
    # The spy above pins the parameter; this pins the effect. Two uncovered stretches
    # 4000 columns apart stay two gaps at the default and become one when --merge-gap
    # spans the distance between them.
    msa = _two_gap_panel(tmp_path)
    gaps = {}
    for label, merge_gap in (("default", None), ("wide", 6000)):
        out = tmp_path / label
        run_recomb(RecombParams(msa=msa, output=out, query="query", plot_format="png",
                                methods=("3seq",), merge_gap=merge_gap), logger)
        # coverage_gaps.tsv opens with a "# best-similarity threshold" comment line.
        lines = [x for x in (out / "coverage_gaps.tsv").read_text().splitlines()
                 if x.strip() and not x.startswith("#")]
        header = lines[0].split("\t")
        gaps[label] = [dict(zip(header, x.split("\t"), strict=False)) for x in lines[1:]]
    assert [(g["msa_start"], g["msa_end"]) for g in gaps["default"]] == [
        ("1600", "6400"), ("8600", "13400")]
    # 2200 columns apart, so a 6000-column --merge-gap spans them and they report as one.
    assert [(g["msa_start"], g["msa_end"]) for g in gaps["wide"]] == [("1600", "13400")]
