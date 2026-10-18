"""Callers that cannot test do not count toward the agreement gate; --merge-gap reaches the
coverage diagnostic.

Audit 4, items B7 (the gate half, handed over by plan F) and B6 (the coverage half).
"""

from __future__ import annotations

import json
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
    assert gate == "1 (requested 2; 1 caller ran)"


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
