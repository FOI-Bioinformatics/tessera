"""End-to-end recombination scan on a synthetic MSA (no external binaries)."""

from __future__ import annotations

import logging
import random
from pathlib import Path

from tessera.recomb.run import RecombParams, run_recomb

_LOG = logging.getLogger("tessera.test")


def _synthetic_msa(path: Path) -> None:
    """Cowpox-like query with a variola segment spliced into the middle."""
    rng = random.Random(1)
    base = "".join(rng.choice("ACGT") for _ in range(6000))

    def mutate(seq: str, frac: float) -> str:
        chars = list(seq)
        for i in range(len(chars)):
            if rng.random() < frac:
                chars[i] = rng.choice("ACGT")
        return "".join(chars)

    cowpox = mutate(base, 0.02)
    variola = mutate(base, 0.06)
    query = list(cowpox)
    query[2000:4000] = list(variola[2000:4000])  # recombinant segment
    records = {
        "query": "".join(query),
        "cowpox": cowpox,
        "variola": variola,
        "other": mutate(base, 0.18),
    }
    with open(path, "w") as fo:
        for name, seq in records.items():
            fo.write(f">{name}\n{seq}\n")


def test_recomb_pipeline_outputs_and_ranking(tmp_path: Path) -> None:
    msa = tmp_path / "msa.fasta"
    _synthetic_msa(msa)
    out = tmp_path / "out"
    params = RecombParams(
        msa=msa, output=out, query="query",
        window_size=1000, window_step=100, plot_format="png",
    )
    run_recomb(params, _LOG)

    # renamed tables, plots and the self-contained report are written
    assert (out / "similarity_stats.tsv").exists()
    assert (out / "window_winners.tsv").exists()
    assert (out / "recombination_regions.tsv").exists()
    assert (out / "report.html").exists()
    assert any(out.glob("similarity_top*.png"))
    assert (out / "similarity_pair.png").exists()

    # the raw per-window matrix carries both coordinate systems
    windows_header = (out / "similarity_windows.tsv").read_text().splitlines()[0].split("\t")
    assert windows_header[:3] == ["msa_position", "query_position", "winner"]

    # the recombinant query is closest to cowpox overall, variola next
    winners = (out / "window_winners.tsv").read_text().splitlines()[1:]
    ranked = [line.split("\t")[0] for line in winners]
    assert ranked[0] == "cowpox"
    assert "variola" in ranked[:2]

    # a recombinant region over the variola insert is called, with query coords
    region_lines = (out / "recombination_regions.tsv").read_text().splitlines()
    assert len(region_lines) >= 2  # header + at least one region
    fields = dict(zip(region_lines[0].split("\t"), region_lines[1].split("\t"), strict=True))
    assert fields["minor_parent"] == "variola"
    assert fields["major_parent"] == "cowpox"
    assert int(fields["query_start"]) < int(fields["query_end"])

    # the report embeds the interactive plot and the region table
    report = (out / "report.html").read_text()
    assert "Recombinant regions" in report
    assert "variola" in report


def _clonal_msa(path: Path, seed: int = 5) -> None:
    """A non-recombinant query against a redundant panel.

    Five near-equidistant relatives (~2 % from a shared ancestor, so ~3 % from each
    other) plus three outgroups. Nothing recombined anywhere, so the scan must report
    no recombinant region -- this is the panel shape in which a caller that only
    counts window votes flips winner between windows by chance.
    """
    rng = random.Random(seed)
    base = "".join(rng.choice("ACGT") for _ in range(8000))

    def mutate(seq: str, frac: float) -> str:
        chars = list(seq)
        for i in range(len(chars)):
            if rng.random() < frac:
                chars[i] = rng.choice("ACGT")
        return "".join(chars)

    records = {"query": mutate(base, 0.02)}
    for i in range(5):
        records[f"sib{i}"] = mutate(base, 0.02)
    for i in range(3):
        records[f"out{i}"] = mutate(base, 0.15)
    with open(path, "w") as fo:
        for name, seq in records.items():
            fo.write(f">{name}\n{seq}\n")


def test_agreement_gate_clears_a_clonal_panel(tmp_path: Path) -> None:
    """Specificity regression for the agreement gate on a redundant panel.

    The panel below is deliberately redundant -- five near-equidistant relatives -- which
    is the regime where a caller that only counts window votes flips winner by chance.
    `--min-methods 2` must clear those: a genuine tract is found by several callers, a
    chance flip by one.

    The gate is *not* the default (see `RecombParams.min_methods`): on the curated panels
    the hybrid harness builds it costs true detections and gains nothing, so it is opted
    into here rather than assumed.

    Coverage-gap rows (``donor_absent``) are excluded -- they are legitimate "no close
    reference here" flags, matching ``run_hybrids.py::_score_neg_pure``.
    """
    msa = tmp_path / "clonal.fasta"
    _clonal_msa(msa)
    out = tmp_path / "out"
    run_recomb(RecombParams(msa=msa, output=out, query="query",
                            window_size=1000, window_step=100, plot_format="png",
                            min_methods=2), _LOG)

    lines = [x for x in (out / "recombination_regions.tsv").read_text().splitlines() if x]
    header = lines[0].split("\t")
    called = [dict(zip(header, x.split("\t"), strict=True)) for x in lines[1:]]
    claimed = [r for r in called if r.get("donor_absent") != "yes"]
    assert claimed == [], f"false positive(s) on clonal data: {claimed}"
