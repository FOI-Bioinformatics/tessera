"""Pure donor re-attribution (no aligner/network)."""

from __future__ import annotations

import csv
import random
from pathlib import Path

import numpy as np

from tessera.recomb.reattribute import reattribute_donors
from tessera.recomb.regions import Region
from tessera.recomb.run import RecombParams, run_recomb

from ..conftest import write_fasta


class _Result:
    """Minimal stand-in for WindowSimilarity: reattribute_donors uses only rows/query."""

    def __init__(self, rows: dict[str, np.ndarray], query: str):
        self.rows = rows
        self.query = query


def _enc(seq: str) -> np.ndarray:
    return np.frombuffer(seq.encode(), dtype=np.uint8)


def _region(minor: str, major: str, lo: int, hi: int) -> Region:
    return Region(minor_parent=minor, major_parent=major, msa_start=lo, msa_end=hi,
                  query_start=lo, query_end=hi, n_windows=1,
                  mean_sim_minor=0.5, mean_sim_major=0.9, margin=0.0)


# Width-20 panel. Over span [5:15) clade C matches the query fully; the current donor
# (clade B) matches 6/10; clade A is the backbone (matches everywhere).
_QUERY = "A" * 20
_A = "A" * 20                                  # backbone member
_B = "AAAAA" + "TTTT" + "A" * 11               # 4 mismatches in the span -> 6/10
_C = "A" * 20                                  # full match in the span


def _panel() -> tuple[_Result, dict[str, str]]:
    rows = {"q": _enc(_QUERY), "a1": _enc(_A), "b1": _enc(_B), "c1": _enc(_C)}
    lineage_map = {"a1": "A", "b1": "B", "c1": "C"}
    return _Result(rows, "q"), lineage_map


def test_reattributes_to_better_clade_above_margin():
    result, lm = _panel()
    regions = [_region("b1", "a1", 5, 15)]
    out = reattribute_donors(regions, result, lm, "A", margin=0.1, min_sites=4)
    assert out[0].minor_parent == "c1"               # re-labelled to clade C's rep
    assert out[0].major_parent == "a1"               # backbone untouched
    assert out[0].mean_sim_minor == 1.0


def test_unchanged_below_margin():
    result, lm = _panel()
    regions = [_region("b1", "a1", 5, 15)]
    out = reattribute_donors(regions, result, lm, "A", margin=0.5, min_sites=4)
    assert out[0].minor_parent == "b1"               # 0.4 gain < 0.5 margin


def test_never_reattributes_to_the_backbone_clade():
    # Only clades A (backbone) and B; A scores best in the span but is excluded, so the
    # current donor B is kept rather than re-labelled to A.
    rows = {"q": _enc(_QUERY), "a1": _enc(_A), "b1": _enc(_B)}
    lm = {"a1": "A", "b1": "B"}
    out = reattribute_donors([_region("b1", "a1", 5, 15)], _Result(rows, "q"), lm,
                             "A", margin=0.1, min_sites=4)
    assert out[0].minor_parent == "b1"


def test_noop_without_lineage_map():
    result, _ = _panel()
    regions = [_region("b1", "a1", 5, 15)]
    assert reattribute_donors(regions, result, None, "A") is regions


def test_noop_with_fewer_than_two_clades():
    rows = {"q": _enc(_QUERY), "a1": _enc(_A)}
    regions = [_region("a1", "a1", 5, 15)]
    out = reattribute_donors(regions, _Result(rows, "q"),
                             {"a1": "A"}, "A", margin=0.1, min_sites=4)
    assert out is regions
    assert out[0].minor_parent == "a1"


def test_skips_span_with_too_few_comparable_sites():
    result, lm = _panel()
    regions = [_region("b1", "a1", 5, 8)]            # 3 columns < min_sites
    out = reattribute_donors(regions, result, lm, "A", margin=0.1, min_sites=4)
    assert out[0].minor_parent == "b1"


def test_keeps_an_untyped_donor_that_cannot_be_scored():
    """A donor absent from the lineage map has no consensus to score, so there is
    nothing to compare a rival clade against and re-attribution must decline.

    Treating the unscorable donor as a similarity of zero would let any clade clear
    the margin -- replacing a better donor with a worse one.
    """
    # 'x1' matches the query over the whole span but is deliberately untyped;
    # clade C also matches fully, clade B only 6/10.
    rows = {"q": _enc(_QUERY), "a1": _enc(_A), "b1": _enc(_B), "c1": _enc(_C),
            "x1": _enc(_QUERY)}
    lm = {"a1": "A", "b1": "B", "c1": "C"}           # x1 absent from the map
    out = reattribute_donors([_region("x1", "a1", 5, 15)], _Result(rows, "q"), lm,
                             "A", margin=0.1, min_sites=4)
    assert out[0].minor_parent == "x1"


# --- the run's other outputs follow the re-attribution ---------------------

def _reattribution_panel(tmp_path: Path) -> tuple[Path, dict[str, str]]:
    """A cowpox backbone with a variola insert at 2000-4000. The insert is copied from
    ``variolaA``, but ``variolaA`` shares clade V with two distant genomes, so clade V's
    consensus matches the insert worse than clade W's (``variolaB`` alone)."""
    rng = random.Random(11)

    def mutate(seq: str, frac: float) -> str:
        chars = list(seq)
        for i in range(len(chars)):
            if rng.random() < frac:
                chars[i] = rng.choice("ACGT")
        return "".join(chars)

    base = "".join(rng.choice("ACGT") for _ in range(6000))
    cowpox, variola = mutate(base, 0.03), mutate(base, 0.08)
    variola_b = mutate(variola, 0.01)
    query = list(cowpox)
    query[2000:4000] = list(variola[2000:4000])
    msa = write_fasta(tmp_path / "panel.fasta", {
        "query": "".join(query), "cowpox": cowpox, "variolaA": variola,
        "variolaB": variola_b, "junk1": mutate(base, 0.15), "junk2": mutate(base, 0.15),
    })
    lineages = {"cowpox": "CPX", "variolaA": "V", "junk1": "V", "junk2": "V", "variolaB": "W"}
    return msa, lineages


def test_method_comparison_names_the_reattributed_donor(tmp_path: Path, logger) -> None:
    """Re-attribution re-labelled the region but not the per-method breakdown, so
    `recombination_methods.tsv` and the report's method table kept the old donor."""
    msa, lineages = _reattribution_panel(tmp_path)
    out = tmp_path / "out"
    run_recomb(
        RecombParams(msa=msa, output=out, query="query", plot_format="png",
                     lineage_map=lineages, reattribute_donors=True, cluster_lineages=False),
        logger,
    )
    regions = list(csv.DictReader((out / "recombination_regions.tsv").open(), delimiter="\t"))
    methods = list(csv.DictReader((out / "recombination_methods.tsv").open(), delimiter="\t"))
    assert [r["minor_parent"] for r in regions] == ["variolaB"]  # re-attributed from variolaA
    assert [m["minor_parent"] for m in methods] == ["variolaB"]
    assert [(m["query_start"], m["query_end"]) for m in methods] == [
        (r["query_start"], r["query_end"]) for r in regions
    ]
