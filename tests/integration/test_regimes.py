"""Pipeline-level regressions for the regimes simulated by validation/run_regimes.py.

One seed each at reduced size, so they run in the fast suite; the harness itself
reports rates over replicates.
"""

from __future__ import annotations

import importlib.util
import logging
import sys
from pathlib import Path

import numpy as np

from tessera.recomb.analyze import analyze
from tessera.recomb.regions import RegionParams, call_regions
from tessera.recomb.similarity import compute_similarity

from ..conftest import write_fasta

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_regimes.py"
_SPEC = importlib.util.spec_from_file_location("run_regimes", _PATH)
rg = importlib.util.module_from_spec(_SPEC)
sys.modules["run_regimes"] = rg
_SPEC.loader.exec_module(rg)

WINDOW, STEP = 1000, 100
_LOG = logging.getLogger("tessera.test")


def _hmm_regions(result, window: int = WINDOW):
    params = RegionParams.with_defaults(window, method="hmm")
    regions, major, _ = call_regions(result, analyze(result), window, params)
    return regions, major


# --- HMM region span -------------------------------------------------------

def test_hmm_region_stops_where_donor_and_major_become_identical(tmp_path) -> None:
    # The query differs from the sibling only over `span`; downstream of it parent_A
    # (the donor relative to the sibling) and the sibling are identical.
    seqs, (lo, hi) = rg.simulate_sibling_tract(5000)
    msa = write_fasta(tmp_path / "sib.fasta", seqs)
    result = compute_similarity(str(msa), rg.QUERY, WINDOW, STEP)
    regions, major = _hmm_regions(result)
    assert major == "sibling"
    assert [r.minor_parent for r in regions] == ["parent_A"]
    region = regions[0]
    assert abs(region.msa_start - lo) <= STEP and abs(region.msa_end - hi) <= STEP
    assert region.length_bp <= (hi - lo) + STEP  # was 11.7 kb for a 600 bp difference


def test_hmm_region_reaches_the_alignment_start_for_a_terminal_tract(tmp_path) -> None:
    rng = np.random.default_rng(1)
    root = rng.integers(0, 4, size=3000)
    parent_a, parent_b, other = (rg._evolve(root, 0.05, rng) for _ in range(3))
    query = parent_a.copy()
    query[:800] = parent_b[:800]  # the donor tract starts at column 0
    seqs = rg._as_text({rg.QUERY: query, "parent_A": parent_a, "parent_B": parent_b,
                        "other": other})
    msa = write_fasta(tmp_path / "term.fasta", seqs)
    result = compute_similarity(str(msa), rg.QUERY, 300, 30)
    regions, _ = _hmm_regions(result, window=300)
    assert [r.minor_parent for r in regions] == ["parent_B"]
    # the first window centre is column 150; the tract must not be reported from there
    assert regions[0].msa_start <= 30
    assert abs(regions[0].msa_end - 800) <= 30


def test_trimmed_region_reports_query_coordinates_across_a_query_gap(tmp_path) -> None:
    # 300 alignment columns the query does not have (a backbone insertion), upstream of
    # the region: MSA and query coordinates then differ by exactly those 300 columns.
    seqs, (lo, hi) = rg.simulate_sibling_tract(5000)
    gap_at, gap = 3000, 300
    padded = {}
    for label, seq in seqs.items():
        insert = "-" * gap if label == rg.QUERY else seqs["parent_A"][gap_at:gap_at + gap]
        padded[label] = seq[:gap_at] + insert + seq[gap_at:]
    msa = write_fasta(tmp_path / "gapped.fasta", padded)
    result = compute_similarity(str(msa), rg.QUERY, WINDOW, STEP)
    regions, _ = _hmm_regions(result)
    assert [r.minor_parent for r in regions] == ["parent_A"]
    region = regions[0]
    assert abs(region.msa_start - (lo + gap)) <= STEP
    assert region.query_start == region.msa_start - gap
    assert region.query_end == region.msa_end - gap
    assert region.length_bp == region.length_msa


def test_trimmed_hmm_region_still_merges_with_the_site_callers(tmp_path) -> None:
    # Agreement is counted on overlapping regions. A trimmed HMM region must still
    # overlap the 3SEQ / MaxChi tract, or --min-methods 2 would drop a real event.
    seqs, (lo, hi) = rg.simulate_sibling_tract(5000)
    rows = rg.scan(seqs, _LOG, min_methods=2)
    hits = [r for r in rows if int(r["query_start"]) < hi and int(r["query_end"]) > lo]
    assert len(hits) == 1
    assert {"hmm", "3seq", "maxchi"} <= set(hits[0]["methods"].split(","))
