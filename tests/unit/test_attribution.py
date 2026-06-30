"""Pure attribution-scoring helpers from the hybrid harness (no aligner/network).

The harness (``validation/run_hybrids.py``) is opt-in and not collected by pytest, so
its pure helpers -- the clade-distance metric, the attribution tiering, and the
source-removed per-clade consensus panel -- are loaded by path and unit-tested here.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_HARNESS = (
    Path(__file__).resolve().parents[2] / "validation" / "run_hybrids.py"
)


def _load_harness():
    spec = importlib.util.spec_from_file_location("run_hybrids", _HARNESS)
    module = importlib.util.module_from_spec(spec)
    # Register before exec: @dataclass with `from __future__ import annotations`
    # resolves cls.__module__ via sys.modules during class processing.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


rh = _load_harness()


def test_shared_clade_depth_counts_shared_prefix() -> None:
    assert rh.shared_clade_depth("A.D.1", "A.D.1.8") == 3
    assert rh.shared_clade_depth("A.D.1.8", "A.D.1.8") == 4
    assert rh.shared_clade_depth("A.3", "A.D.1.8") == 1  # share only top-level A
    assert rh.shared_clade_depth("B", "A.1") == 0
    assert rh.shared_clade_depth("?", "A.1") == 0


def test_shared_clade_depth_ignores_reassortment_suffix() -> None:
    # base_clade strips the /rYYYY re-emergence tag before comparing.
    assert rh.shared_clade_depth("Ebov-2018b/r2021b", "Ebov-2018b") == 1


def test_attribution_tier_exact_sibling_mismatch() -> None:
    assert rh.attribution_tier("A.D.1.8", "A.D.1.8") == "exact"
    assert rh.attribution_tier("A.D", "A.D.1.8") == "exact"      # hierarchical match
    assert rh.attribution_tier("A.3", "A.D.1.8") == "sibling"    # same top-level A
    assert rh.attribution_tier("B.1", "A.D.1.8") == "mismatch"   # different top-level
    assert rh.attribution_tier("?", "A.D.1.8") == "mismatch"


def test_consensus_label_matches_production_naming() -> None:
    assert rh.consensus_label("A.1") == "A.1_consensus"
    assert rh.consensus_label("clade 2") == "clade_2_consensus"


def test_consensus_panel_majority_per_clade(tmp_path) -> None:
    import logging

    reference = "AAAAAA"
    # Two clades; clade A members agree on a C at pos 2, clade B on a G at pos 4.
    # reconstruct_gapped applies <orig><1-based pos><new> mutations to the reference.
    tips = {
        "a1": ("A", ["A3C"]),
        "a2": ("A", ["A3C"]),
        "a3": ("A", ["A3T"]),   # outvoted -> consensus keeps C at pos 3
        "b1": ("B", ["A5G"]),
        "b2": ("B", ["A5G"]),
    }
    members = {"A": ["a1", "a2", "a3"], "B": ["b1", "b2"]}
    paths, label_to_clade = rh.consensus_panel(
        members, reference, tips, tmp_path, logging.getLogger("t")
    )
    assert {p.name for p in paths} == {"A_consensus.fasta", "B_consensus.fasta"}
    assert label_to_clade == {"A_consensus": "A", "B_consensus": "B"}
    seqs = {}
    for p in paths:
        lines = p.read_text().splitlines()
        seqs[lines[0].lstrip(">")] = lines[1]
    assert seqs["A_consensus"] == "AACAAA"   # majority C at pos 3
    assert seqs["B_consensus"] == "AAAAGA"   # G at pos 5


def test_consensus_panel_single_member_clade(tmp_path) -> None:
    import logging

    reference = "GGGG"
    tips = {"x": ("X", ["G2A"])}
    paths, _ = rh.consensus_panel(
        {"X": ["x"]}, reference, tips, tmp_path, logging.getLogger("t")
    )
    assert paths[0].read_text().splitlines()[1] == "GAGG"


def _write_regions(out_dir: Path, rows: list[dict]) -> None:
    header = ["minor_parent", "major_parent", "query_start", "query_end",
              "donor_absent", "methods"]
    lines = ["\t".join(header)]
    for r in rows:
        lines.append("\t".join(str(r.get(h, "")) for h in header))
    (out_dir / "recombination_regions.tsv").write_text("\n".join(lines) + "\n")


def _setup(tmp_path, clade_a="A.1", clade_b="A.D.1.8", divergence=6.6):
    return rh.CaseSetup(
        name="t", out=tmp_path, clade_a=clade_a, clade_b=clade_b, divergence=divergence,
        query=tmp_path / "q.fasta", query_label="q", q_start=2000, q_end=4000,
        window=300, step=30, sel_window=900, aligner="mafft", reference="",
        tips={}, pool=[], members_by_clade={},
    )


def test_score_regions_sibling_donor_is_reported_not_passed(tmp_path) -> None:
    # rsv_a-like: a donor region is recovered but attributed to a sibling A-sub-clade
    # (A.3) of the true donor A.D.1.8. Same top-level A -> donor_match rejects it
    # (no PASS), but the finer tier records it as 'sibling' with the right depth.
    _write_regions(tmp_path, [
        {"minor_parent": "ref_A3", "major_parent": "ref_A1",
         "query_start": 2100, "query_end": 3800, "donor_absent": "no",
         "methods": "hmm,bootscan"},
    ])
    clade = {"ref_A3": "A.3", "ref_A1": "A.1"}
    out = rh._score_regions(tmp_path, lambda lbl: clade.get(lbl, "?"),
                            _setup(tmp_path), n_refs=4, mode="bp", runtime=1.0)
    assert out["detected"] is True
    assert out["donor_ok"] is False           # sibling intra-A is not credited for PASS
    assert out["pass"] is False
    assert out["donor_obs"] == "A.3"
    assert out["donor_tier"] == "sibling"
    assert out["donor_depth"] == 1            # shares only top-level A with A.D.1.8


def test_score_regions_exact_donor_passes_with_agreement(tmp_path) -> None:
    # The donor region is attributed to the true donor clade and called by >1 method.
    _write_regions(tmp_path, [
        {"minor_parent": "ref_AD", "major_parent": "ref_A1",
         "query_start": 2100, "query_end": 3800, "donor_absent": "no",
         "methods": "hmm,3seq,barcode"},
    ])
    clade = {"ref_AD": "A.D.1.8", "ref_A1": "A.1"}
    out = rh._score_regions(tmp_path, lambda lbl: clade.get(lbl, "?"),
                            _setup(tmp_path), n_refs=4, mode="bp", runtime=1.0)
    assert out["donor_ok"] is True and out["backbone_ok"] is True
    assert out["pass"] is True and out["agree"] is True
    assert out["donor_tier"] == "exact" and out["donor_depth"] == 4


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
