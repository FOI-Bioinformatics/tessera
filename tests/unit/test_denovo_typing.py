"""The harness's de-novo ANI tip typing relabels tips by cluster (no binaries)."""

from __future__ import annotations

import importlib.util
import logging
import sys
from pathlib import Path

_PATH = Path(__file__).resolve().parents[2] / "validation" / "run_hybrids.py"
_SPEC = importlib.util.spec_from_file_location("run_hybrids", _PATH)
rh = importlib.util.module_from_spec(_SPEC)
sys.modules["run_hybrids"] = rh
_SPEC.loader.exec_module(rh)

LOG = logging.getLogger("test")


def test_denovo_type_tips_relabels_by_cluster(monkeypatch):
    # Three tips, no clade attribute (empty labels). Stub the ANI clustering to put A+B in
    # one lineage and C in another; the tips must come back relabelled by lineage, keeping
    # their mutation paths.
    tips = {"A.1": ("", ["m1"]), "B.1": ("", ["m2"]), "C.1": ("", ["m3"])}
    genomes = [Path("A.fasta"), Path("B.fasta"), Path("C.fasta"), Path("EX1.fasta")]

    monkeypatch.setattr(rh, "_ani_edges", lambda gs, *, threshold, logger: [])
    monkeypatch.setattr(rh, "cluster_by_ani",
                        lambda gs, edges: {"denovo_1": [Path("A.fasta"), Path("B.fasta")],
                                           "denovo_2": [Path("C.fasta")]})

    typed = rh._denovo_type_tips(genomes, tips, 90.0, LOG)
    assert typed == {"A.1": ("denovo_1", ["m1"]),
                     "B.1": ("denovo_1", ["m2"]),
                     "C.1": ("denovo_2", ["m3"])}


def test_denovo_type_tips_clusters_only_tip_backed_genomes(monkeypatch):
    # The pool also carries example genomes (no tip entry); only tip-backed genomes are
    # passed to the clusterer.
    tips = {"A.1": ("", ["m1"]), "B.1": ("", ["m2"])}
    genomes = [Path("A.fasta"), Path("B.fasta"), Path("EX_example.fasta")]
    seen = {}

    def fake_edges(gs, *, threshold, logger):
        seen["genomes"] = [g.name for g in gs]
        return []

    monkeypatch.setattr(rh, "_ani_edges", fake_edges)
    monkeypatch.setattr(rh, "cluster_by_ani",
                        lambda gs, edges: {"denovo_1": [Path("A.fasta")],
                                           "denovo_2": [Path("B.fasta")]})

    rh._denovo_type_tips(genomes, tips, 90.0, LOG)
    assert seen["genomes"] == ["A.fasta", "B.fasta"]  # EX_example.fasta excluded
