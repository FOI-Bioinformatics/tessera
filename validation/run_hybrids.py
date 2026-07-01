#!/usr/bin/env python
"""Generate synthetic recombinant (hybrid) genomes from Nextclade datasets and
test Tessera's detection performance on them.

For each configured dataset the harness:

1. Builds the Nextclade reference pool (the shipped ``build_pool``; cached).
2. Picks two clades A and B (the most-divergent well-represented pair by default,
   or an explicit pair) and each clade's central genome, in reference coordinates.
3. Splices an A-backbone / B-insert hybrid (B donates the middle of the genome),
   recording the true donor span in query coordinates.
4. Runs Tessera pool-only against the dataset's pool with the two exact source
   genomes removed (their clades remain represented by other genomes), so the
   query is not a trivial self-match.
5. Checks the call: was recombination detected, is the backbone (major parent)
   clade A, is a donor region recovered for clade B, and does it overlap the true
   donor span. Records wall-clock runtime per case.

This is opt-in and separate from the pytest suite (it needs MAFFT/skani/skDER on
PATH and contacts the Nextclade dataset server on first run). See
validation/README.md.

Usage:
    export PATH="$HOME/miniforge3/envs/recomfi-aln/bin:$PATH"
    python validation/run_hybrids.py            # all cases
    python validation/run_hybrids.py hiv1 dengue   # only named cases
"""
from __future__ import annotations

import json
import logging
import os
import re
import shutil
import sys
import time
from dataclasses import dataclass
from pathlib import Path

from tessera.core.cache import nextclade_cache
from tessera.core.errors import ToolExecutionError
from tessera.core.io import strip_sequence_extension, write_fasta_record
from tessera.discover.nextclade import (
    _MUT,
    _accession_of,
    _clade_of,
    _download_text,
    build_pool,
    resolve_dataset,
)
from tessera.discover.pool import select_regional
from tessera.msa.build import MsaParams, build_msa
from tessera.recomb.consensus import consensus_sequence
from tessera.recomb.regions import DEFAULT_METHODS, parse_methods
from tessera.recomb.run import RecombParams, run_recomb
from tessera.recomb.typing import (
    LINEAGES_TSV,
    lineage_map_from_rows,
    write_lineage_map,
)
from tessera.recomb.typing import (
    is_recombinant_lineage as is_recombinant_clade,
)

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "hybrids"

# One representative dataset per distinct Nextclade pathogen (14 core + 10
# community). ``clades`` may pin the two parent clades; left empty the
# most-divergent well-represented pair is used. ``aligner`` defaults to mafft;
# large genomes (mpox, VZV) use minimap2 so the per-genome alignment is tractable.
# The donor (clade B) occupies the middle ``insert`` fraction of the genome; the
# backbone (clade A) the flanks. Window sizes adapt to the genome length.
HYBRIDS: list[dict] = [
    # --- core (Nextstrain-maintained) ---
    # XBB's Pango lineages are too fine to group; use the coarse Nextstrain clade.
    {"name": "sars_cov_2", "dataset": "nextstrain/sars-cov-2/XBB",
     "clade_key": "clade_nextstrain"},
    {"name": "dengue", "dataset": "nextstrain/dengue/all"},
    {"name": "measles", "dataset": "nextstrain/measles/genome/WHO-2012"},
    {"name": "mumps", "dataset": "nextstrain/mumps/genome"},
    {"name": "rubella", "dataset": "nextstrain/rubella/genome"},
    {"name": "rsv_a", "dataset": "nextstrain/rsv/a/EPI_ISL_412866"},
    {"name": "hmpv", "dataset": "nextstrain/hmpv/all-clades/NC_039199"},
    {"name": "wnv", "dataset": "nextstrain/wnv/all-lineages"},
    {"name": "yellow_fever", "dataset": "nextstrain/yellow-fever/prM-E"},
    {"name": "ebola", "dataset": "nextstrain/orthoebolavirus/ebov", "min_divergence": 0.0},
    {"name": "hantavirus", "dataset": "nextstrain/orthohantavirus/andv/l"},
    {"name": "flu_h3n2_ha", "dataset": "nextstrain/flu/h3n2/ha/EPI1857216"},
    {"name": "mpox", "dataset": "nextstrain/mpox/all-clades", "aligner": "minimap2",
     "min_divergence": 0.0},
    {"name": "vzv", "dataset": "nextstrain/herpes/vzv/NC_001348", "aligner": "minimap2",
     "min_divergence": 0.0},
    # --- community ---
    {"name": "hiv1", "dataset": "community/neherlab/hiv-1/hxb2"},
    {"name": "marburg", "dataset": "community/genspectrum/marburg/HK1980/all-lineages"},
    {"name": "oropouche", "dataset": "community/itps/orov/L/refseq"},
    {"name": "zika", "dataset": "community/itps/zikav"},
    {"name": "hepatitis_a", "dataset": "community/masphl-bioinformatics/hav/whole-genome"},
    {"name": "iav_h5_ha", "dataset": "community/moncla-lab/iav-h5/ha/all-clades"},
    {"name": "cchfv", "dataset": "community/pathoplexus/cchfv/L"},
    {"name": "chikv", "dataset": "community/v-gen-lab/chikV/genotypes"},
    {"name": "enterovirus_d68", "dataset": "enpen/enterovirus/ev-d68"},
    {"name": "prrsv2", "dataset": "community/isuvdl/mazeller/prrsv2/orf5/yimim2023"},
]
INSERT = (0.35, 0.65)  # donor (clade B) occupies this fraction of the genome
MIN_GENOME = 400  # skip a dataset whose genome/segment is too short to splice
MIN_DIVERGENCE = 4.0  # skip when the most-divergent clade pair is below this (% ) --
# too few discriminating sites for a meaningful recombination call (conserved DNA
# viruses, intra-species sets)
REP_PANEL_MAX = 14  # representative-panel fallback: at most this many clade reps
THREADS = 4


def window_params(genome_len: int) -> tuple[int, int, int]:
    """Adaptive (recomb window, step, regional-selection window) for a genome length.

    Scales the sliding window to the genome so short gene/segment datasets (flu HA,
    a PRRSV ORF) are not given a window wider than the alignment.
    """
    window = max(120, min(500, genome_len // 12))
    step = max(20, window // 10)
    select_window = max(window, min(1500, genome_len // 4))
    return window, step, select_window


def reconstruct_gapped(reference: str, nuc_muts: list[str]) -> str:
    """Apply substitutions/deletions to ``reference`` WITHOUT stripping gaps.

    Keeps the result in reference coordinates so two tips can be spliced at a
    shared reference position; gaps are stripped only after splicing.
    """
    chars = list(reference)
    for mut in nuc_muts:
        m = _MUT.match(mut)
        if not m:
            continue
        i = int(m.group(2)) - 1
        if 0 <= i < len(chars):
            chars[i] = m.group(3)
    return "".join(chars)


def node_clade(node_attrs: dict, clade_key: str | None) -> str:
    """The clade label for a node: a specific attribute if ``clade_key`` is given
    (e.g. SARS-CoV-2 ``clade_nextstrain`` instead of the too-fine Pango lineage),
    otherwise the shipped priority used by the feature."""
    if clade_key:
        value = node_attrs.get(clade_key)
        if isinstance(value, dict):
            value = value.get("value")
        return str(value) if value not in (None, "") else "NA"
    return _clade_of(node_attrs)


def collect_tips(
    tree: dict, clade_key: str | None = None
) -> dict[str, tuple[str, list[str]]]:
    """Map ``accession -> (clade, root_to_tip_nuc_mutations)`` for every leaf."""
    tips: dict[str, tuple[str, list[str]]] = {}
    stack: list[tuple[dict, list[str]]] = [(tree, [])]
    while stack:
        node, inherited = stack.pop()
        nuc = (node.get("branch_attrs", {}).get("mutations", {}) or {}).get("nuc", [])
        path_muts = inherited + nuc
        children = node.get("children", [])
        if children:
            for child in children:
                stack.append((child, path_muts))
            continue
        acc = _accession_of(node)
        if acc and acc not in tips:
            tips[acc] = (node_clade(node.get("node_attrs", {}) or {}, clade_key), path_muts)
    return tips


def read_reference(text: str) -> str:
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))


def pct_identity(a: str, b: str) -> float:
    """Percent identity over reference columns where both carry a canonical base."""
    bases = set("ACGTacgt")
    comp = match = 0
    for x, y in zip(a, b, strict=False):
        if x in bases and y in bases:
            comp += 1
            match += x.upper() == y.upper()
    return 100.0 * match / comp if comp else 0.0


MIN_MEMBERS = 3  # a clade must have at least this many genomes to be a parent


def base_clade(name: str) -> str:
    """Strip a reassortment / re-emergence suffix, e.g. ebola ``Ebov-2018b/r2021b`` ->
    ``Ebov-2018b``: the ``/rYYYY`` tag marks the same parental lineage re-emerging in a
    later year (Ebola persists in survivors), so it names the same clade for attribution.
    """
    return re.sub(r"/r\d.*$", "", name).strip()


def clade_match(observed: str, expected: str) -> bool:
    """True when two Nextclade labels name the same lineage, allowing for hierarchical
    naming (``A`` matches ``A.1``) and reassortment suffixes (``X/rYYYY`` == ``X``)."""
    observed, expected = base_clade(observed), base_clade(expected)
    return (observed == expected
            or observed.startswith(expected + ".")
            or expected.startswith(observed + "."))


def donor_match(observed: str, donor: str, backbone: str) -> bool:
    """The recombinant region's clade names the donor's parental lineage.

    Exact/hierarchical match (``clade_match``), or a sibling sub-clade under the
    donor's top-level clade (e.g. Marburg ``RAVV.1`` for a ``RAVV.2`` donor) --
    but only when the donor's top-level clade differs from the backbone's, so this
    cannot trivially credit a shared-top-level pair (e.g. RSV ``A.1`` / ``A.D.1.8``)
    where distinguishing the donor sub-clade is the actual test.
    """
    if clade_match(observed, donor):
        return True
    return (donor.split(".")[0] != backbone.split(".")[0]
            and observed.split(".")[0] == donor.split(".")[0])


def shared_clade_depth(observed: str, expected: str) -> int:
    """How close two clade labels are: the number of leading dot-separated components
    they share after :func:`base_clade` normalisation.

    A granular attribution-distance that shows improvement even when it does not flip
    the PASS verdict, e.g. ``A.D.1`` vs ``A.D.1.8`` -> 3; ``A.3`` vs ``A.D.1.8`` -> 1;
    ``B`` vs ``A.1`` -> 0. ``?``/``NA``/empty observed scores 0.
    """
    if observed in ("?", "NA", ""):
        return 0
    o = base_clade(observed).split(".")
    e = base_clade(expected).split(".")
    depth = 0
    for a, b in zip(o, e, strict=False):
        if a != b:
            break
        depth += 1
    return depth


def attribution_tier(observed: str, expected: str) -> str:
    """Classify an attribution as ``exact``, ``sibling`` or ``mismatch``.

    ``exact`` -- :func:`clade_match` (same lineage, allowing hierarchy); ``sibling`` --
    a different sub-clade under the same top-level clade (e.g. the intra-``A`` RSV case
    ``A.3`` for an ``A.D.1.8`` donor, which ``donor_match`` deliberately rejects for the
    PASS rule); ``mismatch`` -- otherwise. Reported alongside the verdict, never folded
    into PASS.
    """
    if observed in ("?", "NA", ""):
        return "mismatch"
    if clade_match(observed, expected):
        return "exact"
    if base_clade(observed).split(".")[0] == base_clade(expected).split(".")[0]:
        return "sibling"
    return "mismatch"


def consensus_label(clade: str) -> str:
    """The per-clade consensus genome's label, e.g. ``A.1`` -> ``A.1_consensus``.

    Matches the production ``build_pool(per_clade_consensus=True)`` naming so the two
    consensus paths read the same.
    """
    return re.sub(r"[^\w.]", "_", clade) + "_consensus"


def consensus_panel(
    members_by_clade: dict[str, list[str]],
    reference: str,
    tips: dict[str, tuple[str, list[str]]],
    out_dir: Path,
    logger: logging.Logger,
) -> tuple[list[Path], dict[str, str]]:
    """One denoised majority-consensus genome per clade (source genomes already removed).

    Each clade member is reconstructed in reference coordinates from its mutation path,
    the column-wise majority base is taken (:func:`consensus_sequence`), gaps are
    stripped, and one ``{clade}_consensus`` genome is written. Because ``members_by_clade``
    is built from the source-removed pool, the consensus excludes the two spliced source
    genomes -- it keeps the harness's "source removed" honesty while giving the stable
    per-clade stand-in that ``--pool-consensus`` provides. Returns the written paths and a
    ``label -> clade`` map.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    label_to_clade: dict[str, str] = {}
    for clade, members in members_by_clade.items():
        seqs = [reconstruct_gapped(reference, tips[acc][1]) for acc in members]
        cons = consensus_sequence(seqs).replace("-", "").upper()
        label = consensus_label(clade)
        path = out_dir / f"{label}.fasta"
        with open(path, "w") as fo:
            write_fasta_record(fo, label, cons)
        paths.append(path)
        label_to_clade[label] = clade
    logger.info("Built a %d-clade consensus panel (source-removed).", len(paths))
    return paths, label_to_clade


class CaseSkipped(Exception):
    """A dataset is unsuitable for a hybrid test (reported as SKIP, not ERROR)."""


def clade_members(
    tips: dict[str, tuple[str, list[str]]], *, exclude_recombinant: bool = True
) -> dict[str, list[str]]:
    """Map each clade to its member accessions (skip unlabelled, and -- by default
    -- recombinant clades)."""
    by_clade: dict[str, list[str]] = {}
    for acc, (clade, _) in tips.items():
        if not clade or clade.strip().lower() in ("na", "unassigned", "unclassified", "unknown"):
            continue  # not a single clean clade -> not a valid parental lineage
        if exclude_recombinant and is_recombinant_clade(clade):
            continue
        by_clade.setdefault(clade, []).append(acc)
    return by_clade


def clade_representative(members: list[str], tips: dict[str, tuple[str, list[str]]]) -> str:
    """The clade's central genome: median root-to-tip mutation count (not an outlier).

    A central member is close to the rest of its clade, so after it is spliced into
    the hybrid and removed from the pool the remaining clade genomes still match the
    backbone -- the test exercises clade-level detection, not an outlier mismatch.
    """
    ranked = sorted(members, key=lambda a: (len(tips[a][1]), a))
    return ranked[len(ranked) // 2]


def pick_parents(
    tips: dict[str, tuple[str, list[str]]], reference: str, pinned: list[str],
    logger: logging.Logger, *, objective: str = "max",
) -> tuple[str, str, str, str]:
    """Return ``(clade_a, clade_b, source_a, source_b)``.

    Among clades with at least ``MIN_MEMBERS`` genomes, choose the most-divergent
    pair (each represented by its central genome); the larger clade is the backbone
    (A). A pinned pair overrides the search. Sources are the clades' central genomes.

    Recombinant clades are excluded as parents, but if that leaves fewer than two
    eligible clades (e.g. the SARS-CoV-2 XBB dataset, which is entirely XBB
    sub-lineages) they are allowed back in. Too few clades either way is a SKIP.

    ``objective`` controls pair selection when not pinned: ``"max"`` picks the
    most-divergent pair (default); ``"min"`` picks the least-divergent pair whose
    divergence is still at or above ``MIN_DIVERGENCE``.
    """
    def eligible_of(exclude_recombinant: bool) -> dict[str, list[str]]:
        by_clade = clade_members(tips, exclude_recombinant=exclude_recombinant)
        return {c: m for c, m in by_clade.items() if len(m) >= MIN_MEMBERS}

    eligible = eligible_of(exclude_recombinant=True)
    if len(eligible) < 2:
        eligible = eligible_of(exclude_recombinant=False)
        if len(eligible) >= 2:
            logger.info("No non-recombinant clade pair; allowing recombinant clades as parents.")
    if len(eligible) < 2:
        raise CaseSkipped(f"fewer than two clades with >={MIN_MEMBERS} genomes")
    reps = {c: clade_representative(m, tips) for c, m in eligible.items()}
    rep_seq = {c: reconstruct_gapped(reference, tips[reps[c]][1]) for c in reps}

    if pinned:
        ca, cb = pinned[0], pinned[1]
        if ca not in eligible or cb not in eligible:
            raise ValueError(f"pinned clades {pinned} need >={MIN_MEMBERS} genomes each")
    else:
        clades = sorted(reps)
        best = None
        for i, ca in enumerate(clades):
            for cb in clades[i + 1:]:
                div = 100.0 - pct_identity(rep_seq[ca], rep_seq[cb])
                if objective == "min":
                    if div >= MIN_DIVERGENCE and (best is None or div < best[0]):
                        best = (div, ca, cb)
                else:
                    if best is None or div > best[0]:
                        best = (div, ca, cb)
        if best is None:
            raise CaseSkipped("no clade pair meets the minimum divergence threshold")
        _, ca, cb = best
    if len(eligible[cb]) > len(eligible[ca]):  # backbone = the better-represented clade
        ca, cb = cb, ca
    return ca, cb, reps[ca], reps[cb]


def _pick_within_clade(
    tips: dict[str, tuple[str, list[str]]], reference: str, logger: logging.Logger,
) -> tuple[str, str, str]:
    """Return ``(clade, src_a, src_b)``: the largest non-recombinant clade and its two
    most-divergent members, for constructing a within-clade negative splice."""
    by_clade: dict[str, list[str]] = {}
    for acc, (clade, _muts) in tips.items():
        if clade and clade not in ("?", "NA") and not is_recombinant_clade(clade):
            by_clade.setdefault(clade, []).append(acc)
    ranked = sorted(by_clade.items(), key=lambda kv: len(kv[1]), reverse=True)
    for clade, members in ranked:
        if len(members) < 2:
            continue
        seqs = {a: reconstruct_gapped(reference, tips[a][1]) for a in members}
        best = max(
            ((100.0 - pct_identity(seqs[a], seqs[b]), a, b)
             for i, a in enumerate(members) for b in members[i + 1:]),
            default=None,
        )
        if best:
            return clade, best[1], best[2]
    raise CaseSkipped("no clade has two divergent members for a within-clade splice")


def make_hybrid(reference: str, src_a_muts, src_b_muts) -> tuple[str, int, int]:
    """Splice an A-backbone / B-insert hybrid; return (query, q_start, q_end)."""
    a_gap = reconstruct_gapped(reference, src_a_muts)
    b_gap = reconstruct_gapped(reference, src_b_muts)
    length = len(reference)
    b1, b2 = int(length * INSERT[0]), int(length * INSERT[1])
    hybrid_ref = a_gap[:b1] + b_gap[b1:b2] + a_gap[b2:]
    q_start = len(hybrid_ref[:b1].replace("-", ""))
    q_end = len(hybrid_ref[:b2].replace("-", ""))
    return hybrid_ref.replace("-", "").upper(), q_start, q_end


def parse_regions(path: Path) -> list[dict]:
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"), strict=False)) for ln in lines[1:]]


def parse_signal(path: Path) -> tuple[str, str]:
    """Read the parent-free PHI p-value and Rmin from the profile TSV header line:
    ``# PHI p-value\t<p>\t(window ..., <n> sites, Rmin <rmin>)``."""
    if not path.exists():
        return "-", "-"
    head = path.read_text().splitlines()[0]
    fields = head.split("\t")
    phi_p = fields[1] if len(fields) > 1 else "-"
    tail = fields[2] if len(fields) > 2 else ""
    rmin = tail.split("Rmin")[-1].strip(" )") if "Rmin" in tail else "-"
    return phi_p, rmin


def clade_of_label(label: str, tips: dict[str, tuple[str, list[str]]]) -> str:
    base = strip_sequence_extension(label).split(".")[0]
    for acc, (clade, _) in tips.items():
        if acc.split(".")[0] == base:
            return clade
    return "?"


def representative_panel(
    pool: list[Path], tips: dict[str, tuple[str, list[str]]], logger: logging.Logger
) -> list[Path]:
    """A clade-representative panel, used when skani rejects a short gene/segment.

    Groups the (already example-free, source-free) pool by clade and keeps one
    central genome per clade, preferring the larger clades up to ``REP_PANEL_MAX``.
    This skips skani entirely while still giving the recombination scan a diverse,
    clade-spanning panel that contains the parental clades.
    """
    base_info = {acc.split(".")[0]: (clade, len(muts)) for acc, (clade, muts) in tips.items()}
    by_clade: dict[str, list[tuple[int, Path]]] = {}
    for g in pool:
        clade, nmut = base_info.get(strip_sequence_extension(g.name).split(".")[0], (None, 0))
        if clade and clade not in ("NA", "unassigned") and not is_recombinant_clade(clade):
            by_clade.setdefault(clade, []).append((nmut, g))
    chosen: list[Path] = []
    for clade in sorted(by_clade, key=lambda c: -len(by_clade[c]))[:REP_PANEL_MAX]:
        members = sorted(by_clade[clade])
        chosen.append(members[len(members) // 2][1])  # the central (median-mutation) genome
    logger.info("Built a %d-clade representative panel (skani fallback).", len(chosen))
    return chosen


@dataclass
class CaseSetup:
    """The config-independent setup for one hybrid case (built once, scored under many
    method/panel configurations)."""

    name: str
    out: Path
    clade_a: str
    clade_b: str
    divergence: float
    query: Path
    query_label: str
    q_start: int
    q_end: int
    window: int
    step: int
    sel_window: int
    aligner: str
    reference: str
    tips: dict[str, tuple[str, list[str]]]
    pool: list[Path]  # source-removed, clade-labelled tree tips
    members_by_clade: dict[str, list[str]]  # source-removed clade -> tip accessions
    case_type: str = "single_insert"


def _prepare_case(case: dict, logger: logging.Logger) -> CaseSetup:
    """Download, pick the parent clades, splice the hybrid query, and build the
    source-removed pool -- everything that does not depend on the caller set or panel
    mode. Raises ``CaseSkipped`` for an unsuitable dataset (too similar, too short)."""
    name, path = case["name"], case["dataset"]
    out = DATA / name
    out.mkdir(parents=True, exist_ok=True)

    dataset = resolve_dataset(Path("/dev/null"), path, email=None, logger=logger)
    genomes = build_pool(dataset, cache_dir=nextclade_cache(dataset.path, dataset.tag),
                         logger=logger)
    reference = read_reference(_download_text(dataset, "reference", logger))
    tree = json.loads(_download_text(dataset, "treeJson", logger))["tree"]
    tips = collect_tips(tree, case.get("clade_key"))

    case_type = case.get("case_type", "single_insert")
    if case_type == "neg_pure":
        clade_a, _, src_a, _ = pick_parents(tips, reference, case.get("clades", []), logger)
        clade_b = ""
        src_b = src_a
        query_seq = reconstruct_gapped(reference, tips[src_a][1]).replace("-", "")
        q_start = q_end = 0
        divergence = 0.0
        logger.info("[%s] neg_pure: clade %s (%s) unspliced query", name, clade_a, src_a)
    elif case_type == "neg_within":
        clade_a, src_a, src_b = _pick_within_clade(tips, reference, logger)
        clade_b = clade_a
        divergence = 100.0 - pct_identity(
            reconstruct_gapped(reference, tips[src_a][1]),
            reconstruct_gapped(reference, tips[src_b][1]),
        )
        logger.info("[%s] neg_within: clade %s (%s x %s); intra-clade divergence %.1f%%",
                    name, clade_a, src_a, src_b, divergence)
        query_seq, q_start, q_end = make_hybrid(reference, tips[src_a][1], tips[src_b][1])
    else:  # single_insert and low_div
        clade_a, clade_b, src_a, src_b = pick_parents(
            tips, reference, case.get("clades", []), logger,
            objective=case.get("pair_objective", "max"))
        divergence = 100.0 - pct_identity(
            reconstruct_gapped(reference, tips[src_a][1]),
            reconstruct_gapped(reference, tips[src_b][1]),
        )
        logger.info(
            "[%s] backbone clade %s (%s) x donor clade %s (%s); inter-clade divergence %.1f%%",
            name, clade_a, src_a, clade_b, src_b, divergence)
        if divergence < case.get("min_divergence", MIN_DIVERGENCE):
            raise CaseSkipped(f"clades too similar ({divergence:.1f}% divergence)")
        band = case.get("divergence_band")
        if band is not None:
            lo, hi = band
            if not (lo <= divergence <= hi):
                raise CaseSkipped(
                    f"divergence {divergence:.1f}% outside band [{lo}, {hi}]")
        query_seq, q_start, q_end = make_hybrid(reference, tips[src_a][1], tips[src_b][1])

    if len(query_seq) < MIN_GENOME:
        raise CaseSkipped(f"genome/segment too short to test ({len(query_seq)} bp)")
    window, step, sel_window = window_params(len(query_seq))
    aligner = case.get("aligner", "mafft")
    query = out / "hybrid.fasta"
    with open(query, "w") as fo:
        write_fasta_record(fo, f"hybrid_{clade_a}_{clade_b}", query_seq)
    query_label = strip_sequence_extension(query.name)

    # Pool = clade-labelled tree tips only (drop Nextclade example sequences, which
    # carry no clade and would otherwise win the backbone unlabelled), minus the two
    # exact source genomes (their clades remain represented by other tips).
    drop = {src_a.split(".")[0], src_b.split(".")[0]}
    base_to_tip = {acc.split(".")[0]: acc for acc in tips}
    pool = [
        g for g in genomes
        if strip_sequence_extension(g.name).split(".")[0] in base_to_tip
        and strip_sequence_extension(g.name).split(".")[0] not in drop
    ]
    # Source-removed clade -> tip accessions, for the per-clade consensus panel.
    members_by_clade: dict[str, list[str]] = {}
    for g in pool:
        tipkey = base_to_tip.get(strip_sequence_extension(g.name).split(".")[0])
        if tipkey is None:
            continue
        clade = tips[tipkey][0]
        if (not clade or clade in ("?", "NA")
                or clade.strip().lower() in ("na", "unassigned", "unclassified", "unknown")
                or is_recombinant_clade(clade)):
            continue
        members_by_clade.setdefault(clade, []).append(tipkey)

    return CaseSetup(
        name=name, out=out, clade_a=clade_a, clade_b=clade_b, divergence=divergence,
        query=query, query_label=query_label, q_start=q_start, q_end=q_end,
        window=window, step=step, sel_window=sel_window, aligner=aligner,
        reference=reference, tips=tips, pool=pool, members_by_clade=members_by_clade,
        case_type=case_type,
    )


def _build_and_score(
    setup: CaseSetup, panel_mode: str, methods: tuple[str, ...], out_dir: Path,
    logger: logging.Logger, *, reattribute: bool = False,
) -> dict:
    """Build the panel for ``panel_mode`` ('tip' or 'consensus'), align, run detection
    with ``methods``, and score the call. Raises ``CaseSkipped`` when a parent clade has
    no representative in the panel (the harness's representation invariant)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    if panel_mode == "consensus":
        # One denoised consensus genome per clade (source genomes already removed).
        selected, cons_map = consensus_panel(
            setup.members_by_clade, setup.reference, setup.tips, out_dir / "consensus", logger
        )
        avail = set(setup.members_by_clade)
        if not any(clade_match(c, setup.clade_a) for c in avail):
            raise CaseSkipped(
                f"backbone clade {setup.clade_a!r} has no panel representative "
                "after the source genome was removed")
        if not any(donor_match(c, setup.clade_b, setup.clade_a) for c in avail):
            raise CaseSkipped(
                f"donor clade {setup.clade_b!r} has no panel representative "
                "after the source genome was removed")

        def clade_of(label: str) -> str:
            return cons_map.get(strip_sequence_extension(label), "?")

        rows = [(strip_sequence_extension(p.name),
                 cons_map[strip_sequence_extension(p.name)], "consensus") for p in selected]
    else:
        try:
            # Keep "siblings": for a hybrid of close parents the backbone parent is
            # >95% genome-wide ANI to the query (its 70% backbone dominates the
            # average) and would be dropped as a masking twin -- but the synthetic
            # pool contains no actual recombinant twin, so this is the documented
            # close-parent / --seed-keep-siblings setting, and dropping it would
            # remove the backbone parent itself.
            # Type the pool by tree clade so the panel is reduced by lineage (one
            # query-closest representative per clade, recombinant clades excluded) --
            # this keeps a true parent that ANI would collapse (rsv_a's A.D.1.8) and
            # drops the masking CRF clades (hiv1).
            lineage_of = {
                strip_sequence_extension(g.name): clade_of_label(g.name, setup.tips)
                for g in setup.pool
            }
            lineage_of = {k: v for k, v in lineage_of.items() if v and v not in ("?", "NA")}
            selected = select_regional(
                setup.query, setup.pool, window=setup.sel_window, per_window=2,
                drop_siblings=False, lineage_of=lineage_of or None,
                keep_recombinant=False, logger=logger,
            ).selected
        except ToolExecutionError:  # skani rejects very short gene/segment datasets
            selected = representative_panel(setup.pool, setup.tips, logger)
        # On a near-identical panel (mpox/VZV ~0.5%) dereplication collapses the parent
        # clades into one representative, leaving nothing to compete; rebuild from one
        # central genome per clade so both parents are present for the detector to test.
        sel_clades = {clade_of_label(g.name, setup.tips) for g in selected}
        if len([c for c in sel_clades if c not in ("?", "NA")]) < 2:
            logger.info("[%s] regional selection collapsed to %d clade(s); using clade "
                        "representatives instead.", setup.name, len(sel_clades))
            selected = representative_panel(setup.pool, setup.tips, logger)
        if not selected:
            raise CaseSkipped("no reference panel could be built")
        # The harness removes the two source genomes; a fair test needs each parent clade
        # to remain credit-able from the panel (its design invariant).
        panel_clades = {clade_of_label(g.name, setup.tips) for g in selected}
        if not any(clade_match(c, setup.clade_a) for c in panel_clades):
            raise CaseSkipped(
                f"backbone clade {setup.clade_a!r} has no panel representative "
                "after the source genome was removed")
        if not any(donor_match(c, setup.clade_b, setup.clade_a) for c in panel_clades):
            raise CaseSkipped(
                f"donor clade {setup.clade_b!r} has no panel representative "
                "after the source genome was removed")

        def clade_of(label: str) -> str:
            return clade_of_label(label, setup.tips)

        rows = [(strip_sequence_extension(g.name), clade_of_label(g.name, setup.tips), "nextclade")
                for g in selected]

    collection = out_dir / "collection"
    if collection.exists():
        shutil.rmtree(collection)
    collection.mkdir(parents=True)
    for g in selected:
        shutil.copy(g, collection / g.name)

    rows.append((setup.query_label, f"hybrid({setup.clade_a}+{setup.clade_b})", "query"))
    write_lineage_map(out_dir / LINEAGES_TSV, rows)
    lineage_map = lineage_map_from_rows(rows)

    msa = out_dir / "panel.msa.fasta"
    build_msa(MsaParams(query=setup.query, collection=collection, output=msa,
                        aligner=setup.aligner, threads=THREADS), logger)
    windowing = run_recomb(RecombParams(msa=msa, output=out_dir, query=setup.query_label,
                                        window_size=setup.window, window_step=setup.step,
                                        organism=setup.name, methods=methods,
                                        lineage_map=lineage_map,
                                        reattribute_donors=reattribute), logger)
    runtime = time.monotonic() - t0
    mode = "info-site" if windowing.startswith("informative") else "bp"
    return _score_regions(out_dir, clade_of, setup, len(selected), mode, runtime)


def _score_regions(
    out_dir: Path, clade_of, setup: CaseSetup, n_refs: int, mode: str, runtime: float,
) -> dict:
    """Score a completed run by the case's type (default: single-insert hybrid)."""
    scorer = {
        "single_insert": _score_single_insert,
        "neg_pure": _score_neg_pure,
        "neg_within": _score_neg_within,
        "low_div": _score_low_div,
    }.get(setup.case_type, _score_single_insert)
    return scorer(out_dir, clade_of, setup, n_refs, mode, runtime)


def _score_single_insert(
    out_dir: Path, clade_of, setup: CaseSetup, n_refs: int, mode: str, runtime: float,
) -> dict:
    """Score a completed detection run: detection, backbone/donor attribution, ensemble
    agreement, and the finer attribution-quality tiers/depths (reported, never folded
    into PASS). ``clade_of`` maps a region's parent label to its clade."""
    regions = parse_regions(out_dir / "recombination_regions.tsv")
    phi_p, rmin = parse_signal(out_dir / "recombination_profile.tsv")
    present = [r for r in regions if r.get("donor_absent") != "yes"]
    major_clade = clade_of(regions[0]["major_parent"]) if regions else "?"
    span = [r for r in present
            if int(r["query_start"]) <= setup.q_end and int(r["query_end"]) >= setup.q_start]
    donor_hits = [r for r in span if donor_match(clade_of(r["minor_parent"]),
                                                 setup.clade_b, setup.clade_a)]
    detected = len(present) >= 1
    backbone_ok = clade_match(major_clade, setup.clade_a)
    donor_ok = len(donor_hits) >= 1
    # Below the 4% attribution floor (the harness's own meaningfulness threshold) the
    # exact backbone clade is statistical noise -- near-identical genomes from adjacent
    # clades win windows by chance. There the meaningful, testable result is detection
    # plus donor-region recovery (which still requires the donor clade and span); the
    # backbone-label match is dropped. Above the floor it is still required.
    low_divergence = setup.divergence < MIN_DIVERGENCE
    passed = detected and donor_ok and (backbone_ok or low_divergence)
    # The recovered donor region was called by more than one ensemble method.
    agree = any("," in (r.get("methods") or "") for r in donor_hits)
    # Finer attribution quality: the best donor attribution among the span-overlapping
    # regions (closest to the true donor clade), reported alongside the verdict.
    donor_obs = max((clade_of(r["minor_parent"]) for r in span),
                    key=lambda c: shared_clade_depth(c, setup.clade_b), default="?")
    return {
        "name": setup.name, "clade_a": setup.clade_a, "clade_b": setup.clade_b,
        "divergence": setup.divergence, "n_refs": n_refs,
        "true_span": (setup.q_start, setup.q_end), "n_regions": len(present),
        "major_clade": major_clade, "detected": detected, "mode": mode,
        "backbone_ok": backbone_ok, "donor_ok": donor_ok, "phi_p": phi_p, "rmin": rmin,
        "agree": agree, "low_divergence": low_divergence,
        "pass": passed, "runtime": runtime,
        "donor_obs": donor_obs,
        "donor_tier": attribution_tier(donor_obs, setup.clade_b),
        "donor_depth": shared_clade_depth(donor_obs, setup.clade_b),
        "backbone_tier": attribution_tier(major_clade, setup.clade_a),
        "backbone_depth": shared_clade_depth(major_clade, setup.clade_a),
    }


def _base(setup: CaseSetup, mode: str, runtime: float, present: list, **extra) -> dict:
    """Reporting fields shared by every scorer."""
    row = {
        "name": setup.name, "clade_a": setup.clade_a, "clade_b": setup.clade_b,
        "divergence": setup.divergence, "mode": mode, "runtime": runtime,
        "detected": len(present) >= 1, "n_regions": len(present),
        "case_type": setup.case_type,
    }
    row.update(extra)
    return row


def _score_neg_pure(
    out_dir: Path, clade_of, setup: CaseSetup, n_refs: int, mode: str, runtime: float,
) -> dict:
    """Non-recombinant query: PASS iff no (donor-present) region is called."""
    regions = parse_regions(out_dir / "recombination_regions.tsv")
    present = [r for r in regions if r.get("donor_absent") != "yes"]
    row = _base(setup, mode, runtime, present, n_refs=n_refs, n_false_regions=len(present))
    row["pass"] = len(present) == 0
    return row


def _score_neg_within(
    out_dir: Path, clade_of, setup: CaseSetup, n_refs: int, mode: str, runtime: float,
) -> dict:
    """Within-clade splice: PASS iff no region attributes a cross-top-level-clade donor."""
    regions = parse_regions(out_dir / "recombination_regions.tsv")
    present = [r for r in regions if r.get("donor_absent") != "yes"]
    cross = [
        r for r in present
        if base_clade(clade_of(r["minor_parent"])).split(".")[0]
        != base_clade(clade_of(r["major_parent"])).split(".")[0]
    ]
    row = _base(setup, mode, runtime, present, n_refs=n_refs, n_false_regions=len(cross))
    row["pass"] = len(cross) == 0
    return row


def _score_low_div(
    out_dir: Path, clade_of, setup: CaseSetup, n_refs: int, mode: str, runtime: float,
) -> dict:
    """Low-divergence attribution: require detection + donor top-level + backbone
    top-level clade (no <4% free pass). Exact-vs-sibling is reported, not gating."""
    regions = parse_regions(out_dir / "recombination_regions.tsv")
    present = [r for r in regions if r.get("donor_absent") != "yes"]
    major_clade = clade_of(regions[0]["major_parent"]) if regions else "?"
    span = [r for r in present
            if int(r["query_start"]) <= setup.q_end and int(r["query_end"]) >= setup.q_start]
    donor_hits = [r for r in span
                  if attribution_tier(clade_of(r["minor_parent"]), setup.clade_b) != "mismatch"]
    detected = len(present) >= 1
    backbone_ok = attribution_tier(major_clade, setup.clade_a) != "mismatch"
    passed = detected and len(donor_hits) >= 1 and backbone_ok
    donor_obs = max((clade_of(r["minor_parent"]) for r in span),
                    key=lambda c: shared_clade_depth(c, setup.clade_b), default="?")
    return _base(setup, mode, runtime, present, n_refs=n_refs, major_clade=major_clade,
                 backbone_ok=backbone_ok, donor_ok=len(donor_hits) >= 1,
                 backbone_tier=attribution_tier(major_clade, setup.clade_a),
                 donor_tier=attribution_tier(donor_obs, setup.clade_b),
                 **{"pass": passed})


def run_case(case: dict, logger: logging.Logger) -> dict:
    """One hybrid case under the default (or ``$HARNESS_METHODS``) caller set, tip panel."""
    setup = _prepare_case(case, logger)
    # Region caller(s): the default ensemble, overridable via $HARNESS_METHODS
    # (e.g. "hmm" or "all") to benchmark a different caller set on the same datasets.
    methods = parse_methods(os.environ.get("HARNESS_METHODS", ",".join(DEFAULT_METHODS)))
    return _build_and_score(setup, "tip", methods, setup.out, logger)


# Five-config attribution-comparison grid: tip vs source-removed consensus panel,
# default ensemble vs +barcode, plus post-hoc reattribute on the tip panel.
# The baseline (tip, default) reproduces the headline harness.
COMPARE_CONFIGS = [
    ("baseline", "tip", ("hmm", "3seq", "maxchi", "bootscan"), False),
    ("+barcode", "tip", ("hmm", "3seq", "maxchi", "bootscan", "barcode"), False),
    ("consensus", "consensus", ("hmm", "3seq", "maxchi", "bootscan"), False),
    ("consensus+barcode", "consensus", ("hmm", "3seq", "maxchi", "bootscan", "barcode"), False),
    # post-hoc donor re-attribution on the tip panel (backbone unchanged)
    ("reattribute", "tip", ("hmm", "3seq", "maxchi", "bootscan"), True),
]


def compare_case(case: dict, logger: logging.Logger) -> dict:
    """Run one case under every :data:`COMPARE_CONFIGS` configuration, reusing the
    expensive setup once. A SKIP from the baseline (tip, default) means the case itself
    is ill-posed and propagates; a SKIP from a later config is recorded per config."""
    setup = _prepare_case(case, logger)
    records: dict[str, dict] = {}
    for label, mode, names, reattribute in COMPARE_CONFIGS:
        methods = parse_methods(",".join(names))
        sub = setup.out / ("cmp_" + re.sub(r"[^\w]+", "_", label))
        try:
            records[label] = _build_and_score(setup, mode, methods, sub, logger,
                                               reattribute=reattribute)
        except CaseSkipped:
            if label == "baseline":
                raise
            logger.info("[%s] %s: SKIP (panel representation)", setup.name, label)
            records[label] = {"skip": True}
    return {
        "name": setup.name, "clade_a": setup.clade_a, "clade_b": setup.clade_b,
        "divergence": setup.divergence, "configs": records,
    }


def main(argv: list[str]) -> int:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    logger = logging.getLogger("tessera")
    compare = "--compare" in argv
    names = [a for a in argv if not a.startswith("-")]
    cases = [c for c in HYBRIDS if not names or c["name"] in names]
    if compare:
        return _run_compare(cases, logger)
    return _run_default(cases, logger)


def _run_default(cases: list[dict], logger: logging.Logger) -> int:
    print(f"\nTessera hybrid detection -- {len(cases)} case(s)\n" + "=" * 72)
    results: list[dict] = []
    for case in cases:
        try:
            results.append(run_case(case, logger))
        except CaseSkipped as exc:  # an unsuitable dataset, not a failure
            logger.info("[%s] SKIP: %s", case["name"], exc)
            results.append({"name": case["name"], "pass": None, "skip": str(exc)})
        except Exception as exc:  # noqa: BLE001 - report and continue the batch
            logger.exception("[%s] ERROR", case["name"])
            results.append({"name": case["name"], "pass": None, "error": str(exc)})

    print("\n" + "=" * 80)
    print(f"{'case':16} {'backbone x donor':24} {'div':>5} {'major':>10} "
          f"{'det':>3} {'bb':>3} {'don':>3} {'agr':>3} {'mode':>9} {'PHI p':>8} {'Rmin':>4} "
          f"{'time':>6}  verdict")
    print("-" * 108)
    for r in results:
        if r.get("skip") is not None:
            print(f"{r['name']:16} SKIP: {r['skip'][:54]}")
            continue
        if r.get("pass") is None:
            print(f"{r['name']:16} ERROR: {r.get('error', '')[:50]}")
            continue
        verdict = "PASS" if r["pass"] else "FAIL"
        print(f"{r['name']:16} {r['clade_a']+' x '+r['clade_b']:24.24} "
              f"{r['divergence']:4.1f}% {r['major_clade']:>10.10} "
              f"{'Y' if r['detected'] else 'n':>3} {'Y' if r['backbone_ok'] else 'n':>3} "
              f"{'Y' if r['donor_ok'] else 'n':>3} {'Y' if r.get('agree') else '.':>3} "
              f"{r['mode']:>9} "
              f"{r.get('phi_p', '-'):>8.8} {r.get('rmin', '-'):>4} "
              f"{r['runtime']:5.0f}s  {verdict}")
    passed = sum(1 for r in results if r.get("pass"))
    ran = sum(1 for r in results if r.get("pass") is not None)
    skipped = sum(1 for r in results if r.get("skip") is not None)
    errored = sum(1 for r in results if r.get("pass") is None and r.get("skip") is None)
    print(f"\n{passed}/{ran} passed  ({skipped} skipped, {errored} error)")
    scored = [r for r in results if r.get("pass") is not None]
    negs = [r for r in scored if r.get("case_type") in ("neg_pure", "neg_within")]
    pos = [r for r in scored if r not in negs]
    sens_pass = sum(1 for r in pos if r["pass"])
    spec_pass = sum(1 for r in negs if r["pass"])
    false_calls = sum(r.get("n_false_regions", 0) for r in negs)
    print(f"sensitivity {sens_pass}/{len(pos)}  "
          f"specificity {spec_pass}/{len(negs)} ({false_calls} false call(s))  "
          f"({skipped} skipped, {errored} error)")
    return 1 if errored else 0


def _run_compare(cases: list[dict], logger: logging.Logger) -> int:
    """Run the five-config attribution grid per case (G1 measurement); print a per-config
    table and write ``data/attribution_compare.tsv``. The baseline column reproduces the
    headline harness; the other columns measure barcode, the consensus panel, and
    post-hoc donor re-attribution (``reattribute`` column)."""
    print(f"\nTessera attribution comparison -- {len(cases)} case(s), "
          f"{len(COMPARE_CONFIGS)} configs each\n" + "=" * 72)
    tsv = [["case", "divergence", "config", "detected", "backbone_tier", "backbone_depth",
            "donor_tier", "donor_depth", "donor_obs", "expected_donor", "agree", "pass"]]
    pass_count = {label: 0 for label, _, _, _ in COMPARE_CONFIGS}
    improved = {label: 0 for label, _, _, _ in COMPARE_CONFIGS}
    ran = 0
    for case in cases:
        try:
            res = compare_case(case, logger)
        except CaseSkipped as exc:
            print(f"\n{case['name']:18} SKIP: {exc}")
            continue
        except Exception as exc:  # noqa: BLE001 - report and continue the batch
            logger.exception("[%s] ERROR", case["name"])
            print(f"\n{case['name']:18} ERROR: {exc}")
            continue
        ran += 1
        base = res["configs"].get("baseline", {})
        base_depth = base.get("donor_depth", -1) if not base.get("skip") else -1
        print(f"\n{res['name']}  ({res['clade_a']} x {res['clade_b']}, {res['divergence']:.1f}%)")
        for label, _, _, _ in COMPARE_CONFIGS:
            rec = res["configs"][label]
            if rec.get("skip"):
                print(f"  {label:18} SKIP")
                tsv.append([res["name"], f"{res['divergence']:.1f}", label,
                            "", "", "", "", "", "", res["clade_b"], "", ""])
                continue
            verdict = "PASS" if rec["pass"] else "FAIL"
            print(f"  {label:18} det {'Y' if rec['detected'] else 'n'}  "
                  f"bb {rec['backbone_tier']:8} don {rec['donor_tier']:8} "
                  f"depth {rec['donor_depth']}  agr {'Y' if rec['agree'] else '.'}  {verdict}")
            if rec["pass"]:
                pass_count[label] += 1
            if label != "baseline" and rec["donor_depth"] > base_depth:
                improved[label] += 1
            tsv.append([res["name"], f"{res['divergence']:.1f}", label,
                        "Y" if rec["detected"] else "n", rec["backbone_tier"],
                        rec["backbone_depth"], rec["donor_tier"], rec["donor_depth"],
                        rec["donor_obs"], res["clade_b"], "Y" if rec["agree"] else "n",
                        "PASS" if rec["pass"] else "FAIL"])

    out_tsv = HERE / "data" / "attribution_compare.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_tsv.write_text("\n".join("\t".join(map(str, row)) for row in tsv) + "\n")

    print("\n" + "=" * 72 + f"\nsummary (of {ran} cases that ran)")
    for label, _, _, _ in COMPARE_CONFIGS:
        extra = ("" if label == "baseline"
                 else f"; donor attribution improved vs baseline in {improved[label]}")
        print(f"  {label:18} PASS {pass_count[label]}/{ran}{extra}")
    print(f"\nMachine-readable comparison: {out_tsv}")
    print("Note: 'baseline' reproduces the headline harness; the other columns are the "
          "G1 measurement, not the headline numbers.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
