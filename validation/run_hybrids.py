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
import re
import shutil
import sys
import time
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
from tessera.recomb.run import RecombParams, run_recomb
from tessera.recomb.typing import LINEAGES_TSV, lineage_map_from_rows, write_lineage_map

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

# Recombinant clade names (HIV CRF/URF, Pango X-lineages): excluded as hybrid
# parents -- a clean recombination test needs non-recombinant parental lineages.
_RECOMBINANT_CLADE = re.compile(r"^(CRF|URF)\d|^X[A-Z]|recombinant", re.IGNORECASE)


def is_recombinant_clade(name: str) -> bool:
    return bool(_RECOMBINANT_CLADE.search(name))


def clade_match(observed: str, expected: str) -> bool:
    """True when two Nextclade labels name the same lineage, allowing for the
    hierarchical naming (e.g. ``A`` matches ``A.1``; a backbone recovered at a
    coarser granularity still counts)."""
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


class CaseSkipped(Exception):
    """A dataset is unsuitable for a hybrid test (reported as SKIP, not ERROR)."""


def clade_members(
    tips: dict[str, tuple[str, list[str]]], *, exclude_recombinant: bool = True
) -> dict[str, list[str]]:
    """Map each clade to its member accessions (skip unlabelled, and -- by default
    -- recombinant clades)."""
    by_clade: dict[str, list[str]] = {}
    for acc, (clade, _) in tips.items():
        if not clade or clade == "NA":
            continue
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
    logger: logging.Logger,
) -> tuple[str, str, str, str]:
    """Return ``(clade_a, clade_b, source_a, source_b)``.

    Among clades with at least ``MIN_MEMBERS`` genomes, choose the most-divergent
    pair (each represented by its central genome); the larger clade is the backbone
    (A). A pinned pair overrides the search. Sources are the clades' central genomes.

    Recombinant clades are excluded as parents, but if that leaves fewer than two
    eligible clades (e.g. the SARS-CoV-2 XBB dataset, which is entirely XBB
    sub-lineages) they are allowed back in. Too few clades either way is a SKIP.
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
                if best is None or div > best[0]:
                    best = (div, ca, cb)
        _, ca, cb = best
    if len(eligible[cb]) > len(eligible[ca]):  # backbone = the better-represented clade
        ca, cb = cb, ca
    return ca, cb, reps[ca], reps[cb]


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


def run_case(case: dict, logger: logging.Logger) -> dict:
    name, path = case["name"], case["dataset"]
    out = DATA / name
    out.mkdir(parents=True, exist_ok=True)

    dataset = resolve_dataset(Path("/dev/null"), path, email=None, logger=logger)
    genomes = build_pool(dataset, cache_dir=nextclade_cache(dataset.path, dataset.tag),
                         logger=logger)
    reference = read_reference(_download_text(dataset, "reference", logger))
    tree = json.loads(_download_text(dataset, "treeJson", logger))["tree"]
    tips = collect_tips(tree, case.get("clade_key"))

    clade_a, clade_b, src_a, src_b = pick_parents(tips, reference, case.get("clades", []), logger)
    divergence = 100.0 - pct_identity(
        reconstruct_gapped(reference, tips[src_a][1]),
        reconstruct_gapped(reference, tips[src_b][1]),
    )
    logger.info("[%s] backbone clade %s (%s) x donor clade %s (%s); inter-clade divergence %.1f%%",
                name, clade_a, src_a, clade_b, src_b, divergence)
    if divergence < case.get("min_divergence", MIN_DIVERGENCE):
        raise CaseSkipped(f"clades too similar ({divergence:.1f}% divergence)")

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
    tip_bases = {acc.split(".")[0] for acc in tips}
    pool = [
        g for g in genomes
        if strip_sequence_extension(g.name).split(".")[0] in tip_bases
        and strip_sequence_extension(g.name).split(".")[0] not in drop
    ]

    t0 = time.monotonic()
    try:
        # Keep "siblings": for a hybrid of close parents the backbone parent is
        # >95% genome-wide ANI to the query (its 70% backbone dominates the
        # average) and would be dropped as a masking twin -- but the synthetic
        # pool contains no actual recombinant twin, so this is the documented
        # close-parent / --seed-keep-siblings setting, and dropping it would
        # remove the backbone parent itself.
        selected = select_regional(query, pool, window=sel_window, per_window=2,
                                   drop_siblings=False, logger=logger).selected
    except ToolExecutionError:  # skani rejects very short gene/segment datasets
        selected = representative_panel(pool, tips, logger)
    # On a near-identical panel (mpox/VZV ~0.5%) dereplication collapses the parent
    # clades into one representative, leaving nothing to compete; rebuild from one
    # central genome per clade so both parents are present for the detector to test.
    sel_clades = {clade_of_label(g.name, tips) for g in selected}
    if len([c for c in sel_clades if c not in ("?", "NA")]) < 2:
        logger.info("[%s] regional selection collapsed to %d clade(s); using clade "
                    "representatives instead.", name, len(sel_clades))
        selected = representative_panel(pool, tips, logger)
    if not selected:
        raise CaseSkipped("no reference panel could be built")
    collection = out / "collection"
    if collection.exists():
        shutil.rmtree(collection)
    collection.mkdir(parents=True)
    for g in selected:
        shutil.copy(g, collection / g.name)

    rows = [(strip_sequence_extension(g.name), clade_of_label(g.name, tips), "nextclade")
            for g in selected]
    rows.append((query_label, f"hybrid({clade_a}+{clade_b})", "query"))
    write_lineage_map(out / LINEAGES_TSV, rows)
    lineage_map = lineage_map_from_rows(rows)

    msa = out / "panel.msa.fasta"
    build_msa(MsaParams(query=query, collection=collection, output=msa,
                        aligner=aligner, threads=THREADS), logger)
    windowing = run_recomb(RecombParams(msa=msa, output=out, query=query_label,
                                        window_size=window, window_step=step,
                                        lineage_map=lineage_map), logger)
    runtime = time.monotonic() - t0
    mode = "info-site" if windowing.startswith("informative") else "bp"

    regions = parse_regions(out / "recombination_regions.tsv")
    phi_p, rmin = parse_signal(out / "recombination_profile.tsv")
    present = [r for r in regions if r.get("donor_absent") != "yes"]
    major_clade = clade_of_label(regions[0]["major_parent"], tips) if regions else "?"
    donor_hits = [
        r for r in present
        if donor_match(clade_of_label(r["minor_parent"], tips), clade_b, clade_a)
        and int(r["query_start"]) <= q_end and int(r["query_end"]) >= q_start
    ]
    detected = len(present) >= 1
    backbone_ok = clade_match(major_clade, clade_a)
    donor_ok = len(donor_hits) >= 1
    # The recovered donor region was called by more than one ensemble method.
    agree = any("," in (r.get("methods") or "") for r in donor_hits)
    return {
        "name": name, "clade_a": clade_a, "clade_b": clade_b,
        "divergence": divergence, "n_refs": len(selected),
        "true_span": (q_start, q_end), "n_regions": len(present),
        "major_clade": major_clade, "detected": detected, "mode": mode,
        "backbone_ok": backbone_ok, "donor_ok": donor_ok, "phi_p": phi_p, "rmin": rmin,
        "agree": agree,
        "pass": detected and backbone_ok and donor_ok, "runtime": runtime,
    }


def main(argv: list[str]) -> int:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    logger = logging.getLogger("tessera")
    cases = [c for c in HYBRIDS if not argv or c["name"] in argv]
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
    return 1 if errored else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
