#!/usr/bin/env python
"""Generate synthetic recombinant (hybrid) genomes from Nextclade datasets and
test RecomFi's detection performance on them.

For each configured dataset the harness:

1. Builds the Nextclade reference pool (the shipped ``build_pool``; cached).
2. Picks two clades A and B (the most-divergent well-represented pair by default,
   or an explicit pair) and each clade's central genome, in reference coordinates.
3. Splices an A-backbone / B-insert hybrid (B donates the middle of the genome),
   recording the true donor span in query coordinates.
4. Runs RecomFi pool-only against the dataset's pool with the two exact source
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

from recomfi.core.cache import nextclade_cache
from recomfi.core.io import strip_sequence_extension, write_fasta_record
from recomfi.discover.nextclade import (
    _MUT,
    _accession_of,
    _clade_of,
    _download_text,
    build_pool,
    resolve_dataset,
)
from recomfi.discover.pool import select_regional
from recomfi.msa.build import MsaParams, build_msa
from recomfi.recomb.run import RecombParams, run_recomb
from recomfi.recomb.typing import LINEAGES_TSV, lineage_map_from_rows, write_lineage_map

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "hybrids"

# One entry per hybrid test case. ``clades`` may pin the two parent clades; left
# empty the most-divergent well-represented pair is used. The donor (clade B)
# occupies the middle ``insert`` fraction of the genome; the backbone (clade A) the
# flanks.
HYBRIDS: list[dict] = [
    {"name": "hiv1", "dataset": "community/neherlab/hiv-1/hxb2"},
    {"name": "dengue", "dataset": "nextstrain/dengue/all"},
    {"name": "measles", "dataset": "nextstrain/measles/genome/WHO-2012"},
    {"name": "rsv_a", "dataset": "nextstrain/rsv/a/EPI_ISL_412866"},
    {"name": "rubella", "dataset": "nextstrain/rubella/genome"},
]
INSERT = (0.35, 0.65)  # donor (clade B) occupies this fraction of the genome
WINDOW, STEP = 500, 50
SELECT_WINDOW = 1500
THREADS = 4


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


def collect_tips(tree: dict) -> dict[str, tuple[str, list[str]]]:
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
            tips[acc] = (_clade_of(node.get("node_attrs", {}) or {}), path_muts)
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


def clade_members(tips: dict[str, tuple[str, list[str]]]) -> dict[str, list[str]]:
    """Map each clade to its member accessions (skip unlabelled/recombinant tips)."""
    by_clade: dict[str, list[str]] = {}
    for acc, (clade, _) in tips.items():
        if clade and clade != "NA" and not is_recombinant_clade(clade):
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
    tips: dict[str, tuple[str, list[str]]], reference: str, pinned: list[str]
) -> tuple[str, str, str, str]:
    """Return ``(clade_a, clade_b, source_a, source_b)``.

    Among clades with at least ``MIN_MEMBERS`` genomes, choose the most-divergent
    pair (each represented by its central genome); the larger clade is the backbone
    (A). A pinned pair overrides the search. Sources are the clades' central genomes.
    """
    by_clade = clade_members(tips)
    eligible = {c: m for c, m in by_clade.items() if len(m) >= MIN_MEMBERS}
    if len(eligible) < 2:
        raise ValueError(f"dataset has fewer than two clades with >={MIN_MEMBERS} genomes")
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


def clade_of_label(label: str, tips: dict[str, tuple[str, list[str]]]) -> str:
    base = strip_sequence_extension(label).split(".")[0]
    for acc, (clade, _) in tips.items():
        if acc.split(".")[0] == base:
            return clade
    return "?"


def run_case(case: dict, logger: logging.Logger) -> dict:
    name, path = case["name"], case["dataset"]
    out = DATA / name
    out.mkdir(parents=True, exist_ok=True)

    dataset = resolve_dataset(Path("/dev/null"), path, email=None, logger=logger)
    genomes = build_pool(dataset, cache_dir=nextclade_cache(dataset.path, dataset.tag),
                         logger=logger)
    reference = read_reference(_download_text(dataset, "reference", logger))
    tree = json.loads(_download_text(dataset, "treeJson", logger))["tree"]
    tips = collect_tips(tree)

    clade_a, clade_b, src_a, src_b = pick_parents(tips, reference, case.get("clades", []))
    divergence = 100.0 - pct_identity(
        reconstruct_gapped(reference, tips[src_a][1]),
        reconstruct_gapped(reference, tips[src_b][1]),
    )
    logger.info("[%s] backbone clade %s (%s) x donor clade %s (%s); inter-clade divergence %.1f%%",
                name, clade_a, src_a, clade_b, src_b, divergence)

    query_seq, q_start, q_end = make_hybrid(reference, tips[src_a][1], tips[src_b][1])
    query = out / "hybrid.fasta"
    with open(query, "w") as fo:
        write_fasta_record(fo, f"hybrid_{clade_a}_{clade_b}", query_seq)
    query_label = strip_sequence_extension(query.name)

    # Pool minus the two exact source genomes (their clades remain represented).
    drop = {src_a.split(".")[0], src_b.split(".")[0]}
    pool = [g for g in genomes if strip_sequence_extension(g.name).split(".")[0] not in drop]

    t0 = time.monotonic()
    selected = select_regional(query, pool, window=SELECT_WINDOW, per_window=2,
                               drop_siblings=True, logger=logger).selected
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
                        aligner="mafft", threads=THREADS), logger)
    run_recomb(RecombParams(msa=msa, output=out, query=query_label,
                            window_size=WINDOW, window_step=STEP,
                            lineage_map=lineage_map), logger)
    runtime = time.monotonic() - t0

    regions = parse_regions(out / "recombination_regions.tsv")
    present = [r for r in regions if r.get("donor_absent") != "yes"]
    major_clade = clade_of_label(regions[0]["major_parent"], tips) if regions else "?"
    donor_hits = [
        r for r in present
        if clade_match(clade_of_label(r["minor_parent"], tips), clade_b)
        and int(r["query_start"]) <= q_end and int(r["query_end"]) >= q_start
    ]
    detected = len(present) >= 1
    backbone_ok = clade_match(major_clade, clade_a)
    donor_ok = len(donor_hits) >= 1
    return {
        "name": name, "clade_a": clade_a, "clade_b": clade_b,
        "divergence": divergence, "n_refs": len(selected),
        "true_span": (q_start, q_end), "n_regions": len(present),
        "major_clade": major_clade, "detected": detected,
        "backbone_ok": backbone_ok, "donor_ok": donor_ok,
        "pass": detected and backbone_ok and donor_ok, "runtime": runtime,
    }


def main(argv: list[str]) -> int:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    logger = logging.getLogger("recomfi")
    cases = [c for c in HYBRIDS if not argv or c["name"] in argv]
    print(f"\nRecomFi hybrid detection -- {len(cases)} case(s)\n" + "=" * 72)
    results: list[dict] = []
    for case in cases:
        try:
            results.append(run_case(case, logger))
        except Exception as exc:  # noqa: BLE001 - report and continue the batch
            logger.exception("[%s] ERROR", case["name"])
            results.append({"name": case["name"], "pass": None, "error": str(exc)})

    print("\n" + "=" * 72)
    print(f"{'case':9} {'backbone x donor':24} {'div':>5} {'major':>10} "
          f"{'det':>3} {'bb':>3} {'don':>3} {'time':>6}  verdict")
    print("-" * 72)
    for r in results:
        if r.get("pass") is None:
            print(f"{r['name']:9} ERROR: {r.get('error', '')[:50]}")
            continue
        verdict = "PASS" if r["pass"] else "FAIL"
        print(f"{r['name']:9} {r['clade_a']+' x '+r['clade_b']:24.24} "
              f"{r['divergence']:4.1f}% {r['major_clade']:>10.10} "
              f"{'Y' if r['detected'] else 'n':>3} {'Y' if r['backbone_ok'] else 'n':>3} "
              f"{'Y' if r['donor_ok'] else 'n':>3} {r['runtime']:5.0f}s  {verdict}")
    passed = sum(1 for r in results if r.get("pass"))
    print(f"\n{passed}/{sum(1 for r in results if r.get('pass') is not None)} passed")
    failed = [r for r in results if r.get("pass") is False or r.get("pass") is None]
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
