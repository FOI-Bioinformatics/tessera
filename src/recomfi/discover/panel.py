"""Curate a reference panel for recombination detection (skani / skDER).

Auto-filling a collection maximises *coverage* by recruiting the query's closest
genomes -- which, for a recombinant query, are its own-lineage relatives. Those
match the query everywhere and **mask** the recombination (the HIV CRF01_AE case).
Detection instead needs a *diverse* panel of distinct parental lineages and must
exclude the query's siblings.

This module curates a downloaded set with two skani-based steps:

- **Sibling exclusion** -- a candidate is a sibling of the query when its
  genome-wide ANI to the query exceeds the current backbone's by a margin AND it
  aligns over most of the query (high query alignment fraction). The fraction is
  the key signal: an env-only parent has high ANI but low query-AF (regional,
  kept); a same-lineage relative has high ANI and high query-AF (whole genome,
  dropped). Being relative to the backbone, it needs no per-organism ANI cutoff.
- **Dereplication** -- skDER collapses near-duplicate references to representatives.

Both tools are optional (bioconda ``skani`` / ``skder``); absence is reported, not
fatal. ANI/AF are percentages (0-100) throughout.
"""

from __future__ import annotations

import html
import logging
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.errors import UserInputError
from ..core.io import strip_sequence_extension
from ..core.plugins import ToolCapabilities
from ..core.process import run_tool
from ..recomb.typing import LineageMap, lineage_of

SKANI = ToolCapabilities(
    name="skani", conda=("bioconda::skani",),
    required_binaries=(BinarySpec("skani", version_args=("-V",)),),
)
SKDER = ToolCapabilities(
    name="skder", conda=("bioconda::skder", "bioconda::skani"),
    required_binaries=(BinarySpec("skder", version_args=("-v",)),),
)

# Defaults, both relative/justified rather than fit to a dataset.
DEFAULT_SIBLING_MARGIN = 3.0  # query-ANI must beat the backbone's by this many %
DEFAULT_AF_MIN = 80.0  # ... over at least this % of the query (whole-genome match)
DEFAULT_DEREP_ANI = 99.0  # skDER: genomes >= this ANI collapse to one representative
# Above this total inline-path length, pass skDER a staging directory instead of
# individual paths (a large pool overflows the OS command-line length limit).
_ARG_LIST_LIMIT = 100_000
# A whole-genome near-twin: when the backbone is a partial-coverage parent (it donated
# only one region, so its ANI is high but its query coverage is not), a genome that
# covers nearly the whole query at comparable ANI is the query's own lineage, even
# though its ANI does not beat the backbone's. Coverage, not ANI, separates them here.
DEFAULT_AF_FULL = 90.0  # ... covers at least this % of the query (the whole genome)
DEFAULT_SIBLING_ANI_TOL = 1.5  # ... at an ANI within this much of the backbone's


def skani_available() -> bool:
    return shutil.which("skani") is not None


def skder_available() -> bool:
    return shutil.which("skder") is not None


@dataclass
class CurationResult:
    kept: list[Path]  # the curated panel (sibling-free, dereplicated)
    siblings: list[Path]  # dropped: same lineage as the query
    redundant: list[Path]  # dropped: near-duplicate of a kept representative
    table: list[dict] = field(default_factory=list)  # per-genome role + ani/af


def skani_query_ani(
    query_fasta: Path, refs: list[Path], logger: logging.Logger
) -> dict[Path, tuple[float, float]]:
    """Genome-wide ``{ref -> (ANI, query_alignment_fraction)}`` via ``skani dist``.

    A reference that does not align to the query at all is absent from skani's
    output; it is reported here as ``(0.0, 0.0)``.
    """
    if not refs:
        return {}
    with tempfile.TemporaryDirectory() as tmp:
        ref_list = Path(tmp) / "refs.txt"
        ref_list.write_text("\n".join(str(r.resolve()) for r in refs) + "\n")
        out = Path(tmp) / "dist.tsv"
        run_tool(
            SKANI,
            ["skani", "dist", "-q", str(query_fasta.resolve()),
             "--rl", str(ref_list), "-o", str(out), "--medium"],
            logger=logger, log_prefix="skani",
        )
        rows = out.read_text().splitlines()
    by_path = {r.resolve(): r for r in refs}
    result: dict[Path, tuple[float, float]] = {r: (0.0, 0.0) for r in refs}
    for line in rows[1:]:  # skip header
        f = line.split("\t")
        if len(f) < 5:
            continue
        ref = Path(f[0]).resolve()
        if ref in by_path:
            result[by_path[ref]] = (float(f[2]), float(f[4]))  # ANI, AF_query
    return result


def pick_backbone(
    query_fasta: Path,
    genomes: list[Path],
    *,
    af_min: float = DEFAULT_AF_MIN,
    logger: logging.Logger,
) -> Path | None:
    """Choose the genome the query matches across its whole length (the backbone).

    The backbone anchors the relative sibling test, so it must be a *full-genome*
    parent of the query's majority lineage, not a regional (env-only) reference.
    Among genomes covering at least ``af_min`` % of the query, the one with the
    highest ANI is the closest whole-genome relative -- the backbone. Pick this
    from a *sibling-free* set (curate each round before siblings accumulate);
    otherwise a sibling, being closer still, would be chosen. Returns ``None`` for
    an empty set; falls back to the highest query-AF genome if none clears
    ``af_min``.
    """
    if not genomes:
        return None
    ani = skani_query_ani(query_fasta, genomes, logger)
    covering = [(g, a, af) for g, (a, af) in ani.items() if af >= af_min]
    if not covering:
        return max(ani.items(), key=lambda kv: kv[1][1])[0]  # highest AF overall
    return max(covering, key=lambda t: t[1])[0]  # highest ANI among full-coverage


def filter_siblings(
    query_fasta: Path,
    backbone: Path,
    candidates: list[Path],
    *,
    ani_margin: float = DEFAULT_SIBLING_MARGIN,
    af_min: float = DEFAULT_AF_MIN,
    af_full: float = DEFAULT_AF_FULL,
    ani_tol: float = DEFAULT_SIBLING_ANI_TOL,
    logger: logging.Logger,
) -> tuple[list[Path], list[Path], dict[Path, tuple[float, float]]]:
    """Split candidates into ``(parents, siblings, ani_map)``.

    A candidate is a *sibling* of the query (its own lineage, which would mask
    recombination) when either:

    - its query-ANI exceeds the backbone's by ``ani_margin`` over at least ``af_min``
      % of the query -- it is closer than the backbone genome-wide; or
    - it covers nearly the whole query (``af_full``) at an ANI within ``ani_tol`` of
      the backbone's -- a whole-genome near-twin. This catches the masking sibling
      when the backbone is a *partial-coverage* parent (e.g. one that donated a single
      ORF), where ANI alone cannot tell the sibling from the parent; coverage does.
    """
    ani = skani_query_ani(query_fasta, [backbone, *candidates], logger)
    backbone_ani = ani.get(backbone, (0.0, 0.0))[0]
    parents, siblings = [], []
    for cand in candidates:
        cand_ani, cand_af = ani[cand]
        closer = cand_ani - backbone_ani >= ani_margin and cand_af >= af_min
        whole_genome_twin = cand_af >= af_full and cand_ani >= backbone_ani - ani_tol
        if closer or whole_genome_twin:
            siblings.append(cand)
        else:
            parents.append(cand)
    return parents, siblings, ani


def dereplicate(
    genomes: list[Path], *, ani: float = DEFAULT_DEREP_ANI, logger: logging.Logger
) -> tuple[list[Path], list[Path]]:
    """Return ``(representatives, redundant)`` via skDER greedy dereplication."""
    if len(genomes) < 2:
        return list(genomes), []
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "skder"
        genome_args = [str(g.resolve()) for g in genomes]
        # skDER's -g accepts file paths or a directory. A large pool (thousands of
        # genomes) passed inline overflows the OS command-line length limit, so stage
        # the genomes into one directory and pass that directory instead.
        if sum(len(a) for a in genome_args) > _ARG_LIST_LIMIT:
            staged = Path(tmp) / "genomes"
            staged.mkdir()
            for i, g in enumerate(genomes):
                (staged / f"g{i}_{g.name}").symlink_to(g.resolve())  # unique, keeps suffix
            genome_args = [str(staged)]
        run_tool(
            SKDER,
            ["skder", "-g", *genome_args, "-o", str(out), "-i", str(ani), "-d", "greedy"],
            logger=logger, log_prefix="skder",
        )
        # Resolve representative paths while the staging symlinks still exist, so a
        # staged symlink dereferences to its real genome rather than a dangling path.
        rep_paths = {
            Path(line.strip()).resolve()
            for line in (out / "skDER_Results.txt").read_text().splitlines()
            if line.strip()
        }
    representatives = [g for g in genomes if g.resolve() in rep_paths]
    redundant = [g for g in genomes if g.resolve() not in rep_paths]
    return representatives, redundant


def curate_panel(
    query_fasta: Path,
    genomes: list[Path],
    backbone: Path,
    *,
    ani_margin: float = DEFAULT_SIBLING_MARGIN,
    af_min: float = DEFAULT_AF_MIN,
    derep_ani: float = DEFAULT_DEREP_ANI,
    logger: logging.Logger,
) -> CurationResult:
    """Drop the query's siblings, then dereplicate the rest into a diverse panel."""
    if not skani_available():
        raise UserInputError(
            "Panel curation needs skani. Install with: conda install -c bioconda skani"
        )
    candidates = [g for g in genomes if g.resolve() != backbone.resolve()]
    parents, siblings, ani = filter_siblings(
        query_fasta, backbone, candidates,
        ani_margin=ani_margin, af_min=af_min, logger=logger,
    )
    if skder_available():
        reps, redundant = dereplicate([backbone, *parents], ani=derep_ani, logger=logger)
    else:
        logger.info("skder not found; skipping dereplication (sibling filter only).")
        reps, redundant = [backbone, *parents], []

    rep_set = {p.resolve() for p in reps}
    sib_set = {p.resolve() for p in siblings}
    table = []
    for g in genomes:
        a, af = ani.get(g, (None, None))
        if g.resolve() == backbone.resolve():
            role = "backbone"
        elif g.resolve() in sib_set:
            role = "sibling-dropped"
        elif g.resolve() in rep_set:
            role = "representative"
        else:
            role = "redundant-dropped"
        table.append({
            "genome": strip_sequence_extension(g.name), "role": role,
            "query_ani": a, "query_af": af,
        })
    logger.info(
        "Panel curation: %d kept, %d sibling(s) dropped, %d redundant dropped.",
        len(reps), len(siblings), len(redundant),
    )
    return CurationResult(kept=reps, siblings=siblings, redundant=redundant, table=table)


def curate_collection_dir(
    query_fasta: Path,
    collection: Path,
    backbone: Path,
    *,
    ani_margin: float = DEFAULT_SIBLING_MARGIN,
    af_min: float = DEFAULT_AF_MIN,
    derep_ani: float = DEFAULT_DEREP_ANI,
    logger: logging.Logger,
) -> CurationResult:
    """Curate the genomes in ``collection`` in place: drop siblings/redundant files.

    Runs :func:`curate_panel`, then deletes the dropped genome files from disk so a
    subsequent MSA rebuild sees only the diverse, sibling-free panel.
    """
    genomes = sorted(p for p in collection.iterdir() if p.is_file())
    result = curate_panel(
        query_fasta, genomes, backbone,
        ani_margin=ani_margin, af_min=af_min, derep_ani=derep_ani, logger=logger,
    )
    keep = {p.resolve() for p in result.kept} | {backbone.resolve()}
    for g in genomes:
        if g.resolve() not in keep:
            g.unlink()
    return result


def _resolve_backbone(genomes: list[Path], reference: str) -> Path:
    """Find the genome named ``reference`` (by label or filename)."""
    for g in genomes:
        if reference in (g.name, strip_sequence_extension(g.name)):
            return g
    raise UserInputError(
        f"Backbone reference {reference!r} is not in the collection. Available: "
        f"{', '.join(sorted(strip_sequence_extension(g.name) for g in genomes))}"
    )


def curate_collection(
    query_fasta: Path,
    collection: Path,
    output: Path,
    *,
    reference: str | None = None,
    ani_margin: float = DEFAULT_SIBLING_MARGIN,
    af_min: float = DEFAULT_AF_MIN,
    derep_ani: float = DEFAULT_DEREP_ANI,
    logger: logging.Logger,
) -> CurationResult:
    """Curate a collection on disk: copy it to ``output``, drop siblings, dereplicate.

    The backbone is the genome given by ``reference``, or -- if omitted -- the
    query's closest whole-genome relative (:func:`pick_backbone`). Writes the
    curated ``collection/`` and ``panel_lineages.tsv`` under ``output``.
    """
    if not skani_available():
        raise UserInputError(
            "curate-panel needs skani. Install with: conda install -c bioconda skani skder"
        )
    output.mkdir(parents=True, exist_ok=True)
    dest = output / "collection"
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(collection, dest)
    genomes = sorted(p for p in dest.iterdir() if p.is_file())
    if not genomes:
        raise UserInputError(f"Collection {collection} has no genome files.")
    if reference:
        backbone = _resolve_backbone(genomes, reference)
    else:
        backbone = pick_backbone(query_fasta, genomes, af_min=af_min, logger=logger)
        if backbone is None:
            raise UserInputError("Could not determine a backbone; pass --reference.")
    logger.info(
        "Backbone (query's whole-genome anchor): %s", strip_sequence_extension(backbone.name)
    )
    result = curate_collection_dir(
        query_fasta, dest, backbone,
        ani_margin=ani_margin, af_min=af_min, derep_ani=derep_ani, logger=logger,
    )
    write_panel_tsv(output / "panel_lineages.tsv", result.table, logger)
    logger.info(
        "Curated panel: %d kept -> %s. Rebuild with 'recomfi msa' then 'recomfi recomb'.",
        len(result.kept), dest,
    )
    return result


def write_panel_tsv(
    path: Path, table: list[dict], logger: logging.Logger, lineage_map: LineageMap | None = None
) -> None:
    """Write the per-genome curation table (genome, role, query-ANI, query-AF).

    With a non-empty ``lineage_map`` a ``lineage`` column carrying each genome's typed
    genotype is inserted; without one the output is the original four columns.
    """
    logger.info("Writing reference panel: %s", path)
    with open(path, "w") as fo:
        if lineage_map:
            fo.write("genome\trole\tlineage\tquery_ani\tquery_af\n")
        else:
            fo.write("genome\trole\tquery_ani\tquery_af\n")
        for r in table:
            ani = "" if r["query_ani"] is None else f"{r['query_ani']:.2f}"
            af = "" if r["query_af"] is None else f"{r['query_af']:.2f}"
            if lineage_map:
                lineage = lineage_of(r["genome"], lineage_map) or ""
                fo.write(f"{r['genome']}\t{r['role']}\t{lineage}\t{ani}\t{af}\n")
            else:
                fo.write(f"{r['genome']}\t{r['role']}\t{ani}\t{af}\n")


_ROLE_LABEL = {
    "backbone": "backbone",
    "representative": "kept (parent)",
    "sibling-dropped": "dropped (sibling)",
    "redundant-dropped": "dropped (redundant)",
}
_ROLE_ORDER = {"backbone": 0, "representative": 1, "sibling-dropped": 2, "redundant-dropped": 3}


def panel_table_html(table: list[dict], lineage_map: LineageMap | None = None) -> str:
    """An HTML block describing the curated panel, for the report.

    A ``Lineage`` column is added when ``lineage_map`` is non-empty.
    """
    rows = ""
    for r in sorted(table, key=lambda x: (_ROLE_ORDER.get(x["role"], 9), -(x["query_ani"] or 0))):
        ani = "&mdash;" if r["query_ani"] is None else f"{r['query_ani']:.1f}%"
        af = "&mdash;" if r["query_af"] is None else f"{r['query_af']:.0f}%"
        lineage_cell = ""
        if lineage_map:
            lineage = lineage_of(r["genome"], lineage_map)
            lineage_cell = f'<td class="lbl">{html.escape(lineage) if lineage else "&mdash;"}</td>'
        rows += (
            "<tr>"
            f'<td class="lbl">{html.escape(r["genome"])}</td>'
            f'<td class="lbl">{_ROLE_LABEL.get(r["role"], r["role"])}</td>'
            f'{lineage_cell}'
            f'<td class="num">{ani}</td>'
            f'<td class="num">{af}</td>'
            "</tr>"
        )
    lineage_head = "<th>Lineage</th>" if lineage_map else ""
    head = (
        f"<tr><th>Reference</th><th>Role</th>{lineage_head}<th>Query ANI</th>"
        "<th>Query coverage</th></tr>"
    )
    cap = (
        '<p class="cap">The panel was curated to keep distinct parental lineages and '
        'exclude the query\'s own-lineage relatives. A sibling matches the query '
        'genome-wide (high ANI and coverage) and would mask recombination; a regional '
        'parent (e.g. an env-only reference) has high ANI but low coverage and is kept. '
        'Redundant near-duplicates were dereplicated to representatives.</p>'
    )
    return f'{cap}<div class="scroll"><table class="table">{head}{rows}</table></div>'
