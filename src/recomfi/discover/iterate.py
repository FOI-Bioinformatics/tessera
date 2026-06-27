"""`recomfi fill-references`: iteratively add missing references until coverage stops improving.

Each round rebuilds the MSA from the (growing) collection, scans it for coverage
gaps, BLASTs the worst gaps against NCBI, and downloads the best new reference per
gap into the collection. The loop stops when the gaps close, when no new reference
can be found, when coverage stops improving (a stubborn residual is reported, not
chased forever), or at ``max_rounds``.

Because every round rebuilds the alignment, this needs an aligner binary and
Entrez Direct, and it contacts NCBI over the network.
"""

from __future__ import annotations

import html
import logging
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from ..core.errors import UserInputError
from ..core.io import read_fasta, strip_sequence_extension
from ..msa.build import MsaParams, build_msa
from ..recomb.coverage import CoverageParams, call_coverage_gaps
from ..recomb.run import RecombParams, run_recomb
from ..recomb.similarity import compute_similarity
from .blast import BlastError, blast_subsequence
from .fetch import efetch_available, efetch_fasta
from .panel import (
    curate_collection_dir,
    panel_table_html,
    pick_backbone,
    skani_available,
    write_panel_tsv,
)
from .run import (
    MIN_SUBSEQ,
    _base_accession,
    _download,
    _existing_labels,
    _is_self_hit,
    collect_candidates,
)


@dataclass
class FillParams:
    query: Path  # query FASTA (not an MSA -- the MSA is rebuilt each round)
    collection: Path | None  # starting collection (a copy grows); None = start fresh from NCBI
    output: Path
    aligner: str = "mafft"
    reference: str | None = None
    max_rounds: int = 3
    seed_mode: str = "windowed"  # fresh-start seeding: "whole" | "windowed" | "parents"
    seed_hits: int = 10  # BLAST hits to keep per seed search (per window, or whole-query)
    seed_window: int = 1500  # window width (bp) for windowed/parents seeding
    window_size: int = 1000
    window_step: int = 100
    coverage_floor: float | None = None
    coverage_rel_drop: float = 0.05
    max_hits: int = 5
    top_gaps: int = 3
    email: str | None = None
    exclude: tuple[str, ...] = ()
    keep_self_hits: bool = False
    threads: int = 4
    min_improvement: float = 0.01  # stop if the worst gap's best-sim gains less than this
    curate: bool = False  # drop the query's siblings and dereplicate each round (skani/skDER)
    sibling_margin: float = 3.0  # query-ANI must beat the backbone's by this many % to be a sibling
    af_min: float = 80.0  # ... over at least this % of the query (whole-genome match)
    derep_ani: float = 99.0  # skDER: collapse references >= this ANI to one representative


@dataclass
class RoundResult:
    round: int
    n_gaps: int
    undercovered_bp: int
    worst_best_sim: float
    threshold: float
    added: list[str] = field(default_factory=list)


def fill_references(params: FillParams, logger: logging.Logger) -> list[RoundResult]:
    """Drive the build -> scan -> find -> download loop; returns the per-round trace."""
    if not efetch_available():
        raise UserInputError(
            "fill-references downloads genomes and needs Entrez Direct. Install with: "
            "conda install -c bioconda entrez-direct"
        )
    if params.curate and not skani_available():
        raise UserInputError(
            "--curate needs skani. Install with: conda install -c bioconda skani skder"
        )
    query_label = strip_sequence_extension(params.query.name)
    params.output.mkdir(parents=True, exist_ok=True)

    # Grow a copy of the collection so the user's input is left untouched. With no
    # starting collection, begin with an empty directory and seed it below.
    collection = params.output / "collection"
    if collection.exists():
        shutil.rmtree(collection)
    if params.collection is not None:
        shutil.copytree(params.collection, collection)
    else:
        collection.mkdir(parents=True)

    exclude = {_base_accession(e) for e in params.exclude}
    # The query's own GenBank record matches itself almost perfectly and would be
    # re-added every round, closing gaps trivially. Its accession is in the query
    # FASTA header (e.g. ">MG572182.1 ..."), so exclude it automatically.
    query_records = read_fasta(str(params.query))
    if query_records:
        own = _base_accession(query_records[0][0])
        exclude.add(own)
        logger.info("Auto-excluding the query's own record '%s' (from its FASTA header).", own)

    # Start fresh: seed an empty collection from a whole-query NCBI search so the
    # first MSA has something to align against.
    if not any(collection.iterdir()):
        _seed_collection(params, collection, query_records, exclude, logger)
        if not any(collection.iterdir()):
            raise UserInputError(
                "Could not seed any reference from NCBI; provide a starting --collection."
            )

    cov = CoverageParams.with_defaults(
        params.window_size, floor=params.coverage_floor, rel_drop=params.coverage_rel_drop,
    )

    trace: list[RoundResult] = []
    # Curated-panel state: a per-genome role table accumulated across rounds (a
    # dropped genome keeps the role it had when removed, even after later rounds).
    panel_rows: dict[str, dict] = {}
    last_msa: Path | None = None
    prev_worst: float | None = None
    for rnd in range(1, params.max_rounds + 1):
        msa = params.output / f"round{rnd}.msa.fasta"
        logger.info("=== Round %d: building MSA from %d reference(s) ===",
                    rnd, sum(1 for _ in collection.iterdir()))
        build_msa(
            MsaParams(
                query=params.query, collection=collection, output=msa,
                aligner=params.aligner, reference=params.reference, threads=params.threads,
            ),
            logger,
        )
        last_msa = msa
        result = compute_similarity(
            str(msa), query_label,
            window_size=params.window_size, window_step=params.window_step,
        )
        gaps, threshold = call_coverage_gaps(result, params.window_size, cov)
        undercovered = sum(max(0, g.query_end - g.query_start) for g in gaps)
        worst = min((g.mean_best for g in gaps), default=1.0)
        rr = RoundResult(rnd, len(gaps), undercovered, round(worst, 4), round(threshold, 4))
        trace.append(rr)
        logger.info(
            "Round %d: %d gap(s), %d bp under-covered, worst best-similarity %.3f.",
            rnd, len(gaps), undercovered, worst,
        )

        if not gaps:
            logger.info("Converged: no coverage gaps remain.")
            break
        if prev_worst is not None and (worst - prev_worst) < params.min_improvement:
            logger.info(
                "Stopping: coverage stopped improving (worst best-sim %.3f vs %.3f last round).",
                worst, prev_worst,
            )
            break
        prev_worst = worst

        records = dict(read_fasta(str(msa)))
        query_row = records[query_label]
        existing = _existing_labels(result, collection, query_label)
        targets = sorted(gaps, key=lambda g: g.length_bp, reverse=True)[: params.top_gaps]
        candidates = collect_candidates(
            targets, query_row, existing,
            max_hits=params.max_hits, email=params.email,
            exclude=exclude, keep_self_hits=params.keep_self_hits, logger=logger,
        )
        # Pick the backbone from the pre-download (curated, sibling-free) collection
        # so a freshly-downloaded sibling cannot be mistaken for it.
        backbone = None
        if params.curate:
            backbone = pick_backbone(
                params.query, [p for p in collection.iterdir() if p.is_file()],
                af_min=params.af_min, logger=logger,
            )
        downloaded = _download(candidates, collection, logger)
        rr.added = [c.hit.accession for c in downloaded]
        if not downloaded:
            logger.info("Stopping: no new references available to add.")
            break
        if params.curate and backbone is not None:
            curation = curate_collection_dir(
                params.query, collection, backbone,
                ani_margin=params.sibling_margin, af_min=params.af_min,
                derep_ani=params.derep_ani, logger=logger,
            )
            for row in curation.table:
                panel_rows[row["genome"]] = row
            dropped = {c.hit.accession for c in downloaded} - {
                strip_sequence_extension(p.name) for p in collection.iterdir() if p.is_file()
            }
            rr.added = [a for a in rr.added if a not in dropped]
    else:
        logger.info("Reached the maximum of %d round(s).", params.max_rounds)

    _write_trace(params.output, trace, logger)
    final_size = sum(1 for _ in collection.iterdir())
    extra_sections = [("Reference recovery", _progress_section(trace, final_size))]
    if params.curate and panel_rows:
        table = list(panel_rows.values())
        write_panel_tsv(params.output / "panel_lineages.tsv", table, logger)
        extra_sections.append(("Reference panel", panel_table_html(table)))
    if last_msa is not None:
        logger.info("Writing the final report for the expanded collection...")
        run_recomb(
            RecombParams(
                msa=last_msa, output=params.output, query=query_label,
                window_size=params.window_size, window_step=params.window_step,
                coverage_floor=params.coverage_floor, coverage_rel_drop=params.coverage_rel_drop,
            ),
            logger,
            extra_sections=extra_sections,
        )
    logger.info("Final collection (%d references): %s", final_size, collection)
    return trace


# Fresh-start seeding. A whole-query search returns the query's closest *whole-genome*
# relatives; for a recombinant query those are its own-lineage siblings, which cover
# the whole query and mask the recombination. The parents that donated each region
# are *regional* best matches, so windowed/parents seeding searches per region instead.
SEED_SIBLING_IDENTITY = 95.0  # a hit this identical over near-full coverage is a sibling
SEED_SIBLING_COVERAGE = 90.0
SEED_PER_WINDOW = 2  # regional candidates to keep per window
SEED_TOTAL_CAP = 25  # bound the number of genomes seeded


def _seed_collection(
    params: FillParams,
    collection: Path,
    query_records: list[tuple[str, str]],
    exclude: set[str],
    logger: logging.Logger,
) -> None:
    """Seed an empty collection by searching NCBI with the query and downloading hits.

    Three strategies (``--seed-mode``):

    - ``whole`` -- one BLAST of the whole query; seeds its closest whole-genome
      relatives. Best for a non-recombinant query.
    - ``windowed`` (default) -- BLAST each window and seed the per-window best hits;
      surfaces regional matches without suppressing anything (safe for any query).
    - ``parents`` -- like windowed, but drop near-identical full-coverage hits (the
      query's siblings) so each region contributes its best *divergent* source. This
      recruits the parental lineages of a recombinant query instead of its siblings.

    The query's own record and any ``--exclude`` accessions are always skipped.
    """
    if not query_records:
        raise UserInputError(f"Query FASTA {params.query} has no sequence to search with.")
    query_seq = query_records[0][1].replace("-", "")
    if params.seed_mode == "whole":
        accessions = _seed_whole(query_seq, params, exclude, logger)
    else:
        accessions = _seed_windowed(
            query_seq, params, exclude, logger,
            drop_siblings=params.seed_mode == "parents",
        )
    if not accessions:
        return  # caller reports that nothing could be seeded
    logger.info("Seeding %d reference(s): %s", len(accessions), ", ".join(accessions))
    seeded = 0
    for acc in accessions:
        try:
            out = efetch_fasta(acc, collection, logger)
            logger.info("  + seed %s -> %s", acc, out.name)
            seeded += 1
        except Exception as exc:  # noqa: BLE001 - report and continue the batch
            logger.warning("  ! failed to download seed %s: %s", acc, exc)
    logger.info("Seeded %d reference(s) to start from.", seeded)


def _blast_or_none(seq: str, params: FillParams, logger: logging.Logger) -> list:
    try:
        return blast_subsequence(seq, max_hits=params.seed_hits, logger=logger, email=params.email)
    except BlastError as exc:
        logger.warning("  BLAST failed for a seed search, skipping: %s", exc)
        return []


def _keep(hit, exclude: set[str], params: FillParams) -> bool:
    return _base_accession(hit.accession) not in exclude and not _is_self_hit(
        hit, params.keep_self_hits
    )


def _seed_whole(
    query_seq: str, params: FillParams, exclude: set[str], logger: logging.Logger
) -> list[str]:
    """Seed from a single whole-query BLAST (the query's closest whole-genome relatives)."""
    logger.info(
        "No starting references: seeding from a whole-query NCBI search (top %d hit(s))...",
        params.seed_hits,
    )
    hits = _blast_or_none(query_seq, params, logger)
    if not hits:
        raise UserInputError(
            "Could not seed the collection (whole-query BLAST returned nothing)."
        )
    out: list[str] = []
    for hit in hits:
        acc = hit.accession
        if _keep(hit, exclude, params) and acc not in out:
            out.append(acc)
    return out[:SEED_TOTAL_CAP]


def _is_sibling_hit(hit) -> bool:
    """A near-identical, near-full-length hit -- the query's own lineage (sibling)."""
    return (
        hit.pct_identity >= SEED_SIBLING_IDENTITY
        and hit.query_coverage >= SEED_SIBLING_COVERAGE
    )


def _seed_windowed(
    query_seq: str,
    params: FillParams,
    exclude: set[str],
    logger: logging.Logger,
    *,
    drop_siblings: bool,
) -> list[str]:
    """Seed from per-window BLAST, optionally suppressing siblings (parents mode)."""
    window = max(params.seed_window, MIN_SUBSEQ)
    chunks = [query_seq[i : i + window] for i in range(0, len(query_seq), window)]
    if len(chunks) >= 2 and len(chunks[-1]) < MIN_SUBSEQ:
        chunks[-2] += chunks.pop()  # fold a tiny tail into the previous window
    label = "regional parents" if drop_siblings else "per-window best hits"
    logger.info(
        "No starting references: seeding from a %d-window NCBI search (%s)...",
        len(chunks), label,
    )
    per_window_hits: list[list] = []
    for chunk in chunks:
        if len(chunk) < MIN_SUBSEQ:
            continue
        kept = [h for h in _blast_or_none(chunk, params, logger) if _keep(h, exclude, params)]
        per_window_hits.append(kept)

    selected: list[str] = []
    dropped_siblings: set[str] = set()

    def add(acc: str) -> None:
        if acc not in selected:
            selected.append(acc)

    for kept in per_window_hits:
        pool = kept
        if drop_siblings:
            pool = [h for h in kept if not _is_sibling_hit(h)]
            dropped_siblings.update(h.accession for h in kept if _is_sibling_hit(h))
        for hit in pool[:SEED_PER_WINDOW]:
            add(hit.accession)

    # Parents mode left nothing (only siblings anywhere, or a clonal query): fall back
    # to the per-window best hits so the collection is never empty.
    if not selected:
        if drop_siblings and dropped_siblings:
            logger.info(
                "Every window's hits were siblings of the query; seeding the best hits "
                "instead (no distinct parental lineage found in NCBI)."
            )
        for kept in per_window_hits:
            if kept:
                add(kept[0].accession)
    elif drop_siblings and dropped_siblings:
        logger.info(
            "Suppressed %d sibling hit(s) (>= %.0f%% identity over near-full coverage) to "
            "recruit parental lineages.", len(dropped_siblings), SEED_SIBLING_IDENTITY,
        )
    return selected[:SEED_TOTAL_CAP]


def _progress_section(trace: list[RoundResult], final_size: int) -> str:
    """An HTML block summarising each fill round, for the report."""
    rows = ""
    for r in trace:
        added = ", ".join(html.escape(a) for a in r.added) if r.added else "&mdash;"
        rows += (
            "<tr>"
            f'<td class="num">{r.round}</td>'
            f'<td class="num">{r.n_gaps}</td>'
            f'<td class="num">{r.undercovered_bp:,} bp</td>'
            f'<td class="num strong">{r.worst_best_sim:.3f}</td>'
            f'<td class="lbl">{added}</td>'
            "</tr>"
        )
    head = (
        "<tr><th>Round</th><th>Gaps</th><th>Under-covered</th>"
        "<th>Worst best sim</th><th>References added</th></tr>"
    )
    cap = (
        '<p class="cap">Each round rebuilt the alignment from the growing collection, '
        'searched the worst gaps on NCBI, and added the best new reference per gap &mdash; '
        f'ending with {final_size} references. The run stops when the gaps close, no new '
        'reference is found, or a round no longer improves the worst gap.</p>'
    )
    return f'{cap}<div class="scroll"><table class="table">{head}{rows}</table></div>'


def _write_trace(output: Path, trace: list[RoundResult], logger: logging.Logger) -> None:
    path = output / "fill_summary.tsv"
    logger.info("Writing fill summary: %s", path)
    header = ["round", "n_gaps", "undercovered_bp", "worst_best_similarity",
              "threshold", "references_added"]
    with open(path, "w") as fo:
        fo.write("\t".join(header) + "\n")
        for r in trace:
            fo.write("\t".join(map(str, [
                r.round, r.n_gaps, r.undercovered_bp, r.worst_best_sim,
                r.threshold, ",".join(r.added) or "-",
            ])) + "\n")
