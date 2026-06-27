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
    seed_hits: int = 10  # when starting fresh, how many whole-query BLAST hits to seed with
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


def _seed_collection(
    params: FillParams,
    collection: Path,
    query_records: list[tuple[str, str]],
    exclude: set[str],
    logger: logging.Logger,
) -> None:
    """Seed an empty collection by BLASTing the whole query and downloading the hits.

    The query's own record (near-identical, full-length) and any ``--exclude``
    accessions are skipped, as in the per-gap search.
    """
    if not query_records:
        raise UserInputError(f"Query FASTA {params.query} has no sequence to search with.")
    query_seq = query_records[0][1].replace("-", "")
    logger.info(
        "No starting references: seeding from a whole-query NCBI search (top %d hit(s))...",
        params.seed_hits,
    )
    try:
        hits = blast_subsequence(
            query_seq, max_hits=params.seed_hits, logger=logger, email=params.email
        )
    except BlastError as exc:
        raise UserInputError(f"Could not seed the collection (BLAST failed): {exc}") from exc
    seeded = 0
    for hit in hits:
        if _base_accession(hit.accession) in exclude or _is_self_hit(hit, params.keep_self_hits):
            continue
        try:
            out = efetch_fasta(hit.accession, collection, logger)
            logger.info("  + seed %s -> %s", hit.accession, out.name)
            seeded += 1
        except Exception as exc:  # noqa: BLE001 - report and continue the batch
            logger.warning("  ! failed to download seed %s: %s", hit.accession, exc)
    logger.info("Seeded %d reference(s) to start from.", seeded)


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
