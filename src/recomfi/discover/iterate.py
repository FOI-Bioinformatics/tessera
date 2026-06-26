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
from .fetch import efetch_available
from .run import (
    _base_accession,
    _download,
    _existing_labels,
    collect_candidates,
)


@dataclass
class FillParams:
    query: Path  # query FASTA (not an MSA -- the MSA is rebuilt each round)
    collection: Path  # starting collection directory (left untouched; a copy grows)
    output: Path
    aligner: str = "mafft"
    reference: str | None = None
    max_rounds: int = 3
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
    query_label = strip_sequence_extension(params.query.name)
    params.output.mkdir(parents=True, exist_ok=True)

    # Grow a copy of the collection so the user's input is left untouched.
    collection = params.output / "collection"
    if collection.exists():
        shutil.rmtree(collection)
    shutil.copytree(params.collection, collection)

    exclude = {_base_accession(e) for e in params.exclude}
    # The query's own GenBank record matches itself almost perfectly and would be
    # re-added every round, closing gaps trivially. Its accession is in the query
    # FASTA header (e.g. ">MG572182.1 ..."), so exclude it automatically.
    query_records = read_fasta(str(params.query))
    if query_records:
        own = _base_accession(query_records[0][0])
        exclude.add(own)
        logger.info("Auto-excluding the query's own record '%s' (from its FASTA header).", own)

    cov = CoverageParams.with_defaults(
        params.window_size, floor=params.coverage_floor, rel_drop=params.coverage_rel_drop,
    )

    trace: list[RoundResult] = []
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
        downloaded = _download(candidates, collection, logger)
        rr.added = [c.hit.accession for c in downloaded]
        if not downloaded:
            logger.info("Stopping: no new references available to add.")
            break
    else:
        logger.info("Reached the maximum of %d round(s).", params.max_rounds)

    _write_trace(params.output, trace, logger)
    final_size = sum(1 for _ in collection.iterdir())
    if last_msa is not None:
        logger.info("Writing the final report for the expanded collection...")
        run_recomb(
            RecombParams(
                msa=last_msa, output=params.output, query=query_label,
                window_size=params.window_size, window_step=params.window_step,
                coverage_floor=params.coverage_floor, coverage_rel_drop=params.coverage_rel_drop,
            ),
            logger,
            extra_sections=[("Reference recovery", _progress_section(trace, final_size))],
        )
    logger.info("Final collection (%d references): %s", final_size, collection)
    return trace


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
