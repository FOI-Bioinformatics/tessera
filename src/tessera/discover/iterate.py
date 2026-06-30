"""`tessera fill-references`: iteratively add missing references until coverage stops improving.

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
from ..recomb.pango import crosscheck_html, expand_recombinant, load_alias_key
from ..recomb.regions import DEFAULT_METHODS
from ..recomb.run import RecombParams, run_recomb
from ..recomb.similarity import compute_similarity
from ..recomb.typing import (
    LINEAGES_TSV,
    build_lineage_map,
    dominant_lineage_token,
    first_header,
    genotype_from_title,
    lineage_map_from_rows,
    organism_from_title,
    titles_from_collection,
    write_lineage_map,
)
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
    seed_source: str = "blast"  # fresh-start source: "blast" | "local" | "ncbi-virus" | "nextclade"
    seed_mode: str = "windowed"  # blast seeding: "whole" | "windowed" | "parents"
    seed_hits: int = 10  # BLAST hits to keep per seed search (per window, or whole-query)
    seed_window: int = 1500  # window width (bp) for windowed/parents/pool seeding
    candidate_pool: Path | None = None  # local genome pool (seed_source="local")
    taxon: str | None = None  # taxon for NCBI Virus ("ncbi-virus"); auto-detected when None
    organism: str | None = None  # species name for the report header; falls back to taxon
    nextclade_dataset: str | None = None  # Nextclade dataset path; None = auto-detect
    source_refseq: bool = True  # NCBI Virus: fetch the RefSeq set (else all complete genomes)
    seed_keep_siblings: bool = False  # keep the query's siblings when selecting from a pool
    # When BLAST seeding finds only siblings (a saturated lineage), switch to recruiting
    # by NCBI Virus taxonomy diversity instead, so the parental lineages are present.
    auto_diversify: bool = True
    negative_lineage: bool = True  # on saturation, re-BLAST excluding the query's lineage
    fetch_limit: int = 2000  # cap a broad NCBI Virus fetch (complete genomes)
    derep_diverse_ani: float = 95.0  # aggressive one-per-lineage dereplication of a pool
    cache_dir: Path | None = None  # where to cache fetched NCBI Virus panels (None = default)
    window_size: int = 1000
    window_step: int = 100
    coverage_floor: float | None = None
    coverage_rel_drop: float = 0.05
    max_hits: int = 5
    top_gaps: int = 3
    subtile: int = 400  # search gaps longer than this in overlapping sub-intervals (0 = off)
    email: str | None = None
    exclude: tuple[str, ...] = ()
    keep_self_hits: bool = False
    threads: int = 4
    min_improvement: float = 0.01  # stop if the worst gap's best-sim gains less than this
    curate: bool = False  # drop the query's siblings and dereplicate each round (skani/skDER)
    sibling_margin: float = 3.0  # query-ANI must beat the backbone's by this many % to be a sibling
    af_min: float = 80.0  # ... over at least this % of the query (whole-genome match)
    derep_ani: float = 99.0  # skDER: collapse references >= this ANI to one representative
    report: bool = True  # call recombination detection after building the panel (False = stop)
    methods: tuple[str, ...] = DEFAULT_METHODS  # region caller(s) for the detection step
    # Nextclade seeding: use one denoised consensus genome per clade as the pool (a stable
    # per-lineage reference) instead of every tree tip. See discover/nextclade.build_pool.
    pool_consensus: bool = False
    lineage_map: Path | None = None  # user TSV (accession<TAB>genotype) to type references

    @classmethod
    def for_detection(
        cls,
        *,
        query: Path,
        output: Path,
        aligner: str = "mafft",
        max_rounds: int = 2,
        window_size: int = 1000,
        window_step: int = 100,
        email: str | None = None,
        threads: int = 4,
        cache_dir: Path | None = None,
        taxon: str | None = None,
        candidate_pool: Path | None = None,
        nextclade: bool = False,
        nextclade_dataset: str | None = None,
        methods: tuple[str, ...] = DEFAULT_METHODS,
        pool_consensus: bool = False,
        organism: str | None = None,
        lineage_map: Path | None = None,
    ) -> FillParams:
        """Build the detection-tuned preset over fill-references.

        Recruit the parental lineages (drop the query's siblings, curate the
        panel) so a recombinant query whose own lineage is common in NCBI is
        not masked. The seed source follows the supplied input: a local pool
        wins, then a Nextclade dataset, then the default BLAST recruitment.
        This is the single home for the preset so ``tessera detect`` does not
        duplicate (and drift from) the field list.
        """
        if candidate_pool:
            seed_source = "local"
        elif nextclade or nextclade_dataset:
            seed_source = "nextclade"
        else:
            seed_source = "blast"
        return cls(
            query=query, collection=None, output=output,
            aligner=aligner, max_rounds=max_rounds,
            window_size=window_size, window_step=window_step,
            email=email, threads=threads, cache_dir=cache_dir,
            seed_source=seed_source, nextclade_dataset=nextclade_dataset,
            candidate_pool=candidate_pool, taxon=taxon,
            seed_mode="parents", curate=True,
            auto_diversify=True, negative_lineage=True,
            methods=methods, pool_consensus=pool_consensus,
            organism=organism, lineage_map=lineage_map,
        )


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
            subtile=params.subtile,
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

    # Type the recruited references: mine a genotype from each genome's header (NCBI
    # lineage note or GenBank title), overlaid by a user-supplied --lineage-map, and
    # persist it next to the panel so the report -- and a later standalone recomb --
    # can name parents by lineage instead of by bare accession.
    coll_files = [p for p in collection.iterdir() if p.is_file()]
    lineage_rows = build_lineage_map(
        user_tsv=params.lineage_map,
        title_by_label=titles_from_collection(coll_files),
        organism=params.taxon,
    )
    # Type the query itself from its own header so the verdict can name its lineage
    # (and a recombinant query can be cross-checked against its designated parents).
    query_lineage = genotype_from_title(first_header(params.query), params.taxon)
    if query_lineage:
        lineage_rows.append((query_label, query_lineage, "query"))
    write_lineage_map(params.output / LINEAGES_TSV, lineage_rows)
    lineage_map = lineage_map_from_rows(lineage_rows)
    if lineage_map:
        logger.info("Typed %d reference(s) with a lineage/genotype name.", len(lineage_map))

    extra_sections = [("Reference recovery", _progress_section(trace, final_size))]
    if params.curate and panel_rows:
        table = list(panel_rows.values())
        write_panel_tsv(params.output / "panel_lineages.tsv", table, logger, lineage_map)
        extra_sections.append(("Reference panel", panel_table_html(table, lineage_map)))
    # For a recombinant SARS-CoV-2 query, list the Pango-designated parents alongside
    # the recruited ones as a cross-check (best-effort; needs the cached alias key).
    if params.report and query_lineage and query_lineage.split(".")[0].startswith("X"):
        parents = expand_recombinant(
            query_lineage, load_alias_key(cache_override=params.cache_dir, logger=logger)
        )
        if parents:
            logger.info("Pango designates %s as a recombinant of %s.",
                        query_lineage, ", ".join(parents))
            extra_sections.append(
                ("Pango cross-check", crosscheck_html(query_lineage, parents))
            )
    # Publish the final alignment under a stable name so a separate detection step
    # (tessera recomb) can consume the panel without guessing the last round number.
    panel_msa = params.output / "panel.msa.fasta"
    if last_msa is not None:
        shutil.copyfile(last_msa, panel_msa)
    if last_msa is not None and params.report:
        logger.info("Writing the final report for the expanded collection...")
        run_recomb(
            RecombParams(
                msa=panel_msa, output=params.output, query=query_label,
                window_size=params.window_size, window_step=params.window_step,
                coverage_floor=params.coverage_floor, coverage_rel_drop=params.coverage_rel_drop,
                methods=params.methods, lineage_map=lineage_map or None,
                organism=params.organism or params.taxon,
            ),
            logger,
            extra_sections=extra_sections,
        )
    elif last_msa is not None:
        logger.info(
            "Panel built (no detection run). To call recombination, run:\n"
            "  tessera recomb -i %s -q %s -o %s",
            panel_msa, query_label, params.output,
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
SEED_MIN_DIVERSE = 3  # below this many RefSeq genomes, broaden to complete genomes


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
    if params.seed_source in ("local", "ncbi-virus", "nextclade"):
        _seed_from_pool(params, collection, logger)
        return
    query_seq = query_records[0][1].replace("-", "")
    if params.seed_mode == "whole":
        accessions, saturated = _seed_whole(query_seq, params, exclude, logger)
    else:
        accessions, saturated = _seed_windowed(
            query_seq, params, exclude, logger,
            drop_siblings=params.seed_mode == "parents",
        )
    from .pool import datasets_available

    # Saturation (only siblings found): first try a negative-lineage search -- exclude
    # the query's own over-represented lineage so the divergent parents surface. Works
    # for mega-taxa (no download). Fall back to NCBI Virus taxonomy diversity.
    if saturated and params.negative_lineage:
        neg = _seed_negative_lineage(query_seq, params, exclude, logger)
        if neg:
            logger.info("Negative-lineage search recruited %d reference(s).", len(neg))
            accessions, saturated = neg, False
    if saturated and params.auto_diversify and datasets_available():
        logger.info(
            "BLAST seeding recruited only siblings of the query (its lineage saturates "
            "NCBI); switching to NCBI Virus taxonomy diversity to recruit the parents..."
        )
        _seed_from_pool(params, collection, logger, force_ncbi=True)
        return
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


def _seed_from_pool(
    params: FillParams, collection: Path, logger: logging.Logger, *, force_ncbi: bool = False
) -> None:
    """Seed by regional selection from a genome pool (local directory or NCBI Virus)."""
    from .pool import iter_pool_genomes

    if params.seed_source == "nextclade" and not force_ncbi:
        genomes = _fetch_nextclade(params, logger)
    elif params.seed_source == "local" and not force_ncbi:
        if params.candidate_pool is None:
            raise UserInputError("--seed-source local needs --candidate-pool <dir>.")
        genomes = iter_pool_genomes(params.candidate_pool)
    else:  # ncbi-virus (or auto-switched into it)
        genomes = _fetch_diverse(params, logger)
    _copy_into(_select_from(params, genomes, logger).selected, collection, logger)


def _fetch_nextclade(params: FillParams, logger: logging.Logger) -> list[Path]:
    """Reconstruct a candidate pool from a Nextclade dataset (cached per path@tag)."""
    from ..core.cache import nextclade_cache
    from .nextclade import build_pool, resolve_dataset

    dataset = resolve_dataset(
        params.query, params.nextclade_dataset, email=params.email, logger=logger
    )
    cache = nextclade_cache(dataset.path, dataset.tag, override=params.cache_dir)
    return build_pool(dataset, cache_dir=cache, logger=logger,
                      per_clade_consensus=params.pool_consensus)


def _fetch_diverse(params: FillParams, logger: logging.Logger) -> list[Path]:
    """Fetch a diverse taxon-scoped set from NCBI Virus (cached per taxon): the RefSeq
    representative set, broadening to the full set of complete genomes when RefSeq is
    too thin.

    The fetched panel is cached on disk by taxon, so a repeat run skips the network.
    The full set is dereplicated locally during regional selection, so a heavily
    sequenced taxon is reduced to diverse representatives rather than truncated to a
    biased sample.
    """
    from ..core.cache import cached_genomes, ncbi_virus_cache
    from .pool import detect_taxon, fetch_ncbi_virus

    taxon = params.taxon or detect_taxon(params.query, email=params.email, logger=logger)
    cache = ncbi_virus_cache(taxon, override=params.cache_dir)
    existing = cached_genomes(cache)
    if existing:
        logger.info(
            "Using the cached NCBI Virus panel for '%s' (%d genome(s)): %s",
            taxon, len(existing), cache,
        )
        return existing
    cache.mkdir(parents=True, exist_ok=True)
    if params.source_refseq:
        fetched = fetch_ncbi_virus(taxon, cache, refseq=True, logger=logger)
        if len(fetched) >= SEED_MIN_DIVERSE:
            return fetched
        logger.info(
            "RefSeq set for '%s' has only %d genome(s); broadening to complete genomes.",
            taxon, len(fetched),
        )
        for g in fetched:
            g.unlink()
    fetched = fetch_ncbi_virus(
        taxon, cache, refseq=False, complete_only=True, logger=logger
    )
    if len(fetched) > params.fetch_limit:
        logger.info(
            "Fetched %d '%s' complete genomes (a heavily sequenced taxon); dereplicating "
            "locally to a diverse panel -- this can take a few minutes.",
            len(fetched), taxon,
        )
    return fetched


def _select_from(params: FillParams, genomes: list[Path], logger: logging.Logger):
    from .pool import select_regional

    return select_regional(
        params.query, genomes, window=params.seed_window, per_window=SEED_PER_WINDOW,
        derep_ani=params.derep_diverse_ani, drop_siblings=not params.seed_keep_siblings,
        logger=logger,
    )


def _copy_into(genomes: list[Path], collection: Path, logger: logging.Logger) -> None:
    if not genomes:
        return
    logger.info("Seeding %d reference(s): %s",
                len(genomes), ", ".join(strip_sequence_extension(g.name) for g in genomes))
    for g in genomes:
        shutil.copy(g, collection / g.name)


def _blast_or_none(
    seq: str, params: FillParams, logger: logging.Logger, *, entrez_query: str | None = None
) -> list:
    try:
        return blast_subsequence(
            seq, max_hits=params.seed_hits, logger=logger, email=params.email,
            entrez_query=entrez_query,
        )
    except BlastError as exc:
        logger.warning("  BLAST failed for a seed search, skipping: %s", exc)
        return []


# Lineage-token extraction lives in recomb.typing (single source of truth). These thin
# shims keep the hit-object signature used by the negative-lineage seeder and its tests.
_organism_from_title = organism_from_title


def _dominant_lineage_token(hits: list, min_frac: float = 0.5) -> str | None:
    """The lineage designation shared by most top hits (the query's own lineage)."""
    return dominant_lineage_token([hit.title for hit in hits], min_frac=min_frac)


def _seed_negative_lineage(
    query_seq: str, params: FillParams, exclude: set[str], logger: logging.Logger
) -> list[str]:
    """Seed by re-BLASTing each window with the query's own lineage excluded.

    Identifies the dominant lineage token from a whole-query BLAST, then searches each
    window with ``<organism> NOT <token>`` so the best *non-lineage* hit per region --
    the regional parent -- surfaces. Returns ``[]`` if no lineage token can be
    identified (the caller then falls back to taxonomy diversity).
    """
    hits = _blast_or_none(query_seq[:1500], params, logger)
    token = _dominant_lineage_token(hits)
    if token is None:
        logger.info("Could not identify the query's lineage from BLAST titles; skipping.")
        return []
    organism = _organism_from_title(hits[0].title) if hits else ""
    entrez_query = (f'"{organism}"[Organism] NOT "{token}"') if organism else f'NOT "{token}"'
    logger.info(
        "Saturation: excluding the query's lineage '%s' to recruit divergent parents...", token,
    )
    window = max(params.seed_window, MIN_SUBSEQ)
    selected: list[str] = []
    for i in range(0, len(query_seq), window):
        chunk = query_seq[i : i + window]
        if len(chunk) < MIN_SUBSEQ:
            continue
        for hit in _blast_or_none(chunk, params, logger, entrez_query=entrez_query):
            if _keep(hit, exclude, params) and not _is_sibling_hit(hit) \
                    and hit.accession not in selected:
                selected.append(hit.accession)
    return selected[:SEED_TOTAL_CAP]


def _keep(hit, exclude: set[str], params: FillParams) -> bool:
    return _base_accession(hit.accession) not in exclude and not _is_self_hit(
        hit, params.keep_self_hits
    )


def _seed_whole(
    query_seq: str, params: FillParams, exclude: set[str], logger: logging.Logger
) -> tuple[list[str], bool]:
    """Seed from a single whole-query BLAST; returns ``(accessions, saturated)``.

    ``saturated`` is True when every kept hit is a near-identical sibling of the query
    (the query's lineage is over-represented, so no parental lineage was recruited).
    """
    logger.info(
        "No starting references: seeding from a whole-query NCBI search (top %d hit(s))...",
        params.seed_hits,
    )
    hits = _blast_or_none(query_seq, params, logger)
    if not hits:
        raise UserInputError(
            "Could not seed the collection (whole-query BLAST returned nothing)."
        )
    kept = [h for h in hits if _keep(h, exclude, params)]
    out: list[str] = []
    for hit in kept:
        if hit.accession not in out:
            out.append(hit.accession)
    saturated = bool(kept) and all(_is_sibling_hit(h) for h in kept)
    return out[:SEED_TOTAL_CAP], saturated


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
) -> tuple[list[str], bool]:
    """Seed from per-window BLAST, optionally suppressing siblings (parents mode).

    Returns ``(accessions, saturated)``; ``saturated`` is True when no window yielded
    a non-sibling hit (the query's lineage saturates NCBI, so no parent was found).
    """
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
    saturated = bool(per_window_hits) and not any(
        not _is_sibling_hit(h) for kept in per_window_hits for h in kept
    )
    return selected[:SEED_TOTAL_CAP], saturated


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
