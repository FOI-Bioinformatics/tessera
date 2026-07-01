"""Build a candidate panel from a genome pool by regional matching (skani).

BLAST against ``nt`` ranks hits by similarity, so an over-represented lineage
saturates every region and the divergent parents never surface (the HIV CRF01_AE
case). This module takes a *finite* pool of genomes instead -- a local directory
the user controls, or a taxon-scoped set fetched from NCBI Virus -- and selects a
diverse, query-relevant panel from it. Because the pool is finite it can be
dereplicated, which collapses an over-represented lineage to a single
representative, so saturation no longer hides the parents.

Selection is *regional*: the query is cut into windows and each window is matched
against the pool with skani, so a genome that donated only one region is recruited
for that region rather than being out-ranked genome-wide. The focused panel then
feeds the normal alignment and HMM caller for the precise call; pair it with
``--curate`` to drop any sibling that survives dereplication.
"""

from __future__ import annotations

import json
import logging
import tempfile
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.errors import UserInputError
from ..core.io import read_fasta, strip_sequence_extension, write_fasta_record
from ..core.plugins import ToolCapabilities
from ..core.process import run_tool
from ..recomb.typing import is_recombinant_lineage
from .blast import BlastError, blast_subsequence
from .panel import dereplicate, skani_available, skani_query_ani

DATASETS = ToolCapabilities(
    name="datasets", conda=("conda-forge::ncbi-datasets-cli",),
    required_binaries=(BinarySpec("datasets", version_args=("--version",)),),
)

DEFAULT_SEED_WINDOW = 1500  # query window width (bp) for regional matching
DEFAULT_PER_WINDOW = 2  # candidates to keep per window
DEFAULT_ANI_FLOOR = 70.0  # drop pool genomes below this genome-wide ANI as unrelated
DEFAULT_DEREP_ANI = 99.0  # collapse pool genomes at or above this ANI to one representative
NCBI_LINEAGES_TSV = "ncbi_lineages.tsv"  # sidecar: accession<TAB>datasets-lineage
# A pool genome this identical to the query genome-wide, over near-full coverage, is a
# sibling (the query's own lineage). It is closest in every region and so out-ranks the
# parents everywhere; dropping it lets the regional parents win. The absolute cutoff
# suits divergent-lineage recombinants (e.g. HIV subtypes); for closely related parents
# (e.g. SARS-CoV-2 sublineages, which have no closer masking sibling) disable the drop.
DEFAULT_SIBLING_IDENTITY = 95.0
DEFAULT_SIBLING_COVERAGE = 90.0


@dataclass
class PoolSelection:
    selected: list[Path]  # the focused, query-relevant panel
    table: list[dict] = field(default_factory=list)  # per-window winners (for reporting)


def datasets_available() -> bool:
    import shutil

    return shutil.which("datasets") is not None


def iter_pool_genomes(pool: Path) -> list[Path]:
    """List the genome FASTA files in a pool directory."""
    if not pool.is_dir():
        raise UserInputError(f"Candidate pool {pool} is not a directory.")
    genomes = sorted(
        p for p in pool.iterdir()
        if p.is_file() and p.suffix.lower() in {".fasta", ".fa", ".fna", ".gz"}
    )
    if not genomes:
        raise UserInputError(f"Candidate pool {pool} has no FASTA genome files.")
    return genomes


def detect_taxon(query_fasta: Path, *, email: str | None, logger: logging.Logger) -> str:
    """Infer the query's taxon name from the top BLAST hit's title (one search)."""
    records = read_fasta(query_fasta)
    if not records:
        raise UserInputError(f"Query FASTA {query_fasta} has no sequence.")
    probe = records[0][1].replace("-", "")[:1500]
    logger.info("Detecting the query's taxon from a short BLAST probe...")
    try:
        hits = blast_subsequence(probe, max_hits=1, logger=logger, email=email)
    except BlastError as exc:
        raise UserInputError(
            f"Could not detect the taxon (BLAST failed): {exc}. Pass --taxon."
        ) from exc
    if not hits:
        raise UserInputError("Could not detect the taxon from BLAST; pass --taxon.")
    # Title looks like 'Human immunodeficiency virus 1 isolate ... , complete genome'.
    taxon = hits[0].title.split(",")[0].split(" isolate ")[0].split(" strain ")[0].strip()
    logger.info("Detected taxon: %s (from %s).", taxon, hits[0].accession)
    return taxon


def _scope_flags(refseq: bool, complete_only: bool, released_after: str | None) -> list[str]:
    flags: list[str] = []
    if refseq:
        flags.append("--refseq")
    if complete_only:
        flags.append("--complete-only")
    if released_after:
        flags += ["--released-after", released_after]
    return flags


def write_ncbi_lineages(dest: Path, labels: dict[str, str]) -> Path | None:
    """Persist an accession -> NCBI-datasets-lineage map as a 2-column sidecar in
    ``dest`` (``accession<TAB>lineage``). Returns the path, or ``None`` when there is
    nothing to write. Lets a later typing pass tag these accessions ``ncbi-datasets``
    after the source ``data_report.jsonl`` has been discarded."""
    if not labels:
        return None
    out = dest / NCBI_LINEAGES_TSV
    with open(out, "w") as fo:
        for acc, lineage in sorted(labels.items()):
            fo.write(f"{acc}\t{lineage}\n")
    return out


def fetch_ncbi_virus(
    taxon: str,
    dest: Path,
    *,
    refseq: bool = True,
    complete_only: bool = False,
    released_after: str | None = None,
    limit: int | None = None,
    logger: logging.Logger,
) -> list[Path]:
    """Download a taxon-scoped genome set from NCBI Virus and split it into ``dest``.

    By default fetches the RefSeq/representative set (small and diverse). Use
    ``complete_only`` (optionally with ``released_after``) for a broader set, then
    dereplicate. ``limit`` caps how many genomes are kept after download (the
    ``datasets`` CLI cannot count-limit a virus-genome download), bounding the
    downstream split/align work; a heavily sequenced taxon still downloads in full,
    so a curated ``--candidate-pool`` is preferable there.
    """
    if not datasets_available():
        raise UserInputError(
            "Fetching from NCBI Virus needs the datasets CLI. Install with: "
            "conda install -c conda-forge ncbi-datasets-cli"
        )
    dest.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        zip_path = Path(tmp) / "package.zip"
        command = ["datasets", "download", "virus", "genome", "taxon", taxon,
                   "--no-progressbar", "--filename", str(zip_path)]
        command += _scope_flags(refseq, complete_only, released_after)
        scope = "RefSeq" if refseq else ("complete" if complete_only else "all")
        logger.info("Fetching %s '%s' genomes from NCBI Virus...", scope, taxon)
        run_tool(DATASETS, command, logger=logger, log_prefix="datasets")
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tmp)
        data_dir = Path(tmp) / "ncbi_dataset" / "data"
        fna = data_dir / "genomic.fna"
        if not fna.exists():
            raise UserInputError(
                f"NCBI Virus returned no genomes for '{taxon}'. Try --taxon or a broader scope."
            )
        labels = _lineage_labels(data_dir / "data_report.jsonl")
        write_ncbi_lineages(dest, labels)
        written = _split_fasta(fna, dest, labels, limit=limit)
    logger.info("Fetched %d genome(s) from NCBI Virus into %s.", len(written), dest)
    return written


def _lineage_labels(report: Path) -> dict[str, str]:
    """Map accession -> lineage from a datasets data_report.jsonl (best-effort)."""
    labels: dict[str, str] = {}
    if not report.exists():
        return labels
    for line in report.read_text().splitlines():
        if not line.strip():
            continue
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        acc = rec.get("accession")
        lineage = (rec.get("isolate") or {}).get("lineage") or (rec.get("virus") or {}).get(
            "pangolinClassification"
        )
        if acc and lineage:
            labels[acc] = str(lineage)
    return labels


def _split_fasta(
    fna: Path, dest: Path, labels: dict[str, str], limit: int | None = None
) -> list[Path]:
    """Split a multi-FASTA into one file per accession under ``dest``.

    With ``limit`` set, at most that many records are written (a downstream cap for a
    heavily sequenced taxon; dereplication then keeps the panel diverse).
    """
    written: list[Path] = []
    for header, seq in read_fasta(fna):
        if limit is not None and len(written) >= limit:
            break
        acc = header.split()[0]
        out = dest / f"{acc}.fasta"
        with open(out, "w") as fo:
            note = f" {labels[acc]}" if acc in labels else ""
            write_fasta_record(fo, f"{acc}{note}", seq)
        written.append(out)
    return written


def _lineage_select(
    genomes: list[Path],
    lineage_of: dict[str, str],
    query_ani: dict[Path, tuple[float, float]],
    *,
    keep_recombinant: bool,
    derep_ani: float,
    logger: logging.Logger,
) -> list[Path]:
    """Reduce a pool by lineage: keep one query-closest representative per lineage.

    Genomes with a lineage label are grouped by lineage and the member with the
    highest whole-genome ANI to the query is kept; recombinant (CRF/URF/X) lineages
    are dropped unless ``keep_recombinant`` is set, since they carry both parents'
    segments and mask the true parents. Genomes with no label fall back to ANI
    dereplication (the pre-lineage behaviour). Returns the union.
    """
    typed_by_lineage: dict[str, list[Path]] = {}
    untyped: list[Path] = []
    for g in genomes:
        lineage = lineage_of.get(strip_sequence_extension(g.name))
        if not lineage:
            untyped.append(g)
        elif keep_recombinant or not is_recombinant_lineage(lineage):
            typed_by_lineage.setdefault(lineage, []).append(g)
    reps = [
        max(members, key=lambda g: query_ani.get(g, (0.0, 0.0))[0])
        for members in typed_by_lineage.values()
    ]
    if len(untyped) > 2:
        kept_untyped, _ = dereplicate(untyped, ani=derep_ani, logger=logger)
    else:
        kept_untyped = untyped
    logger.info(
        "Lineage selection: %d representative(s) from %d lineage(s); "
        "%d untyped genome(s) -> %d after ANI dereplication.",
        len(reps), len(typed_by_lineage), len(untyped), len(kept_untyped),
    )
    return reps + kept_untyped


def select_regional(
    query_fasta: Path,
    genomes: list[Path],
    *,
    window: int = DEFAULT_SEED_WINDOW,
    per_window: int = DEFAULT_PER_WINDOW,
    ani_floor: float = DEFAULT_ANI_FLOOR,
    derep_ani: float = DEFAULT_DEREP_ANI,
    dereplicate_pool: bool = True,
    lineage_of: dict[str, str] | None = None,
    keep_recombinant: bool = False,
    drop_siblings: bool = True,
    sibling_identity: float = DEFAULT_SIBLING_IDENTITY,
    sibling_coverage: float = DEFAULT_SIBLING_COVERAGE,
    logger: logging.Logger,
) -> PoolSelection:
    """Select a diverse, query-relevant panel from a genome pool by regional matching.

    Dereplicates the pool, drops the query's genome-wide siblings (so the regional
    parents are not out-ranked everywhere) and genomes unrelated to the query, then
    cuts the query into windows and keeps each window's best-matching genomes. The
    union is a focused panel covering the query region by region.

    When ``lineage_of`` labels the genomes, the pool is reduced by lineage -- one
    query-closest representative per lineage, recombinant lineages excluded unless
    ``keep_recombinant`` -- instead of by ANI; untyped genomes still dereplicate.
    """
    if not skani_available():
        raise UserInputError(
            "Pool selection needs skani. Install with: conda install -c bioconda skani skder"
        )
    if not genomes:
        raise UserInputError("The candidate pool is empty.")
    records = read_fasta(query_fasta)
    if not records:
        raise UserInputError(f"Query FASTA {query_fasta} has no sequence.")

    # Whole-genome query ANI over the full pool, computed once: it drives the reduction
    # (a lineage's representative is its query-closest member) and the sibling / related
    # filters below. A superset of any reduced pool, so every later ``whole[g]`` resolves.
    whole = skani_query_ani(query_fasta, genomes, logger)
    if lineage_of:
        # Typed references: reduce by lineage, not by clade-blind ANI -- one
        # query-closest representative per lineage, recombinant lineages excluded.
        pool = _lineage_select(
            genomes, lineage_of, whole,
            keep_recombinant=keep_recombinant, derep_ani=derep_ani, logger=logger,
        )
    elif dereplicate_pool and len(genomes) > 2:
        pool, _ = dereplicate(genomes, ani=derep_ani, logger=logger)
        logger.info(
            "Dereplicated pool: %d representative(s) from %d genome(s).",
            len(pool), len(genomes),
        )
    else:
        pool = genomes
    # Drop the query's siblings -- near-identical genome-wide over near-full coverage --
    # so they cannot win every region and mask the recombination.
    if drop_siblings:
        siblings = [
            g for g in pool
            if whole[g][0] >= sibling_identity and whole[g][1] >= sibling_coverage
        ]
        if siblings and len(siblings) < len(pool):
            logger.info(
                "Dropping %d sibling(s) of the query (>= %.0f%% ANI over near-full coverage).",
                len(siblings), sibling_identity,
            )
            pool = [g for g in pool if g not in siblings]
    # Drop genomes unrelated to the query (keep everything if the floor removes all).
    related = [g for g in pool if whole[g][0] >= ani_floor] or pool

    query_seq = records[0][1].replace("-", "")
    n_windows = max(1, len(query_seq) // window)
    selected: list[Path] = []
    table: list[dict] = []
    seen = set()
    with tempfile.TemporaryDirectory() as tmp:
        for w in range(n_windows):
            start = w * window
            end = len(query_seq) if w == n_windows - 1 else start + window
            chunk = query_seq[start:end]
            if len(chunk) < 50:
                continue
            win = Path(tmp) / f"win_{w}.fasta"
            with open(win, "w") as fo:
                write_fasta_record(fo, f"win_{w}", chunk)
            ani = skani_query_ani(win, related, logger)
            ranked = sorted(
                (g for g in related if ani[g][0] > 0), key=lambda g: ani[g][0], reverse=True
            )
            for g in ranked[:per_window]:
                if g.resolve() not in seen:
                    seen.add(g.resolve())
                    selected.append(g)
                table.append({
                    "window": w, "query_start": start, "query_end": end,
                    "genome": strip_sequence_extension(g.name), "ani": ani[g][0],
                })
    logger.info(
        "Regional selection: %d genome(s) cover the query across %d window(s).",
        len(selected), n_windows,
    )
    return PoolSelection(selected=selected, table=table)
