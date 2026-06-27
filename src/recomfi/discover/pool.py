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
    ``complete_only`` (optionally with ``released_after`` / ``limit``) for a broader
    set, then dereplicate; a broad fetch can be large, so ``limit`` caps it.
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
        if limit is not None:
            command += ["--limit", str(limit)]
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
        written = _split_fasta(fna, dest, labels)
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


def _split_fasta(fna: Path, dest: Path, labels: dict[str, str]) -> list[Path]:
    """Split a multi-FASTA into one file per accession under ``dest``."""
    written: list[Path] = []
    for header, seq in read_fasta(fna):
        acc = header.split()[0]
        out = dest / f"{acc}.fasta"
        with open(out, "w") as fo:
            note = f" {labels[acc]}" if acc in labels else ""
            write_fasta_record(fo, f"{acc}{note}", seq)
        written.append(out)
    return written


def select_regional(
    query_fasta: Path,
    genomes: list[Path],
    *,
    window: int = DEFAULT_SEED_WINDOW,
    per_window: int = DEFAULT_PER_WINDOW,
    ani_floor: float = DEFAULT_ANI_FLOOR,
    derep_ani: float = DEFAULT_DEREP_ANI,
    dereplicate_pool: bool = True,
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

    pool = genomes
    if dereplicate_pool and len(genomes) > 2:
        pool, _ = dereplicate(genomes, ani=derep_ani, logger=logger)
        logger.info(
            "Dereplicated pool: %d representative(s) from %d genome(s).", len(pool), len(genomes)
        )

    whole = skani_query_ani(query_fasta, pool, logger)
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
