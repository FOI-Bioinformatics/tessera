"""MAFFT aligner.

Produces a genuine base-level alignment, the canonical input for the
sliding-window similarity method, and the best fit for similar, largely
collinear genomes within a genus. To honour Tessera's reference-anchored
contract each genome is added onto the backbone with
``mafft --addfragments <genome> --keeplength <reference>``: ``--keeplength``
keeps the output in reference coordinates (insertions relative to the backbone
are dropped) and ``--addfragments`` is designed for fragmented assemblies, so a
multi-contig query is handled cleanly. ``--adjustdirection`` lets MAFFT reverse-
complement an added sequence that is on the opposite strand to the backbone (it renames
such a record ``_R_<name>``; the merge below is positional and ignores names). A genome's
contigs are then merged into a single reference-anchored row.

Memory. Without a pair option ``--addfragments`` uses ``--multipair``: a pairwise local
alignment of each added sequence against the backbone, whose memory grows with the product
of the two lengths (1.6-6 GB per job was measured on 8.5-30 kb pairs; over 100 GB was
observed for a 200 kb pair, MAFFT 7.526). MAFFT's own switch to k-mer distances for long
input only applies under ``--auto``, so the adapter applies a length rule itself: when the
backbone or any sequence to add is longer than ``SIXMER_MIN_BACKBONE`` it passes
``--6merpair``. Otherwise the number of concurrent jobs is limited by physical memory (the
host's, as the operating system reports it; a container or scheduler limit is not read).
"""

from __future__ import annotations

import logging
import os
from collections.abc import Sequence
from pathlib import Path

from ..converters.mafft_merge import merge_added_fragments
from ..core.binaries import BinarySpec
from ..core.executors import parallel_map
from ..core.io import normalize_reference, read_fasta, write_fasta_record
from ..core.plugins import ToolCapabilities
from ..core.process import run_tool
from .base import Aligner, AlignParams, AlignResult

# extra-arg key -> MAFFT option taking a value.
_VALUE_TUNING = {
    "maxiterate": "--maxiterate",
    "retree": "--retree",
    "op": "--op",
    "ep": "--ep",
}

# Sequence length (nt; backbone or a genome to add) above which ``--6merpair`` replaces
# the pairwise local alignment.
# An explicit ``sixmerpair=true|false`` aligner argument overrides the rule.
SIXMER_MIN_BACKBONE = 20_000

# Upper estimate of the memory one ``--multipair`` job needs, per cell (backbone length
# times added-sequence length).
# Measured peaks were 7-20 bytes (8.5-30 kb genome pairs, MAFFT 7.526).
_MULTIPAIR_BYTES_PER_CELL = 24
# Share of physical memory the concurrent ``--multipair`` jobs may be estimated to use.
_MEMORY_SHARE = 0.5
_FALLBACK_MEMORY = 8 * 2**30


def _physical_memory() -> int:
    """Physical memory in bytes, or a small fallback where it cannot be read."""
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (AttributeError, OSError, ValueError):
        return _FALLBACK_MEMORY


def _multipair_job_limit(backbone_len: int, added_len: int) -> int:
    """How many ``--multipair`` jobs fit the memory share at once (at least one)."""
    per_job = _MULTIPAIR_BYTES_PER_CELL * backbone_len * added_len
    return max(1, int(_MEMORY_SHARE * _physical_memory()) // max(1, per_job))


class MafftAligner(Aligner):
    capabilities = ToolCapabilities(
        name="mafft",
        conda=("bioconda::mafft",),
        required_binaries=(BinarySpec("mafft", version_args=("--version",)),),
        recommended_max_genomes=1000,
        threads_param="--thread",
    )

    def align(
        self,
        genomes: Sequence[Path],
        reference: Path | None,
        out_dir: Path,
        params: AlignParams,
        logger: logging.Logger,
    ) -> AlignResult:
        genomes, reference = normalize_reference(
            genomes, reference, tool="MAFFT", min_genomes=2, ensure_member=True
        )
        out_dir.mkdir(parents=True, exist_ok=True)

        # MAFFT anchors to a single-sequence backbone; concatenate a multi-contig
        # reference so the output stays one reference row.
        ref_seq = "".join(seq for _, seq in read_fasta(reference))
        ref_fasta = out_dir / "reference.fasta"
        ref_fasta.write_text(f">{reference.stem}\n{ref_seq}\n")

        tuning: list[str] = []
        for key, flag in _VALUE_TUNING.items():
            if key in params.extra:
                tuning += [flag, str(params.extra[key])]
        queries = [g for g in genomes if g != reference]
        # Pairwise memory grows with backbone length times the length of the sequence
        # added, so the rule looks at the longest sequence on either side.
        longest_added = max(
            (len(seq) for g in queries for _, seq in read_fasta(g)), default=0
        )
        longest = max(len(ref_seq), longest_added)
        sixmer = params.flag("sixmerpair", default=longest > SIXMER_MIN_BACKBONE)
        if sixmer:
            tuning.append("--6merpair")
            if "sixmerpair" not in params.extra:
                logger.info(
                    "Longest sequence is %d nt (> %d): using MAFFT k-mer distances "
                    "(--6merpair) to bound memory; override with --aligner-arg "
                    "sixmerpair=false.",
                    longest, SIXMER_MIN_BACKBONE,
                )

        # The thread budget is shared: concurrent jobs each get a slice of it rather than
        # the whole. Pairwise-alignment jobs are also limited by memory.
        budget = max(1, params.threads)
        jobs = max(1, min(budget, len(queries)))
        if not sixmer:
            limit = _multipair_job_limit(len(ref_seq), max(longest_added, 1))
            if limit < jobs:
                logger.info(
                    "Running %d MAFFT job(s) at a time instead of %d to bound memory.",
                    limit, jobs,
                )
                jobs = limit
        threads = str(budget // jobs)

        def add_genome(genome: Path) -> tuple[str, str]:
            aligned = out_dir / f"{genome.stem}.aln.fasta"
            run_tool(
                self.capabilities,
                ["mafft", "--thread", threads, "--keeplength", "--adjustdirection",
                 *tuning, "--addfragments", str(genome.resolve()), str(ref_fasta.resolve())],
                logger=logger,
                log_prefix=f"mafft:{genome.stem}",
                stdout_path=aligned,
            )
            _, merged = merge_added_fragments(aligned)
            return genome.stem, merged

        rows = parallel_map(add_genome, queries, jobs, logger=logger)

        msa = out_dir / "msa.fasta"
        with open(msa, "w") as out:
            write_fasta_record(out, reference.stem, ref_seq)
            for stem, row in rows:
                write_fasta_record(out, stem, row)
        return AlignResult(msa_fasta=msa)
