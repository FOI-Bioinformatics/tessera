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
from ..core.binaries import BinarySpec, check_binaries
from ..core.executors import parallel_map
from ..core.io import normalize_reference, read_fasta, write_fasta_record
from ..core.plugins import ToolCapabilities
from ..core.process import run_tool
from .base import (
    Aligner,
    AlignParams,
    AlignResult,
    split_thread_budget,
    warn_unplaced,
)

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

# MAFFT 7.525 segfaults inside `pairlocalalign` whenever `--addfragments` is given more
# than one record and the pairwise-local path is taken. Measured against 7.525 here: it
# crashes on two 300 bp fragments as readily as on three 2 kb ones, with and without
# `--adjustdirection` and `--keeplength`, so neither the strand handling nor the length
# rule above has anything to do with it; a single-record file is fine, and 7.526 aligns
# every one of those inputs cleanly. Since SIXMER_MIN_BACKBONE only engages above 20 kb,
# an unpatched 7.525 could not align *any* draft assembly of a smaller genome.
#
# `--6merpair` derives its guide distances from k-mer counts and never enters the
# crashing routine, so a multi-record genome added by a MAFFT older than this takes that
# path. It is the same escape the memory rule above already uses, which is what makes it
# an acceptable substitution rather than a new mode: on anything over 20 kb -- most real
# genomes -- Tessera is already asking for it.
MULTIFRAGMENT_PAIRWISE_FIXED_IN = (7, 526)

# Upper estimate of the memory one ``--multipair`` job needs, per cell (backbone length
# times added-sequence length).
# Measured peaks were 7-20 bytes (8.5-30 kb genome pairs, MAFFT 7.526).
_MULTIPAIR_BYTES_PER_CELL = 24
# Share of physical memory the concurrent ``--multipair`` jobs may be estimated to use.
_MEMORY_SHARE = 0.5
_FALLBACK_MEMORY = 8 * 2**30


def _mafft_at_least(version: str, minimum: tuple[int, ...]) -> bool | None:
    """Is the reported MAFFT ``version`` at least ``minimum``? ``None`` if unreadable."""
    parts: list[int] = []
    for field_text in version.split(".")[: len(minimum)]:
        digits = "".join(c for c in field_text if c.isdigit())
        if not digits:
            return None
        parts.append(int(digits))
    return tuple(parts) >= minimum if len(parts) == len(minimum) else None


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
        ref_fasta.write_text(f">{reference.stem}\n{ref_seq}\n", encoding="utf-8")

        tuning: list[str] = []
        for key, flag in _VALUE_TUNING.items():
            if key in params.extra:
                tuning += [flag, str(params.extra[key])]
        queries = [g for g in genomes if g != reference]
        # Pairwise memory grows with backbone length times the length of the sequence
        # added, so the rule looks at the longest sequence on either side. The record
        # count is kept from the same read, for the 7.525 workaround below.
        fragments: dict[Path, int] = {}
        longest_added = 0
        for genome in queries:
            records = read_fasta(genome)
            fragments[genome] = len(records)
            longest_added = max(
                longest_added, max((len(seq) for _, seq in records), default=0)
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

        # A MAFFT too old to add several fragments on the pairwise path (see
        # MULTIFRAGMENT_PAIRWISE_FIXED_IN) would crash on a draft assembly, so those
        # genomes take the k-mer path instead. An unreadable version is left alone: a
        # clear failure from the tool is better than quietly coarsening every alignment
        # on the strength of a version string nobody could parse.
        multifragment_sixmer: set[Path] = set()
        if not sixmer and any(count > 1 for count in fragments.values()):
            reported = check_binaries(self.capabilities.required_binaries).get("mafft", "")
            if _mafft_at_least(reported, MULTIFRAGMENT_PAIRWISE_FIXED_IN) is False:
                multifragment_sixmer = {g for g, count in fragments.items() if count > 1}
                logger.warning(
                    "MAFFT %s segfaults when it adds more than one fragment on the "
                    "pairwise path (fixed in %s); using k-mer distances (--6merpair) "
                    "for %d multi-contig genome(s). Upgrade MAFFT to keep the pairwise "
                    "alignment for them.",
                    reported, ".".join(map(str, MULTIFRAGMENT_PAIRWISE_FIXED_IN)),
                    len(multifragment_sixmer),
                )

        # The thread budget is shared: concurrent jobs each get a slice of it rather than
        # the whole. Pairwise-alignment jobs are also limited by memory.
        jobs, per_job = split_thread_budget(params.threads, len(queries))
        if not sixmer:
            limit = _multipair_job_limit(len(ref_seq), max(longest_added, 1))
            if limit < jobs:
                logger.info(
                    "Running %d MAFFT job(s) at a time instead of %d to bound memory.",
                    limit, jobs,
                )
                jobs, per_job = split_thread_budget(params.threads, len(queries), limit)
        threads = str(per_job)

        def add_genome(genome: Path) -> tuple[str, str]:
            aligned = out_dir / f"{genome.stem}.aln.fasta"
            per_genome = [*tuning]
            if genome in multifragment_sixmer:
                per_genome.append("--6merpair")
            run_tool(
                self.capabilities,
                ["mafft", "--thread", threads, "--keeplength", "--adjustdirection",
                 *per_genome, "--addfragments", str(genome.resolve()),
                 str(ref_fasta.resolve())],
                logger=logger,
                log_prefix=f"mafft:{genome.stem}",
                stdout_path=aligned,
            )
            _, merged = merge_added_fragments(aligned)
            return genome.stem, merged

        rows = parallel_map(add_genome, queries, jobs, logger=logger)

        warn_unplaced(rows, reference.stem, logger)
        msa = out_dir / "msa.fasta"
        with open(msa, "w", encoding="utf-8") as out:
            write_fasta_record(out, reference.stem, ref_seq)
            for stem, row in rows:
                write_fasta_record(out, stem, row)
        return AlignResult(msa_fasta=msa)
