"""minimap2 aligner.

Aligns every genome (including a fragmented assembly query) to the backbone with
an assembly preset and projects the SAM onto reference coordinates. This is the
fastest backend and the closest match to Tessera's reference-anchored model --
minimap2 is the standard assembly-to-reference aligner, so a multi-contig query
maps cleanly onto the backbone. Best for similar genomes; the ``preset``
extra-arg (default ``asm20``) sets the expected divergence.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..converters.sam_to_fasta import project_sam_to_row
from ..core.binaries import BinarySpec
from ..core.executors import parallel_map
from ..core.io import normalize_reference, read_fasta, write_fasta_record
from ..core.plugins import ToolCapabilities
from ..core.process import run_tool
from .base import Aligner, AlignParams, AlignResult


class Minimap2Aligner(Aligner):
    capabilities = ToolCapabilities(
        name="minimap2",
        conda=("bioconda::minimap2",),
        required_binaries=(BinarySpec("minimap2", version_args=("--version",)),),
        recommended_max_genomes=5000,
        threads_param="-t",
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
            genomes, reference, tool="minimap2", min_genomes=2, ensure_member=True
        )
        out_dir.mkdir(parents=True, exist_ok=True)

        # Concatenate the reference contigs into one row and record each contig's
        # offset, so SAM positions (per reference contig) map to a global column.
        ref_records = read_fasta(reference)
        offsets: dict[str, int] = {}
        cursor = 0
        ref_parts: list[str] = []
        for name, seq in ref_records:
            offsets[name] = cursor
            cursor += len(seq)
            ref_parts.append(seq)
        ref_seq = "".join(ref_parts)
        total_len = cursor

        preset = str(params.extra.get("preset", "asm20"))
        threads = str(max(1, params.threads))
        queries = [g for g in genomes if g != reference]

        def map_genome(genome: Path) -> tuple[str, str]:
            sam = out_dir / f"{genome.stem}.sam"
            run_tool(
                self.capabilities,
                ["minimap2", "-a", "-t", threads, "-x", preset,
                 str(reference.resolve()), str(genome.resolve())],
                logger=logger,
                log_prefix=f"minimap2:{genome.stem}",
                stdout_path=sam,
            )
            return genome.stem, project_sam_to_row(sam, offsets, total_len)

        rows = parallel_map(map_genome, queries, params.threads, logger=logger)

        msa = out_dir / "msa.fasta"
        with open(msa, "w") as out:
            write_fasta_record(out, reference.stem, ref_seq)
            for stem, row in rows:
                write_fasta_record(out, stem, row)
        return AlignResult(msa_fasta=msa)
