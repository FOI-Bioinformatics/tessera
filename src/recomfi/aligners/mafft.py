"""MAFFT aligner.

Produces a genuine base-level alignment, the canonical input for the
sliding-window similarity method, and the best fit for similar, largely
collinear genomes within a genus. To honour RecomFi's reference-anchored
contract each genome is added onto the backbone with
``mafft --addfragments <genome> --keeplength <reference>``: ``--keeplength``
keeps the output in reference coordinates (insertions relative to the backbone
are dropped) and ``--addfragments`` is designed for fragmented assemblies, so a
multi-contig query is handled cleanly. A genome's contigs are then merged into a
single reference-anchored row.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..converters.mafft_merge import merge_added_fragments, read_fasta
from ..core.binaries import BinarySpec
from ..core.errors import UserInputError
from ..core.executors import parallel_map
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
        genomes = list(genomes)
        if reference is None:
            reference = genomes[0]
        if reference not in genomes:
            genomes = [reference, *genomes]
        if len(genomes) < 2:
            raise UserInputError("MAFFT alignment needs at least 2 genomes.")

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
        if params.flag("sixmerpair"):
            tuning.append("--6merpair")

        threads = str(max(1, params.threads))
        queries = [g for g in genomes if g != reference]

        def add_genome(genome: Path) -> tuple[str, str]:
            aligned = out_dir / f"{genome.stem}.aln.fasta"
            run_tool(
                self.capabilities,
                ["mafft", "--thread", threads, "--keeplength", *tuning,
                 "--addfragments", str(genome.resolve()), str(ref_fasta.resolve())],
                logger=logger,
                log_prefix=f"mafft:{genome.stem}",
                stdout_path=aligned,
            )
            _, merged = merge_added_fragments(aligned)
            return genome.stem, merged

        rows = parallel_map(add_genome, queries, params.threads, logger=logger)

        msa = out_dir / "msa.fasta"
        with open(msa, "w") as out:
            _write_row(out, reference.stem, ref_seq)
            for stem, row in rows:
                _write_row(out, stem, row)
        return AlignResult(msa_fasta=msa)


def _write_row(out, name: str, seq: str, width: int = 80) -> None:
    out.write(f">{name}\n")
    for pos in range(0, len(seq), width):
        out.write(seq[pos : pos + width] + "\n")
