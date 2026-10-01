"""progressiveMauve aligner.

Reproduces the legacy ``phylo.py`` accurate path: align every genome to a single
reference with progressiveMauve, project each pairwise XMFA onto reference
coordinates (:mod:`tessera.converters.xmfa_to_fasta`), then concatenate into one
reference-anchored MSA-FASTA. Aligning to a single reference keeps the work
linear in the number of genomes and yields the reference-anchored coordinates
the recombination scan assumes. It tolerates large rearrangements but is slow
and heavy; ``sibeliaz`` is Tessera's default backend.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..converters.xmfa_to_fasta import xmfa_to_fasta
from ..core.binaries import BinarySpec
from ..core.errors import OutputError
from ..core.executors import parallel_map
from ..core.io import normalize_reference, read_fasta, write_fasta_record
from ..core.plugins import ToolCapabilities
from ..core.process import run_tool
from .base import Aligner, AlignParams, AlignResult


class ProgressiveMauveAligner(Aligner):
    capabilities = ToolCapabilities(
        name="progressivemauve",
        conda=("bioconda::mauve", "conda-forge::boost-cpp=1.74.0"),
        required_binaries=(BinarySpec("progressiveMauve", version_args=()),),
        recommended_max_genomes=500,
        threads_param=None,  # progressiveMauve is single-threaded per alignment
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
            genomes, reference, tool="progressiveMauve", min_genomes=3, ensure_member=True
        )
        out_dir.mkdir(parents=True, exist_ok=True)
        xmfa_dir = out_dir / "xmfa"
        xmfa_dir.mkdir(exist_ok=True)
        ref_arg = str(reference.resolve())

        # Pin every per-query projection to the full reference length so the
        # rows concatenate into a rectangular MSA regardless of how far each
        # query aligns.
        ref_length = sum(len(seq) for _, seq in read_fasta(reference))

        queries = [g for g in genomes if g != reference]

        # A lower --seed-weight raises sensitivity for divergent genomes (more,
        # shorter anchor seeds); unset keeps progressiveMauve's default.
        seed_opt: list[str] = []
        if "seed_weight" in params.extra:
            seed_opt = ["--seed-weight", str(params.extra["seed_weight"])]

        # progressiveMauve is single-threaded per alignment; run one process per
        # thread budget, each aligning an independent query to the reference.
        # Setting threads to 1 (or extra single=true) serialises them, which has
        # resolved a yet-unexplained progressiveMauve error on some systems.
        workers = 1 if params.flag("single") else params.threads

        def align_query(query: Path) -> tuple[str, Path]:
            stem = query.stem
            xmfa = xmfa_dir / f"{stem}.xmfa"
            fa = xmfa_dir / f"{stem}.fa"
            run_tool(
                self.capabilities,
                ["progressiveMauve", "--output", xmfa, *seed_opt, ref_arg, str(query.resolve())],
                logger=logger,
                log_prefix=f"progressivemauve:{stem}",
            )
            xmfa_to_fasta(xmfa, ref_arg, 0, fa, reference_length=ref_length)
            return stem, fa

        per_query = parallel_map(align_query, queries, workers, logger=logger)

        msa = out_dir / "msa.fasta"
        _concatenate(per_query, reference.stem, msa)
        return AlignResult(msa_fasta=msa)


def _concatenate(
    per_query: list[tuple[str, Path]], reference_label: str, out_path: Path
) -> None:
    """Write the reference row once, then one row per query, named by staged label.

    Each per-query FASTA holds exactly two records, in a fixed order: the reference
    projection, then the query's. The names inside those files are the paths
    progressiveMauve echoed from its command line -- resolved symlink targets, which can
    carry an unrecognised extension, whitespace, or a basename shared with another
    genome -- so rows are identified by position and labelled from the staged file,
    never from those names.
    """
    with open(out_path, "w") as out:
        for i, (label, fa) in enumerate(per_query):
            records = read_fasta(fa)
            if len(records) != 2:
                raise OutputError(
                    f"Expected a reference row and one query row in {fa} (the "
                    f"projection of '{label}' onto '{reference_label}'), found "
                    f"{len(records)} record(s). progressiveMauve's output is not a "
                    "pairwise alignment."
                )
            if i == 0:
                write_fasta_record(out, reference_label, records[0][1])
            write_fasta_record(out, label, records[1][1])
