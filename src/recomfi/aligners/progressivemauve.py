"""progressiveMauve aligner.

Reproduces the legacy ``phylo.py`` accurate path: align every genome to a single
reference with progressiveMauve, project each pairwise XMFA onto reference
coordinates (:mod:`recomfi.converters.xmfa_to_fasta`), then concatenate into one
reference-anchored MSA-FASTA. Aligning to a single reference keeps the work
linear in the number of genomes and yields the reference-anchored coordinates
the recombination scan assumes. It tolerates large rearrangements but is slow
and heavy; ``sibeliaz`` is RecomFi's default backend.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..converters.xmfa_to_fasta import xmfa_to_fasta
from ..core.binaries import BinarySpec
from ..core.errors import UserInputError
from ..core.executors import parallel_map
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
        genomes = list(genomes)
        if reference is None:
            reference = genomes[0]
        if reference not in genomes:
            genomes = [reference, *genomes]
        if len(genomes) < 3:
            raise UserInputError(
                "progressiveMauve alignment needs at least 3 genomes "
                "(query plus a reference collection)."
            )

        out_dir.mkdir(parents=True, exist_ok=True)
        xmfa_dir = out_dir / "xmfa"
        xmfa_dir.mkdir(exist_ok=True)
        ref_arg = str(reference.resolve())

        # Pin every per-query projection to the full reference length so the
        # rows concatenate into a rectangular MSA regardless of how far each
        # query aligns.
        ref_length = sum(len(seq) for _, seq in _read_fasta(reference))

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

        def align_query(query: Path) -> Path:
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
            return fa

        per_query_fastas = parallel_map(align_query, queries, workers, logger=logger)

        msa = out_dir / "msa.fasta"
        _concatenate(per_query_fastas, reference, msa)
        return AlignResult(msa_fasta=msa)


def _concatenate(per_query_fastas: list[Path], reference: Path, out_path: Path) -> None:
    """Write the reference row once, then each query row; leaf names are stems."""
    ref_stem = reference.stem
    written_ref = False
    with open(out_path, "w") as out:
        for fa in per_query_fastas:
            for name, seq in _read_fasta(fa):
                leaf = Path(name).stem
                if leaf == ref_stem:
                    if written_ref:
                        continue
                    written_ref = True
                out.write(f">{leaf}\n")
                for pos in range(0, len(seq), 80):
                    out.write(seq[pos : pos + 80] + "\n")


def _read_fasta(path: Path):
    name = None
    seq: list[str] = []
    with open(path) as fo:
        for line in fo:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)
