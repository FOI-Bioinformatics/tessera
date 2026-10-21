"""SibeliaZ aligner.

SibeliaZ builds locally collinear blocks and emits a MAF for moderately
divergent genomes (up to ~0.09 substitutions/site). The MAF is projected onto
reference coordinates by :mod:`tessera.converters.maf_to_fasta`.
"""

from __future__ import annotations

import logging
import shutil
import sys
from collections.abc import Sequence
from pathlib import Path

from ..converters.maf_to_fasta import maf_to_fasta
from ..core.binaries import BinarySpec
from ..core.errors import OutputError, UserInputError
from ..core.io import normalize_reference, read_fasta
from ..core.plugins import ToolCapabilities
from ..core.process import run_tool
from .base import Aligner, AlignParams, AlignResult

# SibeliaZ's bash wrapper uses GNU/Linux-only constructs in its alignment step
# (`free`, GNU `find -printf`, `stat -c`, `mktemp --suffix`). On macOS/BSD these
# fail silently, so blocks are found but no MAF is written. We run a patched copy
# of the wrapper on those platforms.
_BSD_PATCHES: tuple[tuple[str, str], ...] = (
    (
        "free -g -w | head -2 | tail -1 | awk '{print $2}'",
        "echo $(( $(sysctl -n hw.memsize) / 1073741824 ))",
    ),
    (
        "free -k -w | head -2 | tail -1 | awk '{print $2}'",
        "echo $(( $(sysctl -n hw.memsize) / 1024 ))",
    ),
    ('find $outdir -name "*.tmp" -printf "%p\\n"', 'find $outdir -name "*.tmp"'),
    ('stat -c "%s" $i', "stat -f%z $i"),
    # BSD mktemp has no --suffix; append .fa after the call (spoa needs a known
    # FASTA extension). Patch the whole $(...) so the result still ends in .fa.
    ("$(mktemp --suffix=.fa $outdir/block.XXXXX)", "$(mktemp $outdir/block.XXXXX).fa"),
    ("ulimit $memory_min", ":"),
)


class SibeliazAligner(Aligner):
    capabilities = ToolCapabilities(
        name="sibeliaz",
        conda=("bioconda::sibeliaz",),
        required_binaries=(BinarySpec("sibeliaz", version_args=("-v",)),),
        recommended_max_genomes=2000,
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
        genomes, reference = normalize_reference(genomes, reference, tool="SibeliaZ")
        out_dir.mkdir(parents=True, exist_ok=True)

        # Default mode runs the block alignment step, producing alignment.maf.
        # (The -n flag would emit only block coordinates, no MAF.)
        sibeliaz_cmd = _sibeliaz_invocation(out_dir, logger)
        # Pass -f (twopaco bloom-filter memory, GB) explicitly so the wrapper
        # skips its system-memory probe, which uses Linux-only `free`/`stat`.
        filtermemory = params.extra.get("filtermemory", 8)
        # Divergence/sensitivity tuning: a smaller -k (k-mer) finds more anchors
        # across divergent genomes; -a (vertex abundance) and -b (max bubble size)
        # trade sensitivity for speed/memory. Unset keys keep SibeliaZ's defaults.
        tuning: list[str] = []
        for key, flag in (("kmer", "-k"), ("abundance", "-a"), ("bubble", "-b")):
            if key in params.extra:
                tuning += [flag, str(params.extra[key])]
        run_tool(
            self.capabilities,
            [
                *sibeliaz_cmd,
                *tuning,
                "-f", str(filtermemory),
                "-t", str(params.threads),
                "-o", out_dir,
                *[str(g.resolve()) for g in genomes],
            ],
            logger=logger,
            log_prefix="sibeliaz",
        )
        maf = out_dir / "alignment.maf"
        if not maf.exists():
            candidates = sorted(out_dir.rglob("*.maf"))
            if not candidates:
                raise OutputError("SibeliaZ did not produce a MAF file")
            maf = candidates[0]

        # SibeliaZ may find blocks yet write a MAF with no alignment rows: its
        # wrapper runs spoa per block with stderr suppressed, so an OOM-killed
        # spoa is silent and yields an empty MAF.
        if not any(line.startswith("s") for line in maf.read_text().splitlines()):
            raise OutputError(
                f"SibeliaZ wrote an empty MAF (no alignment rows) at {maf}. Its "
                "blocks were found but the per-block spoa alignment produced "
                "nothing -- most often spoa was OOM-killed on a large collinear "
                "block. Give the machine more memory or use the progressivemauve "
                "backend."
            )

        # SibeliaZ's MAF uses sequence/contig IDs (FASTA header first token), not
        # genome filenames; build the seqid -> genome-stem map for the converter.
        name_map = _build_seqid_map(genomes)
        msa = out_dir / "msa.fasta"
        # The MAF names only the backbone contigs that fall in a block, in no particular
        # order, and only the genomes it placed. Hand the converter the backbone's own
        # contig order and the full genome list so neither is inferred from the MAF.
        maf_to_fasta(
            maf, reference.stem, msa, name_map=name_map,
            ref_contigs=[(seqid, len(seq)) for seqid, seq in read_fasta(reference)],
            expected=[g.stem for g in genomes],
            logger=logger,
        )
        return AlignResult(msa_fasta=msa, native_format=maf)


def _sibeliaz_invocation(out_dir: Path, logger: logging.Logger) -> list[str]:
    """Return the command prefix to run SibeliaZ.

    On non-macOS just ``["sibeliaz"]``. On macOS/BSD, write a compatible copy of
    the wrapper (its alignment step otherwise fails silently) and run it via bash.
    """
    if sys.platform != "darwin":
        return ["sibeliaz"]

    real = shutil.which("sibeliaz")
    if real is None:
        return ["sibeliaz"]
    script = Path(real).read_text()
    if "-printf" not in script:  # already BSD-friendly / unexpected layout
        return ["sibeliaz"]

    patched = script
    for old, new in _BSD_PATCHES:
        if old in patched:
            patched = patched.replace(old, new)
    patched_path = out_dir / "sibeliaz_bsd.sh"
    patched_path.write_text(patched)
    logger.info("Using BSD-compatible SibeliaZ wrapper on macOS: %s", patched_path)
    return ["bash", str(patched_path)]


def _build_seqid_map(genomes) -> dict[str, str]:
    """Map each FASTA sequence ID (header first token) to its genome stem."""
    name_map: dict[str, str] = {}
    for genome in genomes:
        stem = genome.stem
        with open(genome) as fo:
            for lineno, line in enumerate(fo, start=1):
                if line.startswith(">"):
                    tokens = line[1:].split()
                    if not tokens:
                        raise UserInputError(
                            f"{genome} has a record with no sequence ID (line {lineno}). "
                            "The sibeliaz backend identifies genomes by sequence ID, so "
                            "every record needs a name after '>'."
                        )
                    seqid = tokens[0]
                    owner = name_map.setdefault(seqid, stem)
                    if owner != stem:
                        # SibeliaZ names alignment rows by sequence ID alone, so two
                        # genomes sharing one would be merged under a single label
                        # and the other would vanish from the MSA without a word.
                        raise UserInputError(
                            f"Sequence ID '{seqid}' occurs in both '{owner}' and "
                            f"'{stem}'. The sibeliaz backend identifies genomes by "
                            "sequence ID, so IDs must be unique across the query and "
                            "collection -- rename one, or use another --aligner."
                        )
    return name_map
