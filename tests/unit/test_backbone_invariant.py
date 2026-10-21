"""The backbone row of every backend's alignment is the backbone's own sequence.

The backbone is both the coordinate anchor and a donor candidate, so a backbone row
that lost columns is a wrong alignment that still looks well-formed: full width,
rectangular, no ragged-row guard fires. Each backend is driven by a fake of its tool
that aligns two queries to opposite halves of a two-contig backbone, so no single
query covers the whole backbone. The fakes cover every backbone column, so for the
MAF-based backends (sibeliaz, cactus) equality is the right check here too.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from tessera.aligners import cactus as cactus_mod
from tessera.aligners import mafft as mafft_mod
from tessera.aligners import minimap2 as mm_mod
from tessera.aligners import progressivemauve as pm_mod
from tessera.aligners import sibeliaz as sz_mod
from tessera.aligners.base import AlignParams

_LOG = logging.getLogger("tessera.test")

# Two backbone contigs of 10 bp. Query A aligns only to c1, query B only to c2.
C1 = "ACGTACGTAA"
C2 = "TTGCATGCCG"
BACKBONE = C1 + C2


def _read_msa(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            out[name] = ""
        elif name is not None:
            out[name] += line
    return out


def _genomes(tmp_path: Path) -> list[Path]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    ref = tmp_path / "ref.fasta"
    ref.write_text(f">c1\n{C1}\n>c2\n{C2}\n")
    qa = tmp_path / "qryA.fasta"
    qa.write_text(f">a1\n{C1}\n")
    qb = tmp_path / "qryB.fasta"
    qb.write_text(f">b1\n{C2}\n")
    return [ref, qa, qb]


def _run_progressivemauve(monkeypatch, genomes: list[Path], out: Path) -> Path:
    def fake(caps, cmd, **kw):
        cmd = [str(c) for c in cmd]
        xmfa = Path(cmd[cmd.index("--output") + 1])
        ref_path, query_path = cmd[-2], cmd[-1]
        # progressiveMauve reads a multi-record FASTA as one concatenated sequence.
        lo, hi = (1, 10) if "qryA" in query_path else (11, 20)
        bases = BACKBONE[lo - 1:hi]
        xmfa.write_text(
            f"#Sequence1File\t{ref_path}\n#Sequence2File\t{query_path}\n"
            f"> 1:{lo}-{hi} + {ref_path}\n{bases}\n> 2:1-10 + {query_path}\n{bases}\n=\n"
        )
        return ""

    monkeypatch.setattr(pm_mod, "run_tool", fake)
    return pm_mod.ProgressiveMauveAligner().align(
        genomes, genomes[0], out, AlignParams(threads=1), _LOG
    ).msa_fasta


def _run_mafft(monkeypatch, genomes: list[Path], out: Path) -> Path:
    def fake(caps, cmd, **kw):
        added = str(cmd[cmd.index("--addfragments") + 1])
        row = C1 + "-" * 10 if "qryA" in added else "-" * 10 + C2
        kw["stdout_path"].write_text(f">ref\n{BACKBONE}\n>contig\n{row}\n")
        return ""

    monkeypatch.setattr(mafft_mod, "run_tool", fake)
    return mafft_mod.MafftAligner().align(
        genomes, genomes[0], out, AlignParams(threads=1), _LOG
    ).msa_fasta


def _run_minimap2(monkeypatch, genomes: list[Path], out: Path) -> Path:
    def fake(caps, cmd, **kw):
        query = str(cmd[-1])
        rname, seq = ("c1", C1) if "qryA" in query else ("c2", C2)
        kw["stdout_path"].write_text(
            "@SQ\tSN:c1\tLN:10\n@SQ\tSN:c2\tLN:10\n"
            f"q\t0\t{rname}\t1\t60\t10M\t*\t0\t0\t{seq}\t*\n"
        )
        return ""

    monkeypatch.setattr(mm_mod, "run_tool", fake)
    return mm_mod.Minimap2Aligner().align(
        genomes, genomes[0], out, AlignParams(threads=1), _LOG
    ).msa_fasta


def _run_sibeliaz(monkeypatch, genomes: list[Path], out: Path) -> Path:
    def fake(caps, cmd, **kw):
        out_dir = Path(cmd[cmd.index("-o") + 1])
        (out_dir / "alignment.maf").write_text(
            f"a\ns c1 0 10 + 10 {C1}\ns a1 0 10 + 10 {C1}\n\n"
            f"a\ns c2 0 10 + 10 {C2}\ns b1 0 10 + 10 {C2}\n"
        )
        return ""

    monkeypatch.setattr(sz_mod, "run_tool", fake)
    monkeypatch.setattr(sz_mod, "_sibeliaz_invocation", lambda out_dir, logger: ["sibeliaz"])
    return sz_mod.SibeliazAligner().align(
        genomes, genomes[0], out, AlignParams(threads=1), _LOG
    ).msa_fasta


def _run_cactus(monkeypatch, genomes: list[Path], out: Path) -> Path:
    def fake(caps, cmd, **kw):
        results = Path(cmd[cmd.index("--outDir") + 1])
        results.mkdir(parents=True, exist_ok=True)
        (results / "pangenome.full.hal").write_text("")
        return ""

    def fake_hal_to_maf(hal, ref_name, maf, logger, caps=None):
        # hal2maf names rows <sample>.<contig>.
        Path(maf).write_text(
            f"a\ns ref.c1 0 10 + 10 {C1}\ns qryA.a1 0 10 + 10 {C1}\n\n"
            f"a\ns ref.c2 0 10 + 10 {C2}\ns qryB.b1 0 10 + 10 {C2}\n"
        )
        return Path(maf)

    monkeypatch.setattr(cactus_mod, "run_tool", fake)
    monkeypatch.setattr(cactus_mod, "hal_to_maf", fake_hal_to_maf)
    return cactus_mod.CactusAligner().align(
        genomes, genomes[0], out, AlignParams(threads=1), _LOG
    ).msa_fasta


_BACKENDS = {
    "progressivemauve": _run_progressivemauve,
    "mafft": _run_mafft,
    "minimap2": _run_minimap2,
    "sibeliaz": _run_sibeliaz,
    "cactus": _run_cactus,
}


@pytest.mark.parametrize("backend", sorted(_BACKENDS))
@pytest.mark.parametrize("order", ["AB", "BA"])
def test_backbone_row_is_the_backbone_sequence(
    monkeypatch, tmp_path: Path, backend: str, order: str
) -> None:
    ref, qa, qb = _genomes(tmp_path / "in")
    genomes = [ref, qa, qb] if order == "AB" else [ref, qb, qa]
    msa = _read_msa(_BACKENDS[backend](monkeypatch, genomes, tmp_path / "out"))
    assert msa["ref"] == BACKBONE
    assert msa["qryA"] == C1 + "-" * 10
    assert msa["qryB"] == "-" * 10 + C2


def test_progressivemauve_rejects_a_projection_wider_than_the_backbone(
    monkeypatch, tmp_path: Path
) -> None:
    # The backbone row is now the reference itself, so a projection that reaches past
    # its end would make the rows ragged. It is reported against the file instead.
    from tessera.core.errors import OutputError

    def fake(caps, cmd, **kw):
        cmd = [str(c) for c in cmd]
        xmfa = Path(cmd[cmd.index("--output") + 1])
        ref_path, query_path = cmd[-2], cmd[-1]
        bases = BACKBONE + "ACGT"
        xmfa.write_text(
            f"#Sequence1File\t{ref_path}\n#Sequence2File\t{query_path}\n"
            f"> 1:1-24 + {ref_path}\n{bases}\n> 2:1-24 + {query_path}\n{bases}\n=\n"
        )
        return ""

    monkeypatch.setattr(pm_mod, "run_tool", fake)
    genomes = _genomes(tmp_path / "in")
    with pytest.raises(OutputError, match=r"qryA.*24 column.*20 bp"):
        pm_mod.ProgressiveMauveAligner().align(
            genomes, genomes[0], tmp_path / "out", AlignParams(threads=1), _LOG
        )


def test_cactus_refuses_labels_that_collide_once_dots_are_replaced(
    monkeypatch, tmp_path: Path
) -> None:
    # Cactus sample names may not contain '.', so the adapter maps it to '_'. Two
    # accessions differing only there would then share one sample name and one row.
    from tessera.core.errors import UserInputError

    def must_not_run(caps, cmd, **kw):
        raise AssertionError("cactus-pangenome must not be started")

    monkeypatch.setattr(cactus_mod, "run_tool", must_not_run)
    genomes = []
    for label in ("ref", "MN908947.3", "MN908947_3"):
        path = tmp_path / f"{label}.fasta"
        path.write_text(f">{label}\nACGT\n")
        genomes.append(path)
    with pytest.raises(UserInputError, match=r"MN908947\.3.*MN908947_3.*MN908947_3"):
        cactus_mod.CactusAligner().align(
            genomes, genomes[0], tmp_path / "out", AlignParams(threads=1), _LOG
        )
