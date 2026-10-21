"""MAFFT and minimap2 backends: converters and command construction."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from tessera.aligners import mafft as mafft_mod
from tessera.aligners import minimap2 as mm_mod
from tessera.aligners import progressivemauve as pm_mod
from tessera.aligners.base import AlignParams
from tessera.converters.mafft_merge import merge_added_fragments
from tessera.converters.sam_to_fasta import project_sam_to_row

_LOG = logging.getLogger("tessera.test")


# --- MAFFT keeplength merge ------------------------------------------------
def test_merge_added_fragments_picks_non_gap(tmp_path: Path) -> None:
    # reference first, then two contig rows that each cover part of the reference
    aligned = tmp_path / "a.fasta"
    aligned.write_text(
        ">ref\nACGTACGT\n"
        ">contig1\nAC------\n"
        ">contig2\n----ACGT\n"
    )
    ref_row, merged = merge_added_fragments(aligned)
    assert ref_row == "ACGTACGT"
    assert merged == "AC--ACGT"  # contig1 covers 0-1, contig2 covers 4-7


# --- SAM -> reference projection -------------------------------------------
def _sam(*records: str) -> str:
    header = "@HD\tVN:1.6\n@SQ\tSN:ref\tLN:10\n"
    return header + "".join(r + "\n" for r in records)


def _record(qname, flag, rname, pos, cigar, seq) -> str:
    return f"{qname}\t{flag}\t{rname}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t*"


def test_project_sam_match_insertion_deletion(tmp_path: Path) -> None:
    sam = tmp_path / "a.sam"
    # 2M1I2M: AA at ref0-1, C inserted (dropped), GT at ref2-3
    sam.write_text(_sam(_record("q", 0, "ref", 1, "2M1I2M", "AACGT")))
    row = project_sam_to_row(str(sam), {"ref": 0}, 10)
    assert row == "AAGT------"

    sam.write_text(_sam(_record("q", 0, "ref", 1, "2M1D2M", "AAGT")))
    row = project_sam_to_row(str(sam), {"ref": 0}, 10)
    assert row == "AA-GT-----"  # 1D leaves a gap at ref2


def test_project_sam_skips_secondary_and_unmapped(tmp_path: Path) -> None:
    sam = tmp_path / "a.sam"
    sam.write_text(_sam(
        _record("q", 0x4, "ref", 1, "4M", "AAAA"),     # unmapped: skipped
        _record("q", 0x100, "ref", 1, "4M", "TTTT"),   # secondary: skipped
        _record("q", 0, "ref", 5, "3M", "GGG"),        # primary: placed at 4-6
    ))
    row = project_sam_to_row(str(sam), {"ref": 0}, 10)
    assert row == "----GGG---"


def test_project_sam_multicontig_reference_offsets(tmp_path: Path) -> None:
    sam = tmp_path / "a.sam"
    sam.write_text(
        "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:5\n@SQ\tSN:c2\tLN:5\n"
        + _record("q", 0, "c2", 1, "3M", "TTT") + "\n"
    )
    row = project_sam_to_row(str(sam), {"c1": 0, "c2": 5}, 10)
    assert row == "-----TTT--"  # c2 starts at global offset 5


# --- command construction (monkeypatched run_tool) -------------------------
def _two_genomes(tmp_path: Path) -> list[Path]:
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref\nACGTACGTACGT\n")
    qry = tmp_path / "qry.fasta"
    qry.write_text(">qry\nACGTACGTACGT\n")
    return [ref, qry]


def test_mafft_forwards_keeplength_and_tuning(monkeypatch, tmp_path: Path) -> None:
    captured: dict[str, list[str]] = {}

    def fake_run(caps, cmd, **kw):
        captured["cmd"] = [str(c) for c in cmd]
        raise RuntimeError("stop after capture")

    monkeypatch.setattr(mafft_mod, "run_tool", fake_run)
    ref, qry = _two_genomes(tmp_path)
    params = AlignParams(threads=2, extra={"maxiterate": 1000})
    with pytest.raises(RuntimeError):
        mafft_mod.MafftAligner().align([ref, qry], ref, tmp_path / "out", params, _LOG)

    cmd = captured["cmd"]
    assert "--keeplength" in cmd
    assert "--addfragments" in cmd
    assert cmd[cmd.index("--thread") + 1] == "2"
    assert cmd[cmd.index("--maxiterate") + 1] == "1000"


def test_minimap2_forwards_preset_and_threads(monkeypatch, tmp_path: Path) -> None:
    captured: dict[str, list[str]] = {}

    def fake_run(caps, cmd, **kw):
        captured["cmd"] = [str(c) for c in cmd]
        raise RuntimeError("stop after capture")

    monkeypatch.setattr(mm_mod, "run_tool", fake_run)
    ref, qry = _two_genomes(tmp_path)
    params = AlignParams(threads=3, extra={"preset": "asm10"})
    with pytest.raises(RuntimeError):
        mm_mod.Minimap2Aligner().align([ref, qry], ref, tmp_path / "out", params, _LOG)

    cmd = captured["cmd"]
    assert "-a" in cmd
    assert cmd[cmd.index("-x") + 1] == "asm10"
    assert cmd[cmd.index("-t") + 1] == "3"


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


def test_mafft_assembles_reference_anchored_msa(monkeypatch, tmp_path: Path) -> None:
    def fake_run(caps, cmd, **kw):
        # emulate `mafft --keeplength`: reference row + one aligned contig row
        kw["stdout_path"].write_text(">ref\nACGTACGTACGT\n>qry_contig\nACGTACGTACGT\n")
        return ""

    monkeypatch.setattr(mafft_mod, "run_tool", fake_run)
    ref, qry = _two_genomes(tmp_path)
    result = mafft_mod.MafftAligner().align([ref, qry], ref, tmp_path / "out",
                                            AlignParams(threads=1), _LOG)
    msa = _read_msa(result.msa_fasta)
    assert set(msa) == {"ref", "qry"}
    assert len(msa["ref"]) == len(msa["qry"]) == 12


def test_progressivemauve_concatenation_is_rectangular(monkeypatch, tmp_path: Path) -> None:
    # The reference is 8 bp. Query A aligns to only the first 4 reference bases,
    # query B to all 8. Each per-query XMFA must still project onto the full
    # reference length so the concatenated MSA has equal-length rows.
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref\nACGTACGT\n")
    qry_a = tmp_path / "qryA.fasta"
    qry_a.write_text(">qryA\nACGT\n")
    qry_b = tmp_path / "qryB.fasta"
    qry_b.write_text(">qryB\nACGTACGT\n")

    def fake_run(caps, cmd, **kw):
        cmd = [str(c) for c in cmd]
        xmfa = Path(cmd[cmd.index("--output") + 1])
        ref_path, query_path = cmd[-2], cmd[-1]
        extent = 4 if "qryA" in Path(query_path).stem else 8
        bases = "ACGTACGT"[:extent]
        xmfa.write_text(
            f"#Sequence1File\t{ref_path}\n"
            f"#Sequence2File\t{query_path}\n"
            f"> 1:1-{extent} + {ref_path}\n{bases}\n"
            f"> 2:1-{extent} + {query_path}\n{bases}\n"
            "=\n"
        )
        return ""

    monkeypatch.setattr(pm_mod, "run_tool", fake_run)
    result = pm_mod.ProgressiveMauveAligner().align(
        [ref, qry_a, qry_b], ref, tmp_path / "out", AlignParams(threads=1), _LOG
    )
    msa = _read_msa(result.msa_fasta)
    assert set(msa) == {"ref", "qryA", "qryB"}
    widths = {len(seq) for seq in msa.values()}
    assert widths == {8}, f"rows are ragged: {[(k, len(v)) for k, v in msa.items()]}"


def test_minimap2_assembles_reference_anchored_msa(monkeypatch, tmp_path: Path) -> None:
    def fake_run(caps, cmd, **kw):
        kw["stdout_path"].write_text(
            "@SQ\tSN:ref\tLN:12\n"
            "q\t0\tref\t1\t60\t12M\t*\t0\t0\tACGTACGTACGT\t*\n"
        )
        return ""

    monkeypatch.setattr(mm_mod, "run_tool", fake_run)
    ref, qry = _two_genomes(tmp_path)
    result = mm_mod.Minimap2Aligner().align([ref, qry], ref, tmp_path / "out",
                                            AlignParams(threads=1), _LOG)
    msa = _read_msa(result.msa_fasta)
    assert set(msa) == {"ref", "qry"}
    assert msa["ref"] == "ACGTACGTACGT"
    assert msa["qry"] == "ACGTACGTACGT"


def test_mafft_adjusts_direction_of_added_sequences(monkeypatch, tmp_path: Path) -> None:
    """A draft contig or a whole genome may be on the opposite strand to the backbone.
    Without --adjustdirection MAFFT aligns it as given and the row is noise."""
    captured: dict[str, list[str]] = {}

    def fake_run(caps, cmd, **kw):
        captured["cmd"] = [str(c) for c in cmd]
        raise RuntimeError("stop after capture")

    monkeypatch.setattr(mafft_mod, "run_tool", fake_run)
    ref, qry = _two_genomes(tmp_path)
    with pytest.raises(RuntimeError):
        mafft_mod.MafftAligner().align([ref, qry], ref, tmp_path / "out",
                                       AlignParams(threads=1), _LOG)

    cmd = captured["cmd"]
    assert "--adjustdirection" in cmd
    # The option must precede --addfragments, whose two arguments are positional.
    assert cmd.index("--adjustdirection") < cmd.index("--addfragments")


def _mafft_commands(monkeypatch, tmp_path: Path, *, ref_len: int, n_queries: int,
                    params: AlignParams, memory: int | None = None) -> list[list[str]]:
    """Run the MAFFT adapter against a fake `mafft` and return the commands it issued."""
    commands: list[list[str]] = []
    seq = "ACGT" * (ref_len // 4)

    def fake_run(caps, cmd, **kw):
        commands.append([str(c) for c in cmd])
        kw["stdout_path"].write_text(f">ref\n{seq}\n>added\n{seq}\n")
        return ""

    monkeypatch.setattr(mafft_mod, "run_tool", fake_run)
    if memory is not None:
        monkeypatch.setattr(mafft_mod, "_physical_memory", lambda: memory)
    tmp_path.mkdir(parents=True, exist_ok=True)
    ref = tmp_path / "ref.fasta"
    ref.write_text(f">ref\n{seq}\n")
    genomes = [ref]
    for i in range(n_queries):
        qry = tmp_path / f"qry{i}.fasta"
        qry.write_text(f">qry{i}\n{seq}\n")
        genomes.append(qry)
    mafft_mod.MafftAligner().align(genomes, ref, tmp_path / "out", params, _LOG)
    return commands


def test_mafft_uses_kmer_distances_on_a_long_backbone(monkeypatch, tmp_path: Path) -> None:
    """`--addfragments` otherwise runs a pairwise local alignment whose memory grows with
    the product of the two lengths (over 100 GB for a 200 kb genome pair)."""
    (cmd,) = _mafft_commands(monkeypatch, tmp_path, ref_len=20_004, n_queries=1,
                             params=AlignParams(threads=1))
    assert "--6merpair" in cmd
    assert cmd.index("--6merpair") < cmd.index("--addfragments")


def test_mafft_keeps_pairwise_distances_on_a_short_backbone(monkeypatch, tmp_path: Path) -> None:
    (cmd,) = _mafft_commands(monkeypatch, tmp_path, ref_len=20_000, n_queries=1,
                             params=AlignParams(threads=1))
    assert "--6merpair" not in cmd


def test_mafft_explicit_sixmerpair_overrides_the_length_rule(monkeypatch, tmp_path: Path) -> None:
    (long_cmd,) = _mafft_commands(monkeypatch, tmp_path / "long", ref_len=20_004, n_queries=1,
                                  params=AlignParams(threads=1, extra={"sixmerpair": "false"}))
    assert "--6merpair" not in long_cmd
    (short_cmd,) = _mafft_commands(monkeypatch, tmp_path / "short", ref_len=400, n_queries=1,
                                   params=AlignParams(threads=1, extra={"sixmerpair": "true"}))
    assert "--6merpair" in short_cmd


def test_mafft_splits_threads_across_concurrent_jobs(monkeypatch, tmp_path: Path) -> None:
    """Four jobs at once, each with `--thread 4`, would use sixteen threads."""
    commands = _mafft_commands(monkeypatch, tmp_path, ref_len=400, n_queries=5,
                               params=AlignParams(threads=4), memory=64 * 2**30)
    assert len(commands) == 5
    assert {cmd[cmd.index("--thread") + 1] for cmd in commands} == {"1"}


def test_mafft_gives_a_single_job_all_threads(monkeypatch, tmp_path: Path) -> None:
    (cmd,) = _mafft_commands(monkeypatch, tmp_path, ref_len=400, n_queries=1,
                             params=AlignParams(threads=4), memory=64 * 2**30)
    assert cmd[cmd.index("--thread") + 1] == "4"


def test_mafft_limits_concurrent_pairwise_jobs_by_memory(monkeypatch, tmp_path: Path) -> None:
    """Below the length threshold each pairwise job still needs several GB; with little
    memory the jobs run one at a time (seen as one job holding every thread)."""
    workers: list[int] = []
    real_map = mafft_mod.parallel_map

    def spy_map(fn, items, n, **kw):
        workers.append(n)
        return real_map(fn, items, n, **kw)

    monkeypatch.setattr(mafft_mod, "parallel_map", spy_map)
    commands = _mafft_commands(monkeypatch, tmp_path, ref_len=16_000, n_queries=5,
                               params=AlignParams(threads=4), memory=8 * 2**30)
    assert workers == [1]
    assert {cmd[cmd.index("--thread") + 1] for cmd in commands} == {"4"}


def test_mafft_kmer_jobs_are_not_limited_by_memory(monkeypatch, tmp_path: Path) -> None:
    commands = _mafft_commands(monkeypatch, tmp_path, ref_len=40_000, n_queries=5,
                               params=AlignParams(threads=4), memory=8 * 2**30)
    assert {cmd[cmd.index("--thread") + 1] for cmd in commands} == {"1"}


def _mafft_long_query_commands(monkeypatch, tmp_path: Path, *, memory: int) -> list[list[str]]:
    """A short backbone with two genomes to add, one of them 40 kb."""
    commands: list[list[str]] = []

    def fake_run(caps, cmd, **kw):
        commands.append([str(c) for c in cmd])
        kw["stdout_path"].write_text(">ref\n" + "ACGT" * 2000 + "\n>added\n" + "ACGT" * 2000 + "\n")
        return ""

    monkeypatch.setattr(mafft_mod, "run_tool", fake_run)
    monkeypatch.setattr(mafft_mod, "_physical_memory", lambda: memory)
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref\n" + "ACGT" * 2000 + "\n")          # 8 kb
    long_qry = tmp_path / "long.fasta"
    long_qry.write_text(">long\n" + "ACGT" * 10_000 + "\n")  # 40 kb
    short_qry = tmp_path / "short.fasta"
    short_qry.write_text(">short\n" + "ACGT" * 2000 + "\n")
    mafft_mod.MafftAligner().align([ref, long_qry, short_qry], ref, tmp_path / "out",
                                   AlignParams(threads=2), _LOG)
    return commands


def test_mafft_length_rule_counts_the_genomes_being_added(monkeypatch, tmp_path: Path) -> None:
    """Pairwise memory grows with backbone length times added-sequence length, so a long
    genome added to a short backbone needs k-mer distances as well."""
    commands = _mafft_long_query_commands(monkeypatch, tmp_path, memory=64 * 2**30)
    assert len(commands) == 2
    assert all("--6merpair" in cmd for cmd in commands)


def test_mafft_memory_estimate_uses_the_longest_added_sequence(
    monkeypatch, tmp_path: Path
) -> None:
    # 8 kb backbone, 18 kb added genome: pairwise mode, estimated at 24 * 8e3 * 18e3
    # = 3.5 GB per job. Half of 6 GB holds one such job; an 8 kb x 8 kb estimate
    # (1.5 GB) would have allowed two.
    commands: list[list[str]] = []

    def fake_run(caps, cmd, **kw):
        commands.append([str(c) for c in cmd])
        kw["stdout_path"].write_text(">ref\n" + "ACGT" * 2000 + "\n>added\n" + "ACGT" * 2000 + "\n")
        return ""

    monkeypatch.setattr(mafft_mod, "run_tool", fake_run)
    monkeypatch.setattr(mafft_mod, "_physical_memory", lambda: 6 * 2**30)
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref\n" + "ACGT" * 2000 + "\n")
    genomes = [ref]
    for name, n in (("a", 4500), ("b", 2000)):
        qry = tmp_path / f"{name}.fasta"
        qry.write_text(f">{name}\n" + "ACGT" * n + "\n")
        genomes.append(qry)
    mafft_mod.MafftAligner().align(genomes, ref, tmp_path / "out", AlignParams(threads=2), _LOG)
    assert all("--6merpair" not in cmd for cmd in commands)
    assert {cmd[cmd.index("--thread") + 1] for cmd in commands} == {"2"}  # one job at a time


def test_merge_added_fragments_ignores_reversed_name_prefix(tmp_path: Path) -> None:
    # `mafft --adjustdirection` renames a record it reverse-complemented to `_R_<name>`.
    # The merge is positional (first record = reference, the rest = contigs), so the
    # renamed contig must still be merged.
    aligned = tmp_path / "a.fasta"
    aligned.write_text(
        ">ref\nACGTACGT\n"
        ">contig1\nAC------\n"
        ">_R_contig2\n----ACGT\n"
    )
    _, merged = merge_added_fragments(aligned)
    assert merged == "AC--ACGT"


# --- progressiveMauve row names come from the staged labels ----------------------
def _fake_mauve(caps, cmd, **kw):
    """Write the XMFA the way progressiveMauve does: sequence names are the paths it
    was given on the command line (which the adapter passes resolved)."""
    cmd = [str(c) for c in cmd]
    xmfa = Path(cmd[cmd.index("--output") + 1])
    ref_path, query_path = cmd[-2], cmd[-1]
    ref_seq = "".join(
        line for line in Path(ref_path).read_text().splitlines() if not line.startswith(">")
    )
    qry_seq = "".join(
        line for line in Path(query_path).read_text().splitlines() if not line.startswith(">")
    )
    xmfa.write_text(
        f"#Sequence1File\t{ref_path}\n#Sequence2File\t{query_path}\n"
        f"> 1:1-8 + {ref_path}\n{ref_seq}\n> 2:1-8 + {query_path}\n{qry_seq}\n=\n"
    )
    return ""


def _mauve_rows(monkeypatch, staged: list[Path], out_dir: Path) -> list[tuple[str, str]]:
    monkeypatch.setattr(pm_mod, "run_tool", _fake_mauve)
    result = pm_mod.ProgressiveMauveAligner().align(
        staged, staged[0], out_dir, AlignParams(threads=1), _LOG
    )
    rows: list[tuple[str, str]] = []
    for line in result.msa_fasta.read_text().splitlines():
        if line.startswith(">"):
            rows.append((line[1:], ""))
        else:
            rows[-1] = (rows[-1][0], rows[-1][1] + line)
    return rows


def test_progressivemauve_names_rows_by_staged_label_not_symlink_target(
    monkeypatch, tmp_path: Path
) -> None:
    # Staged genomes are symlinks named <label>.fasta. Two of them point at files that
    # are both called genome.fna, and one target is named like the backbone's label.
    store = tmp_path / "store"
    for sub in ("x", "y", "z"):
        (store / sub).mkdir(parents=True)
    (store / "x" / "genome.fna").write_text(">r\nACGTACGT\n")
    (store / "y" / "genome.fna").write_text(">b\nACGTACGA\n")
    (store / "z" / "aRef.fna").write_text(">c\nTTTTTTTT\n")
    stage = tmp_path / "stage"
    stage.mkdir()
    staged = []
    for label, target in (("aRef", "x/genome.fna"), ("panelB", "y/genome.fna"),
                          ("panelC", "z/aRef.fna")):
        link = stage / f"{label}.fasta"
        link.symlink_to((store / target).resolve())
        staged.append(link)

    rows = _mauve_rows(monkeypatch, staged, tmp_path / "out")

    assert rows == [("aRef", "ACGTACGT"), ("panelB", "ACGTACGA"), ("panelC", "TTTTTTTT")]


def test_progressivemauve_tolerates_whitespace_in_the_path(monkeypatch, tmp_path: Path) -> None:
    stage = tmp_path / "dir with space"
    stage.mkdir()
    staged = []
    for label, seq in (("ref", "ACGTACGT"), ("qryA", "ACGTACGA"), ("qryB", "ACGTACGC")):
        path = stage / f"{label}.fasta"
        path.write_text(f">{label}\n{seq}\n")
        staged.append(path)

    rows = _mauve_rows(monkeypatch, staged, tmp_path / "out")

    assert [name for name, _ in rows] == ["ref", "qryA", "qryB"]


def test_progressivemauve_rejects_a_projection_that_is_not_pairwise(
    monkeypatch, tmp_path: Path
) -> None:
    from tessera.core.errors import OutputError

    def one_sequence_xmfa(caps, cmd, **kw):
        cmd = [str(c) for c in cmd]
        xmfa = Path(cmd[cmd.index("--output") + 1])
        ref_path = cmd[-2]
        xmfa.write_text(f"#Sequence1File\t{ref_path}\n> 1:1-8 + {ref_path}\nACGTACGT\n=\n")
        return ""

    monkeypatch.setattr(pm_mod, "run_tool", one_sequence_xmfa)
    genomes = []
    for label in ("ref", "qryA", "qryB"):
        path = tmp_path / f"{label}.fasta"
        path.write_text(f">{label}\nACGTACGT\n")
        genomes.append(path)
    with pytest.raises(OutputError, match=r"qryA\.fa.*found 1 record"):
        pm_mod.ProgressiveMauveAligner().align(
            genomes, genomes[0], tmp_path / "out", AlignParams(threads=1), _LOG
        )


# --- minimap2 shares its thread budget across concurrent jobs ----------------------
def _minimap2_runs(monkeypatch, tmp_path: Path, *, threads: int, n_queries: int
                   ) -> tuple[list[int], list[list[str]]]:
    """Run the minimap2 adapter against a fake `minimap2`; return (pool sizes, commands)."""
    commands: list[list[str]] = []
    workers: list[int] = []
    real_map = mm_mod.parallel_map

    def spy_map(fn, items, n, **kw):
        workers.append(n)
        return real_map(fn, items, n, **kw)

    def fake_run(caps, cmd, **kw):
        commands.append([str(c) for c in cmd])
        kw["stdout_path"].write_text(
            "@SQ\tSN:ref\tLN:12\nq\t0\tref\t1\t60\t12M\t*\t0\t0\tACGTACGTACGT\t*\n"
        )
        return ""

    monkeypatch.setattr(mm_mod, "parallel_map", spy_map)
    monkeypatch.setattr(mm_mod, "run_tool", fake_run)
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref\nACGTACGTACGT\n")
    genomes = [ref]
    for i in range(n_queries):
        qry = tmp_path / f"qry{i}.fasta"
        qry.write_text(f">qry{i}\nACGTACGTACGT\n")
        genomes.append(qry)
    mm_mod.Minimap2Aligner().align(genomes, ref, tmp_path / "out",
                                   AlignParams(threads=threads), _LOG)
    return workers, commands


def test_minimap2_splits_threads_across_concurrent_jobs(monkeypatch, tmp_path: Path) -> None:
    """--threads 4 used to start four minimap2 processes of four threads each."""
    workers, commands = _minimap2_runs(monkeypatch, tmp_path, threads=4, n_queries=5)
    assert workers == [4]
    assert len(commands) == 5
    assert {cmd[cmd.index("-t") + 1] for cmd in commands} == {"1"}


def test_minimap2_gives_a_single_job_all_threads(monkeypatch, tmp_path: Path) -> None:
    workers, commands = _minimap2_runs(monkeypatch, tmp_path, threads=4, n_queries=1)
    assert workers == [1]
    assert commands[0][commands[0].index("-t") + 1] == "4"


def test_minimap2_never_starts_more_threads_than_the_budget(monkeypatch, tmp_path: Path) -> None:
    # A budget that does not divide evenly: five threads over two queries. Two jobs of
    # two threads stay within it; two jobs of five, or of three, would not.
    workers, commands = _minimap2_runs(monkeypatch, tmp_path, threads=5, n_queries=2)
    per_job = {int(cmd[cmd.index("-t") + 1]) for cmd in commands}
    assert workers == [2]
    assert per_job == {2}
    assert workers[0] * max(per_job) <= 5


# --- a genome the aligner did not place is named -----------------------------------
# Not under the "tessera" logger: the CLI tests switch its propagation off, and caplog
# only sees records that reach the root logger.
_WARN_LOG = logging.getLogger("aligner_unplaced_test")


def _ref_and_queries(tmp_path: Path) -> list[Path]:
    paths = []
    for name, seq in (("ref", "ACGTACGTACGT"), ("near", "ACGTACGTACGT"), ("far", "TTTTTTTTTTTT")):
        path = tmp_path / f"{name}.fasta"
        path.write_text(f">{name}\n{seq}\n")
        paths.append(path)
    return paths


def test_minimap2_names_a_genome_it_could_not_map(monkeypatch, tmp_path: Path, caplog) -> None:
    def fake_run(caps, cmd, **kw):
        if "far" in str(cmd[-1]):
            body = "far\t4\t*\t0\t0\t*\t*\t0\t0\tTTTTTTTTTTTT\t*\n"  # unmapped
        else:
            body = "near\t0\tref\t1\t60\t12M\t*\t0\t0\tACGTACGTACGT\t*\n"
        kw["stdout_path"].write_text("@SQ\tSN:ref\tLN:12\n" + body)
        return ""

    monkeypatch.setattr(mm_mod, "run_tool", fake_run)
    genomes = _ref_and_queries(tmp_path)
    with caplog.at_level(logging.WARNING, logger=_WARN_LOG.name):
        result = mm_mod.Minimap2Aligner().align(genomes, genomes[0], tmp_path / "out",
                                                AlignParams(threads=1), _WARN_LOG)
    assert _read_msa(result.msa_fasta)["far"] == "-" * 12
    assert "1 genome(s)" in caplog.text
    assert "far" in caplog.text and "near" not in caplog.text


def test_mafft_names_a_genome_it_could_not_place(monkeypatch, tmp_path: Path, caplog) -> None:
    def fake_run(caps, cmd, **kw):
        row = "-" * 12 if "far" in str(cmd[cmd.index("--addfragments") + 1]) else "ACGTACGTACGT"
        kw["stdout_path"].write_text(f">ref\nACGTACGTACGT\n>added\n{row}\n")
        return ""

    monkeypatch.setattr(mafft_mod, "run_tool", fake_run)
    genomes = _ref_and_queries(tmp_path)
    with caplog.at_level(logging.WARNING, logger=_WARN_LOG.name):
        mafft_mod.MafftAligner().align(genomes, genomes[0], tmp_path / "out",
                                       AlignParams(threads=1), _WARN_LOG)
    assert "far" in caplog.text and "near" not in caplog.text


def test_progressivemauve_names_a_genome_in_no_block(monkeypatch, tmp_path: Path, caplog) -> None:
    def fake_run(caps, cmd, **kw):
        cmd = [str(c) for c in cmd]
        xmfa = Path(cmd[cmd.index("--output") + 1])
        ref_path, query_path = cmd[-2], cmd[-1]
        header = f"#Sequence1File\t{ref_path}\n#Sequence2File\t{query_path}\n"
        if "far" in query_path:  # the reference aligned to nothing
            xmfa.write_text(header + f"> 1:1-12 + {ref_path}\nACGTACGTACGT\n=\n")
        else:
            xmfa.write_text(header + f"> 1:1-12 + {ref_path}\nACGTACGTACGT\n"
                            f"> 2:1-12 + {query_path}\nACGTACGTACGT\n=\n")
        return ""

    monkeypatch.setattr(pm_mod, "run_tool", fake_run)
    genomes = _ref_and_queries(tmp_path)
    with caplog.at_level(logging.WARNING, logger=_WARN_LOG.name):
        pm_mod.ProgressiveMauveAligner().align(genomes, genomes[0], tmp_path / "out",
                                               AlignParams(threads=1), _WARN_LOG)
    assert "far" in caplog.text and "near" not in caplog.text
