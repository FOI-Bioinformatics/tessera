"""Inputs the CLI must reject before any work starts.

Each case used to be accepted: the option was silently ignored, or the run failed
later under "Unexpected error". None of these tests may reach the network -- the
panel-building entry points are replaced with a function that fails the test if the
command gets that far.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

from tessera.cli.main import app

runner = CliRunner()

EXAMPLE = "example_data/divergent.msa.fasta"


def _must_not_run(*args, **kwargs):
    raise AssertionError("the command reached the pipeline; it should have been rejected")


@pytest.fixture
def no_pipeline(monkeypatch):
    """Fail the test if a command gets past validation into panel building or typing."""
    monkeypatch.setattr("tessera.discover.iterate.fill_references", _must_not_run)
    monkeypatch.setattr("tessera.discover.lineage_assign.assign_lineages", _must_not_run)
    monkeypatch.setattr("tessera.discover.run.find_references", _must_not_run)


def _query(tmp_path: Path) -> Path:
    path = tmp_path / "query.fasta"
    path.write_text(">query\nACGTACGTACGT\n")
    return path


def _collection(tmp_path: Path) -> Path:
    directory = tmp_path / "collection"
    directory.mkdir()
    (directory / "ref.fasta").write_text(">ref\nACGTACGTACGT\n")
    return directory


def _rejected(result, needle: str) -> None:
    assert result.exit_code == 1, result.output
    assert "Unexpected error" not in result.output
    assert needle in " ".join(result.output.split())


# --- --lineage-map must exist ---------------------------------------------

def test_recomb_rejects_a_missing_lineage_map(tmp_path: Path) -> None:
    """A mistyped path used to be ignored: the run exited 0 with an untyped report, and
    the barcode caller and donor re-attribution silently did nothing."""
    result = runner.invoke(app, [
        "recomb", "--msa", EXAMPLE, "--query", "query",
        "--output", str(tmp_path / "out"), "--lineage-map", str(tmp_path / "nope.tsv"),
    ])
    _rejected(result, "--lineage-map file not found")
    assert not (tmp_path / "out" / "report.html").exists()


def test_recomb_rejects_a_lineage_map_that_is_a_directory(tmp_path: Path) -> None:
    result = runner.invoke(app, [
        "recomb", "--msa", EXAMPLE, "--query", "query",
        "--output", str(tmp_path / "out"), "--lineage-map", str(tmp_path),
    ])
    _rejected(result, "is a directory, not a file")


def test_type_lineages_rejects_a_missing_lineage_map(tmp_path: Path, no_pipeline) -> None:
    result = runner.invoke(app, [
        "type-lineages", "--collection", str(_collection(tmp_path)),
        "--output", str(tmp_path / "out"), "--lineage-map", str(tmp_path / "nope.tsv"),
    ])
    _rejected(result, "--lineage-map file not found")


@pytest.mark.parametrize("command", ["detect", "fill-references", "build-panel"])
def test_panel_commands_reject_a_missing_lineage_map(
    tmp_path: Path, no_pipeline, command: str
) -> None:
    result = runner.invoke(app, [
        command, "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        "--lineage-map", str(tmp_path / "nope.tsv"),
    ])
    _rejected(result, "--lineage-map file not found")


# --- paths ------------------------------------------------------------------

def test_find_references_rejects_a_missing_msa(tmp_path: Path, no_pipeline) -> None:
    result = runner.invoke(app, [
        "find-references", "--msa", str(tmp_path / "nope.fasta"), "--query", "query",
        "--output", str(tmp_path / "out"),
    ])
    _rejected(result, "MSA file not found")


def test_recomb_rejects_an_output_path_that_is_a_file(tmp_path: Path) -> None:
    occupied = tmp_path / "out"
    occupied.write_text("not a directory\n")
    result = runner.invoke(app, [
        "recomb", "--msa", EXAMPLE, "--query", "query", "--output", str(occupied),
    ])
    _rejected(result, "is not a directory")
    assert occupied.read_text() == "not a directory\n"  # left untouched


# --- numeric options on the panel-building commands ------------------------

@pytest.mark.parametrize("command", ["detect", "fill-references", "build-panel"])
@pytest.mark.parametrize(
    ("option", "value"),
    [("--max-rounds", "0"), ("--window-size", "0"), ("--window-step", "0")],
)
def test_panel_commands_reject_out_of_range_options(
    tmp_path: Path, no_pipeline, command: str, option: str, value: str
) -> None:
    """`fill-references --max-rounds 0` exited 0 having built no alignment and no
    report; a bad window reached the scan only after seeding and alignment."""
    result = runner.invoke(app, [
        command, "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        option, value,
    ])
    _rejected(result, f"Invalid {option}")


@pytest.mark.parametrize(("option", "value"), [("--window-size", "0"), ("--window-step", "0")])
def test_find_references_rejects_out_of_range_windows(
    tmp_path: Path, no_pipeline, option: str, value: str
) -> None:
    result = runner.invoke(app, [
        "find-references", "--msa", EXAMPLE, "--query", "query",
        "--output", str(tmp_path / "out"), option, value,
    ])
    _rejected(result, f"Invalid {option}")


# --- reassort ---------------------------------------------------------------

def _segments(tmp_path: Path) -> Path:
    path = tmp_path / "segments.fasta"
    path.write_text(">HA\nACGTACGT\n>NA\nACGTACGT\n")
    return path


@pytest.fixture
def no_reassort(monkeypatch):
    monkeypatch.setattr("tessera.cli.cmd_reassort.assign_segments", _must_not_run)


@pytest.mark.parametrize(
    ("option", "value", "needle"),
    [("--ani-floor", "500", "Invalid --ani-floor"),
     ("--ani-floor", "-1", "Invalid --ani-floor"),
     ("--margin", "-3", "Invalid --margin")],
)
def test_reassort_rejects_out_of_range_options(
    tmp_path: Path, no_reassort, option: str, value: str, needle: str
) -> None:
    """A negative margin makes every near-best set empty (a clonal pair then reads as
    undetermined); an ANI floor above 100 leaves every segment unassigned. Both exited 0."""
    result = runner.invoke(app, [
        "reassort", "--query", str(_segments(tmp_path)), "--output", str(tmp_path / "out"),
        option, value,
    ])
    _rejected(result, needle)


def test_reassort_rejects_a_dataset_override_for_an_unknown_segment(
    tmp_path: Path, no_reassort
) -> None:
    """A typo in the segment name silently fell back to dataset auto-detection."""
    result = runner.invoke(app, [
        "reassort", "--query", str(_segments(tmp_path)), "--output", str(tmp_path / "out"),
        "--dataset", "HAA=nextstrain/flu/h3n2/ha",
    ])
    _rejected(result, "HAA")
    assert "HA, NA" in " ".join(result.output.split())  # names the segments that do exist


def test_segment_scan_tsv_keeps_one_row_per_segment(tmp_path: Path, monkeypatch) -> None:
    """An aligner error message spans several lines and may hold tabs; written verbatim
    into the note column it broke the table."""
    from tessera.reassort.assign import ReassortmentResult, SegmentAssignment
    from tessera.reassort.scan import SegmentScan

    def fake_assign(query, **kwargs):
        return ReassortmentResult(
            segments=[SegmentAssignment("HA", "dataset", None, None, 0.0, "unassigned")],
            verdict="undetermined", groups=[], pair_notes=[],
            scans=[SegmentScan("HA", False, False, 0,
                               "scan failed: mafft failed:\nline1\tx\nline2")],
        )

    monkeypatch.setattr("tessera.cli.cmd_reassort.assign_segments", fake_assign)
    out = tmp_path / "out"
    result = runner.invoke(app, [
        "reassort", "--query", str(_segments(tmp_path)), "--output", str(out),
    ])
    assert result.exit_code == 0, result.output
    lines = (out / "segment_scan.tsv").read_text().splitlines()
    assert len(lines) == 2  # header + one segment
    assert all(len(line.split("\t")) == 4 for line in lines)
    assert lines[1].split("\t")[3] == "scan failed: mafft failed: line1 x line2"


# --- type-lineages reads the same collection the other commands do ---------

def test_type_lineages_accepts_what_a_collection_may_hold(tmp_path: Path, monkeypatch) -> None:
    """Every other command treats each file in the collection as a genome, whatever its
    extension and gzip-compressed or not; type-lineages filtered on three suffixes and
    reported "No FASTA genomes found"."""
    import gzip

    collection = tmp_path / "collection"
    collection.mkdir()
    with gzip.open(collection / "a.fasta.gz", "wt") as fo:
        fo.write(">a\nACGTACGT\n")
    (collection / "b.fas").write_text(">b\nACGTACGT\n")
    (collection / ".DS_Store").write_text("not a genome")
    seen: list[str] = []

    def fake_assign(genomes, **kwargs):
        seen.extend(p.name for p in genomes)
        return [("a", "L1", "denovo"), ("b.fas", "L1", "denovo")]

    monkeypatch.setattr("tessera.discover.lineage_assign.assign_lineages", fake_assign)
    out = tmp_path / "out"
    result = runner.invoke(app, [
        "type-lineages", "--collection", str(collection), "--output", str(out),
    ])
    assert result.exit_code == 0, result.output
    assert seen == ["a.fasta.gz", "b.fas"]  # hidden files are not genomes
    assert (out / "lineages.tsv").exists()


def test_type_lineages_rejects_a_file_that_is_not_fasta(tmp_path: Path, no_pipeline) -> None:
    collection = _collection(tmp_path)
    (collection / "notes.txt").write_text("these are my notes\n")
    result = runner.invoke(app, [
        "type-lineages", "--collection", str(collection), "--output", str(tmp_path / "out"),
    ])
    _rejected(result, "does not look like a FASTA file")


@pytest.mark.parametrize("value", ["-0.1", "1.5"])
def test_fill_references_rejects_a_min_improvement_outside_zero_one(
    tmp_path: Path, no_pipeline, value: str
) -> None:
    """--min-improvement is a fraction of the previous round's coverage shortfall."""
    result = runner.invoke(app, [
        "fill-references", "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        "--min-improvement", value,
    ])
    _rejected(result, "Invalid --min-improvement")


# --- seed-source options that cannot act together (audit 4, E5) ------------------

@pytest.mark.parametrize("command", ["detect", "build-panel"])
@pytest.mark.parametrize(
    "nextclade_args",
    [["--nextclade"], ["--nextclade-dataset", "nextstrain/sars-cov-2/XBB"]],
)
def test_a_candidate_pool_and_a_nextclade_request_are_refused_together(
    tmp_path: Path, no_pipeline, command: str, nextclade_args: list[str]
) -> None:
    """Both choose where the panel is recruited from. The pool used to win in silence,
    so a run asked to use Nextclade clades used the pool and did not say so."""
    pool_dir = tmp_path / "pool"
    pool_dir.mkdir()
    result = runner.invoke(app, [
        command, "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        "--candidate-pool", str(pool_dir), *nextclade_args,
    ])
    _rejected(result, "--candidate-pool and --nextclade/--nextclade-dataset")


@pytest.mark.parametrize(
    ("command", "args"),
    [("detect", []),
     ("detect", ["--candidate-pool", "POOL"]),
     ("fill-references", ["--seed-source", "local", "--candidate-pool", "POOL"]),
     ("fill-references", ["--seed-source", "ncbi-virus"])],
)
def test_pool_consensus_without_a_nextclade_pool_is_refused(
    tmp_path: Path, no_pipeline, command: str, args: list[str]
) -> None:
    """--pool-consensus acts only on a Nextclade pool; under any other seed source it
    did nothing and the run looked like one that had used it."""
    pool_dir = tmp_path / "pool"
    pool_dir.mkdir()
    args = [str(pool_dir) if a == "POOL" else a for a in args]
    result = runner.invoke(app, [
        command, "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        "--pool-consensus", *args,
    ])
    _rejected(result, "--pool-consensus")


def test_fill_references_refuses_a_candidate_pool_under_another_seed_source(
    tmp_path: Path, no_pipeline
) -> None:
    """fill-references reads --candidate-pool only with --seed-source local; without it
    the default BLAST seeding ran and the pool was ignored."""
    pool_dir = tmp_path / "pool"
    pool_dir.mkdir()
    result = runner.invoke(app, [
        "fill-references", "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        "--candidate-pool", str(pool_dir),
    ])
    _rejected(result, "--seed-source local")


# --- discover commands: coverage fractions and search counts (audit 4, E5) ---

@pytest.mark.parametrize(
    ("option", "value"),
    [("--coverage-floor", "50"), ("--coverage-rel-drop", "5"), ("--top-gaps", "0"),
     ("--max-hits", "0"), ("--subtile", "-1")],
)
@pytest.mark.parametrize("command", ["fill-references", "find-references"])
def test_discover_commands_reject_out_of_range_search_options(
    tmp_path: Path, no_pipeline, command: str, option: str, value: str
) -> None:
    """recomb range-checked the coverage fractions; here --coverage-floor 50 made every
    window a gap, and --top-gaps 0 searched nothing and exited cleanly."""
    head = (["--msa", EXAMPLE, "--query", "query"] if command == "find-references"
            else ["--query", str(_query(tmp_path))])
    result = runner.invoke(app, [command, *head, "--output", str(tmp_path / "out"),
                                 option, value])
    _rejected(result, f"Invalid {option}")


@pytest.mark.parametrize(("option", "value"), [("--seed-hits", "0"), ("--seed-window", "0")])
def test_fill_references_rejects_empty_seed_searches(
    tmp_path: Path, no_pipeline, option: str, value: str
) -> None:
    result = runner.invoke(app, [
        "fill-references", "--query", str(_query(tmp_path)),
        "--output", str(tmp_path / "out"), option, value,
    ])
    _rejected(result, f"Invalid {option}")


@pytest.mark.parametrize("command", ["detect", "fill-references", "build-panel"])
def test_panel_commands_reject_zero_threads(tmp_path: Path, no_pipeline, command: str) -> None:
    result = runner.invoke(app, [
        command, "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        "--threads", "0",
    ])
    _rejected(result, "Invalid --threads")


def test_find_references_keeps_the_query_label_verbatim(
    tmp_path: Path, monkeypatch
) -> None:
    """-q is an alignment label, but was declared a path: 'sample/' became 'sample'."""
    seen = {}

    def capture(params, logger):
        seen["query"] = params.query

    monkeypatch.setattr("tessera.discover.run.find_references", capture)
    result = runner.invoke(app, [
        "find-references", "--msa", EXAMPLE, "--query", "sample/",
        "--output", str(tmp_path / "out"),
    ])
    assert result.exit_code == 0, result.output
    assert seen["query"] == "sample/"


def test_build_panel_passes_the_panel_shaping_options(tmp_path: Path, monkeypatch) -> None:
    """build-panel could not set three options that shape the panel it builds."""
    seen = {}

    def capture(params, logger):
        seen.update(vars(params))

    monkeypatch.setattr("tessera.discover.iterate.fill_references", capture)
    result = runner.invoke(app, [
        "build-panel", "--query", str(_query(tmp_path)), "--output", str(tmp_path / "out"),
        # --nextclade: plan E Task 11 refuses --pool-consensus under any other seed
        # source, so this is the one combination valid before and after plan E.
        "--nextclade",
        "--pool-consensus", "--keep-recombinant-lineages", "--deep-typing",
    ])
    assert result.exit_code == 0, result.output
    assert seen["pool_consensus"] and seen["keep_recombinant"] and seen["deep_typing"]


def test_seed_source_help_names_nextclade() -> None:
    import typer

    cmd = typer.main.get_command(app).commands["fill-references"]
    param = next(p for p in cmd.params if "--seed-source" in getattr(p, "opts", []))
    assert "'nextclade'" in (param.help or "")


# --- msa: backbone choice and overwrite (audit 4, E5) ------------------------

def test_msa_rejects_reference_with_query_as_backbone(tmp_path: Path, monkeypatch) -> None:
    """--reference silently won, so the query was not the backbone the user asked for."""
    monkeypatch.setattr("tessera.msa.build.build_msa", _must_not_run)
    result = runner.invoke(app, [
        "msa", "--query", str(_query(tmp_path)), "--collection", str(_collection(tmp_path)),
        "--output", str(tmp_path / "aln.fasta"), "--aligner", "mafft",
        "--reference", "ref", "--query-as-backbone",
    ])
    _rejected(result, "give one")


def test_msa_says_when_it_overwrites_an_alignment(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr("tessera.msa.build.build_msa", lambda params, logger: None)
    existing = tmp_path / "aln.fasta"
    existing.write_text(">old\nACGT\n")
    result = runner.invoke(app, [
        "msa", "--query", str(_query(tmp_path)), "--collection", str(_collection(tmp_path)),
        "--output", str(existing), "--aligner", "mafft",
    ])
    assert result.exit_code == 0, result.output
    assert "Overwriting the existing alignment" in result.output
