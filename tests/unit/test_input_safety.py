"""Regressions from the second production-readiness review.

Each test pins a failure that ran to completion (or crashed as "Unexpected error")
on input a user could plausibly supply: a collection inside the output directory, a
stray non-genome file, a gzipped query, a stale cache directory, a missing binary.
"""

from __future__ import annotations

import gzip
import logging
import sys
from pathlib import Path

import pytest
from typer.testing import CliRunner

from tessera.aligners import sibeliaz
from tessera.cli.main import app
from tessera.core import cache
from tessera.core.errors import MissingBinaryError, OutputError, UserInputError
from tessera.core.io import collection_genomes, copy_collection, read_fasta, stage_genomes
from tessera.core.plugins import ToolCapabilities
from tessera.core.process import run_tool
from tessera.discover import fetch, iterate, panel
from tessera.discover import lineage_assign as la
from tessera.discover import run as discover_run
from tessera.discover.blast import Hit
from tessera.discover.iterate import FillParams, fill_references
from tessera.reassort import assign

_LOG = logging.getLogger("tessera.test")
runner = CliRunner()


def _collection(directory: Path, names: tuple[str, ...] = ("refA", "refB")) -> Path:
    directory.mkdir(parents=True)
    for name in names:
        (directory / f"{name}.fasta").write_text(f">{name}\nACGTACGT\n")
    return directory


# --- the working copy must never be the input ------------------------------

@pytest.mark.parametrize("relative", ["out/collection", "out", "out/collection/inner"])
def test_copy_collection_refuses_an_overlapping_source(tmp_path: Path, relative: str) -> None:
    source = _collection(tmp_path / relative)
    with pytest.raises(UserInputError, match="working copy"):
        copy_collection(source, tmp_path / "out" / "collection")
    assert sorted(p.name for p in source.iterdir()) == ["refA.fasta", "refB.fasta"]


def test_copy_collection_replaces_a_previous_working_copy(tmp_path: Path) -> None:
    stale = _collection(tmp_path / "earlier", names=("stale",))
    source = _collection(tmp_path / "coll")
    dest = tmp_path / "out" / "collection"
    copy_collection(stale, dest)  # an earlier run's working copy
    copy_collection(source, dest)
    assert sorted(p.name for p in dest.iterdir()) == ["refA.fasta", "refB.fasta"]


def test_curate_panel_keeps_a_collection_that_is_the_output_copy(
    monkeypatch, tmp_path: Path
) -> None:
    """`curate-panel -c out/collection -o out` used to delete the collection."""
    monkeypatch.setattr(panel, "skani_available", lambda: True)
    out = tmp_path / "out"
    collection = _collection(out / "collection")
    query = tmp_path / "q.fasta"
    query.write_text(">q\nACGT\n")
    with pytest.raises(UserInputError, match="working copy"):
        panel.curate_collection(query, collection, out, logger=_LOG)
    assert len(list(collection.iterdir())) == 2


def test_fill_references_keeps_a_collection_that_is_the_output_copy(
    monkeypatch, tmp_path: Path
) -> None:
    # Without this the run stops earlier, on a machine that has no Entrez Direct.
    monkeypatch.setattr(iterate, "efetch_available", lambda: True)
    out = tmp_path / "out"
    collection = _collection(out / "collection")
    query = tmp_path / "q.fasta"
    query.write_text(">q\nACGT\n")
    with pytest.raises(UserInputError, match="working copy"):
        fill_references(FillParams(query=query, collection=collection, output=out), _LOG)
    assert len(list(collection.iterdir())) == 2


# --- what counts as a genome -----------------------------------------------

def test_hidden_files_are_not_genomes(tmp_path: Path) -> None:
    collection = _collection(tmp_path / "coll")
    (collection / ".DS_Store").write_bytes(b"\x00\x00\x00\x01Bud1")
    assert [p.name for p in collection_genomes(collection)] == ["refA.fasta", "refB.fasta"]

    query = tmp_path / "q.fasta"
    query.write_text(">q\nACGT\n")
    staged, _ = stage_genomes(query, collection, tmp_path / "stage", _LOG)
    assert sorted(p.name for p in staged) == ["q.fasta", "refA.fasta", "refB.fasta"]


def test_a_non_fasta_file_in_the_collection_is_rejected_by_name(tmp_path: Path) -> None:
    collection = _collection(tmp_path / "coll")
    (collection / "README.txt").write_text("notes about this panel\n")
    query = tmp_path / "q.fasta"
    query.write_text(">q\nACGT\n")
    with pytest.raises(UserInputError, match=r"README\.txt does not look like a FASTA"):
        stage_genomes(query, collection, tmp_path / "stage", _LOG)


def test_read_fasta_reads_gzip_and_crlf(tmp_path: Path) -> None:
    plain = tmp_path / "q.fasta"
    plain.write_bytes(b">q desc\r\nACGT\r\nTTAA\r\n")
    assert read_fasta(plain) == [("q", "ACGTTTAA")]
    zipped = tmp_path / "q.fasta.gz"
    with gzip.open(zipped, "wt") as fo:
        fo.write(">q desc\nACGT\nTTAA\n")
    assert read_fasta(zipped) == [("q", "ACGTTTAA")]


def test_sibeliaz_rejects_a_sequence_id_shared_by_two_genomes(tmp_path: Path) -> None:
    genomes = []
    for stem in ("g1", "g2"):
        path = tmp_path / f"{stem}.fasta"
        path.write_text(">contig_1\nACGT\n")
        genomes.append(path)
    with pytest.raises(UserInputError, match="contig_1.*g1.*g2"):
        sibeliaz._build_seqid_map(genomes)


def test_segment_name_cannot_leave_the_temp_directory(monkeypatch, tmp_path: Path) -> None:
    def no_dataset(*args, **kwargs):
        raise UserInputError("no dataset")

    monkeypatch.setattr(assign, "resolve_dataset", no_dataset)
    work = tmp_path / "a" / "b" / "work"
    work.mkdir(parents=True)
    for name in ("../../escaped", "A/California/07/2009|HA"):
        assign._type_segment(name, "ACGT", {}, 0.8, 0.5, None, None, str(work), _LOG)
    written = sorted(p.relative_to(tmp_path).as_posix() for p in tmp_path.rglob("*.fasta"))
    assert len(written) == 2
    assert all(path.startswith("a/b/work/") for path in written)


# --- caches -----------------------------------------------------------------

def test_a_manifestless_cache_directory_is_replaced_not_obeyed(tmp_path: Path) -> None:
    """A stale directory at the cache path used to discard every fresh build."""
    target = tmp_path / "ncbi" / "taxon"
    target.mkdir(parents=True)
    (target / "old.fasta").write_text(">old\nACGT\n")
    with cache.atomic_cache_dir(target) as staging:
        (staging / "new.fasta").write_text(">new\nACGT\n")
        cache.write_cache_manifest(staging, taxon="t")
    names = [p.name for p in cache.cached_genomes(target, manifest_required=True)]
    assert names == ["new.fasta"]
    assert len(list(target.parent.iterdir())) == 1  # no staging directory left behind


def test_a_complete_cache_from_a_concurrent_run_is_kept(tmp_path: Path) -> None:
    target = tmp_path / "ncbi" / "taxon"
    with cache.atomic_cache_dir(target) as staging:
        (staging / "theirs.fasta").write_text(">theirs\nACGT\n")
        cache.write_cache_manifest(staging, taxon="t")
    with cache.atomic_cache_dir(target) as staging:
        (staging / "ours.fasta").write_text(">ours\nACGT\n")
        cache.write_cache_manifest(staging, taxon="t")
    names = [p.name for p in cache.cached_genomes(target, manifest_required=True)]
    assert names == ["theirs.fasta"]
    assert len(list(target.parent.iterdir())) == 1


def test_reference_typing_builds_its_pool_under_the_cache_root(
    monkeypatch, tmp_path: Path
) -> None:
    """`--cache-dir` is the cache root; the pool must go in a per-dataset directory."""

    class _Dataset:
        path = "org/virus"
        tag = "2026-01-01"

    seen: dict[str, Path] = {}

    def fake_build_pool(dataset, *, cache_dir, logger):
        seen["cache_dir"] = cache_dir
        return []

    monkeypatch.setattr(la, "resolve_dataset", lambda *a, **k: _Dataset())
    monkeypatch.setattr(la, "build_pool", fake_build_pool)
    root = tmp_path / "cache"
    la._reference_tips(
        query=tmp_path / "q.fasta", nextclade_dataset=None, email=None,
        cache_dir=root, logger=_LOG,
    )
    assert seen["cache_dir"] == cache.nextclade_cache("org/virus", "2026-01-01", override=root)
    assert seen["cache_dir"] != root


# --- external tools ---------------------------------------------------------

def test_a_failed_download_leaves_no_file_in_the_collection(
    monkeypatch, tmp_path: Path
) -> None:
    def empty_output(caps, command, *, stdout_path, **kwargs):
        Path(stdout_path).write_text("")
        return ""

    monkeypatch.setattr(fetch, "run_tool", empty_output)
    monkeypatch.setattr(fetch.eutils_throttle, "wait", lambda logger: None)
    with pytest.raises(OutputError):
        fetch.efetch_fasta("XX000000.1", tmp_path, _LOG)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("to_file", [False, True])
def test_a_missing_binary_is_reported_as_one(tmp_path: Path, to_file: bool) -> None:
    caps = ToolCapabilities(name="nope")
    kwargs = {"stdout_path": tmp_path / "out.txt"} if to_file else {}
    with pytest.raises(MissingBinaryError, match="tessera-no-such-binary"):
        run_tool(caps, ["tessera-no-such-binary", "--version"], logger=_LOG, **kwargs)
    assert list(tmp_path.iterdir()) == []


def test_a_missing_input_path_is_not_mistaken_for_a_missing_binary(tmp_path: Path) -> None:
    caps = ToolCapabilities(name="python")
    with pytest.raises(FileNotFoundError):
        run_tool(
            caps, [sys.executable, "-c", "pass"], logger=_LOG,
            stdout_path=tmp_path / "missing-dir" / "out.txt",
        )


# --- recruitment bookkeeping ------------------------------------------------

def test_collection_membership_ignores_the_accession_version(tmp_path: Path) -> None:
    collection = _collection(tmp_path / "coll", names=("NC_045512.2",))

    class _Result:
        similarities = {"NC_045512.2": [1.0]}

    existing = discover_run._existing_labels(_Result(), collection, "query")
    assert "NC_045512" in existing
    assert "NC_045512.2" in existing


def test_seeding_with_no_hits_at_all_is_not_called_saturated(monkeypatch) -> None:
    """Every search failing is not evidence that only siblings exist."""
    monkeypatch.setattr(iterate, "_blast_or_none", lambda *a, **k: [])
    params = FillParams(query=Path("q.fasta"), collection=None, output=Path("out"))
    selected, saturated = iterate._seed_windowed(
        "ACGT" * 2000, params, set(), _LOG, drop_siblings=True
    )
    assert selected == []
    assert saturated is False


def test_seeding_that_finds_only_siblings_is_saturated(monkeypatch) -> None:
    sibling = Hit("SIB1", "sibling", 97.0, 99.0, 0.0)  # a sibling, below the self-hit bar
    monkeypatch.setattr(iterate, "_blast_or_none", lambda *a, **k: [sibling])
    params = FillParams(query=Path("q.fasta"), collection=None, output=Path("out"))
    _, saturated = iterate._seed_windowed(
        "ACGT" * 2000, params, set(), _LOG, drop_siblings=True
    )
    assert saturated is True


# --- CLI boundary -----------------------------------------------------------

@pytest.mark.parametrize(
    "args",
    [
        ["detect", "-q", "nope.fasta", "-o", "out"],
        ["build-panel", "-q", "nope.fasta", "-o", "out"],
        ["fill-references", "-q", "nope.fasta", "-o", "out"],
    ],
)
def test_a_missing_query_is_reported_as_missing(tmp_path: Path, monkeypatch, args) -> None:
    monkeypatch.chdir(tmp_path)
    result = runner.invoke(app, args)
    assert result.exit_code == 1
    assert "Query file not found" in result.output
    assert "Unexpected error" not in result.output


def test_a_malformed_dataset_override_is_a_user_error(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    (tmp_path / "q.fasta").write_text(">HA\nACGT\n")
    result = runner.invoke(
        app, ["reassort", "-q", "q.fasta", "-o", "out", "--dataset", "no-equals-sign"]
    )
    assert result.exit_code == 1
    assert "--dataset must be SEGMENT=path" in result.output
    assert "Unexpected error" not in result.output


# --- provenance travels with the published panel -----------------------------

def test_the_published_panel_alignment_keeps_its_provenance(monkeypatch, tmp_path: Path) -> None:
    """`panel.msa.fasta` is a renamed copy; its sidecar must be copied with it."""
    from tessera.msa.build import provenance_path

    class _Result:
        similarities = {"refA": [1.0]}

    def fake_build_msa(params, logger):
        params.output.write_text(">q\nACGT\n")
        provenance_path(params.output).write_text('{"aligner": "mafft"}\n')

    monkeypatch.setattr(iterate, "efetch_available", lambda: True)
    monkeypatch.setattr(iterate, "build_msa", fake_build_msa)
    monkeypatch.setattr(iterate, "compute_similarity", lambda *a, **k: _Result())
    monkeypatch.setattr(iterate, "run_recomb", lambda p, logger, **kw: None)
    monkeypatch.setattr(iterate, "call_coverage_gaps", lambda *a, **k: ([], 0.95))

    collection = _collection(tmp_path / "coll")
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "ACGT" * 100 + "\n")
    out = tmp_path / "out"
    fill_references(FillParams(query=query, collection=collection, output=out), _LOG)

    assert (out / "panel.msa.fasta").exists()
    assert provenance_path(out / "panel.msa.fasta").read_text() == '{"aligner": "mafft"}\n'


# --- malformed records and tool output ---------------------------------------

def test_read_fasta_drops_whitespace_inside_sequence_lines(tmp_path: Path) -> None:
    """A trailing space or tab is not a base. Kept, it widened the backbone row and made
    the mafft MSA ragged (reference 6100 columns, query 6000)."""
    path = tmp_path / "ws.fasta"
    path.write_text(">ref desc\nACGT \nAC GT\t\n\nTTAA\n")
    assert read_fasta(path) == [("ref", "ACGTACGTTTAA")]


def test_read_fasta_accepts_a_header_with_no_name(tmp_path: Path) -> None:
    # ">" alone already read as an unnamed record; "> " raised IndexError instead.
    path = tmp_path / "blank.fasta"
    path.write_text("> \nACGT\n>\nTTAA\n")
    assert read_fasta(path) == [("", "ACGT"), ("", "TTAA")]


def test_sibeliaz_rejects_a_record_with_no_sequence_id(tmp_path: Path) -> None:
    path = tmp_path / "g1.fasta"
    path.write_text(">ok\nACGT\n> \nACGT\n")
    with pytest.raises(UserInputError, match=r"g1\.fasta.*line 3"):
        sibeliaz._build_seqid_map([path])


def test_truncated_maf_row_is_reported_against_the_file(tmp_path: Path) -> None:
    from tessera.converters.maf_to_fasta import maf_to_fasta

    maf = tmp_path / "cut.maf"
    maf.write_text("a\ns ref.c 0 8 + 8 ACGTACGT\ns qry.x 0 8 + 8 ACGT")
    with pytest.raises(OutputError, match=r"cut\.maf.*qry\.x"):
        maf_to_fasta(maf, "ref", tmp_path / "msa.fasta")


def test_xmfa_without_the_reference_is_reported_against_the_file(tmp_path: Path) -> None:
    from tessera.converters.xmfa_to_fasta import xmfa_to_fasta

    xmfa = tmp_path / "other.xmfa"
    xmfa.write_text(
        "#Sequence1File\t/data/a.fasta\n#Sequence2File\t/data/b.fasta\n"
        "> 1:1-4 + /data/a.fasta\nACGT\n> 2:1-4 + /data/b.fasta\nACGT\n=\n"
    )
    with pytest.raises(UserInputError, match=r"other\.xmfa.*/data/ref\.fasta"):
        xmfa_to_fasta(xmfa, "/data/ref.fasta", 0, tmp_path / "out.fasta")


def test_empty_mafft_output_is_a_tessera_error(tmp_path: Path) -> None:
    from tessera.converters.mafft_merge import merge_added_fragments

    empty = tmp_path / "empty.aln.fasta"
    empty.write_text("")
    with pytest.raises(OutputError, match=r"empty\.aln\.fasta"):
        merge_added_fragments(empty)


# --- the working copy is only cleared when Tessera made it ---------------------------

def test_copy_collection_refuses_to_replace_a_directory_it_did_not_create(tmp_path: Path) -> None:
    source = tmp_path / "src"
    source.mkdir()
    (source / "a.fasta").write_text(">a\nACGT\n")
    dest = tmp_path / "project" / "collection"
    dest.mkdir(parents=True)
    (dest / "mine.fasta").write_text(">mine\nACGT\n")

    with pytest.raises(UserInputError, match="was not created by Tessera"):
        copy_collection(source, dest)

    assert (dest / "mine.fasta").exists()


def test_copy_collection_replaces_the_working_copy_of_an_older_release(tmp_path: Path) -> None:
    # Output directories written before the marker existed are recognised by the files a
    # run leaves beside the working copy, so re-running into one still works.
    source = tmp_path / "src"
    source.mkdir()
    (source / "a.fasta").write_text(">a\nACGT\n")
    out = tmp_path / "out"
    (out / "collection").mkdir(parents=True)
    (out / "collection" / "old.fasta").write_text(">old\nACGT\n")
    (out / "round1.msa.fasta").write_text(">q\nACGT\n")

    copy_collection(source, out / "collection")

    assert [p.name for p in collection_genomes(out / "collection")] == ["a.fasta"]


def test_copy_collection_accepts_an_empty_existing_directory(tmp_path: Path) -> None:
    source = tmp_path / "src"
    source.mkdir()
    (source / "a.fasta").write_text(">a\nACGT\n")
    dest = tmp_path / "out" / "collection"
    dest.mkdir(parents=True)

    copy_collection(source, dest)

    assert [p.name for p in collection_genomes(dest)] == ["a.fasta"]
