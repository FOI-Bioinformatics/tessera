"""Pool sourcing: regional selection, sibling drop, and NCBI-Virus package parsing."""

from __future__ import annotations

from pathlib import Path

from tessera.discover import pool


def _genomes(tmp_path: Path, names: list[str]) -> list[Path]:
    out = []
    for n in names:
        p = tmp_path / f"{n}.fasta"
        p.write_text(f">{n}\nACGT\n")
        out.append(p)
    return out


def _fake_ani(by_stem_whole, by_stem_window):
    """skani mock: dispatch on whether the query path is a window file."""

    def ani(qpath, refs, logger):
        name = Path(qpath).name
        which = by_stem_window if name.startswith("win_") else by_stem_whole
        win = name.split(".")[0] if name.startswith("win_") else None
        out = {}
        for r in refs:
            val = which[r.stem]
            out[r] = val[win] if isinstance(val, dict) else val
        return out

    return ani


def test_select_regional_drops_siblings_and_picks_regional_parents(monkeypatch, tmp_path, logger):
    sib, par_a, par_b = _genomes(tmp_path, ["SIB", "PAR_A", "PAR_B"])
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "A" * 3000 + "\n")  # 2 windows at window=1500

    whole = {"SIB": (98.0, 95.0), "PAR_A": (90.0, 85.0), "PAR_B": (88.0, 84.0)}
    # PAR_A wins window 0, PAR_B wins window 1.
    window = {
        "PAR_A": {"win_0": (96.0, 80.0), "win_1": (85.0, 80.0)},
        "PAR_B": {"win_0": (85.0, 80.0), "win_1": (96.0, 80.0)},
    }
    monkeypatch.setattr(pool, "skani_available", lambda: True)
    monkeypatch.setattr(pool, "dereplicate", lambda g, ani, logger: (list(g), []))
    monkeypatch.setattr(pool, "skani_query_ani", _fake_ani(whole, window))

    sel = pool.select_regional(query, [sib, par_a, par_b], window=1500, per_window=1, logger=logger)
    assert {p.stem for p in sel.selected} == {"PAR_A", "PAR_B"}  # sibling dropped


def test_select_regional_keeps_siblings_when_disabled(monkeypatch, tmp_path, logger):
    sib, par_a = _genomes(tmp_path, ["SIB", "PAR_A"])
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "A" * 1500 + "\n")
    whole = {"SIB": (98.0, 95.0), "PAR_A": (90.0, 85.0)}
    window = {"SIB": {"win_0": (98.0, 80.0)}, "PAR_A": {"win_0": (90.0, 80.0)}}
    monkeypatch.setattr(pool, "skani_available", lambda: True)
    monkeypatch.setattr(pool, "dereplicate", lambda g, ani, logger: (list(g), []))
    monkeypatch.setattr(pool, "skani_query_ani", _fake_ani(whole, window))

    sel = pool.select_regional(
        query, [sib, par_a], window=1500, per_window=2, drop_siblings=False, logger=logger
    )
    assert "SIB" in {p.stem for p in sel.selected}


def test_split_fasta_writes_one_file_per_accession_with_lineage(tmp_path):
    fna = tmp_path / "genomic.fna"
    fna.write_text(">ACC1.1 some virus\nACGT\n>ACC2.1 other\nTTTT\n")
    dest = tmp_path / "out"
    dest.mkdir()
    written = pool._split_fasta(fna, dest, {"ACC1.1": "GII.4", "ACC2.1": "GII.17"})
    assert {p.name for p in written} == {"ACC1.1.fasta", "ACC2.1.fasta"}
    assert "GII.4" in (dest / "ACC1.1.fasta").read_text()


def test_scope_flags_builds_datasets_filters():
    assert pool._scope_flags(refseq=True, complete_only=False, released_after=None) == ["--refseq"]
    assert pool._scope_flags(False, True, "2024-01-01") == [
        "--complete-only", "--released-after", "2024-01-01"
    ]


def test_lineage_labels_parses_isolate_and_pangolin(tmp_path):
    report = tmp_path / "data_report.jsonl"
    report.write_text(
        '{"accession": "A1", "isolate": {"lineage": "GII.4"}}\n'
        '{"accession": "A2", "virus": {"pangolinClassification": "XBB.1"}}\n'
        "\n"  # blank line tolerated
        '{"accession": "A3"}\n'  # no lineage -> omitted
    )
    labels = pool._lineage_labels(report)
    assert labels == {"A1": "GII.4", "A2": "XBB.1"}


def _patch_selection(monkeypatch, whole, window):
    monkeypatch.setattr(pool, "skani_available", lambda: True)
    monkeypatch.setattr(pool, "dereplicate", lambda g, ani, logger: (list(g), []))
    monkeypatch.setattr(pool, "skani_query_ani", _fake_ani(whole, window))


# --- audit 4: an all-sibling pool is reported (D7) ---------------------------------

def test_select_regional_warns_when_every_genome_is_a_sibling(monkeypatch, tmp_path, caplog):
    """Four genomes at 99 % ANI over 99 % coverage: the sibling drop would empty the pool,
    so it is skipped. That used to happen in silence and yield a normal-looking panel."""
    import logging

    names = ["S1", "S2", "S3", "S4"]
    genomes = _genomes(tmp_path, names)
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "A" * 3000 + "\n")
    _patch_selection(
        monkeypatch,
        {n: (99.0, 99.0) for n in names},
        {n: {"win_0": (99.0, 80.0), "win_1": (99.0, 80.0)} for n in names},
    )
    log = logging.getLogger("pool_test")
    with caplog.at_level(logging.WARNING, logger="pool_test"):
        sel = pool.select_regional(query, genomes, window=1500, per_window=2, logger=log)

    assert sel.selected  # still a panel: the siblings are all there is
    assert "Every pool genome is a sibling of the query" in caplog.text
    assert "no distinct parental lineage found in the pool" in caplog.text


# --- audit 4: a short query is still matched region by region (D9) -----------------

def test_select_regional_cuts_a_short_query_into_two_windows(monkeypatch, tmp_path, logger):
    """A 2.9 kb query at the default 1500 bp window was one window over the whole query,
    so a genome that matches only one half could not be recruited for that half."""
    par_a, par_b = _genomes(tmp_path, ["PAR_A", "PAR_B"])
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "A" * 2900 + "\n")
    _patch_selection(
        monkeypatch,
        {"PAR_A": (90.0, 85.0), "PAR_B": (88.0, 84.0)},
        {"PAR_A": {"win_0": (96.0, 80.0), "win_1": (85.0, 80.0)},
         "PAR_B": {"win_0": (85.0, 80.0), "win_1": (96.0, 80.0)}},
    )

    sel = pool.select_regional(query, [par_a, par_b], window=1500, per_window=1, logger=logger)

    assert {p.stem for p in sel.selected} == {"PAR_A", "PAR_B"}
    assert [(r["query_start"], r["query_end"]) for r in sel.table] == [(0, 1450), (1450, 2900)]


def test_select_regional_keeps_one_window_for_a_query_no_longer_than_it(
    monkeypatch, tmp_path, logger
):
    """Guard (passes before and after): nothing to split at or below one window."""
    par_a, par_b = _genomes(tmp_path, ["PAR_A", "PAR_B"])
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "A" * 1500 + "\n")
    _patch_selection(
        monkeypatch,
        {"PAR_A": (90.0, 85.0), "PAR_B": (88.0, 84.0)},
        {"PAR_A": {"win_0": (96.0, 80.0)}, "PAR_B": {"win_0": (85.0, 80.0)}},
    )

    sel = pool.select_regional(query, [par_a, par_b], window=1500, per_window=1, logger=logger)

    assert [(r["query_start"], r["query_end"]) for r in sel.table] == [(0, 1500)]


def test_select_regional_falls_back_to_the_whole_query_when_no_half_matches(
    monkeypatch, tmp_path, logger
):
    """skani may report nothing for a sub-window chunk (it is a genome-scale tool). The
    split must not then leave a short query with an empty panel where the single
    whole-query window used to select one: the whole-genome ranking is used instead."""
    par_a, par_b = _genomes(tmp_path, ["PAR_A", "PAR_B"])
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "A" * 2900 + "\n")
    _patch_selection(
        monkeypatch,
        {"PAR_A": (90.0, 85.0), "PAR_B": (88.0, 84.0)},
        {"PAR_A": {"win_0": (0.0, 0.0), "win_1": (0.0, 0.0)},
         "PAR_B": {"win_0": (0.0, 0.0), "win_1": (0.0, 0.0)}},
    )

    sel = pool.select_regional(query, [par_a, par_b], window=1500, per_window=1, logger=logger)

    assert [p.stem for p in sel.selected] == ["PAR_A"]


# --- audit 4: a candidate pool is read like a collection (E5) ----------------------

def test_candidate_pool_reads_the_same_files_as_a_collection(tmp_path):
    """The pool reader accepted four extensions where every collection reader takes each
    non-hidden file, so a .fas or .fsa genome was dropped from the pool in silence."""
    import gzip

    from tessera.core.io import collection_genomes

    for name in ("a.fasta", "b.fas", "c.fsa", "d"):
        (tmp_path / name).write_text(f">{name}\nACGT\n")
    with gzip.open(tmp_path / "e.fa.gz", "wt") as fo:
        fo.write(">e\nACGT\n")
    (tmp_path / ".DS_Store").write_bytes(b"\x00\x01")

    genomes = pool.iter_pool_genomes(tmp_path)

    assert genomes == collection_genomes(tmp_path)
    assert [p.name for p in genomes] == ["a.fasta", "b.fas", "c.fsa", "d", "e.fa.gz"]


def test_candidate_pool_refuses_a_file_that_is_not_fasta(tmp_path):
    import pytest

    from tessera.core.errors import UserInputError

    (tmp_path / "a.fasta").write_text(">a\nACGT\n")
    (tmp_path / "README.txt").write_text("pool assembled 2026-09\n")

    with pytest.raises(UserInputError, match="does not look like a FASTA file"):
        pool.iter_pool_genomes(tmp_path)
