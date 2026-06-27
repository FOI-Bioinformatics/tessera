"""Pool sourcing: regional selection, sibling drop, and NCBI-Virus package parsing."""

from __future__ import annotations

from pathlib import Path

from recomfi.discover import pool


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
