"""Panel curation: sibling exclusion, backbone choice, dereplication (skani/skDER mocked)."""

from __future__ import annotations

from pathlib import Path

from recomfi.discover import panel


def _genomes(tmp_path: Path, names: list[str]) -> list[Path]:
    paths = []
    for n in names:
        p = tmp_path / f"{n}.fasta"
        p.write_text(f">{n}\nACGT\n")
        paths.append(p)
    return paths


# --- skani output parsing ------------------------------------------------------

def test_skani_query_ani_parses_and_marks_absent(monkeypatch, tmp_path, logger):
    backbone, cand, absent = _genomes(tmp_path, ["A1", "C", "AE"])

    def fake_run_tool(caps, command, *, logger, log_prefix="", **kw):
        out = Path(command[command.index("-o") + 1])
        # Ref_file, Query_file, ANI, AF_ref, AF_query, ...  (absent ref omitted)
        lines = ["Ref_file\tQuery_file\tANI\tAF_ref\tAF_query\tnames"]
        lines.append(f"{backbone.resolve()}\tq\t92.5\t83\t81.0\tx")
        lines.append(f"{cand.resolve()}\tq\t89.7\t95\t93.6\tx")
        out.write_text("\n".join(lines) + "\n")
        return ""

    monkeypatch.setattr(panel, "run_tool", fake_run_tool)
    ani = panel.skani_query_ani(tmp_path / "q.fasta", [backbone, cand, absent], logger)
    assert ani[backbone] == (92.5, 81.0)
    assert ani[cand] == (89.7, 93.6)
    assert ani[absent] == (0.0, 0.0)  # not in skani output -> reported absent


# --- sibling exclusion (the masking fix) --------------------------------------

def test_filter_siblings_drops_relatives_keeps_env_and_parent(monkeypatch, tmp_path, logger):
    backbone, sibling, env, parent = _genomes(tmp_path, ["A1", "AE_rel", "AE_env", "C"])
    ani = {
        backbone: (92.0, 81.0),   # the query's whole-genome anchor
        sibling: (97.0, 94.0),    # high ANI AND high coverage -> sibling (drop)
        env: (96.0, 25.0),        # high ANI but regional (low coverage) -> keep
        parent: (89.0, 94.0),     # lower ANI than backbone -> keep
    }
    monkeypatch.setattr(panel, "skani_query_ani", lambda *a, **k: ani)
    parents, siblings, _ = panel.filter_siblings(
        tmp_path / "q.fasta", backbone, [sibling, env, parent], logger=logger,
    )
    assert siblings == [sibling]
    assert set(parents) == {env, parent}


def test_filter_siblings_catches_whole_genome_twin_of_partial_backbone(
    monkeypatch, tmp_path, logger
):
    # Backbone is a partial-coverage parent (high ANI, but covers only one ORF). The
    # sibling has ~equal ANI yet covers the whole query -> caught by coverage, not ANI.
    backbone, twin, minor = _genomes(tmp_path, ["ORF1_parent", "recombinant_twin", "ORF2_parent"])
    ani = {
        backbone: (98.5, 69.0),   # high ANI over its donated ORF only
        twin: (98.0, 96.0),       # ~equal ANI, whole-query coverage -> sibling
        minor: (84.0, 27.0),      # the other parent (partial coverage) -> keep
    }
    monkeypatch.setattr(panel, "skani_query_ani", lambda *a, **k: ani)
    parents, siblings, _ = panel.filter_siblings(
        tmp_path / "q.fasta", backbone, [twin, minor], logger=logger,
    )
    assert siblings == [twin]
    assert parents == [minor]


def test_pick_backbone_prefers_highest_ani_among_full_coverage(monkeypatch, tmp_path, logger):
    a1, c, env = _genomes(tmp_path, ["A1", "C", "AE_env"])
    ani = {
        a1: (92.6, 81.0),   # highest ANI among full-coverage -> backbone
        c: (89.7, 94.0),
        env: (96.2, 25.0),  # higher ANI but regional (AF below floor) -> not backbone
    }
    monkeypatch.setattr(panel, "skani_query_ani", lambda *a, **k: ani)
    assert panel.pick_backbone(tmp_path / "q.fasta", [a1, c, env], logger=logger) == a1


# --- dereplication -------------------------------------------------------------

def test_dereplicate_splits_representatives_and_redundant(monkeypatch, tmp_path, logger):
    g1, g2, g3 = _genomes(tmp_path, ["g1", "g2", "g3"])

    def fake_run_tool(caps, command, *, logger, log_prefix="", **kw):
        out = Path(command[command.index("-o") + 1])
        out.mkdir(parents=True, exist_ok=True)
        # g3 is a near-duplicate of g1; skDER keeps g1 and g2.
        (out / "skDER_Results.txt").write_text(f"{g1.resolve()}\n{g2.resolve()}\n")
        return ""

    monkeypatch.setattr(panel, "run_tool", fake_run_tool)
    reps, redundant = panel.dereplicate([g1, g2, g3], logger=logger)
    assert set(reps) == {g1, g2}
    assert redundant == [g3]


def test_dereplicate_noop_for_single_genome(tmp_path, logger):
    (g1,) = _genomes(tmp_path, ["g1"])
    reps, redundant = panel.dereplicate([g1], logger=logger)
    assert reps == [g1]
    assert redundant == []


# --- end-to-end curation table -------------------------------------------------

def test_curate_panel_table_roles(monkeypatch, tmp_path, logger):
    backbone, sibling, parent = _genomes(tmp_path, ["A1", "AE_rel", "C"])
    ani = {backbone: (92.0, 81.0), sibling: (97.0, 94.0), parent: (89.0, 94.0)}
    monkeypatch.setattr(panel, "skani_available", lambda: True)
    monkeypatch.setattr(panel, "skder_available", lambda: False)  # sibling filter only
    monkeypatch.setattr(panel, "skani_query_ani", lambda *a, **k: ani)

    result = panel.curate_panel(
        tmp_path / "q.fasta", [backbone, sibling, parent], backbone, logger=logger,
    )
    assert set(result.kept) == {backbone, parent}
    assert result.siblings == [sibling]
    roles = {r["genome"]: r["role"] for r in result.table}
    assert roles == {"A1": "backbone", "AE_rel": "sibling-dropped", "C": "representative"}


# --- typed Lineage column (conditional) ----------------------------------------

def _panel_row(genome: str, role: str = "representative") -> dict:
    return {"genome": genome, "role": role, "query_ani": 95.0, "query_af": 90.0}


def test_panel_tsv_no_map_keeps_four_columns(tmp_path: Path, logger):
    path = tmp_path / "panel.tsv"
    panel.write_panel_tsv(path, [_panel_row("MK573073")], logger)
    header = path.read_text().splitlines()[0]
    assert header == "genome\trole\tquery_ani\tquery_af"


def test_panel_tsv_with_map_adds_lineage_column(tmp_path: Path, logger):
    path = tmp_path / "panel.tsv"
    panel.write_panel_tsv(path, [_panel_row("MK573073")], logger, {"MK573073": "GII.P16-GII.4"})
    lines = path.read_text().splitlines()
    assert lines[0] == "genome\trole\tlineage\tquery_ani\tquery_af"
    assert lines[1].split("\t")[2] == "GII.P16-GII.4"


def test_panel_html_lineage_column_conditional():
    rows = [_panel_row("MK573073")]
    assert "Lineage" not in panel.panel_table_html(rows)
    html = panel.panel_table_html(rows, {"MK573073": "GII.P16-GII.4"})
    assert "Lineage" in html and "GII.P16-GII.4" in html
