"""Iterative fill-references loop: rounds, convergence, and the run summary."""

from __future__ import annotations

from pathlib import Path

from recomfi.discover import iterate
from recomfi.discover.blast import Hit
from recomfi.discover.iterate import FillParams, fill_references
from recomfi.discover.run import Candidate
from recomfi.recomb.coverage import CoverageGap


class _StubResult:
    similarities = {"refA": [1.0]}


def _gap(mean_best: float) -> CoverageGap:
    return CoverageGap(
        msa_start=100, msa_end=600, query_start=100, query_end=600,
        length_bp=500, n_windows=5, best_label="refA", mean_best=mean_best, kind="divergent",
    )


def _setup(tmp_path: Path) -> tuple[Path, Path, Path]:
    coll = tmp_path / "coll"
    coll.mkdir()
    (coll / "refA.fasta").write_text(">refA\nACGT\n")
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "ACGT" * 100 + "\n")
    return query, coll, tmp_path / "out"


def _common_mocks(monkeypatch, coverage_returns):
    monkeypatch.setattr(iterate, "efetch_available", lambda: True)
    monkeypatch.setattr(iterate, "build_msa", lambda p, logger: p.output.write_text(">q\nACGT\n"))
    monkeypatch.setattr(iterate, "compute_similarity", lambda *a, **k: _StubResult())
    monkeypatch.setattr(iterate, "read_fasta", lambda path: [("q", "ACGT" * 100)])
    monkeypatch.setattr(iterate, "run_recomb", lambda p, logger, **kw: None)
    returns = list(coverage_returns)
    monkeypatch.setattr(iterate, "call_coverage_gaps", lambda *a, **k: returns.pop(0))


def test_loop_converges_when_gaps_close(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    # round 1 has a gap; after adding a reference, round 2 has none.
    _common_mocks(monkeypatch, [([_gap(0.84)], 0.94), ([], 0.95)])
    monkeypatch.setattr(
        iterate, "collect_candidates",
        lambda *a, **k: [Candidate(_gap(0.84), Hit("NEW1", "new", 90.0, 95.0, 1e-9), False)],
    )

    def fake_download(cands, dest, logger):
        (dest / "NEW1.fasta").write_text(">NEW1\nACGT\n")
        return cands

    monkeypatch.setattr(iterate, "_download", fake_download)

    trace = fill_references(FillParams(query=query, collection=coll, output=out), logger)

    assert [r.round for r in trace] == [1, 2]
    assert trace[0].added == ["NEW1"]
    assert trace[1].n_gaps == 0  # converged
    assert (out / "collection" / "NEW1.fasta").exists()  # grown copy, original untouched
    assert not (coll / "NEW1.fasta").exists()
    summary = (out / "fill_summary.tsv").read_text().splitlines()
    assert summary[0].split("\t")[0] == "round"
    assert len(summary) == 3  # header + 2 rounds


def test_query_own_accession_is_auto_excluded(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    # query FASTA header carries its own accession
    query.write_text(">MG572182.1 Norovirus GII\n" + "ACGT" * 100 + "\n")
    _common_mocks(monkeypatch, [([_gap(0.84)], 0.94), ([], 0.95)])

    # the query FASTA reports its accession header; the MSA reports the query label
    def fake_read(path):
        if str(path).endswith("q.fasta"):
            return [("MG572182.1", "ACGT" * 100)]
        return [("q", "ACGT" * 100)]

    monkeypatch.setattr(iterate, "read_fasta", fake_read)

    seen: dict[str, set] = {}

    def capture(targets, query_row, existing, *, exclude, **kw):
        seen["exclude"] = exclude
        return []  # no candidates -> loop stops after round 1

    monkeypatch.setattr(iterate, "collect_candidates", capture)
    monkeypatch.setattr(iterate, "_download", lambda c, d, logger: [])

    fill_references(FillParams(query=query, collection=coll, output=out), logger)
    assert "MG572182" in seen["exclude"]  # auto-excluded from its header


def test_fresh_start_seeds_collection_from_whole_query_blast(monkeypatch, tmp_path, logger):
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "ACGT" * 100 + "\n")
    out = tmp_path / "out"
    # No starting collection -> seed from a whole-query BLAST, then converge.
    _common_mocks(monkeypatch, [([_gap(0.84)], 0.94), ([], 0.95)])
    monkeypatch.setattr(
        iterate, "collect_candidates",
        lambda *a, **k: [Candidate(_gap(0.84), Hit("NEW1", "new", 90.0, 95.0, 1e-9), False)],
    )
    monkeypatch.setattr(iterate, "_download", lambda c, d, logger: [])

    blasted: dict[str, str] = {}

    def fake_blast(seq, *, max_hits, logger, email=None, entrez_query=None):
        blasted["seq"] = seq
        blasted["max_hits"] = max_hits
        return [
            Hit("SELF", "the query itself", 99.9, 99.0, 0.0),  # auto-skipped as self-hit
            Hit("SEED1", "a relative", 92.0, 95.0, 1e-30),
            Hit("SEED2", "another relative", 90.0, 90.0, 1e-20),
        ]

    def fake_efetch(accession, dest, logger):
        path = dest / f"{accession}.fasta"
        path.write_text(f">{accession}\nACGT\n")
        return path

    monkeypatch.setattr(iterate, "blast_subsequence", fake_blast)
    monkeypatch.setattr(iterate, "efetch_fasta", fake_efetch)

    fill_references(FillParams(query=query, collection=None, output=out, seed_hits=7), logger)

    assert blasted["max_hits"] == 7
    assert "ACGTACGT" in blasted["seq"]  # the whole (de-gapped) query
    coll = out / "collection"
    assert (coll / "SEED1.fasta").exists()
    assert (coll / "SEED2.fasta").exists()
    assert not (coll / "SELF.fasta").exists()  # near-identical self-hit not seeded


def _seed_run(monkeypatch, tmp_path, logger, *, seed_mode, fake_blast):
    """Run a fresh-start fill that converges immediately, capturing the seeded files."""
    query = tmp_path / "q.fasta"
    # three distinct 100 bp segments -> three seed windows at seed_window=100
    query.write_text(">q\n" + "A" * 100 + "C" * 100 + "G" * 100 + "\n")
    out = tmp_path / "out"
    _common_mocks(monkeypatch, [([], 0.95)])  # round 1 finds no gaps -> converge
    monkeypatch.setattr(iterate, "read_fasta", lambda p: [("q", "A" * 100 + "C" * 100 + "G" * 100)])
    monkeypatch.setattr(iterate, "blast_subsequence", fake_blast)
    monkeypatch.setattr(
        iterate, "efetch_fasta",
        lambda acc, dest, logger: (dest / f"{acc}.fasta").write_text(f">{acc}\nA\n")
        or (dest / f"{acc}.fasta"),
    )
    fill_references(
        FillParams(
            query=query, collection=None, output=out,
            seed_mode=seed_mode, seed_window=100, seed_hits=5, auto_diversify=False,
        ),
        logger,
    )
    return {p.stem for p in (out / "collection").glob("*.fasta")}


def test_saturation_auto_switches_to_ncbi_virus_diversity(monkeypatch, tmp_path, logger):
    # All BLAST hits are siblings -> with auto_diversify, seeding switches to the
    # NCBI Virus taxonomy-diversity path instead of seeding the siblings.
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "A" * 300 + "\n")
    out = tmp_path / "out"
    _common_mocks(monkeypatch, [([], 0.95)])
    monkeypatch.setattr(iterate, "read_fasta", lambda p: [("q", "A" * 300)])
    monkeypatch.setattr(
        iterate, "blast_subsequence",
        lambda seq, **k: [Hit("SIB", "sibling", 98.0, 99.0, 0.0)],
    )
    from recomfi.discover import pool as pool_mod
    monkeypatch.setattr(pool_mod, "datasets_available", lambda: True)

    called = {}

    def fake_from_pool(params, collection, logger, *, force_ncbi=False):
        called["force_ncbi"] = force_ncbi
        (collection / "DIVERSE.fasta").write_text(">DIVERSE\nA\n")

    monkeypatch.setattr(iterate, "_seed_from_pool", fake_from_pool)
    fill_references(
        FillParams(query=query, collection=None, output=out, seed_window=300, auto_diversify=True),
        logger,
    )
    assert called.get("force_ncbi") is True
    assert (out / "collection" / "DIVERSE.fasta").exists()


def _regional_blast(seq, *, max_hits, logger, email=None, entrez_query=None):
    # one sibling present in every window (near-identical, full coverage) + a distinct
    # regional parent per window (lower identity).
    sib = Hit("SIB", "sibling", 98.0, 99.0, 0.0)
    region = {"A": Hit("PAR_A", "parent A", 91.0, 80.0, 1e-9),
              "C": Hit("PAR_C", "parent C", 90.0, 82.0, 1e-9),
              "G": Hit("PAR_G", "parent G", 89.0, 81.0, 1e-9)}[seq[0]]
    return [sib, region]


def test_parents_mode_suppresses_siblings_and_seeds_parents(monkeypatch, tmp_path, logger):
    seeded = _seed_run(
        monkeypatch, tmp_path, logger, seed_mode="parents", fake_blast=_regional_blast
    )
    assert seeded == {"PAR_A", "PAR_C", "PAR_G"}  # the sibling is dropped


def test_windowed_mode_keeps_per_window_best_including_siblings(monkeypatch, tmp_path, logger):
    seeded = _seed_run(
        monkeypatch, tmp_path, logger, seed_mode="windowed", fake_blast=_regional_blast
    )
    assert "SIB" in seeded and {"PAR_A", "PAR_C", "PAR_G"} <= seeded


def test_parents_mode_falls_back_when_only_siblings(monkeypatch, tmp_path, logger):
    # every window returns only a sibling -> nothing to suppress down to -> seed the best.
    def only_siblings(seq, *, max_hits, logger, email=None, entrez_query=None):
        return [Hit("SIB", "sibling", 98.0, 99.0, 0.0)]

    seeded = _seed_run(monkeypatch, tmp_path, logger, seed_mode="parents", fake_blast=only_siblings)
    assert seeded == {"SIB"}  # fallback keeps the collection non-empty


def test_fetch_diverse_broadens_caps_and_caches(monkeypatch, tmp_path, logger):
    # RefSeq too thin -> broaden to a --limit-capped complete fetch; result is cached.
    from recomfi.discover import pool as pool_mod

    calls = {"fetch": 0, "limit": None}

    def fake_fetch(taxon, d, *, refseq=True, complete_only=False, released_after=None,
                   limit=None, logger):
        calls["fetch"] += 1
        if refseq:
            p = d / "NC_1.fasta"  # 1 < SEED_MIN_DIVERSE -> broaden
            p.write_text(">x\nA\n")
            return [p]
        calls["limit"] = limit
        out = []
        for i in range(limit):  # cap hit -> warning path
            p = d / f"G{i}.fasta"
            p.write_text(">x\nA\n")
            out.append(p)
        return out

    monkeypatch.setattr(pool_mod, "fetch_ncbi_virus", fake_fetch)
    params = FillParams(
        query=tmp_path / "q.fasta", collection=None, output=tmp_path / "o",
        taxon="SARS-CoV-2", fetch_limit=5, cache_dir=tmp_path / "cache",
    )
    result = iterate._fetch_diverse(params, logger)
    assert calls["limit"] == 5
    assert len(result) == 5  # broadened to the capped complete set

    # A second call for the same taxon hits the cache -- no further fetch.
    before = calls["fetch"]
    cached = iterate._fetch_diverse(params, logger)
    assert calls["fetch"] == before  # network skipped
    assert len(cached) == 5


def test_dominant_lineage_token_extracted_from_titles():
    from recomfi.discover.iterate import _dominant_lineage_token

    hits = [
        Hit("A", "Norovirus GII isolate Hu/GII.P16-GII.1/RUS/NS18", 98.0, 99.0, 0.0),
        Hit("B", "Norovirus GII isolate Hu/GII.P16-GII.1/JP/Yuzawa", 98.0, 99.0, 0.0),
        Hit("C", "Norovirus GII strain GII.P16-GII.1 clone X", 98.0, 99.0, 0.0),
    ]
    assert _dominant_lineage_token(hits) == "GII.P16-GII.1"
    # No digit-bearing token shared across hits -> nothing to exclude.
    assert _dominant_lineage_token(
        [Hit("A", "Some virus strain ABC", 98.0, 99.0, 0.0),
         Hit("B", "Other isolate from host", 98.0, 99.0, 0.0)]
    ) is None


def test_negative_lineage_seeding_recruits_parents(monkeypatch, tmp_path, logger):
    from recomfi.discover import iterate as it

    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "ACGT" * 200 + "\n")
    captured = {}

    def fake_blast(seq, *, max_hits, logger, email=None, entrez_query=None):
        if entrez_query is None:  # the whole-query probe -> the saturating lineage
            return [Hit("SIB1", "Norovirus GII isolate Hu/GII.P16-GII.1/A", 98.0, 99.0, 0.0),
                    Hit("SIB2", "Norovirus GII isolate Hu/GII.P16-GII.1/B", 98.0, 99.0, 0.0)]
        captured["entrez_query"] = entrez_query  # per-region negative search
        return [Hit("PARENT", "Norovirus GII.P16 polymerase", 91.0, 80.0, 1e-9),
                Hit("SIBX", "Norovirus GII.P16-GII.1 again", 98.0, 99.0, 0.0)]

    monkeypatch.setattr(it, "blast_subsequence", fake_blast)
    params = FillParams(query=query, collection=None, output=tmp_path / "o", seed_window=400)
    seeds = it._seed_negative_lineage("ACGT" * 200, params, set(), logger)

    assert "NOT \"GII.P16-GII.1\"" in captured["entrez_query"]
    assert "Norovirus GII"[:10] in captured["entrez_query"]  # organism restriction
    assert seeds == ["PARENT"]  # the divergent parent kept; the residual sibling dropped


def test_loop_stops_when_coverage_stalls(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    # the best reference stays at 0.84 both rounds -> no improvement -> stop
    _common_mocks(monkeypatch, [([_gap(0.84)], 0.94), ([_gap(0.84)], 0.94)])
    monkeypatch.setattr(
        iterate, "collect_candidates",
        lambda *a, **k: [Candidate(_gap(0.84), Hit("NEW1", "new", 90.0, 95.0, 1e-9), False)],
    )

    def fake_download(cands, dest, logger):
        (dest / "NEW1.fasta").write_text(">NEW1\nACGT\n")
        return cands

    monkeypatch.setattr(iterate, "_download", fake_download)

    trace = fill_references(
        FillParams(query=query, collection=coll, output=out, min_improvement=0.01), logger
    )
    assert [r.round for r in trace] == [1, 2]
    assert trace[1].added == []  # round 2 stalled, no second download


def test_no_report_skips_detection_but_writes_panel(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    _common_mocks(monkeypatch, [([], 0.95)])  # converge immediately, no downloads
    calls = []
    monkeypatch.setattr(iterate, "run_recomb", lambda p, logger, **kw: calls.append(p))

    fill_references(FillParams(query=query, collection=coll, output=out, report=False), logger)

    assert calls == []  # detection not run with --no-report
    assert (out / "panel.msa.fasta").exists()  # stable panel alignment still published


def test_report_runs_detection_on_stable_panel(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    _common_mocks(monkeypatch, [([], 0.95)])
    calls = []
    monkeypatch.setattr(iterate, "run_recomb", lambda p, logger, **kw: calls.append(p))

    fill_references(FillParams(query=query, collection=coll, output=out), logger)

    assert len(calls) == 1  # detection runs by default
    assert calls[0].msa == out / "panel.msa.fasta"  # consumes the stable copy, not round{N}


def test_capture_writes_typed_lineage_sidecar(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    # a collection genome whose header carries a genotype to be mined
    (coll / "MK573073.fasta").write_text(
        ">MK573073.1 Norovirus GII.P16-GII.4 isolate Hu\nACGT\n"
    )
    _common_mocks(monkeypatch, [([], 0.95)])  # converge immediately, no downloads

    fill_references(FillParams(query=query, collection=coll, output=out, report=False), logger)

    sidecar = (out / "lineages.tsv").read_text()
    assert "MK573073" in sidecar and "GII.P16-GII.4" in sidecar


def test_query_self_typing_writes_query_row(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    # the query's own header carries a genotype that is not its file name
    query.write_text(">MG572182.1 Norovirus GII.P16-GII.1 isolate Hu\n" + "ACGT" * 100 + "\n")
    _common_mocks(monkeypatch, [([], 0.95)])

    fill_references(FillParams(query=query, collection=coll, output=out, report=False), logger)

    rows = (out / "lineages.tsv").read_text().splitlines()
    query_label = query.stem  # "q"
    assert any(r.startswith(f"{query_label}\tGII.P16-GII.1\tquery") for r in rows)


def test_pango_crosscheck_section_for_recombinant_query(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    query.write_text(
        ">OM.1 Severe acute respiratory syndrome coronavirus 2 XBB.1.5\n" + "ACGT" * 100 + "\n"
    )
    _common_mocks(monkeypatch, [([], 0.95)])
    monkeypatch.setattr(iterate, "load_alias_key", lambda **k: {"XBB": ["BJ.1", "CJ.1"]})
    captured = {}
    monkeypatch.setattr(
        iterate, "run_recomb",
        lambda p, logger, **kw: captured.update(kw),
    )

    fill_references(
        FillParams(query=query, collection=coll, output=out, taxon="SARS-CoV-2"), logger
    )

    titles = [t for t, _ in captured["extra_sections"]]
    assert "Pango cross-check" in titles
    body = dict(captured["extra_sections"])["Pango cross-check"]
    assert "BJ.1" in body and "XBB.1.5" in body
