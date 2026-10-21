"""Iterative fill-references loop: rounds, convergence, and the run summary."""

from __future__ import annotations

from pathlib import Path

import pytest

from tessera.core.errors import UserInputError
from tessera.discover import iterate
from tessera.discover.blast import Hit
from tessera.discover.iterate import FillParams, fill_references
from tessera.discover.run import Candidate
from tessera.recomb.coverage import CoverageGap


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

    def fake_blast(seq, *, max_hits, logger, email=None, entrez_query=None, cache_dir=None):
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
    from tessera.discover import pool as pool_mod
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


def _regional_blast(seq, *, max_hits, logger, email=None, entrez_query=None, cache_dir=None):
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
    def only_siblings(seq, *, max_hits, logger, email=None, entrez_query=None, cache_dir=None):
        return [Hit("SIB", "sibling", 98.0, 99.0, 0.0)]

    seeded = _seed_run(monkeypatch, tmp_path, logger, seed_mode="parents", fake_blast=only_siblings)
    assert seeded == {"SIB"}  # fallback keeps the collection non-empty


def test_fetch_diverse_broadens_downloads_all_and_caches(monkeypatch, tmp_path, logger):
    # RefSeq too thin -> broaden to the FULL complete set (no pre-cap; dereplicated later).
    from tessera.discover import pool as pool_mod

    calls = {"fetch": 0, "limit": "unset"}

    def fake_fetch(taxon, d, *, refseq=True, complete_only=False, released_after=None,
                   limit=None, logger):
        calls["fetch"] += 1
        if refseq:
            p = d / "NC_1.fasta"  # 1 < SEED_MIN_DIVERSE -> broaden
            p.write_text(">x\nA\n")
            return [p]
        calls["limit"] = limit
        out = []
        for i in range(8):  # the full complete set (more than fetch_limit), not truncated
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
    assert calls["limit"] is None  # the full set is downloaded, not pre-capped
    assert len(result) == 8  # all genomes returned for local dereplication

    # A second call for the same taxon hits the cache -- no further fetch.
    before = calls["fetch"]
    cached = iterate._fetch_diverse(params, logger)
    assert calls["fetch"] == before  # network skipped
    assert len(cached) == 8  # the full cached set is returned (dereplicated downstream)


def test_interrupted_fetch_leaves_no_reusable_cache(monkeypatch, tmp_path, logger):
    """A fetch that dies part way must not leave a panel the next run trusts.

    Panel composition decides which donors are findable, so silently reusing a
    truncated set narrows detection without saying so.
    """
    from tessera.discover import iterate
    from tessera.discover import pool as pool_mod
    from tessera.discover.iterate import FillParams

    calls = {"fetch": 0}

    def dying_fetch(taxon, d, *, refseq=True, complete_only=False, released_after=None,
                    limit=None, logger):
        calls["fetch"] += 1
        for i in range(3):  # some genomes land on disk ...
            (d / f"G{i}.fasta").write_text(">x\nA\n")
        raise KeyboardInterrupt("user pressed Ctrl-C mid-fetch")

    monkeypatch.setattr(pool_mod, "fetch_ncbi_virus", dying_fetch)
    params = FillParams(
        query=tmp_path / "q.fasta", collection=None, output=tmp_path / "o",
        taxon="SARS-CoV-2", source_refseq=False, cache_dir=tmp_path / "cache",
    )
    with pytest.raises(KeyboardInterrupt):
        iterate._fetch_diverse(params, logger)

    # ... but nothing is installed, and a retry goes back to the network.
    cache_root = tmp_path / "cache" / "ncbi_virus"
    assert not cache_root.exists() or not any(p.is_dir() for p in cache_root.iterdir())

    def good_fetch(taxon, d, *, refseq=True, complete_only=False, released_after=None,
                   limit=None, logger):
        calls["fetch"] += 1
        (d / "G0.fasta").write_text(">x\nA\n")
        return [d / "G0.fasta"]

    monkeypatch.setattr(pool_mod, "fetch_ncbi_virus", good_fetch)
    (tmp_path / "o").mkdir(exist_ok=True)
    assert len(iterate._fetch_diverse(params, logger)) == 1
    assert calls["fetch"] == 2  # the interrupted attempt was not reused


def test_fetch_scope_does_not_share_a_cache_slot(monkeypatch, tmp_path, logger):
    """--source-refseq and the default ask for different genome sets."""
    from tessera.discover import iterate
    from tessera.discover import pool as pool_mod
    from tessera.discover.iterate import FillParams

    def fake_fetch(taxon, d, *, refseq=True, complete_only=False, released_after=None,
                   limit=None, logger):
        # Enough RefSeq genomes to clear SEED_MIN_DIVERSE, so the two scopes stay distinct
        # instead of both broadening to the complete set.
        name = "refseq" if refseq else "complete"
        written = []
        for i in range(4):
            p = d / f"{name}{i}.fasta"
            p.write_text(">x\nA\n")
            written.append(p)
        return written

    monkeypatch.setattr(pool_mod, "fetch_ncbi_virus", fake_fetch)
    out = tmp_path / "o"
    out.mkdir()
    common = {"query": tmp_path / "q.fasta", "collection": None, "output": out,
              "taxon": "SARS-CoV-2", "cache_dir": tmp_path / "cache"}

    complete = iterate._fetch_diverse(FillParams(source_refseq=False, **common), logger)
    refseq = iterate._fetch_diverse(FillParams(source_refseq=True, **common), logger)

    assert all(p.name.startswith("complete") for p in complete)
    # The RefSeq request gets its own slot, not the complete set cached a moment ago.
    assert all(p.name.startswith("refseq") for p in refseq)


def test_dominant_lineage_token_extracted_from_titles():
    from tessera.discover.iterate import _dominant_lineage_token

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
    from tessera.discover import iterate as it

    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "ACGT" * 200 + "\n")
    captured = {}

    def fake_blast(seq, *, max_hits, logger, email=None, entrez_query=None, cache_dir=None):
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


def test_seed_source_nextclade_routes_through_pool_selection(monkeypatch, tmp_path, logger):
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "ACGT" * 100 + "\n")
    out = tmp_path / "out"
    _common_mocks(monkeypatch, [([], 0.95)])  # converge immediately, no NCBI rounds
    monkeypatch.setattr(iterate, "read_fasta", lambda p: [("q", "ACGT" * 100)])

    pool_dir = tmp_path / "pool"
    pool_dir.mkdir()
    (pool_dir / "REF1.fasta").write_text(">REF1 A1\nACGT\n")

    captured = {}

    def fake_fetch(params, logger):
        captured["dataset"] = params.nextclade_dataset
        return [pool_dir / "REF1.fasta"]

    def fake_select(params, genomes, logger):
        from tessera.discover.pool import PoolSelection
        return PoolSelection(selected=list(genomes))

    monkeypatch.setattr(iterate, "_fetch_nextclade", fake_fetch)
    monkeypatch.setattr(iterate, "_select_from", fake_select)

    fill_references(
        FillParams(query=query, collection=None, output=out, seed_source="nextclade",
                   nextclade_dataset="nextstrain/sars-cov-2/XBB"),
        logger,
    )
    assert captured["dataset"] == "nextstrain/sars-cov-2/XBB"
    assert (out / "collection" / "REF1.fasta").exists()


# --- the round structure: what is aligned, and when the panel is curated ----------

def _recording_build(monkeypatch, builds):
    """Replace build_msa with a stub that records which genomes each MSA was built from
    and writes them into the file, so the published panel can be inspected."""
    def fake_build(p, logger):
        members = sorted(f.name for f in p.collection.iterdir())
        builds.append((p.output.name, members))
        p.output.write_text("".join(f">{m}\nACGT\n" for m in ["q", *members]))
    monkeypatch.setattr(iterate, "build_msa", fake_build)


def _one_new_hit_per_round(monkeypatch):
    counter = iter(range(1, 50))

    def collect(*a, **k):
        return [Candidate(_gap(0.8), Hit(f"NEW{next(counter)}", "t", 90.0, 95.0, 1e-9), False)]

    def download(cands, dest, logger):
        for c in cands:
            (dest / f"{c.hit.accession}.fasta").write_text(f">{c.hit.accession}\nACGT\n")
        return cands

    monkeypatch.setattr(iterate, "collect_candidates", collect)
    monkeypatch.setattr(iterate, "_download", download)


def _fake_skani(monkeypatch, table):
    """Stub skani at the panel layer: ``table`` maps a genome label to (ANI %, AF %)."""
    from tessera.discover import panel

    def fake_ani(query, refs, logger):
        return {r: table[r.name.split(".")[0]] for r in refs}

    monkeypatch.setattr(panel, "skani_query_ani", fake_ani)
    monkeypatch.setattr(panel, "skani_available", lambda: True)
    monkeypatch.setattr(panel, "skder_available", lambda: False)
    monkeypatch.setattr(iterate, "skani_available", lambda: True)


def test_last_round_downloads_are_aligned(monkeypatch, tmp_path, logger):
    """On a max_rounds exit the final round's downloads used to sit in collection/ and be
    counted in the summary without ever reaching the published alignment."""
    query, coll, out = _setup(tmp_path)
    _common_mocks(monkeypatch, [([_gap(0.80)], 0.94), ([_gap(0.90)], 0.94)])
    builds: list[tuple[str, list[str]]] = []
    _recording_build(monkeypatch, builds)
    _one_new_hit_per_round(monkeypatch)

    trace = fill_references(
        FillParams(query=query, collection=coll, output=out, max_rounds=2), logger
    )

    assert [(r.round, r.added) for r in trace] == [(1, ["NEW1"]), (2, ["NEW2"])]
    final_members = ["NEW1.fasta", "NEW2.fasta", "refA.fasta"]
    assert sorted(p.name for p in (out / "collection").iterdir()) == final_members
    assert builds == [
        ("round1.msa.fasta", ["refA.fasta"]),
        ("round2.msa.fasta", ["NEW1.fasta", "refA.fasta"]),
        ("final.msa.fasta", final_members),
    ]
    published = (out / "panel.msa.fasta").read_text()
    assert ">NEW2.fasta" in published


def test_no_extra_build_when_the_last_round_adds_nothing(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    _common_mocks(monkeypatch, [([_gap(0.84)], 0.94), ([], 0.95)])
    builds: list[tuple[str, list[str]]] = []
    _recording_build(monkeypatch, builds)
    _one_new_hit_per_round(monkeypatch)

    fill_references(FillParams(query=query, collection=coll, output=out, max_rounds=2), logger)

    assert [name for name, _ in builds] == ["round1.msa.fasta", "round2.msa.fasta"]
    assert not (out / "final.msa.fasta").exists()


def test_curate_runs_on_a_supplied_collection_before_the_first_build(
    monkeypatch, tmp_path, logger
):
    """A whole-genome sibling in the starting collection leaves no coverage gap, so the
    loop converges in round 1. Curation used to sit after the download step and never ran
    -- in exactly the case it exists for."""
    query, coll, out = _setup(tmp_path)
    (coll / "SIBLING.fasta").write_text(">SIBLING\nACGT\n")
    (coll / "PARENT.fasta").write_text(">PARENT\nACGT\n")
    _common_mocks(monkeypatch, [([], 0.95)])
    builds: list[tuple[str, list[str]]] = []
    _recording_build(monkeypatch, builds)
    # refA is picked as backbone only if SIBLING does not outrank it; make the sibling a
    # whole-genome twin of refA (comparable ANI, full coverage) and PARENT a regional donor.
    _fake_skani(monkeypatch, {
        "refA": (97.0, 95.0), "SIBLING": (96.5, 95.0), "PARENT": (88.0, 30.0),
    })
    sections = {}
    monkeypatch.setattr(
        iterate, "run_recomb", lambda p, logger, **kw: sections.update(kw),
    )

    fill_references(FillParams(query=query, collection=coll, output=out, curate=True), logger)

    assert builds == [("round1.msa.fasta", ["PARENT.fasta", "refA.fasta"])]
    assert (coll / "SIBLING.fasta").exists()  # the user's own collection is untouched
    panel_tsv = (out / "panel_lineages.tsv").read_text()
    assert "SIBLING\tsibling-dropped" in panel_tsv
    assert "Reference panel" in [title for title, _ in sections["extra_sections"]]


def test_curate_does_not_run_on_a_freshly_seeded_collection(monkeypatch, tmp_path, logger):
    """Seeding already filters siblings (parents mode, pool selection). Curating the seed
    again before round 1 would change what `detect` recruits, so it is not done."""
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "ACGT" * 100 + "\n")
    out = tmp_path / "out"
    _common_mocks(monkeypatch, [([], 0.95)])
    monkeypatch.setattr(iterate, "skani_available", lambda: True)

    def seed(params, collection, query_records, exclude, logger):
        for name in ("S1", "S2"):
            (collection / f"{name}.fasta").write_text(f">{name}\nACGT\n")

    monkeypatch.setattr(iterate, "_seed_collection", seed)
    calls = []
    monkeypatch.setattr(iterate, "curate_collection_dir",
                        lambda *a, **k: calls.append("curate"))
    monkeypatch.setattr(iterate, "pick_backbone", lambda *a, **k: calls.append("backbone"))

    fill_references(FillParams(query=query, collection=None, output=out, curate=True), logger)

    assert calls == []


def test_curate_never_removes_the_user_reference(monkeypatch, tmp_path, logger):
    """`--curate --reference X` used to auto-pick another backbone, drop X as its twin,
    and fail the next round with "Reference 'X' not found among the staged genomes".
    X is the alignment's coordinate reference, not the sibling test's anchor: the anchor
    stays the query's closest relative and X is kept whatever the comparison says."""
    from tessera.core.io import select_reference

    query, coll, out = _setup(tmp_path)
    (coll / "MYREF.fasta").write_text(">MYREF\nACGT\n")
    _common_mocks(monkeypatch, [([_gap(0.80)], 0.94), ([_gap(0.90)], 0.94)])
    _one_new_hit_per_round(monkeypatch)
    # Auto-picking would choose refA (highest ANI); MYREF is then its whole-genome twin.
    _fake_skani(monkeypatch, {
        "refA": (97.0, 95.0), "MYREF": (96.5, 95.0),
        "NEW1": (88.0, 30.0), "NEW2": (87.0, 30.0),
    })
    builds: list[str] = []

    def fake_build(p, logger):  # resolve the reference the way the real build_msa does
        select_reference(sorted(p.collection.iterdir()), p.query, False, p.reference)
        builds.append(p.output.name)
        p.output.write_text(">q\nACGT\n")

    monkeypatch.setattr(iterate, "build_msa", fake_build)

    fill_references(
        FillParams(query=query, collection=coll, output=out, curate=True,
                   reference="MYREF", max_rounds=2),
        logger,
    )

    assert builds == ["round1.msa.fasta", "round2.msa.fasta", "final.msa.fasta"]
    assert (out / "collection" / "MYREF.fasta").exists()
    rows = dict(
        line.split("\t")[:2]
        for line in (out / "panel_lineages.tsv").read_text().splitlines()[1:]
    )
    assert rows["refA"] == "backbone"
    assert rows["MYREF"] == "sibling-kept"


def test_last_round_downloads_are_curated_before_the_final_build(monkeypatch, tmp_path, logger):
    """The final build sees a curated collection. Here the only download is a sibling, so
    after curation the collection is what round 1 was built from and no rebuild is needed."""
    query, coll, out = _setup(tmp_path)
    _common_mocks(monkeypatch, [([_gap(0.80)], 0.94)])
    builds: list[tuple[str, list[str]]] = []
    _recording_build(monkeypatch, builds)
    _one_new_hit_per_round(monkeypatch)
    # The one download (NEW1) is a sibling: closer than the backbone, whole-genome.
    _fake_skani(monkeypatch, {"refA": (92.0, 95.0), "NEW1": (99.0, 98.0)})

    trace = fill_references(
        FillParams(query=query, collection=coll, output=out, curate=True, max_rounds=1),
        logger,
    )

    assert builds == [("round1.msa.fasta", ["refA.fasta"])]
    assert not (out / "final.msa.fasta").exists()
    assert sorted(p.name for p in (out / "collection").iterdir()) == ["refA.fasta"]
    assert trace[0].added == []  # the sibling was downloaded, then curated away


def test_curate_accepts_a_reference_given_with_its_extension(monkeypatch, tmp_path, logger):
    query, coll, out = _setup(tmp_path)
    (coll / "MYREF.fasta").write_text(">MYREF\nACGT\n")
    _common_mocks(monkeypatch, [([], 0.95)])
    _fake_skani(monkeypatch, {"refA": (97.0, 95.0), "MYREF": (96.5, 95.0)})

    fill_references(
        FillParams(query=query, collection=coll, output=out, curate=True,
                   reference="MYREF.fasta"),
        logger,
    )

    rows = dict(
        line.split("\t")[:2]
        for line in (out / "panel_lineages.tsv").read_text().splitlines()[1:]
    )
    assert rows == {"MYREF": "sibling-kept", "refA": "backbone"}


def test_curate_warns_when_it_leaves_too_few_references(monkeypatch, tmp_path, caplog):
    """On a panel as close to the query as its own lineage, every genome but the backbone
    is classed as a sibling. Detection cannot run on one reference; say why."""
    import logging

    query, coll, out = _setup(tmp_path)
    for name in ("refB", "refC"):
        (coll / f"{name}.fasta").write_text(f">{name}\nACGT\n")
    _common_mocks(monkeypatch, [([], 0.95)])
    _fake_skani(monkeypatch, {
        "refA": (99.6, 99.0), "refB": (99.5, 99.0), "refC": (99.4, 99.0),
    })
    # Not under the "tessera" logger: the CLI tests switch its propagation off, and
    # caplog only sees records that reach the root logger.
    log = logging.getLogger("fill_loop_test")
    with caplog.at_level(logging.WARNING, logger="fill_loop_test"):
        fill_references(FillParams(query=query, collection=coll, output=out, curate=True), log)

    assert sorted(p.name for p in (out / "collection").iterdir()) == ["refA.fasta"]
    assert "Curation left 1 reference" in caplog.text
    assert "without --curate" in caplog.text


def test_curate_with_a_distant_reference_keeps_the_closer_genomes(monkeypatch, tmp_path, logger):
    """The sibling test is relative to its anchor. Anchored on a distant coordinate
    reference, every genome closer to the query than that reference looks like a sibling
    and the panel loses its parental lineages. The anchor is the closest relative."""
    query, coll, out = _setup(tmp_path)
    for name in ("MYREF", "P1", "P2", "P3"):
        (coll / f"{name}.fasta").write_text(f">{name}\nACGT\n")
    _common_mocks(monkeypatch, [([], 0.95)])
    builds: list[tuple[str, list[str]]] = []
    _recording_build(monkeypatch, builds)
    _fake_skani(monkeypatch, {
        "MYREF": (85.0, 95.0), "refA": (95.0, 95.0),
        "P1": (91.0, 95.0), "P2": (89.0, 95.0), "P3": (84.0, 95.0),
    })

    fill_references(
        FillParams(query=query, collection=coll, output=out, curate=True, reference="MYREF"),
        logger,
    )

    assert builds == [(
        "round1.msa.fasta",
        ["MYREF.fasta", "P1.fasta", "P2.fasta", "P3.fasta", "refA.fasta"],
    )]


def test_fresh_seed_refuses_to_clear_a_collection_directory_it_did_not_create(
    monkeypatch, tmp_path, logger
):
    """`detect -o project/` clears `project/collection/` to seed it. If that directory
    holds the user's own genomes and no earlier run made it, refuse, do not delete."""
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "ACGT" * 100 + "\n")
    out = tmp_path / "project"
    (out / "collection").mkdir(parents=True)
    precious = out / "collection" / "my_genome.fasta"
    precious.write_text(">mine\nACGT\n")
    _common_mocks(monkeypatch, [([], 0.95)])
    monkeypatch.setattr(iterate, "_seed_collection", lambda *a, **k: None)

    with pytest.raises(UserInputError, match="was not created by Tessera"):
        fill_references(FillParams(query=query, collection=None, output=out), logger)

    assert precious.exists()


def test_fresh_seed_clears_its_own_working_copy_on_a_rerun(monkeypatch, tmp_path, logger):
    query = tmp_path / "q.fasta"
    query.write_text(">q\n" + "ACGT" * 100 + "\n")
    out = tmp_path / "out"
    _common_mocks(monkeypatch, [([], 0.95), ([], 0.95)])
    names = iter(["S1", "S2"])

    def seed(params, collection, query_records, exclude, logger):
        name = next(names)
        (collection / f"{name}.fasta").write_text(f">{name}\nACGT\n")
        (collection / "other.fasta").write_text(">other\nACGT\n")

    monkeypatch.setattr(iterate, "_seed_collection", seed)
    for _ in range(2):
        fill_references(FillParams(query=query, collection=None, output=out), logger)

    assert sorted(p.name for p in (out / "collection").iterdir()) == ["S2.fasta", "other.fasta"]
