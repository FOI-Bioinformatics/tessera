"""Stubbed assign_segments orchestration (no skani / Nextclade / network)."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from tessera.reassort import assign
from tessera.reassort.assign import assign_segments

LOG = logging.getLogger("test")


class _DS:
    def __init__(self, path):
        self.path = path
        self.tag = "tag"


def _patch(monkeypatch, *, resolve, tips_by_path, ani_by_path):
    monkeypatch.setattr(assign, "skani_available", lambda: True)
    monkeypatch.setattr(assign, "resolve_dataset", resolve)
    monkeypatch.setattr(assign, "nextclade_cache", lambda path, tag, override=None: Path("/x"))
    monkeypatch.setattr(assign, "build_pool",
                        lambda ds, *, cache_dir, logger: tips_by_path[ds.path])
    monkeypatch.setattr(assign, "skani_query_ani",
                        lambda q, refs, logger: ani_by_path[refs[0].parent.name])
    monkeypatch.setattr(assign, "_clade_of_tip", lambda tip: "cladeX")


def _write_query(tmp_path, records):
    q = tmp_path / "q.fasta"
    q.write_text("".join(f">{n}\n{s}\n" for n, s in records))
    return q


def test_low_af_tip_is_dropped_from_ranking(tmp_path, monkeypatch):
    # One HA tip at high ANI but AF below MIN_AF, one NA tip clean. HA has no eligible tip
    # -> unassigned; NA assigned. Verdict undetermined (only one assigned).
    ha_tip = tmp_path / "HA_pool" / "strainA.fasta"
    na_tip = tmp_path / "NA_pool" / "strainB.fasta"
    for t in (ha_tip, na_tip):
        t.parent.mkdir(parents=True, exist_ok=True)
        t.write_text(">x\nACGT\n")

    def resolve(fasta, override, *, email, logger):
        return _DS("HA_ds") if "HA" in fasta.read_text() else _DS("NA_ds")

    _patch(monkeypatch,
           resolve=resolve,
           tips_by_path={"HA_ds": [ha_tip], "NA_ds": [na_tip]},
           ani_by_path={"HA_pool": {ha_tip: (99.0, 0.10)},   # AF 0.10 < MIN_AF -> dropped
                        "NA_pool": {na_tip: (99.0, 0.99)}})
    q = _write_query(tmp_path, [("HA", "HAxx"), ("NA", "NAyy")])
    result = assign_segments(q, logger=LOG)
    status = {s.segment: s.status for s in result.segments}
    assert status == {"HA": "unassigned", "NA": "assigned"}
    assert result.verdict == "undetermined"


def test_build_pool_error_propagates(tmp_path, monkeypatch):
    # A build_pool failure is a download/environment error, not a per-segment biological skip:
    # it must surface (propagate), not be silently reported as "unassigned" (audit #10).
    na_tip = tmp_path / "NA_pool" / "strainB.fasta"
    na_tip.parent.mkdir(parents=True, exist_ok=True)
    na_tip.write_text(">x\nACGT\n")

    def build_pool(ds, *, cache_dir, logger):
        raise RuntimeError("pool build blew up")

    monkeypatch.setattr(assign, "skani_available", lambda: True)
    monkeypatch.setattr(assign, "resolve_dataset",
                        lambda fasta, override, *, email, logger: _DS("HA_ds"))
    monkeypatch.setattr(assign, "nextclade_cache", lambda p, t, override=None: Path("/x"))
    monkeypatch.setattr(assign, "build_pool", build_pool)

    q = _write_query(tmp_path, [("HA", "HAxx")])
    with pytest.raises(RuntimeError, match="pool build blew up"):
        assign_segments(q, logger=LOG)


def test_build_pool_download_error_propagates(tmp_path, monkeypatch):
    # The real download-failure mode is UserInputError from _download_text; it must surface,
    # not be caught by the ToolExecutionError (short-segment) handler (audit #10).
    from tessera.core.errors import UserInputError

    def build_pool(ds, *, cache_dir, logger):
        raise UserInputError("Could not download https://.../tree.json: timed out")

    monkeypatch.setattr(assign, "skani_available", lambda: True)
    monkeypatch.setattr(assign, "resolve_dataset",
                        lambda fasta, override, *, email, logger: _DS("HA_ds"))
    monkeypatch.setattr(assign, "nextclade_cache", lambda p, t, override=None: Path("/x"))
    monkeypatch.setattr(assign, "build_pool", build_pool)

    q = _write_query(tmp_path, [("HA", "HAxx")])
    with pytest.raises(UserInputError, match="Could not download"):
        assign_segments(q, logger=LOG)


def test_skani_short_segment_is_non_fatal(tmp_path, monkeypatch):
    # skani rejects a very short / odd segment with ToolExecutionError -> a genuine per-segment
    # skip: HA unassigned, NA still assigned, run continues.
    from tessera.core.errors import ToolExecutionError
    ha_tip = tmp_path / "HA_pool" / "strainA.fasta"
    na_tip = tmp_path / "NA_pool" / "strainB.fasta"
    for t in (ha_tip, na_tip):
        t.parent.mkdir(parents=True, exist_ok=True)
        t.write_text(">x\nACGT\n")

    def resolve(fasta, override, *, email, logger):
        return _DS("HA_ds") if "HA" in fasta.read_text() else _DS("NA_ds")

    def skani(q, refs, logger):
        if refs[0].parent.name == "HA_pool":
            raise ToolExecutionError(["skani", "dist"], 1, "sequence too short")
        return {na_tip: (99.0, 0.99)}

    def build_pool(ds, *, cache_dir, logger):
        return [ha_tip] if ds.path == "HA_ds" else [na_tip]

    monkeypatch.setattr(assign, "skani_available", lambda: True)
    monkeypatch.setattr(assign, "resolve_dataset", resolve)
    monkeypatch.setattr(assign, "nextclade_cache", lambda p, t, override=None: Path("/x"))
    monkeypatch.setattr(assign, "build_pool", build_pool)
    monkeypatch.setattr(assign, "skani_query_ani", skani)
    monkeypatch.setattr(assign, "_clade_of_tip", lambda tip: "cladeX")

    q = _write_query(tmp_path, [("HA", "HAxx"), ("NA", "NAyy")])
    result = assign_segments(q, logger=LOG)
    status = {s.segment: s.status for s in result.segments}
    assert status["HA"] == "unassigned"
    assert status["NA"] == "assigned"


def test_resolve_dataset_failure_is_non_fatal(tmp_path, monkeypatch):
    # resolve_dataset raises UserInputError for HA (no dataset maps); HA -> unassigned, NA assigned.
    from tessera.core.errors import UserInputError
    na_tip = tmp_path / "NA_pool" / "strainB.fasta"
    na_tip.parent.mkdir(parents=True, exist_ok=True)
    na_tip.write_text(">x\nACGT\n")

    def resolve(fasta, override, *, email, logger):
        if "HA" in fasta.read_text():
            raise UserInputError("no dataset maps to this segment")
        return _DS("NA_ds")

    monkeypatch.setattr(assign, "skani_available", lambda: True)
    monkeypatch.setattr(assign, "resolve_dataset", resolve)
    monkeypatch.setattr(assign, "nextclade_cache", lambda p, t, override=None: Path("/x"))
    monkeypatch.setattr(assign, "build_pool", lambda ds, *, cache_dir, logger: [na_tip])
    monkeypatch.setattr(assign, "skani_query_ani",
                        lambda q, refs, logger: {na_tip: (99.0, 0.99)})
    monkeypatch.setattr(assign, "_clade_of_tip", lambda tip: "cladeX")

    q = _write_query(tmp_path, [("HA", "HAxx"), ("NA", "NAyy")])
    result = assign_segments(q, logger=LOG)
    status = {s.segment: s.status for s in result.segments}
    assert status["HA"] == "unassigned"
    assert status["NA"] == "assigned"


def test_clade_of_tip_ignores_non_clade_markers(tmp_path):
    # A tree tip's 2nd header token is its clade; an example genome carries "example" and an
    # untyped tip "NA" -> neither is a real clade, so _clade_of_tip returns None (audit #9).
    def _tip(text):
        p = tmp_path / f"{abs(hash(text))}.fasta"
        p.write_text(text)
        return p

    assert assign._clade_of_tip(_tip(">A_Darwin_6_2021 J.2.2\nACGT\n")) == "J.2.2"
    assert assign._clade_of_tip(_tip(">EPI_1 example\nACGT\n")) is None
    assert assign._clade_of_tip(_tip(">tipx NA\nACGT\n")) is None
    assert assign._clade_of_tip(_tip(">tipx\nACGT\n")) is None


def test_duplicate_segment_names_raise(tmp_path, monkeypatch):
    # Two records both named HA would collide on the temp file and candidate key; reject them
    # up front rather than silently drop one (which could manufacture a false reassortant).
    from tessera.core.errors import UserInputError
    monkeypatch.setattr(assign, "skani_available", lambda: True)
    q = _write_query(tmp_path, [("HA", "AAAA"), ("HA", "CCCC")])
    with pytest.raises(UserInputError, match="duplicate segment name"):
        assign_segments(q, logger=LOG)


def test_missing_skani_raises(tmp_path, monkeypatch):
    from tessera.core.errors import UserInputError
    monkeypatch.setattr(assign, "skani_available", lambda: False)
    q = _write_query(tmp_path, [("HA", "AAAA")])
    with pytest.raises(UserInputError):
        assign_segments(q, logger=LOG)


def test_scan_segments_scans_assigned_only(tmp_path, monkeypatch):
    # HA and NA both assign; scan_segments=True -> each assigned segment gets a SegmentScan;
    # scan_segment is stubbed so no aligner is needed. require_aligner is stubbed to a no-op.
    ha_tip = tmp_path / "HA_pool" / "strainA.fasta"
    na_tip = tmp_path / "NA_pool" / "strainB.fasta"
    for t in (ha_tip, na_tip):
        t.parent.mkdir(parents=True, exist_ok=True)
        t.write_text(">x\nACGT\n")

    def resolve(fasta, override, *, email, logger):
        return _DS("HA_ds") if "HA" in fasta.read_text() else _DS("NA_ds")

    _patch(monkeypatch, resolve=resolve,
           tips_by_path={"HA_ds": [ha_tip], "NA_ds": [na_tip]},
           ani_by_path={"HA_pool": {ha_tip: (99.0, 0.99)},
                        "NA_pool": {na_tip: (99.0, 0.99)}})
    monkeypatch.setattr(assign, "require_aligner", lambda aligner: None)
    seen = []

    def fake_scan(segment, seq, dataset, out_dir, *, aligner, cache_dir, logger):
        from tessera.reassort.scan import SegmentScan
        seen.append(segment)
        return SegmentScan(segment, True, segment == "HA", 1 if segment == "HA" else 0,
                           "1 region(s)" if segment == "HA" else "none")
    monkeypatch.setattr(assign, "scan_segment", fake_scan)

    q = _write_query(tmp_path, [("HA", "HAxx"), ("NA", "NAyy")])
    result = assign_segments(q, output=tmp_path / "out", scan_segments=True, logger=LOG)
    assert seen == ["HA", "NA"]
    scans = {s.segment: s for s in result.scans}
    assert scans["HA"].recombinant is True and scans["NA"].recombinant is False


def test_scan_segments_marks_unassigned(tmp_path, monkeypatch):
    # HA unassigned (no eligible tip); scan_segments=True -> HA gets a not-scanned row, no scan.
    ha_tip = tmp_path / "HA_pool" / "strainA.fasta"
    na_tip = tmp_path / "NA_pool" / "strainB.fasta"
    for t in (ha_tip, na_tip):
        t.parent.mkdir(parents=True, exist_ok=True)
        t.write_text(">x\nACGT\n")

    def resolve(fasta, override, *, email, logger):
        return _DS("HA_ds") if "HA" in fasta.read_text() else _DS("NA_ds")

    _patch(monkeypatch, resolve=resolve,
           tips_by_path={"HA_ds": [ha_tip], "NA_ds": [na_tip]},
           ani_by_path={"HA_pool": {ha_tip: (10.0, 0.99)},   # below ani_floor -> unassigned
                        "NA_pool": {na_tip: (99.0, 0.99)}})
    monkeypatch.setattr(assign, "require_aligner", lambda aligner: None)

    def fake_scan(segment, seq, dataset, out_dir, *, aligner, cache_dir, logger):
        from tessera.reassort.scan import SegmentScan
        return SegmentScan(segment, True, False, 0, "none")
    monkeypatch.setattr(assign, "scan_segment", fake_scan)

    q = _write_query(tmp_path, [("HA", "HAxx"), ("NA", "NAyy")])
    result = assign_segments(q, output=tmp_path / "out", scan_segments=True, logger=LOG)
    scans = {s.segment: s for s in result.scans}
    assert scans["HA"].scanned is False and scans["HA"].note == "unassigned"
    assert scans["NA"].scanned is True


def test_scan_segments_needs_output(tmp_path, monkeypatch):
    from tessera.core.errors import UserInputError
    monkeypatch.setattr(assign, "skani_available", lambda: True)
    monkeypatch.setattr(assign, "require_aligner", lambda aligner: None)
    q = _write_query(tmp_path, [("HA", "AAAA")])
    with pytest.raises(UserInputError, match="output"):
        assign_segments(q, scan_segments=True, logger=LOG)
