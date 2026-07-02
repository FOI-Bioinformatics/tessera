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


def test_build_pool_failure_is_non_fatal(tmp_path, monkeypatch):
    # build_pool raises for HA; the run must continue and mark HA unassigned, not abort.
    na_tip = tmp_path / "NA_pool" / "strainB.fasta"
    na_tip.parent.mkdir(parents=True, exist_ok=True)
    na_tip.write_text(">x\nACGT\n")

    def resolve(fasta, override, *, email, logger):
        return _DS("HA_ds") if "HA" in fasta.read_text() else _DS("NA_ds")

    def build_pool(ds, *, cache_dir, logger):
        if ds.path == "HA_ds":
            raise RuntimeError("pool build blew up")
        return [na_tip]

    monkeypatch.setattr(assign, "skani_available", lambda: True)
    monkeypatch.setattr(assign, "resolve_dataset", resolve)
    monkeypatch.setattr(assign, "nextclade_cache", lambda p, t, override=None: Path("/x"))
    monkeypatch.setattr(assign, "build_pool", build_pool)
    monkeypatch.setattr(assign, "skani_query_ani",
                        lambda q, refs, logger: {na_tip: (99.0, 0.99)})
    monkeypatch.setattr(assign, "_clade_of_tip", lambda tip: "cladeX")

    q = _write_query(tmp_path, [("HA", "HAxx"), ("NA", "NAyy")])
    result = assign_segments(q, logger=LOG)
    status = {s.segment: s.status for s in result.segments}
    assert status["HA"] == "unassigned"
    assert status["NA"] == "assigned"


def test_missing_skani_raises(tmp_path, monkeypatch):
    from tessera.core.errors import UserInputError
    monkeypatch.setattr(assign, "skani_available", lambda: False)
    q = _write_query(tmp_path, [("HA", "AAAA")])
    with pytest.raises(UserInputError):
        assign_segments(q, logger=LOG)
