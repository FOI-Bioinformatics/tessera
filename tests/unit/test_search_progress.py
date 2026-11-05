"""A gap searched in many pieces says how many searches that is (issue #129).

Round 1 of a run on a query far from every reference planned 52 sequential remote searches,
about 17 minutes, and said only "in 52 sub-intervals" before going quiet.
"""

from __future__ import annotations

import logging

import pytest

from tessera.discover import run as discover_run
from tessera.discover.blast import Hit
from tessera.discover.run import collect_candidates
from tessera.recomb.coverage import CoverageGap


def _gap(length: int) -> CoverageGap:
    return CoverageGap(msa_start=0, msa_end=length, query_start=0, query_end=length,
                       length_bp=length, n_windows=length // 100, best_label="ref",
                       mean_best=0.6, kind="divergent")


def _collect(monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture, length: int):
    monkeypatch.setattr(discover_run, "blast_subsequence",
                        lambda seq, **k: [Hit("X1", "x", 90.0, 95.0, 1e-9)])
    logger = logging.getLogger("tessera.test.progress")
    with caplog.at_level(logging.INFO, logger=logger.name):
        collect_candidates([_gap(length)], "ACGT" * (length // 4), set(), max_hits=5,
                           email=None, exclude=set(), keep_self_hits=True, logger=logger,
                           subtile=400)
    return caplog.text


def test_a_many_piece_gap_announces_the_number_of_searches_and_a_time(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    text = _collect(monkeypatch, caplog, 10000)
    assert "remote searches" in text
    assert "minute(s)" in text
    assert "cache cost nothing" in text


def test_each_search_is_numbered(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    text = _collect(monkeypatch, caplog, 10000)
    assert "sub-interval 1 of " in text
    assert "sub-interval 2 of " in text


def test_a_small_gap_is_not_announced(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    text = _collect(monkeypatch, caplog, 800)
    assert "remote searches" not in text
