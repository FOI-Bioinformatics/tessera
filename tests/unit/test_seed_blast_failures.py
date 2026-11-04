"""Seeding tells a BLAST search that failed from one that found nothing (issue #118).

During a DNS outage two of five seed-window searches failed; seeding then logged "no
distinct parental lineage found in NCBI" for a conclusion it had not established. A prior
fix already stopped an all-failed run from counting as saturation; these tests cover the
partial-failure case and the two neighbouring places with the same conflation.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from tessera.core.errors import UserInputError
from tessera.discover import iterate
from tessera.discover.blast import BlastError, Hit
from tessera.discover.iterate import FillParams

SIBLING = Hit("SIB1", "sibling", 98.0, 99.0, 0.0)
QUERY = "A" * 100 + "C" * 100 + "G" * 100 + "T" * 100 + "A" * 100  # five 100 bp windows


class _Records(logging.Handler):
    def __init__(self) -> None:
        super().__init__(logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)

    def text(self, level: int | None = None) -> str:
        return "\n".join(
            r.getMessage() for r in self.records if level is None or r.levelno == level
        )


@pytest.fixture
def seed_log():
    log = logging.getLogger("tessera.test.seed118")
    log.setLevel(logging.DEBUG)
    handler = _Records()
    log.addHandler(handler)
    yield log, handler
    log.removeHandler(handler)


def _params(tmp_path: Path) -> FillParams:
    return FillParams(
        query=tmp_path / "q.fasta", collection=None, output=tmp_path / "o", seed_window=100
    )


def _blast(failing_first_bases: str, hits: list[Hit]):
    """A stub for ``blast_subsequence``: a window whose first base is in
    ``failing_first_bases`` raises ``BlastError``; the rest return ``hits``."""

    def fake(seq, *, max_hits, logger, email=None, entrez_query=None, cache_dir=None):
        if seq[0] in failing_first_bases:
            raise BlastError("NCBI BLAST request failed: name resolution")
        return list(hits)

    return fake


# --- the helper --------------------------------------------------------------------

def test_blast_or_none_records_the_failure_and_still_returns_nothing(
    monkeypatch, tmp_path, seed_log
):
    log, _ = seed_log
    monkeypatch.setattr(iterate, "blast_subsequence", _blast("A", []))
    failures: list[BlastError] = []

    assert iterate._blast_or_none("A" * 60, _params(tmp_path), log, failures=failures) == []
    assert len(failures) == 1
    # without a recorder the old behaviour is unchanged
    assert iterate._blast_or_none("A" * 60, _params(tmp_path), log) == []


# --- windowed seeding --------------------------------------------------------------

def test_partial_failure_is_reported_and_not_called_absence(monkeypatch, tmp_path, seed_log):
    """Windows 2 and 3 fail; the other three return only siblings."""
    log, records = seed_log
    monkeypatch.setattr(iterate, "blast_subsequence", _blast("CG", [SIBLING]))

    selected, saturated = iterate._seed_windowed(
        QUERY, _params(tmp_path), set(), log, drop_siblings=True
    )

    assert saturated is True  # the protective switch is kept
    assert selected == ["SIB1"]
    warnings = records.text(logging.WARNING)
    assert "2 of 5 seed windows could not be searched" in warnings
    assert "no distinct parental lineage" not in records.text()
    assert "3 of 5" in records.text()  # the windows the conclusion rests on


def test_without_failures_the_original_message_is_kept(monkeypatch, tmp_path, seed_log):
    log, records = seed_log
    monkeypatch.setattr(iterate, "blast_subsequence", _blast("", [SIBLING]))

    _, saturated = iterate._seed_windowed(
        QUERY, _params(tmp_path), set(), log, drop_siblings=True
    )

    assert saturated is True
    assert "no distinct parental lineage found in NCBI" in records.text()
    assert "could not be searched" not in records.text()


def test_every_window_failing_stops_with_a_clear_error(monkeypatch, tmp_path, seed_log):
    log, _ = seed_log
    monkeypatch.setattr(iterate, "blast_subsequence", _blast("ACGT", []))

    with pytest.raises(BlastError, match="every seed search failed"):
        iterate._seed_windowed(QUERY, _params(tmp_path), set(), log, drop_siblings=True)


def test_windows_that_return_nothing_are_not_failures(monkeypatch, tmp_path, seed_log):
    """A search that ran and found nothing is a different outcome from one that failed."""
    log, records = seed_log
    monkeypatch.setattr(iterate, "blast_subsequence", _blast("", []))

    selected, saturated = iterate._seed_windowed(
        QUERY, _params(tmp_path), set(), log, drop_siblings=True
    )

    assert (selected, saturated) == ([], False)
    assert "could not be searched" not in records.text()


# --- whole-query seeding -----------------------------------------------------------

def test_whole_query_failure_is_not_reported_as_an_empty_result(monkeypatch, tmp_path, seed_log):
    log, _ = seed_log
    monkeypatch.setattr(iterate, "blast_subsequence", _blast("A", []))

    with pytest.raises(BlastError, match="whole-query BLAST failed"):
        iterate._seed_whole(QUERY, _params(tmp_path), set(), log)


def test_whole_query_empty_result_keeps_its_own_message(monkeypatch, tmp_path, seed_log):
    log, _ = seed_log
    monkeypatch.setattr(iterate, "blast_subsequence", _blast("", []))

    with pytest.raises(UserInputError, match="returned nothing"):
        iterate._seed_whole(QUERY, _params(tmp_path), set(), log)


# --- negative-lineage seeding ------------------------------------------------------

_LINEAGE_HITS = [
    Hit("SIB1", "Norovirus GII isolate Hu/GII.P16-GII.1/A", 98.0, 99.0, 0.0),
    Hit("SIB2", "Norovirus GII isolate Hu/GII.P16-GII.1/B", 98.0, 99.0, 0.0),
]


def test_negative_lineage_probe_failure_is_not_a_missing_lineage(monkeypatch, tmp_path, seed_log):
    log, records = seed_log

    def fake(seq, *, max_hits, logger, email=None, entrez_query=None, cache_dir=None):
        raise BlastError("NCBI BLAST request failed: name resolution")

    monkeypatch.setattr(iterate, "blast_subsequence", fake)

    assert iterate._seed_negative_lineage("ACGT" * 200, _params(tmp_path), set(), log) == []
    text = records.text()
    assert "lineage probe failed" in text
    assert "Could not identify the query's lineage" not in text


def test_negative_lineage_missing_token_keeps_its_own_message(monkeypatch, tmp_path, seed_log):
    log, records = seed_log
    monkeypatch.setattr(iterate, "blast_subsequence", _blast("", []))

    assert iterate._seed_negative_lineage("ACGT" * 200, _params(tmp_path), set(), log) == []
    assert "Could not identify the query's lineage" in records.text()
    assert "lineage probe failed" not in records.text()


def test_negative_lineage_window_failures_are_counted(monkeypatch, tmp_path, seed_log):
    log, records = seed_log

    def fake(seq, *, max_hits, logger, email=None, entrez_query=None, cache_dir=None):
        if entrez_query is None:  # the whole-query lineage probe succeeds
            return list(_LINEAGE_HITS)
        raise BlastError("NCBI BLAST request failed: name resolution")

    monkeypatch.setattr(iterate, "blast_subsequence", fake)

    assert iterate._seed_negative_lineage("ACGT" * 200, _params(tmp_path), set(), log) == []
    assert "8 of 8 negative-lineage windows could not be searched" in records.text(
        logging.WARNING
    )
