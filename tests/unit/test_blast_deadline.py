"""A BLAST search that never finishes ends in a BlastError (issue #127).

``NCBIWWW.qblast`` polls NCBI until the job is ready and takes no timeout. A search that
stalled left ``fill-references`` waiting for more than half an hour with nothing logged.
"""

from __future__ import annotations

import logging
import threading
import time

import pytest

from tessera.core import ncbi
from tessera.discover import blast
from tessera.discover.blast import BlastError, _within_deadline

LOGGER = logging.getLogger("tessera.test")


def test_a_call_that_finishes_returns_its_value() -> None:
    assert _within_deadline(lambda: 7, 5.0) == 7


def test_a_call_that_fails_raises_its_own_error() -> None:
    def boom() -> None:
        raise ValueError("nope")

    with pytest.raises(ValueError, match="nope"):
        _within_deadline(boom, 5.0)


def test_a_call_that_does_not_finish_raises_blast_error() -> None:
    release = threading.Event()
    started = time.monotonic()
    try:
        with pytest.raises(BlastError, match="did not finish within 0.2 s"):
            _within_deadline(lambda: release.wait(30), 0.2)
    finally:
        release.set()
    assert time.monotonic() - started < 5


def test_no_limit_waits_for_the_call() -> None:
    assert _within_deadline(lambda: "done", None) == "done"


def test_the_message_names_the_setting() -> None:
    release = threading.Event()
    try:
        with pytest.raises(BlastError, match="TESSERA_BLAST_TIMEOUT"):
            _within_deadline(lambda: release.wait(30), 0.1)
    finally:
        release.set()


def test_the_default_timeout_is_fifteen_minutes(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("TESSERA_BLAST_TIMEOUT", raising=False)
    assert ncbi.blast_timeout() == 900.0


@pytest.mark.parametrize(("raw", "expected"), [("60", 60.0), ("0", None), ("-5", None),
                                               ("soon", 900.0)])
def test_the_timeout_can_be_set_from_the_environment(
    monkeypatch: pytest.MonkeyPatch, raw: str, expected: float | None
) -> None:
    monkeypatch.setenv("TESSERA_BLAST_TIMEOUT", raw)
    assert ncbi.blast_timeout() == expected


def test_a_stuck_search_is_a_blast_error_not_a_hang(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    """End to end through ``blast_subsequence``: the search never returns, the call does."""
    release = threading.Event()
    monkeypatch.setenv("TESSERA_BLAST_TIMEOUT", "0.3")
    monkeypatch.setenv("TESSERA_BLAST_INTERVAL", "0")
    import Bio.Blast.NCBIWWW as www

    monkeypatch.setattr(www, "qblast", lambda *a, **k: release.wait(30))
    try:
        with pytest.raises(BlastError, match="did not finish"):
            blast.blast_subsequence("ACGT" * 50, max_hits=3, logger=LOGGER, cache_dir=tmp_path)
    finally:
        release.set()
