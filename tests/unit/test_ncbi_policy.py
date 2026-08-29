"""NCBI politeness: credentials, pacing between requests, and bounded retries.

Tessera queries services other people share, and it runs on every machine it is
installed on. These check that it identifies itself where it can, spaces its
requests, and backs off rather than hammering a service that is already unwell.
"""

from __future__ import annotations

import logging

import pytest

from tessera.core.ncbi import (
    BLAST_INTERVAL,
    EUTILS_INTERVAL,
    EUTILS_INTERVAL_WITH_KEY,
    Throttle,
    api_key,
    blast_interval,
    eutils_env,
    resolve_email,
    with_retries,
)

LOG = logging.getLogger("tessera")


# --- credentials ----------------------------------------------------------

def test_api_key_read_from_the_environment(monkeypatch) -> None:
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    assert api_key() is None
    monkeypatch.setenv("NCBI_API_KEY", "  abc123  ")
    assert api_key() == "abc123"
    monkeypatch.setenv("NCBI_API_KEY", "   ")
    assert api_key() is None  # blank is not a key


def test_email_prefers_explicit_then_ncbi_email_then_email(monkeypatch) -> None:
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    monkeypatch.delenv("EMAIL", raising=False)
    assert resolve_email(None) is None
    assert resolve_email("me@example.org") == "me@example.org"

    monkeypatch.setenv("EMAIL", "shell@example.org")
    assert resolve_email(None) == "shell@example.org"
    monkeypatch.setenv("NCBI_EMAIL", "ncbi@example.org")
    assert resolve_email(None) == "ncbi@example.org"
    assert resolve_email("explicit@example.org") == "explicit@example.org"


def test_eutils_env_carries_credentials_to_edirect(monkeypatch) -> None:
    # EDirect takes the API key from the environment, not as a flag, so this is the
    # only route by which efetch can be authenticated.
    monkeypatch.setenv("NCBI_API_KEY", "abc123")
    monkeypatch.setenv("NCBI_EMAIL", "me@example.org")
    assert eutils_env() == {"NCBI_API_KEY": "abc123", "EMAIL": "me@example.org"}

    monkeypatch.delenv("NCBI_API_KEY")
    monkeypatch.delenv("NCBI_EMAIL")
    monkeypatch.delenv("EMAIL", raising=False)
    assert eutils_env() == {}


# --- pacing ---------------------------------------------------------------

def test_blast_interval_defaults_and_overrides(monkeypatch) -> None:
    monkeypatch.delenv("TESSERA_BLAST_INTERVAL", raising=False)
    assert blast_interval() == BLAST_INTERVAL  # NCBI's published BLAST guidance

    monkeypatch.setenv("TESSERA_BLAST_INTERVAL", "0")
    assert blast_interval() == 0.0  # opt out, e.g. against a local BLAST

    monkeypatch.setenv("TESSERA_BLAST_INTERVAL", "nonsense")
    assert blast_interval() == BLAST_INTERVAL  # a typo must not remove the limit


def test_throttle_does_not_delay_the_first_request() -> None:
    t = Throttle(30.0, name="test")
    t.wait(LOG)  # would hang the suite if the first call waited


def test_throttle_spaces_consecutive_requests(monkeypatch) -> None:
    slept: list[float] = []
    monkeypatch.setattr("tessera.core.ncbi.time.sleep", slept.append)

    # first wait() reads the clock once; the second reads it before and after sleeping
    clock = iter([100.0, 100.5, 100.5])
    monkeypatch.setattr("tessera.core.ncbi.time.monotonic", lambda: next(clock))

    t = Throttle(2.0, name="test")
    t.wait(LOG)  # first: no wait
    t.wait(LOG)  # 0.5s later, so it must wait the remaining 1.5s
    assert slept == [pytest.approx(1.5)]


def test_eutils_rate_follows_whether_a_key_is_set(monkeypatch) -> None:
    # 3/s without a key, 10/s with one -- NCBI's published E-utilities limits.
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    assert EUTILS_INTERVAL > EUTILS_INTERVAL_WITH_KEY
    t = Throttle(lambda: EUTILS_INTERVAL_WITH_KEY if api_key() else EUTILS_INTERVAL)
    assert t._interval() == EUTILS_INTERVAL
    monkeypatch.setenv("NCBI_API_KEY", "abc123")
    assert t._interval() == EUTILS_INTERVAL_WITH_KEY


# --- retries --------------------------------------------------------------

def test_retries_recover_a_transient_failure(monkeypatch) -> None:
    monkeypatch.setattr("tessera.core.ncbi.time.sleep", lambda _s: None)
    attempts = {"n": 0}

    def flaky() -> str:
        attempts["n"] += 1
        if attempts["n"] < 3:
            raise OSError("connection reset")
        return "hits"

    assert with_retries(flaky, what="test call", logger=LOG) == "hits"
    assert attempts["n"] == 3


def test_a_persistent_failure_still_surfaces(monkeypatch) -> None:
    # Retrying must not turn a real outage into a silent empty result: a lost gap and
    # a gap with no hits are not the same thing.
    monkeypatch.setattr("tessera.core.ncbi.time.sleep", lambda _s: None)
    calls = {"n": 0}

    def broken() -> str:
        calls["n"] += 1
        raise OSError("nope")

    with pytest.raises(OSError, match="nope"):
        with_retries(broken, what="test call", logger=LOG, attempts=3)
    assert calls["n"] == 3


def test_a_call_that_succeeds_is_not_retried(monkeypatch) -> None:
    calls = {"n": 0}

    def fine() -> str:
        calls["n"] += 1
        return "ok"

    assert with_retries(fine, what="test call", logger=LOG) == "ok"
    assert calls["n"] == 1


# --- BLAST result cache ---------------------------------------------------

def test_blast_cache_days_defaults_and_overrides(monkeypatch) -> None:
    from tessera.core.ncbi import BLAST_CACHE_DAYS, blast_cache_days

    monkeypatch.delenv("TESSERA_BLAST_CACHE_DAYS", raising=False)
    assert blast_cache_days() == BLAST_CACHE_DAYS

    monkeypatch.setenv("TESSERA_BLAST_CACHE_DAYS", "1.5")
    assert blast_cache_days() == 1.5

    monkeypatch.setenv("TESSERA_BLAST_CACHE_DAYS", "0")
    assert blast_cache_days() == 0.0  # disables reuse

    monkeypatch.setenv("TESSERA_BLAST_CACHE_DAYS", "nonsense")
    assert blast_cache_days() == BLAST_CACHE_DAYS  # a typo must not change the policy


# --- the cache actually spares the network --------------------------------

def test_a_repeated_search_does_not_hit_the_network(tmp_path, monkeypatch, caplog):
    """The point of the whole exercise: the second identical search is free.

    Remote BLAST is minutes per search and paced at one submission per 10 seconds,
    so a re-run or a second round asking about the same window used to pay again.
    """
    import logging

    from tessera.discover import blast as blast_mod

    calls = {"n": 0}

    class _Aln:
        accession, hit_def = "NC_1", "a relative"
        class _H:
            identities, align_length, expect = 95, 100, 1e-30
        hsps = [_H()]

    class _Record:
        query_length = 100
        alignments = [_Aln()]

    def fake_qblast(*a, **k):
        calls["n"] += 1
        return object()

    monkeypatch.setattr(blast_mod, "blast_throttle", type("T", (), {"wait": lambda *a: None})())
    monkeypatch.setitem(
        __import__("sys").modules,
        "Bio.Blast",
        type("M", (), {"NCBIWWW": type("W", (), {"qblast": staticmethod(fake_qblast)}),
                       "NCBIXML": type("X", (), {"read": staticmethod(lambda h: _Record())})}),
    )

    log = logging.getLogger("tessera")
    kwargs = dict(max_hits=5, logger=log, cache_dir=tmp_path)

    first = blast_mod.blast_subsequence("ACGTACGTAC", **kwargs)
    assert calls["n"] == 1
    assert [h.accession for h in first] == ["NC_1"]

    second = blast_mod.blast_subsequence("ACGTACGTAC", **kwargs)
    assert calls["n"] == 1, "the second identical search should have been served from cache"
    assert second == first  # and identical, not merely cheap

    # a different sequence is a different question, so it does search
    blast_mod.blast_subsequence("TTTTTTTTTT", **kwargs)
    assert calls["n"] == 2
