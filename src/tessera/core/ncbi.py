"""Politeness policy for NCBI network access: credentials, pacing, and retries.

Tessera queries public NCBI services that other people share. NCBI publishes usage
limits for them, and a tool distributed to many machines has to observe those limits
itself -- an installed copy that ignores them is indistinguishable from abuse, and the
practical consequence is a blocked address for whoever is running it.

Three things live here:

* **Credentials.** ``$NCBI_API_KEY`` raises the E-utilities rate limit from 3 to 10
  requests per second, and an e-mail address lets NCBI contact a heavy user before
  blocking them. Note the asymmetry: E-utilities (``efetch``) accept both, but the
  BLAST URL API that :mod:`tessera.discover.blast` uses does not -- Biopython's
  ``qblast`` has no parameter for either -- so for BLAST the only lever is pacing.
* **Pacing.** :class:`Throttle` enforces a minimum interval between requests to one
  service. It is process-wide, not per-call-site, because the limit applies to the
  address rather than to a code path.
* **Retries.** A transient network failure previously lost the whole gap or seed
  window it was serving. :func:`with_retries` re-attempts with exponential backoff,
  which also backs off rather than hammering when NCBI is the thing that is unwell.
"""

from __future__ import annotations

import logging
import os
import random
import threading
import time
from collections.abc import Callable
from typing import TypeVar

T = TypeVar("T")

# NCBI's BLAST URL API guidance is not to contact the server more than once every
# 10 seconds. Applied to submissions; Biopython does its own polling within a search.
BLAST_INTERVAL = 10.0
# E-utilities: 3 requests/second without an API key, 10 with one.
EUTILS_INTERVAL = 1.0 / 3.0
EUTILS_INTERVAL_WITH_KEY = 1.0 / 10.0

_INTERVAL_ENV = "TESSERA_BLAST_INTERVAL"
# How long a cached BLAST result stays usable. Long enough that a re-run, a resumed
# run or a second round costs nothing, short enough that a hit list does not quietly
# stop reflecting a database that grows daily.
BLAST_CACHE_DAYS = 14.0
_CACHE_DAYS_ENV = "TESSERA_BLAST_CACHE_DAYS"
_RETRIES = 3
_BACKOFF_BASE = 2.0


def api_key() -> str | None:
    """``$NCBI_API_KEY`` if set and non-empty."""
    return os.environ.get("NCBI_API_KEY", "").strip() or None


def resolve_email(explicit: str | None = None) -> str | None:
    """The contact address: ``explicit``, else ``$NCBI_EMAIL``, else ``$EMAIL``."""
    for candidate in (explicit, os.environ.get("NCBI_EMAIL"), os.environ.get("EMAIL")):
        if candidate and candidate.strip():
            return candidate.strip()
    return None


def eutils_env() -> dict[str, str]:
    """Environment additions that let EDirect tools identify the caller.

    EDirect reads ``NCBI_API_KEY`` from the environment rather than taking a flag,
    so credentials reach ``efetch`` this way or not at all.
    """
    env: dict[str, str] = {}
    key = api_key()
    if key:
        env["NCBI_API_KEY"] = key
    email = resolve_email()
    if email:
        env["EMAIL"] = email
    return env


def blast_interval() -> float:
    """Minimum seconds between BLAST submissions (``$TESSERA_BLAST_INTERVAL``)."""
    raw = os.environ.get(_INTERVAL_ENV)
    if raw is None:
        return BLAST_INTERVAL
    try:
        seconds = float(raw.strip())
    except ValueError:
        return BLAST_INTERVAL
    return max(0.0, seconds)


def blast_cache_days() -> float:
    """Cached-BLAST lifetime in days (``$TESSERA_BLAST_CACHE_DAYS``; 0 disables reuse)."""
    raw = os.environ.get(_CACHE_DAYS_ENV)
    if raw is None:
        return BLAST_CACHE_DAYS
    try:
        days = float(raw.strip())
    except ValueError:
        return BLAST_CACHE_DAYS
    return max(0.0, days)


class Throttle:
    """Enforce a minimum interval between requests to one service.

    Process-wide and thread-safe: the limit NCBI applies is per address, so pacing
    each call site independently would not observe it. Cheap when idle -- the first
    call never waits.
    """

    def __init__(self, interval: float, name: str = "") -> None:
        self._interval = interval
        self._name = name
        self._lock = threading.Lock()
        self._last: float | None = None

    def wait(self, logger: logging.Logger | None = None) -> None:
        with self._lock:
            interval = self._interval() if callable(self._interval) else self._interval
            now = time.monotonic()
            if self._last is not None:
                delay = interval - (now - self._last)
                if delay > 0:
                    if logger is not None and delay >= 1.0:
                        logger.info(
                            "Pausing %.0fs before the next %s request (NCBI rate limit).",
                            delay, self._name or "NCBI",
                        )
                    time.sleep(delay)
                    now = time.monotonic()
            self._last = now


# One throttle per service, shared by every call site in the process.
blast_throttle = Throttle(blast_interval, name="BLAST")
eutils_throttle = Throttle(
    lambda: EUTILS_INTERVAL_WITH_KEY if api_key() else EUTILS_INTERVAL, name="E-utilities"
)


def with_retries(
    call: Callable[[], T],
    *,
    what: str,
    logger: logging.Logger,
    attempts: int = _RETRIES,
    retry_on: tuple[type[BaseException], ...] = (Exception,),
) -> T:
    """Run ``call``, retrying a transient failure with exponential backoff.

    Without this one dropped connection lost whatever the call was serving -- a
    coverage gap's candidates, or a seed window -- and the run continued as though
    that region simply had no hits, which is not the same thing. The last failure is
    re-raised so a persistent problem still surfaces.
    """
    for attempt in range(1, attempts + 1):
        try:
            return call()
        except retry_on as exc:
            if attempt == attempts:
                raise
            # Jitter so parallel workers do not retry in lockstep.
            delay = _BACKOFF_BASE ** (attempt - 1) + random.random()  # noqa: S311
            logger.warning(
                "%s failed (attempt %d/%d: %s); retrying in %.1fs.",
                what, attempt, attempts, exc, delay,
            )
            time.sleep(delay)
    raise AssertionError("unreachable")  # pragma: no cover
