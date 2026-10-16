"""Logging configuration: a reconfigured run log releases the file it held."""

from __future__ import annotations

import logging
from pathlib import Path

from tessera.core.logging import configure_logging


def _run_log_handler(logger: logging.Logger) -> logging.FileHandler:
    (handler,) = [h for h in logger.handlers if isinstance(h, logging.FileHandler)]
    return handler


def test_reconfiguring_closes_the_previous_run_log(tmp_path: Path) -> None:
    # A process that configures logging more than once (the test suite, a notebook)
    # replaced the file handler without closing it, leaving one open file per call.
    first = _run_log_handler(configure_logging(tmp_path / "a"))
    try:
        configure_logging(tmp_path / "b")
        assert first.stream is None  # FileHandler.close() releases the stream
    finally:
        configure_logging(None)
    assert (tmp_path / "a" / "tessera.log").exists()
