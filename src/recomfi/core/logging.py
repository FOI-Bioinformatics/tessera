"""Logging configuration.

A single named logger (``recomfi``) writes to the console; when an output
directory is given, a DEBUG-level run log is also written there so failures keep
a full traceback for diagnosis.
"""

from __future__ import annotations

import logging
from pathlib import Path

LOGGER_NAME = "recomfi"


def configure_logging(output_dir: Path | None = None, level: int = logging.INFO) -> logging.Logger:
    """Configure and return the ``recomfi`` logger.

    Idempotent: re-configuring replaces existing handlers so repeated calls (e.g.
    in tests) do not duplicate log lines.
    """
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    logger.propagate = False

    console = logging.StreamHandler()
    console.setLevel(level)
    console.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    logger.addHandler(console)

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        run_log = logging.FileHandler(output_dir / "recomfi.log")
        run_log.setLevel(logging.DEBUG)
        run_log.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
        )
        logger.addHandler(run_log)

    return logger
