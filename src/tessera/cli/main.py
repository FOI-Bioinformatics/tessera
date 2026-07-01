"""Shared CLI app, callback and helpers, plus subcommand registration.

The Typer ``app`` and common helpers live here so the per-command modules
(``cmd_*.py``) register against a single app without circular imports. Importing
those modules has the side effect of registering their commands.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import typer

from .. import __version__
from ..core.errors import TesseraError, UserInputError
from ..core.logging import configure_logging

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Tessera: identify recombination events in a query against a reference collection.",
)

# Top-level run state shared by subcommands (set in the callback).
_RUN_STATE: dict[str, object] = {"log_level": logging.INFO}

DEFAULT_THREADS = 1


def _require_choice(value: str, choices: set[str], label: str) -> None:
    if value not in choices:
        raise UserInputError(
            f"Invalid {label} {value!r}. Choose from: {', '.join(sorted(choices))}."
        )


def _parse_key_values(items: list[str], label: str) -> dict[str, str]:
    """Parse repeated ``key=value`` options into a dict (used for tool extras)."""
    out: dict[str, str] = {}
    for item in items:
        if "=" not in item:
            raise UserInputError(f"{label} must be key=value, got '{item}'.")
        key, value = item.split("=", 1)
        key = key.strip()
        if not key:
            raise UserInputError(f"{label} has an empty key in '{item}'.")
        out[key] = value.strip()
    return out


def get_logger(output_dir: Path | None = None) -> logging.Logger:
    """Return the configured logger; with ``output_dir`` also persist a run log there."""
    return configure_logging(output_dir=output_dir, level=_RUN_STATE["log_level"])  # type: ignore[arg-type]


@contextmanager
def stage_errors(logger: logging.Logger) -> Iterator[None]:
    """Turn errors into clean CLI exits instead of raw tracebacks.

    A :class:`TesseraError` (expected, user-facing) is logged concisely. Any
    other exception is unexpected: a concise message goes to the console and the
    full traceback is logged at DEBUG (shown with ``--verbose``). Both exit
    non-zero.
    """
    try:
        yield
    except typer.Exit:
        raise
    except TesseraError as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=1) from exc
    except Exception as exc:
        logger.error("Unexpected error: %s", exc)
        logger.debug("Full traceback:", exc_info=True)
        logger.error("Re-run with --verbose for the full traceback.")
        raise typer.Exit(code=1) from exc


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"tessera {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        False, "--version", callback=_version_callback, is_eager=True,
        help="Show version and exit.",
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose (DEBUG) logging."),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Only warnings and errors."),
) -> None:
    """Tessera top-level entry point."""
    if quiet:
        _RUN_STATE["log_level"] = logging.WARNING
    elif verbose:
        _RUN_STATE["log_level"] = logging.DEBUG
    else:
        _RUN_STATE["log_level"] = logging.INFO


# Register subcommands (import for side effects).
from . import (  # noqa: E402,F401
    cmd_build_panel,
    cmd_curate_panel,
    cmd_detect,
    cmd_fill_references,
    cmd_find_references,
    cmd_msa,
    cmd_recomb,
    cmd_type_lineages,
)


if __name__ == "__main__":
    app()
