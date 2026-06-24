"""Exception hierarchy for RecomFi.

Stages and plugins raise these instead of calling ``sys.exit``. The CLI layer
(:mod:`recomfi.cli.main`) catches :class:`RecomFiError`, logs it cleanly and
sets a non-zero exit code, so library callers and tests can handle failures
without a process exit.
"""

from __future__ import annotations


class RecomFiError(Exception):
    """Base class for all expected, user-facing RecomFi failures."""


class UserInputError(RecomFiError):
    """Invalid or missing user input (bad arguments, missing target, etc.)."""


class MissingBinaryError(RecomFiError):
    """A required external tool was not found on PATH or failed its version check."""


class OutputError(RecomFiError):
    """An external tool ran but did not produce the expected output."""


class ToolExecutionError(RecomFiError):
    """An external tool exited with a non-zero status."""

    def __init__(self, command: list[str], returncode: int, output: str | None = None):
        self.command = command
        self.returncode = returncode
        self.output = output
        rendered = " ".join(command)
        msg = f"command failed (exit {returncode}): {rendered}"
        if output:
            msg += f"\n--- output tail ---\n{output}"
        super().__init__(msg)


class PluginError(RecomFiError):
    """A requested plugin/tool could not be found or loaded."""
