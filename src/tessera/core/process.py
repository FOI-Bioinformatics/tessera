"""Single choke point for running external tools.

Every adapter runs its binary through :func:`run_tool` rather than calling
``subprocess`` directly. Commands are passed as argument lists and executed
without a shell, so user-supplied paths and filenames cannot inject shell
metacharacters. A non-zero exit raises :class:`ToolExecutionError`. Standard input
is ``/dev/null``: no tool here is fed on stdin, and one that waited on an inherited
terminal would hang the run, as alignment has no default timeout.

Two properties matter for unattended runs. A tool can be given a ``timeout``, so a
wedged child fails the run instead of hanging it forever -- relevant because these
calls are dispatched through a thread pool, where one stuck child holds its slot and
blocks interpreter shutdown. And output redirected to a file is written through a
temporary path and moved into place on success, so a tool killed part way through
leaves no truncated file for the next stage to read as if it were complete.
"""

from __future__ import annotations

import logging
import os
import signal
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import IO

from .errors import MissingBinaryError, ToolExecutionError
from .plugins import ToolCapabilities

# How many trailing lines of captured output to attach to a failure message.
_OUTPUT_TAIL_LINES = 40

# Timeouts for tools that are bounded operations: a download or a query that should
# finish in minutes. Alignment is deliberately not on this list -- Cactus or SibeliaZ
# on a large panel legitimately runs for hours, and a default that killed those would
# be worse than no default at all. Override either with $TESSERA_TOOL_TIMEOUT (seconds,
# "0" or "none" to disable).
FETCH_TIMEOUT = 1800.0
QUERY_TIMEOUT = 300.0
_TIMEOUT_ENV = "TESSERA_TOOL_TIMEOUT"


def default_timeout(fallback: float | None) -> float | None:
    """``$TESSERA_TOOL_TIMEOUT`` if set and valid, else ``fallback``."""
    raw = os.environ.get(_TIMEOUT_ENV)
    if raw is None:
        return fallback
    text = raw.strip().lower()
    if text in {"", "0", "none", "off"}:
        return None
    try:
        seconds = float(text)
    except ValueError:
        return fallback
    return seconds if seconds > 0 else None


def run_tool(
    capabilities: ToolCapabilities,
    command: Sequence[str | Path],
    *,
    logger: logging.Logger,
    log_prefix: str = "",
    extra_mounts: Sequence[str] | None = None,
    stdout_path: Path | None = None,
    timeout: float | None = None,
    extra_env: dict[str, str] | None = None,
) -> str:
    """Run ``command`` (an argument list) without a shell, returning its output.

    ``capabilities`` is accepted for interface parity with the aligner adapters
    (and future container backends); ``extra_mounts`` is likewise accepted and
    ignored when running natively. Raises :class:`ToolExecutionError` on a
    non-zero exit.

    When ``stdout_path`` is given the tool's standard output is written there
    (not captured), and only standard error is captured for diagnostics. Use this
    for tools whose primary result is emitted on stdout (e.g. MAFFT alignment,
    minimap2 SAM), so verbose progress on stderr does not corrupt the output. The
    write goes through a ``.partial`` sibling and is renamed on success, so a tool
    that dies mid-write does not leave a truncated file at the final path.

    ``timeout`` bounds the run in seconds (``None`` for no limit); exceeding it kills
    the tool and raises :class:`ToolExecutionError` rather than hanging the run.
    ``extra_env`` adds variables to the child's environment -- some tools (EDirect)
    take credentials that way rather than as flags.
    """
    cmd = [str(part) for part in command]
    prefix = f"[{log_prefix}] " if log_prefix else ""
    logger.debug("%srunning: %s", prefix, " ".join(cmd))
    env = {**os.environ, **extra_env} if extra_env else None

    if stdout_path is not None:
        # Write beside the target, then rename: a partial file at the final path would
        # be read by the next stage as though the tool had finished.
        partial = Path(stdout_path).with_name(Path(stdout_path).name + ".partial")
        try:
            with open(partial, "w") as sink:
                code, _, err = _run(cmd, stdout=sink, timeout=timeout, env=env)
        except subprocess.TimeoutExpired as exc:
            partial.unlink(missing_ok=True)
            raise _timeout_error(cmd, timeout) from exc
        except FileNotFoundError as exc:
            partial.unlink(missing_ok=True)
            raise _missing_binary(cmd, exc) from exc
        except BaseException:
            partial.unlink(missing_ok=True)
            raise
        output = err or ""
        if code == 0:
            os.replace(partial, stdout_path)
        else:
            partial.unlink(missing_ok=True)
    else:
        try:
            code, out, err = _run(cmd, stdout=subprocess.PIPE, timeout=timeout, env=env)
        except subprocess.TimeoutExpired as exc:
            raise _timeout_error(cmd, timeout) from exc
        except FileNotFoundError as exc:
            raise _missing_binary(cmd, exc) from exc
        output = (out or "") + (err or "")
    if code != 0:
        tail = "\n".join(output.strip().splitlines()[-_OUTPUT_TAIL_LINES:])
        raise ToolExecutionError(cmd, code, tail or None)
    return output


def _run(
    cmd: list[str], *, stdout: IO[str] | int, timeout: float | None,
    env: dict[str, str] | None,
) -> tuple[int, str | None, str | None]:
    """Run ``cmd`` in a process group of its own; on a timeout stop the whole group.

    Several of the tools here are wrappers that do their work in children of their own
    (the SibeliaZ shell script, Cactus through Toil, the EDirect perl scripts).
    ``subprocess.run``'s timeout signals only the process it started, so those children
    outlived the run that had given up on them -- reparented to init, still holding CPU
    and scratch space, which on a shared machine is somebody else's problem. Starting a
    new session makes the tool a process-group leader, so one ``killpg`` reaches
    everything it spawned.
    """
    with subprocess.Popen(
        cmd, stdin=subprocess.DEVNULL, stdout=stdout, stderr=subprocess.PIPE,
        text=True, env=env, start_new_session=True,
    ) as proc:
        try:
            out, err = proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            _kill_group(proc)
            proc.communicate()  # reap the group, so no zombie is left behind
            raise
        return proc.returncode, out, err


def _kill_group(proc: subprocess.Popen[str]) -> None:
    """SIGKILL the tool's whole process group, falling back to the tool alone."""
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except (ProcessLookupError, PermissionError, OSError):
        proc.kill()


def _missing_binary(cmd: list[str], exc: FileNotFoundError) -> Exception:
    """A tool that is not installed, reported as such.

    Most call sites check for their binary first, but not all can (a tool reached
    only on one branch, a wrapper that shells out to a second program), and the bare
    ``FileNotFoundError`` surfaced as "Unexpected error: [Errno 2]" -- which reads as
    a missing input file. Only the executable itself is treated this way; any other
    missing path is re-raised unchanged.
    """
    if exc.filename is not None and str(exc.filename) != cmd[0]:
        return exc
    return MissingBinaryError(
        f"'{cmd[0]}' was not found on PATH. Install it (see docs/aligners.md and "
        "environment.yml) or activate the environment that provides it."
    )


def _timeout_error(cmd: list[str], timeout: float | None) -> ToolExecutionError:
    return ToolExecutionError(
        cmd, -1,
        f"Still running after {timeout:g}s; the tool was stopped, with every process "
        f"it had started. Raise or remove the "
        f"limit with ${_TIMEOUT_ENV} (seconds, or 'none') if it needs longer.",
    )
