"""Single choke point for running external tools.

Every adapter runs its binary through :func:`run_tool` rather than calling
``subprocess`` directly. Commands are passed as argument lists and executed
without a shell, so user-supplied paths and filenames cannot inject shell
metacharacters. A non-zero exit raises :class:`ToolExecutionError`.

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
import subprocess
from collections.abc import Sequence
from pathlib import Path

from .errors import ToolExecutionError
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
            with open(partial, "w") as out:
                proc = subprocess.run(
                    cmd, stdout=out, stderr=subprocess.PIPE, text=True,
                    timeout=timeout, env=env,
                )
        except subprocess.TimeoutExpired as exc:
            partial.unlink(missing_ok=True)
            raise _timeout_error(cmd, timeout) from exc
        except BaseException:
            partial.unlink(missing_ok=True)
            raise
        output = proc.stderr or ""
        if proc.returncode == 0:
            os.replace(partial, stdout_path)
        else:
            partial.unlink(missing_ok=True)
    else:
        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=timeout, env=env
            )
        except subprocess.TimeoutExpired as exc:
            raise _timeout_error(cmd, timeout) from exc
        output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        tail = "\n".join(output.strip().splitlines()[-_OUTPUT_TAIL_LINES:])
        raise ToolExecutionError(cmd, proc.returncode, tail or None)
    return output


def _timeout_error(cmd: list[str], timeout: float | None) -> ToolExecutionError:
    return ToolExecutionError(
        cmd, -1,
        f"No output within {timeout:g}s; the tool was stopped. Raise or remove the "
        f"limit with ${_TIMEOUT_ENV} (seconds, or 'none') if it needs longer.",
    )
