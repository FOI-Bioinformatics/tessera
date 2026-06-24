"""Single choke point for running external tools.

Every adapter runs its binary through :func:`run_tool` rather than calling
``subprocess`` directly. Commands are passed as argument lists and executed
without a shell, so user-supplied paths and filenames cannot inject shell
metacharacters. A non-zero exit raises :class:`ToolExecutionError`.
"""

from __future__ import annotations

import logging
import subprocess
from collections.abc import Sequence
from pathlib import Path

from .errors import ToolExecutionError
from .plugins import ToolCapabilities

# How many trailing lines of captured output to attach to a failure message.
_OUTPUT_TAIL_LINES = 40


def run_tool(
    capabilities: ToolCapabilities,
    command: Sequence[str | Path],
    *,
    logger: logging.Logger,
    log_prefix: str = "",
    extra_mounts: Sequence[str] | None = None,
) -> str:
    """Run ``command`` (an argument list) without a shell, returning its output.

    ``capabilities`` is accepted for interface parity with the aligner adapters
    (and future container backends); ``extra_mounts`` is likewise accepted and
    ignored when running natively. Raises :class:`ToolExecutionError` on a
    non-zero exit.
    """
    cmd = [str(part) for part in command]
    prefix = f"[{log_prefix}] " if log_prefix else ""
    logger.debug("%srunning: %s", prefix, " ".join(cmd))

    proc = subprocess.run(cmd, capture_output=True, text=True)
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        tail = "\n".join(output.strip().splitlines()[-_OUTPUT_TAIL_LINES:])
        raise ToolExecutionError(cmd, proc.returncode, tail or None)
    return output
