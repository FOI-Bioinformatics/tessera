"""External binary discovery and version preflight.

Every adapter declares the binaries it needs via :class:`BinarySpec` in its
``ToolCapabilities``. Before a tool runs, :func:`check_binaries` confirms each
binary is on PATH and, when a minimum version is declared, that it is new
enough. Resolved versions are returned so a caller can record them for
provenance.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass, field

from .errors import MissingBinaryError

_VERSION_RE = re.compile(r"(\d+)\.(\d+)(?:\.(\d+))?")


@dataclass(frozen=True)
class BinarySpec:
    """A required external executable.

    ``version_args`` is the argument vector that prints a version (e.g.
    ``("--version",)``). ``min_version`` is an optional dotted-string requirement.
    """

    name: str
    version_args: tuple[str, ...] = field(default=("--version",))
    min_version: str | None = None


def _parse_version(text: str) -> tuple[int, int, int] | None:
    match = _VERSION_RE.search(text)
    if not match:
        return None
    major, minor, patch = match.group(1), match.group(2), match.group(3)
    return (int(major), int(minor), int(patch or 0))


def _query_version(name: str, version_args: tuple[str, ...]) -> str | None:
    try:
        proc = subprocess.run(
            [name, *version_args],
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    blob = (proc.stdout or "") + (proc.stderr or "")
    parsed = _parse_version(blob)
    if parsed:
        return ".".join(map(str, parsed))
    return blob.strip().splitlines()[0] if blob.strip() else None


def check_binaries(specs: tuple[BinarySpec, ...]) -> dict[str, str]:
    """Confirm all ``specs`` are present (and new enough). Return name -> version.

    Raises :class:`MissingBinaryError` listing every missing or too-old binary.
    """
    versions: dict[str, str] = {}
    problems: list[str] = []

    for spec in specs:
        if shutil.which(spec.name) is None:
            problems.append(f"{spec.name}: not found on PATH")
            continue

        reported = _query_version(spec.name, spec.version_args)
        versions[spec.name] = reported or "unknown"

        if spec.min_version is not None:
            have = _parse_version(reported or "")
            want = _parse_version(spec.min_version)
            if have is not None and want is not None and have < want:
                problems.append(
                    f"{spec.name}: version {reported} < required {spec.min_version}"
                )

    if problems:
        raise MissingBinaryError(
            "Required external tools are missing or outdated:\n  " + "\n  ".join(problems)
        )
    return versions
