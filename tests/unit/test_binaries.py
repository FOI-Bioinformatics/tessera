"""External-binary preflight: presence, version parsing, and the failure message.

This is the code that decides whether a user's install works at all, and it had no
tests. The binaries are faked with small executable scripts on a temporary PATH, so
the tests exercise the real `shutil.which` / `subprocess` path rather than mocks.
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from tessera.core.binaries import BinarySpec, _parse_version, check_binaries
from tessera.core.errors import MissingBinaryError


def fake_binary(directory: Path, name: str, output: str, *, to_stderr: bool = False) -> Path:
    """A tiny executable that prints ``output`` and exits 0."""
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    stream = "2" if to_stderr else "1"
    path.write_text(f'#!/bin/sh\nprintf %s "{output}" >&{stream}\n')
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


@pytest.fixture
def on_path(tmp_path, monkeypatch):
    """A directory that is the entire PATH, so only the fakes are visible."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    monkeypatch.setenv("PATH", str(bindir))
    return bindir


# --- version parsing ------------------------------------------------------

@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("v7.526 (2024/Apr/26)", (7, 526, 0)),
        ("mafft 7.526.1", (7, 526, 1)),
        ("skani 0.2.2", (0, 2, 2)),
        # progressiveMauve reports a build date, not a dotted version: no match, and the
        # caller falls back to reporting the raw string.
        ("progressiveMauve build date Feb 13 2015 at 09:47:04", None),
        ("no digits here", None),
        ("", None),
    ],
)
def test_version_parsing(text, expected) -> None:
    assert _parse_version(text) == expected


def test_progressivemauve_build_date_is_reported_verbatim(on_path) -> None:
    """A tool without a dotted version is still usable; provenance keeps what it said."""
    fake_binary(on_path, "progressiveMauve", "progressiveMauve build date Feb 13 2015")
    versions = check_binaries((BinarySpec("progressiveMauve", version_args=("--version",)),))
    assert versions["progressiveMauve"] == "progressiveMauve build date Feb 13 2015"


# --- presence -------------------------------------------------------------

def test_missing_binary_is_reported_by_name(on_path) -> None:
    with pytest.raises(MissingBinaryError, match="nosuchtool: not found on PATH"):
        check_binaries((BinarySpec("nosuchtool"),))


def test_every_problem_is_reported_at_once(on_path) -> None:
    """One message listing all of them, so a user fixes their environment in one pass
    rather than one binary per run."""
    with pytest.raises(MissingBinaryError) as exc:
        check_binaries((BinarySpec("toolA"), BinarySpec("toolB"), BinarySpec("toolC")))
    message = str(exc.value)
    assert "toolA" in message and "toolB" in message and "toolC" in message


def test_present_binary_returns_its_version(on_path) -> None:
    fake_binary(on_path, "mafft", "v7.526 (2024/Apr/26)")
    assert check_binaries((BinarySpec("mafft"),)) == {"mafft": "7.526.0"}


def test_version_read_from_stderr(on_path) -> None:
    # Several aligners print their version to stderr.
    fake_binary(on_path, "noisy", "noisy 2.4.1", to_stderr=True)
    assert check_binaries((BinarySpec("noisy"),))["noisy"] == "2.4.1"


def test_unparseable_version_is_recorded_not_fatal(on_path) -> None:
    # A tool that is present but does not report a parseable version should still run;
    # provenance records what it said.
    fake_binary(on_path, "odd", "some build, no version")
    assert check_binaries((BinarySpec("odd"),))["odd"] == "some build, no version"


def test_silent_binary_is_recorded_as_unknown(on_path) -> None:
    fake_binary(on_path, "quiet", "")
    assert check_binaries((BinarySpec("quiet"),))["quiet"] == "unknown"


# --- minimum versions -----------------------------------------------------

def test_version_below_the_minimum_is_rejected(on_path) -> None:
    fake_binary(on_path, "skani", "skani 0.1.0")
    with pytest.raises(MissingBinaryError, match=r"skani: version 0.1.0 < required 0.2.0"):
        check_binaries((BinarySpec("skani", min_version="0.2.0"),))


def test_version_at_or_above_the_minimum_is_accepted(on_path) -> None:
    for reported in ("skani 0.2.0", "skani 0.2.5", "skani 1.0.0"):
        fake_binary(on_path, "skani", reported)
        versions = check_binaries((BinarySpec("skani", min_version="0.2.0"),))
        assert versions["skani"]


def test_unparseable_version_does_not_fail_a_minimum_check(on_path) -> None:
    # We cannot compare what we cannot parse; refusing to run would be worse than
    # trusting a tool the user has installed deliberately.
    fake_binary(on_path, "skani", "custom build")
    assert check_binaries((BinarySpec("skani", min_version="0.2.0"),))["skani"] == "custom build"


def test_versions_are_returned_for_provenance(on_path) -> None:
    # The point of returning these: run_provenance.json records which tools produced a
    # result, which is what makes it reproducible by someone else.
    fake_binary(on_path, "mafft", "v7.526")
    fake_binary(on_path, "skani", "skani 0.2.2")
    assert check_binaries((BinarySpec("mafft"), BinarySpec("skani"))) == {
        "mafft": "7.526.0",
        "skani": "0.2.2",
    }


def test_no_specs_is_a_no_op(on_path) -> None:
    assert check_binaries(()) == {}


def test_a_binary_that_cannot_execute_is_not_a_crash(on_path) -> None:
    # Present on PATH but not executable: `which` may still find it on some systems.
    broken = on_path / "broken"
    broken.write_text("#!/bin/sh\nexit 0\n")
    broken.chmod(broken.stat().st_mode & ~stat.S_IEXEC & ~stat.S_IXGRP & ~stat.S_IXOTH)
    if os.access(broken, os.X_OK):  # pragma: no cover - platform dependent
        pytest.skip("file is still executable on this platform")
    with pytest.raises(MissingBinaryError):
        check_binaries((BinarySpec("broken"),))
