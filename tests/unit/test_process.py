"""The external-tool choke point: timeouts, atomic stdout, and failure reporting.

`core/process.py` is the most environment-fragile code in the package and had no
tests. These use real short-lived subprocesses rather than mocks, because what is
being checked is the interaction with the process itself.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import pytest

from tessera.core.errors import ToolExecutionError
from tessera.core.plugins import ToolCapabilities
from tessera.core.process import default_timeout, run_tool

CAPS = ToolCapabilities(name="test")
LOG = logging.getLogger("tessera")


def _python(*code: str) -> list[str]:
    """A command running inline Python, so the tests need no external binary."""
    return [sys.executable, "-c", "\n".join(code)]


# --- timeouts -------------------------------------------------------------

def test_timeout_stops_a_wedged_tool_instead_of_hanging() -> None:
    with pytest.raises(ToolExecutionError) as exc:
        run_tool(CAPS, _python("import time", "time.sleep(30)"), logger=LOG, timeout=0.4)
    assert "was stopped" in str(exc.value)
    assert "TESSERA_TOOL_TIMEOUT" in str(exc.value)  # says how to raise the limit


def test_no_timeout_by_default_lets_a_slow_tool_finish() -> None:
    # Alignment legitimately runs for hours; an unrequested default would be worse
    # than none.
    assert run_tool(CAPS, _python("import time", "time.sleep(0.2)", "print('ok')"),
                    logger=LOG).strip() == "ok"


def test_default_timeout_reads_the_environment(monkeypatch) -> None:
    monkeypatch.delenv("TESSERA_TOOL_TIMEOUT", raising=False)
    assert default_timeout(300.0) == 300.0

    monkeypatch.setenv("TESSERA_TOOL_TIMEOUT", "45")
    assert default_timeout(300.0) == 45.0

    for disabling in ("0", "none", "off", ""):
        monkeypatch.setenv("TESSERA_TOOL_TIMEOUT", disabling)
        assert default_timeout(300.0) is None

    monkeypatch.setenv("TESSERA_TOOL_TIMEOUT", "not-a-number")
    assert default_timeout(300.0) == 300.0  # a typo must not silently disable the limit


# --- stdout redirection is atomic ----------------------------------------

def test_stdout_redirect_writes_the_output(tmp_path: Path) -> None:
    out = tmp_path / "result.txt"
    run_tool(CAPS, _python("print('>seq')", "print('ACGT')"), logger=LOG, stdout_path=out)
    assert out.read_text() == ">seq\nACGT\n"
    assert not list(tmp_path.glob("*.partial"))  # staging file cleaned up


def test_no_truncated_file_when_the_tool_dies_mid_write(tmp_path: Path) -> None:
    """The failure that motivated this: a killed aligner left a partial alignment at
    the final path, which the next stage read as though it were complete."""
    out = tmp_path / "alignment.fasta"
    with pytest.raises(ToolExecutionError):
        run_tool(
            CAPS,
            _python("import sys", "print('>partial')", "sys.stdout.flush()",
                    "sys.exit(1)"),
            logger=LOG, stdout_path=out,
        )
    assert not out.exists()  # nothing for a later stage to mistake for a finished file
    assert not list(tmp_path.glob("*.partial"))


def test_no_truncated_file_when_the_tool_times_out(tmp_path: Path) -> None:
    out = tmp_path / "alignment.fasta"
    with pytest.raises(ToolExecutionError):
        run_tool(
            CAPS,
            _python("import sys, time", "print('>partial')", "sys.stdout.flush()",
                    "time.sleep(30)"),
            logger=LOG, stdout_path=out, timeout=0.4,
        )
    assert not out.exists()
    assert not list(tmp_path.glob("*.partial"))


def test_existing_output_survives_a_failed_rerun(tmp_path: Path) -> None:
    out = tmp_path / "alignment.fasta"
    out.write_text(">good\nACGT\n")
    with pytest.raises(ToolExecutionError):
        run_tool(CAPS, _python("import sys", "sys.exit(2)"), logger=LOG, stdout_path=out)
    assert out.read_text() == ">good\nACGT\n"  # the previous result is not clobbered


# --- failure reporting ----------------------------------------------------

def test_nonzero_exit_reports_the_command_and_output_tail() -> None:
    with pytest.raises(ToolExecutionError) as exc:
        run_tool(
            CAPS,
            _python("import sys", "print('boom: missing input', file=sys.stderr)",
                    "sys.exit(3)"),
            logger=LOG,
        )
    assert "boom: missing input" in str(exc.value)


def test_arguments_are_not_shell_interpreted(tmp_path: Path) -> None:
    """Commands run without a shell, so a filename with metacharacters is inert."""
    weird = tmp_path / "a b; touch pwned.txt $(whoami).fasta"
    weird.write_text(">x\nACGT\n")
    # the path is a separate argv entry, not part of the -c source
    cmd = [*_python("import sys", "print(sys.argv[1])"), str(weird)]
    assert run_tool(CAPS, cmd, logger=LOG).strip() == str(weird)
    assert not (tmp_path / "pwned.txt").exists()
