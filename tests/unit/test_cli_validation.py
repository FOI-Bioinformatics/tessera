"""CLI argument validation and the ``python -m`` entry point.

An out-of-range option used to reach the statistics and surface as an internal
exception under "Unexpected error", which reads like a defect in Tessera rather than
a bad argument.
"""

from __future__ import annotations

import subprocess
import sys

import pytest
from typer.testing import CliRunner

from tessera.cli.main import _require_range, app
from tessera.core.errors import UserInputError

runner = CliRunner()


# --- the range helper -----------------------------------------------------

def test_require_range_accepts_values_inside_the_bounds() -> None:
    _require_range(0.05, "--alpha", lo=0.0, hi=1.0, lo_open=True, hi_open=True)
    _require_range(1, "--top-n", lo=1)  # inclusive lower bound
    _require_range(0.0, "--margin", lo=0.0, hi=1.0)


@pytest.mark.parametrize(
    ("value", "kwargs", "expected"),
    [
        (0.0, {"lo": 0.0, "hi": 1.0, "lo_open": True, "hi_open": True}, "> 0 and < 1"),
        (1.0, {"lo": 0.0, "hi": 1.0, "lo_open": True, "hi_open": True}, "> 0 and < 1"),
        (0, {"lo": 1}, ">= 1"),
        (1.5, {"lo": 0.0, "hi": 1.0}, ">= 0 and <= 1"),
    ],
)
def test_require_range_rejects_and_states_the_bounds(value, kwargs, expected) -> None:
    with pytest.raises(UserInputError, match="Invalid --opt"):
        _require_range(value, "--opt", **kwargs)
    try:
        _require_range(value, "--opt", **kwargs)
    except UserInputError as exc:
        assert expected in str(exc)


# --- through the CLI ------------------------------------------------------

@pytest.mark.parametrize(
    ("option", "value"),
    [("--alpha", "0"), ("--alpha", "1"), ("--jump-rate", "1.5"),
     ("--top-n", "0"), ("--min-methods", "0"), ("--phi-window", "0")],
)
def test_out_of_range_options_exit_cleanly(tmp_path, option, value) -> None:
    result = runner.invoke(app, [
        "recomb", "--msa", "example_data/divergent.msa.fasta", "--query", "query",
        "--output", str(tmp_path / "out"), option, value,
    ])
    assert result.exit_code == 1  # a clean user error, not a crash
    assert "Unexpected error" not in result.output


# --- entry points ---------------------------------------------------------

def test_python_m_tessera_registers_subcommands() -> None:
    """`python -m tessera.cli.main` bound a second copy of the module, so every
    subcommand registered on a different app than the one `-m` invoked."""
    done = subprocess.run(
        [sys.executable, "-m", "tessera", "curate-panel", "--help"],
        capture_output=True, text=True, timeout=60,
    )
    assert done.returncode == 0, done.stderr
    assert "No such command" not in done.stdout + done.stderr
    assert "--derep-ani" in done.stdout


def test_help_text_renders_percent_signs_literally() -> None:
    # Click does not %-format help text, so "%%" reached the terminal verbatim.
    result = runner.invoke(app, ["curate-panel", "--help"])
    assert result.exit_code == 0
    assert "%%" not in result.output
    assert "% of the query" in result.output
