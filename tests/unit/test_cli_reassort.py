"""The reassort command is registered and exposes its query option (no binaries)."""

from __future__ import annotations

from typer.testing import CliRunner

from tessera.cli.main import app

runner = CliRunner()


def test_reassort_listed_in_help():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "reassort" in result.stdout


def test_reassort_help_shows_query_option():
    result = runner.invoke(app, ["reassort", "--help"])
    assert result.exit_code == 0
    assert "--query" in result.stdout
