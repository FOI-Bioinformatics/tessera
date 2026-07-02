"""The reassort command is registered and exposes its query option (no binaries).

Introspects the command's parameters rather than parsing rendered ``--help`` text,
which wraps at the terminal width and varies across typer/rich/click versions.
"""

from __future__ import annotations

import typer

from tessera.cli.main import app


def _reassort_command():
    return typer.main.get_command(app).commands["reassort"]


def test_reassort_command_registered():
    assert "reassort" in typer.main.get_command(app).commands


def test_reassort_exposes_query_option():
    names = {opt for p in _reassort_command().params for opt in getattr(p, "opts", [])}
    assert "--query" in names
    assert "-q" in names


def test_reassort_exposes_margin_option():
    names = {opt for p in _reassort_command().params for opt in getattr(p, "opts", [])}
    assert "--margin" in names


def test_reassort_exposes_scan_segments_option():
    names = {opt for p in _reassort_command().params for opt in getattr(p, "opts", [])}
    assert "--scan-segments" in names


def test_reassort_exposes_aligner_option():
    names = {opt for p in _reassort_command().params for opt in getattr(p, "opts", [])}
    assert "--aligner" in names
