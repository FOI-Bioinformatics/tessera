"""Options read by only some callers say so in their help.

Audit 4, items B2, B6 and A4: ``--margin`` and ``--merge-gap`` are read by the heuristic
caller only, and ``--min-region`` bounds every reported region only when given, but their
help did not say so. A user who set one and saw
no change could not tell a no-op from a null result. Introspects the parameters rather
than rendered ``--help``, which wraps differently across typer/rich/click versions.
"""

from __future__ import annotations

import typer

from tessera.cli.main import app


def _help(command: str, option: str) -> str:
    cmd = typer.main.get_command(app).commands[command]
    param = next(p for p in cmd.params if option in getattr(p, "opts", []))
    return " ".join((param.help or "").split())


def test_min_region_names_its_callers_and_what_it_bounds() -> None:
    text = _help("recomb", "--min-region")
    assert "whichever caller" in text


def test_margin_and_merge_gap_are_marked_heuristic_only() -> None:
    assert "(heuristic only;" in _help("recomb", "--margin")
    assert "a fraction, 0-1" in _help("recomb", "--margin")
    assert "(heuristic only;" in _help("recomb", "--merge-gap")
