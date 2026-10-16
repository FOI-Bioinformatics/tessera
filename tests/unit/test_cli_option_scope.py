"""Options read by only some callers say so in their help.

Audit 4, items B2, B6 and A4: ``--min-region`` is read by the HMM (on the segment, before
the span is trimmed) and the heuristic caller only, and ``--margin`` and ``--merge-gap``
by the heuristic caller only, but their help did not say so. A user who set one and saw
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
    assert "(hmm, heuristic;" in text
    assert "A reported HMM region can be shorter" in text
    assert "3seq, maxchi, bootscan, geneconv and barcode do not apply it" in text


def test_margin_and_merge_gap_are_marked_heuristic_only() -> None:
    assert "(heuristic only;" in _help("recomb", "--margin")
    assert "a fraction, 0-1" in _help("recomb", "--margin")
    assert "(heuristic only;" in _help("recomb", "--merge-gap")
