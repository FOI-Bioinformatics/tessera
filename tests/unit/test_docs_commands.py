"""The development commands the documentation gives are the ones that work (audit 4, F2-F4).

Prepending the aligner environment to PATH also puts that environment's Python first,
so the harness ran under an interpreter without Tessera installed. The README lint line
left out ``validation``, which CI lints, and CONTRIBUTING.md said mypy was not blocking.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
DOCS = ["CLAUDE.md", "CONTRIBUTING.md", "README.md", "validation/README.md"]
PREPEND = re.compile(r'export PATH="\$HOME/[^"]*/bin:\$PATH"')


@pytest.mark.parametrize("name", DOCS)
def test_no_document_prepends_an_environment_to_path(name: str) -> None:
    text = (ROOT / name).read_text()
    assert not PREPEND.search(text), PREPEND.search(text).group(0)


def test_readme_lints_what_ci_lints() -> None:
    assert "ruff check src tests validation" in (ROOT / "README.md").read_text()


def test_contributing_states_that_mypy_is_blocking() -> None:
    text = (ROOT / "CONTRIBUTING.md").read_text()
    assert "not yet blocking" not in text
