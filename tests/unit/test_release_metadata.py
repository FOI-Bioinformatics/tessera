"""The version is stated in three places, and they must agree.

``CITATION.cff`` was left at 1.2.0 through the 1.3.0 release: nothing compared it with
the package version. The release workflow now checks the tag against it as well, and this
test catches the drift earlier, at the commit that bumps ``__version__``.
"""

from __future__ import annotations

import re
from pathlib import Path

import tessera

ROOT = Path(__file__).resolve().parents[2]


def _cff_version() -> str:
    match = re.search(r"^version:\s*(\S+)\s*$", (ROOT / "CITATION.cff").read_text(), re.M)
    assert match, "CITATION.cff has no version line"
    return match.group(1)


def test_citation_version_matches_the_package() -> None:
    assert _cff_version() == tessera.__version__


def test_citation_records_a_release_date() -> None:
    text = (ROOT / "CITATION.cff").read_text()
    assert re.search(r"^date-released:\s*\d{4}-\d{2}-\d{2}\s*$", text, re.M)


def test_the_changelog_has_a_section_for_the_current_version() -> None:
    text = (ROOT / "CHANGELOG.md").read_text()
    heading = rf"^## \[{re.escape(tessera.__version__)}\] - \d{{4}}-\d{{2}}-\d{{2}}"
    assert re.search(heading, text, re.M)
