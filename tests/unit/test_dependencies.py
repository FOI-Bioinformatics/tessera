"""Every declared runtime dependency is one the package imports.

Tessera is dependency-light by design, so a declared dependency nothing imports is a
cost with no benefit: it is installed for every user and audited in CI.
"""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

# Distribution name -> import name, where they differ.
IMPORT_NAME = {"biopython": "Bio"}


def test_every_runtime_dependency_is_imported() -> None:
    project = tomllib.loads((REPO / "pyproject.toml").read_text())["project"]
    names = [re.split(r"[<>=!~\[ ;]", dep, maxsplit=1)[0] for dep in project["dependencies"]]
    source = "\n".join(p.read_text() for p in (REPO / "src" / "tessera").rglob("*.py"))
    unused = [
        name for name in names
        if not re.search(rf"^\s*(?:import|from)\s+{IMPORT_NAME.get(name, name)}\b",
                         source, flags=re.M)
    ]
    assert unused == []
