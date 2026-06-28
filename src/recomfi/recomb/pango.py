"""Pango lineage aliases: the designated parents of a SARS-CoV-2 recombinant lineage.

The cov-lineages/pango-designation project publishes ``alias_key.json``, which maps
each aliased lineage to its full form -- and each recombinant (``X``-prefixed) lineage
to the *list* of parental lineages it descends from (e.g. ``XBB -> [BJ.1, CJ.1]``).
RecomFi uses it as a ground-truth parent map to cross-check the parents it recruited
for a recombinant query against the designated ones.

The file is fetched once into the on-disk cache and reused; a network failure is
non-fatal (the cross-check is simply skipped, with a warning).
"""

from __future__ import annotations

import html
import json
import logging
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen

from ..core.cache import pango_alias_path

ALIAS_KEY_URL = (
    "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/"
    "pango_designation/alias_key.json"
)


def load_alias_key(
    *, cache_override: str | Path | None = None, logger: logging.Logger, offline: bool = False
) -> dict[str, object]:
    """Load ``alias_key.json`` from the cache, fetching it once if absent.

    Returns ``{}`` when the file is unavailable (offline with no cache, or a network
    or parse error) so callers can treat the cross-check as best-effort.
    """
    path = pango_alias_path(cache_override)
    if path.exists():
        try:
            return json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            logger.warning("Cached Pango alias key is unreadable (%s); ignoring.", exc)
            return {}
    if offline:
        return {}
    try:
        with urlopen(ALIAS_KEY_URL, timeout=30) as resp:  # noqa: S310 - fixed https URL
            data = resp.read()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        logger.info("Fetched the Pango alias key into the cache: %s", path)
        return json.loads(data)
    except (URLError, OSError, json.JSONDecodeError) as exc:
        logger.warning("Could not fetch the Pango alias key: %s", exc)
        return {}


def expand_recombinant(lineage: str, alias_key: dict[str, object]) -> list[str]:
    """The parental lineages of a recombinant lineage, or ``[]`` if it is not one.

    In ``alias_key.json`` a recombinant maps to a *list* of parents (``XBB`` ->
    ``[BJ.1, CJ.1]``) while an ordinary alias maps to a single string. A sublineage
    (``XBB.1.5``) resolves via its ``X`` root.
    """
    if not lineage:
        return []
    root = lineage.split(".")[0]
    if not root.startswith("X"):
        return []
    value = alias_key.get(root)
    return [str(p) for p in value] if isinstance(value, list) and value else []


def crosscheck_html(query_lineage: str, parents: list[str]) -> str:
    """A report block stating the Pango-designated parents of a recombinant query."""
    parent_list = ", ".join(f"<strong>{html.escape(p)}</strong>" for p in parents)
    return (
        f'<p class="cap">Pango designates the query lineage '
        f'<strong>{html.escape(query_lineage)}</strong> as a recombinant of {parent_list}. '
        f'Compare these with the donor regions called above: a match supports the call, a '
        f'mismatch suggests a recruitment gap or a different breakpoint structure.</p>'
    )
