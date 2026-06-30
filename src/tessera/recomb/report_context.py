"""The optional metadata and presentation inputs shared by the report writers.

``write_reports`` and ``write_html_report`` both take a backbone of core scan
results (the similarity matrix, the analysis, the called regions, provenance)
plus a long tail of optional metadata -- coverage, the parent-free signal, the
lineage map, the ensemble breakdown, the organism name. That tail kept growing
one threaded argument at a time. ``ReportContext`` bundles it so a new report
field is added here once, not in every writer's signature and call site.
"""

from __future__ import annotations

from dataclasses import dataclass

from .coverage import CoverageGap
from .diagnostics import RecombinationSignal
from .typing import LineageMap


@dataclass
class ReportContext:
    """Optional inputs for the report writers, beyond the core scan results."""

    coverage_gaps: list[CoverageGap] | None = None
    coverage_threshold: float = 0.0
    extra_sections: list[tuple[str, str]] | None = None
    lineage_map: LineageMap | None = None
    query_lineage: str | None = None
    signal: RecombinationSignal | None = None
    organism: str | None = None
    methods_run: tuple[str, ...] = ()
    method_breakdown: list[dict] | None = None
    per_major: dict[str, str] | None = None

    @property
    def gaps(self) -> list[CoverageGap]:
        """``coverage_gaps`` with ``None`` normalised to an empty list."""
        return self.coverage_gaps or []
