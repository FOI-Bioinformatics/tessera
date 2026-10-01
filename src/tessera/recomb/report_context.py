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

from .analyze import AnalysisResult
from .coverage import BREAKPOINT_KIND, CoverageGap
from .diagnostics import RecombinationSignal
from .similarity import WindowSimilarity
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
    # Selected callers that could not run (barcode on an untyped panel). They stay in
    # ``methods_run`` so the method table keeps a column for them, marked "not run".
    methods_not_run: tuple[str, ...] = ()
    method_breakdown: list[dict] | None = None
    per_major: dict[str, str] | None = None
    # Set only when the scan used informative-site windowing: the per-window identity
    # at polymorphic columns that the HMM segmented, and its analysis. Reported as its
    # own track and file, never under the name "similarity" -- it is on a different
    # scale from identity over all columns.
    site_result: WindowSimilarity | None = None
    site_analysis: AnalysisResult | None = None

    @property
    def gaps(self) -> list[CoverageGap]:
        """``coverage_gaps`` with ``None`` normalised to an empty list."""
        return self.coverage_gaps or []

    @property
    def caveat_gaps(self) -> list[CoverageGap]:
        """The gaps that may mean a missing reference: every kind except ``breakpoint``.

        A breakpoint gap is a window straddling a called breakpoint. It is listed in the
        coverage table and ``coverage_gaps.tsv`` under its own kind, but it is not a
        poorly covered stretch, so the headline caveat, the mosaic and the plots leave
        it out.
        """
        return [g for g in self.gaps if g.kind != BREAKPOINT_KIND]
