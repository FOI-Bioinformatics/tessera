"""The optional metadata and presentation inputs shared by the report writers.

``write_reports`` and ``write_html_report`` both take a backbone of core scan
results (the similarity matrix, the analysis, the called regions, provenance)
plus a long tail of optional metadata -- coverage, the parent-free signal, the
lineage map, the ensemble breakdown, the organism name. That tail kept growing
one threaded argument at a time. ``ReportContext`` bundles it so a new report
field is added here once, not in every writer's signature and call site.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .analyze import AnalysisResult
from .coverage import (
    BREAKPOINT_KIND,
    DIVERGENT_KIND,
    LOCAL_DIP_KIND,
    LOW_INFORMATION_KIND,
    CoverageGap,
)
from .diagnostics import RecombinationSignal
from .regions import Region
from .similarity import WindowSimilarity
from .typing import LineageMap

# The key under which the parent-free PHI test is named in ``ReportContext.not_run``.
PHI = "phi"


@dataclass
class ReportContext:
    """Optional inputs for the report writers, beyond the core scan results."""

    coverage_gaps: list[CoverageGap] | None = None
    coverage_threshold: float = 0.0
    extra_sections: list[tuple[str, str]] | None = None
    lineage_map: LineageMap | None = None
    query_lineage: str | None = None
    signal: RecombinationSignal | None = None
    # The run's significance level, so the report judges the PHI p-value at the same
    # alpha the per-region corroboration used.
    alpha: float = 0.05
    organism: str | None = None
    methods_run: tuple[str, ...] = ()
    # Callers and tests that produced no result, by name, with the reason in the run's
    # own words (the text the log gives). A caller named here stays in ``methods_run``,
    # so the method table keeps a column for it, marked "not run"; ``PHI`` is the
    # parent-free test. Whatever is named here is reported as "not run", never as a
    # negative.
    not_run: dict[str, str] = field(default_factory=dict)
    # Regions set aside because the backbone is locally gap-degraded over them, so the
    # comparison that produced them was not between comparably aligned genomes
    # (issue #81). Reported rather than silently dropped.
    withheld_regions: list[Region] = field(default_factory=list)
    inverted_regions: list[Region] = field(default_factory=list)
    # The parent-free diagnostic ran and found that no recombination event is
    # required to explain the alignment, while regions were nonetheless called:
    # a contradiction rather than an absence of corroboration (issue #84).
    parent_free_contradicts: bool = False
    method_breakdown: list[dict] | None = None
    per_major: dict[str, str] | None = None
    # The run's reconciled major parent -- the one the provenance records and the
    # region rows are relabelled to. The headline names it; ``None`` when no caller
    # reported one, and the report then falls back to the leading dataset.
    major_parent: str | None = None
    # Set only when the scan used informative-site windowing: the per-window identity
    # at polymorphic columns that the HMM segmented, and its analysis. Reported as its
    # own track and file, never under the name "similarity" -- it is on a different
    # scale from identity over all columns.
    site_result: WindowSimilarity | None = None
    site_analysis: AnalysisResult | None = None

    @property
    def methods_not_run(self) -> tuple[str, ...]:
        """The selected callers that could not run, in ``methods_run`` order."""
        return tuple(m for m in self.methods_run if m in self.not_run)

    @property
    def phi_not_run(self) -> str | None:
        """Why the parent-free test gave no result, or ``None`` when it did."""
        return self.not_run.get(PHI)

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

    @property
    def missing_reference_gaps(self) -> list[CoverageGap]:
        """The gaps that may mean a missing reference in the headline: ``divergent`` only.

        A ``low_information`` gap has too few comparable bases to judge -- an unsequenced
        or masked stretch of the query, say -- so it says nothing about the panel and is
        given its own note instead of the advice to search for a reference.
        """
        return [g for g in self.gaps if g.kind == DIVERGENT_KIND]

    @property
    def low_information_gaps(self) -> list[CoverageGap]:
        """The gaps with too few comparable bases to judge."""
        return [g for g in self.gaps if g.kind == LOW_INFORMATION_KIND]

    @property
    def local_dip_gaps(self) -> list[CoverageGap]:
        """The stretches the closest reference matches less well than elsewhere.

        It still matches them closely in absolute terms, so these are neither a gap
        in the panel nor a stretch too uncertain to judge: usually a locally poor
        alignment. They get their own note (issue #82).
        """
        return [g for g in self.gaps if g.kind == LOCAL_DIP_KIND]
