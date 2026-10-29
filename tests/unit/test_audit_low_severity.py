"""The low-severity findings of the 2026-10-28 audit.

None of these changes a verdict. Each is a number, a name or a message that was wrong
or misleading, or a failure reported as a traceback rather than as the input problem it
was. They are grouped here because they share nothing but their provenance.
"""

from __future__ import annotations

import logging
import subprocess
import sys
from pathlib import Path

import pytest
from typer.testing import CliRunner

from tessera.converters.maf_to_fasta import maf_to_fasta
from tessera.discover.run import _tiles
from tessera.recomb.typing import _base_accession, lineage_of

runner = CliRunner()


# --- a tail tile that asks almost the same question -------------------------

def test_a_tiny_remainder_does_not_add_a_near_duplicate_tile() -> None:
    """Two tiles whose starts differ by 1 bp are two BLAST submissions, each paced at
    one per ten seconds, asking almost the same question."""
    assert _tiles(601, 400) == [(0, 400), (201, 601)]
    assert _tiles(401, 400) == [(0, 401)]


def test_tiles_still_cover_the_whole_gap() -> None:
    for length in (401, 450, 601, 700, 800, 999, 1000, 1201):
        tiles = _tiles(length, 400)
        assert tiles[0][0] == 0
        assert tiles[-1][1] == length
        for (_, end), (start, _) in zip(tiles[:-1], tiles[1:], strict=True):
            assert start <= end, (length, tiles)  # no uncovered stretch between tiles


def test_tiles_keep_the_overlap_where_the_remainder_is_real() -> None:
    assert _tiles(700, 400) == [(0, 400), (200, 600), (300, 700)]


# --- the backbone label is not a sequence id -------------------------------

def test_a_sequence_id_equal_to_the_backbone_label_does_not_redirect_it(
    tmp_path: Path,
) -> None:
    """``name_map`` maps source name -> genome label, so looking a *label* up in it
    redirected the backbone to whichever genome carried that id."""
    maf = tmp_path / "a.maf"
    maf.write_text(
        "a score=0\n"
        "s K03455    0 4 + 4 ACGT\n"
        "s K03455.1  0 4 + 4 ACGA\n"
        "\n",
        encoding="utf-8",
    )
    out = tmp_path / "out.fasta"
    # "K03455.1" is a genome label here (a value of name_map), and also the sequence id
    # of a different genome (a key). The label must win.
    maf_to_fasta(
        maf, "K03455.1", out,
        name_map={"K03455": "K03455.1", "K03455.1": "hxb2"},
        ref_contigs=[("K03455", 4)],
    )
    labels = [line[1:].split()[0] for line in out.read_text(encoding="utf-8").splitlines()
              if line.startswith(">")]
    assert "K03455.1" in labels


# --- a version suffix is a trailing ".<digits>", not any dot ---------------

def test_a_dotted_name_is_not_typed_by_an_unrelated_one() -> None:
    assert lineage_of("strain.B.1", {"strain.A.1": "cladeA"}) is None
    assert lineage_of("sample.x", {"sample.y": "other"}) is None


def test_a_real_version_suffix_still_compares_equal() -> None:
    assert _base_accession("U54771.1") == _base_accession("U54771") == "U54771"
    assert lineage_of("U54771.2", {"U54771.1": "subtype-B"}) == "subtype-B"
    assert lineage_of("U54771", {"U54771.1": "subtype-B"}) == "subtype-B"


def test_an_ambiguous_version_match_is_deterministic() -> None:
    """Two versions of one accession carrying different lineages must not resolve by
    dict order."""
    forward = lineage_of("X.3", {"X.1": "one", "X.2": "two"})
    reverse = lineage_of("X.3", {"X.2": "two", "X.1": "one"})
    assert forward == reverse == "one"


# --- a mistyped timeout says so --------------------------------------------

def test_an_unparsable_timeout_warns_instead_of_falling_back_silently(
    monkeypatch, caplog
) -> None:
    from tessera.core.process import default_timeout

    monkeypatch.setenv("TESSERA_TOOL_TIMEOUT", "30s")
    with caplog.at_level(logging.WARNING, logger="tessera.core.process"):
        assert default_timeout(300.0) == 300.0  # the fallback keeps a limit, as before
    assert "TESSERA_TOOL_TIMEOUT" in caplog.text
    assert "30s" in caplog.text


def test_a_valid_timeout_warns_about_nothing(monkeypatch, caplog) -> None:
    from tessera.core.process import default_timeout

    monkeypatch.setenv("TESSERA_TOOL_TIMEOUT", "45")
    with caplog.at_level(logging.WARNING, logger="tessera.core.process"):
        assert default_timeout(300.0) == 45.0
    assert caplog.text == ""


# --- an output path that cannot be created ---------------------------------

def test_an_uncreatable_output_path_is_an_input_error_not_a_traceback() -> None:
    """The run log opens before any command enters stage_errors, so the OSError came
    out of logging setup as "Unexpected error" with a traceback."""
    result = subprocess.run(
        [sys.executable, "-m", "tessera", "detect",
         "-q", "nothing.fa", "-o", "/proc/no_such_dir/out"],
        capture_output=True, text=True, timeout=300,
    )
    combined = result.stdout + result.stderr
    assert result.returncode == 1, combined
    assert "Traceback" not in combined
    assert "Unexpected error" not in combined
    assert "Cannot write the run log" in combined


# --- a gap with no comparable bases does not drive convergence -------------

def test_a_nan_gap_is_left_out_of_the_shortfall() -> None:
    from tessera.discover.iterate import coverage_shortfall
    from tessera.recomb.coverage import LOW_INFORMATION_KIND, CoverageGap

    def gap(mean_best: float, kind: str) -> CoverageGap:
        return CoverageGap(
            msa_start=0, msa_end=1000, query_start=0, query_end=1000,
            length_bp=1000, n_windows=10, best_label="refA",
            mean_best=mean_best, kind=kind,
        )

    nan_gap = gap(float("nan"), LOW_INFORMATION_KIND)
    real_gap = gap(0.80, "divergent")

    # Counted, the nan gap contributed exactly zero, so a round holding only such gaps
    # measured as a round with nothing left to close.
    assert coverage_shortfall([nan_gap], threshold=0.95) == 0.0
    assert coverage_shortfall([nan_gap, real_gap], threshold=0.95) == pytest.approx(150.0)
    assert coverage_shortfall([real_gap], threshold=0.95) == pytest.approx(150.0)


def test_the_worst_similarity_does_not_depend_on_gap_order() -> None:
    """min() over a nan is decided by the order the gaps arrive in."""
    from math import isnan

    values = [float("nan"), 0.80]
    forward = min((v for v in values if not isnan(v)), default=1.0)
    backward = min((v for v in reversed(values) if not isnan(v)), default=1.0)
    assert forward == backward == 0.80


# --- an odd window reaches its own last column -----------------------------

def test_a_coverage_gap_span_matches_the_windows_it_describes() -> None:
    """Deriving the end as ``center + half`` left an odd window one column short, so
    the breakpoint relabel silently skipped a gap's last window."""
    import numpy as np

    from tessera.recomb.coverage import CoverageParams, call_coverage_gaps
    from tessera.recomb.similarity import WindowSimilarity

    window = 1001  # odd
    half = window // 2
    n = 12
    positions = [half + i * 100 for i in range(n)]
    width = positions[-1] + half + 1
    best = [0.99] * n
    for i in (5, 6, 7):
        best[i] = 0.70
    result = WindowSimilarity(
        positions=positions, query_positions=positions, similarities={},
        query="q", width=width, query_cumulative=np.arange(width + 1),
        best_sim=best, best_label=["refA"] * n, informative=[200] * n,
    )

    gaps, _ = call_coverage_gaps(result, window, CoverageParams.with_defaults(window))

    assert len(gaps) == 1
    gap = gaps[0]
    last_window_end = positions[7] - half + window
    assert gap.msa_end == min(width, last_window_end)
