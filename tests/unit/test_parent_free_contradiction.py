"""When the parent-free diagnostic positively contradicts the callers (issue #84).

On a clonal mpox tree tip the callers reported two multi-kilobase regions while the
parent-free diagnostic on the same alignment returned ``PHI p=1, Rmin=0``: the maximum
p-value, and a Hudson-Kaplan minimum of zero, which says no recombination event is
*needed* to explain the data. The region rows already carried
``parent_free_support = no``, but that reads the same whether the diagnostic found
nothing or was never able to look -- and those are different statements.

Why this and not the margin the issue proposed. Measured across the 43 hybrid
validation cases, a true mpox region carries ``margin = 0.0003``, the same as the false
call's, and a four-caller vzv region 0.0013; any margin bar that removes the false call
removes those too. ``Rmin == 0`` separates them cleanly instead: no true positive in
those cases has it, the minimum being 6 where the diagnostic ran at all.
"""

from __future__ import annotations

from tessera.recomb.diagnostics import RecombinationSignal, contradicts_regions
from tessera.recomb.regions import Region


def _signal(phi_p: float | None, rmin: int, n_informative: int = 17) -> RecombinationSignal:
    return RecombinationSignal(
        n_informative=n_informative, phi_p=phi_p, phi_observed=0.5,
        phi_window=100, rmin=rmin,
    )


def test_the_issue_84_signal_is_a_contradiction() -> None:
    """PHI at its maximum and no event required: the diagnostic says there is nothing."""
    assert contradicts_regions(_signal(phi_p=1.0, rmin=0), alpha=0.05) is True


def test_a_significant_phi_is_not_a_contradiction() -> None:
    """The diagnostic agrees that the alignment carries recombination."""
    assert contradicts_regions(_signal(phi_p=0.000999, rmin=8), alpha=0.05) is False


def test_a_nonsignificant_phi_with_events_required_is_not_a_contradiction() -> None:
    """Rmin above zero means some event is needed, whatever PHI's power was.

    Two of the hybrid positives sit here (PHI p=0.015 and p=0.045 with Rmin 149 and
    76), so this must not count as a contradiction.
    """
    assert contradicts_regions(_signal(phi_p=0.20, rmin=76), alpha=0.05) is False


def test_a_diagnostic_that_could_not_run_is_not_a_contradiction() -> None:
    """``phi_p is None`` means the test could not have rejected, not that it did not.

    Three hybrid positives have no PHI result at all (too few informative sites) and
    are correctly detected; treating silence as contradiction would caveat them.
    """
    assert contradicts_regions(_signal(phi_p=None, rmin=0), alpha=0.05) is False


def test_the_boundary_follows_alpha() -> None:
    """Judged against the run's own alpha rather than a second, hidden threshold."""
    assert contradicts_regions(_signal(phi_p=0.04, rmin=0), alpha=0.05) is False
    assert contradicts_regions(_signal(phi_p=0.06, rmin=0), alpha=0.05) is True


def _region() -> Region:
    return Region(
        minor_parent="donor", major_parent="backbone",
        msa_start=0, msa_end=1000, query_start=0, query_end=1000,
        n_windows=5, mean_sim_minor=0.9998, mean_sim_major=0.9995, margin=0.0003,
        qvalue=0.012, support=0.883, methods=("hmm", "bootscan"),
    )


def test_the_report_caveats_the_regions_rather_than_hiding_them() -> None:
    """The rows stay; what changes is that the disagreement is stated.

    Not suppression: mpox does recombine, so a tree tip can carry a real event, and the
    diagnostic is weak on a panel with few informative sites. Removing the rows on its
    word would trade a false positive for a possible false negative.
    """
    from tessera.recomb.report_html import _contradiction_html

    note = _contradiction_html(True, 2)

    assert "does not support" in note
    assert "no" in note and "required" in note
    assert "disagreement" in note


def test_there_is_no_note_without_a_contradiction_or_without_regions() -> None:
    from tessera.recomb.report_html import _contradiction_html

    assert _contradiction_html(False, 2) == ""
    assert _contradiction_html(True, 0) == ""


def test_the_contradiction_caps_the_confidence_verdict() -> None:
    """Whatever the per-region statistics say, the call is not high-confidence."""
    from tessera.recomb.report_html import _confidence

    region = _region()
    region.methods = ("hmm", "3seq", "maxchi", "bootscan")
    region.qvalue = 1e-30

    assert _confidence([region], [], False) == "high"
    assert _confidence([region], [], True) == "low"
