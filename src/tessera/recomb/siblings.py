"""Sibling-aware reference selection for the recombination caller.

A masking **sibling** of the query carries the same recombination mosaic, so it is
the closest reference in (almost) every window: it explains the whole query alone
and, made the major parent, hides the event -- no minor region can beat a reference
that already matches the query everywhere. A true **parent** donated only one
region: it leads there but is *contiguously beaten* elsewhere by the co-parent or
backbone.

The distinction is **structural and relative**, not an absolute identity cutoff, so
it holds for HIV subtypes (~12% apart) and SARS-CoV-2 sublineages (<1% apart)
alike. A whole-genome sibling:

- **covers** nearly the whole query (high informative fraction) -- this excludes a
  partial, region-only parent such as an env-only reference, which leads its own
  windows but is not a whole-genome explainer;
- is significantly **beaten** in only a negligible fraction of windows, and never
  over a region-sized contiguous run. This is the defining signature: nothing is
  closer to the query than its own twin, whereas a true parent is contiguously
  beaten over the co-parent's donor tract. (A sibling need not uniquely *lead* every
  window -- it may merely tie other references where they share one of its regions,
  e.g. a recombinant twin and a pure-capsid parent both matching the capsid.)

Siblings are removed from the HMM's competition only when doing so reveals a
significant mosaic (the existing discordant-site sign test); otherwise the full
panel is kept, so a genuinely non-recombinant query is unaffected.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .hmm import segment_query
from .similarity import WindowSimilarity, discordant_counts
from .stats import proportion_diff_significant, sign_test_pvalue

_MAX_PEEL = 3  # at most this many whole-genome siblings removed (bounds the cost)


@dataclass
class SiblingEvidence:
    label: str
    coverage: float  # fraction of windows where the reference is comparable
    lead_frac: float  # fraction of its informative windows where it significantly leads
    beaten_frac: float  # fraction where another reference significantly beats it
    max_beaten_run: int  # longest contiguous run (windows) of being beaten


def _profiles(
    result: WindowSimilarity, labels: list[str], alpha: float
) -> dict[str, SiblingEvidence]:
    """Per-label leadership/beaten profile, computed among ``labels`` only."""
    n_win = len(result.positions)
    num, den = result.numerators, result.denominators
    leads = dict.fromkeys(labels, 0)
    beaten = dict.fromkeys(labels, 0)
    informative = dict.fromkeys(labels, 0)
    run = dict.fromkeys(labels, 0)
    max_run = dict.fromkeys(labels, 0)
    for i in range(n_win):
        comp = [(label, num[label][i], den[label][i]) for label in labels if den[label][i] > 0]
        for label, _, _ in comp:
            informative[label] += 1
        beaten_now: set[str] = set()
        if comp:
            comp.sort(key=lambda t: t[1] / t[2], reverse=True)
            top, k1, n1 = comp[0]
            # A unique significant leader (used only to order which sibling to peel).
            if len(comp) >= 2:
                k2, n2 = comp[1][1], comp[1][2]
                if proportion_diff_significant(k1, n1, k2, n2, alpha):
                    leads[top] += 1
            # A reference is beaten when the window's best significantly beats it --
            # independent of whether the best also beats the runner-up. (Otherwise two
            # tied top references would leave a distant third looking "never beaten".)
            for label, k, n in comp[1:]:
                if proportion_diff_significant(k1, n1, k, n, alpha):
                    beaten[label] += 1
                    beaten_now.add(label)
        for label in labels:
            if label in beaten_now:
                run[label] += 1
                max_run[label] = max(max_run[label], run[label])
            else:
                run[label] = 0
    out: dict[str, SiblingEvidence] = {}
    for label in labels:
        inf = informative[label] or 1
        out[label] = SiblingEvidence(
            label=label,
            coverage=informative[label] / n_win if n_win else 0.0,
            lead_frac=leads[label] / inf,
            beaten_frac=beaten[label] / inf,
            max_beaten_run=max_run[label],
        )
    return out


def _window_step(result: WindowSimilarity) -> int:
    pos = result.positions
    return (pos[1] - pos[0]) if len(pos) > 1 else 1


def _is_sibling(ev: SiblingEvidence, params, min_run_windows: int) -> bool:
    """A whole-genome sibling covers the query and is never significantly beaten.

    The defining signature is *never being beaten over a region-sized run*: nothing
    is closer to the query than its own twin. A true parent is contiguously beaten
    over the co-parent's donated tract, so it fails ``max_beaten_run``. Leadership is
    not required -- a sibling may merely tie other references where they share one of
    its regions (so ``lead_frac`` only orders which explainer to peel first).
    """
    return (
        ev.coverage >= params.sibling_coverage
        and ev.beaten_frac <= params.sibling_beaten
        and ev.max_beaten_run < min_run_windows
    )


def _reveals_mosaic(result: WindowSimilarity, states: list[str], params) -> bool:
    """True if the HMM on ``states`` yields a significant non-major donor region."""
    segments, major = segment_query(
        result, identity=params.identity, jump_rate=params.jump_rate, states=states
    )
    if major is None:
        return False
    for seg in segments:
        if seg.state == major or seg.msa_end - seg.msa_start < params.min_region:
            continue
        favor_minor, favor_major = discordant_counts(
            result.rows, result.query, major, seg.state, seg.msa_start, seg.msa_end
        )
        if favor_minor > favor_major and sign_test_pvalue(favor_minor, favor_major) <= params.alpha:
            return True
    return False


def sibling_aware_states(
    result: WindowSimilarity, params
) -> tuple[list[str], list[SiblingEvidence], dict[str, SiblingEvidence] | None]:
    """Return ``(states, dropped, profiles)``: the reference set to compete after
    removing whole-genome siblings of the query.

    A sibling is removed only if its removal exposes a sign-test-significant mosaic
    among the survivors; otherwise the full label set is returned unchanged (the
    non-recombinant guard). At most :data:`_MAX_PEEL` siblings are peeled.
    """
    labels = list(result.similarities)
    if len(labels) < 3:
        return labels, [], None
    min_run_windows = max(1, math.ceil(params.min_region / _window_step(result)))

    survivors = list(labels)
    dropped: list[SiblingEvidence] = []
    profiles = _profiles(result, survivors, params.sibling_alpha)
    for _ in range(_MAX_PEEL):
        candidates = [
            profiles[label] for label in survivors
            if _is_sibling(profiles[label], params, min_run_windows)
        ]
        if not candidates:
            break
        # Remove the strongest whole-genome explainer first.
        candidates.sort(key=lambda ev: ev.lead_frac, reverse=True)
        target = candidates[0].label
        trial = [label for label in survivors if label != target]
        if len(trial) < 2 or not _reveals_mosaic(result, trial, params):
            break  # removal does not expose a real mosaic -> keep the panel as is
        survivors = trial
        dropped.append(profiles[target])
        profiles = _profiles(result, survivors, params.sibling_alpha)
    return survivors, dropped, profiles
