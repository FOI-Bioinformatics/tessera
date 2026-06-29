"""HMM segmentation of the query against the reference panel (jpHMM-style).

The query is modelled as a hidden chain of "which reference is the local source".
Each window emits, for every reference state, the per-site copying log-likelihood
(:func:`tessera.recomb.stats.emission_loglik`) of its binomial counts; a small
per-window jump rate penalises switching reference. Viterbi gives the maximum a
posteriori segmentation (the parental mosaic) and forward-backward gives the
per-window posterior, from which breakpoint uncertainty intervals are derived.

This replaces the greedy run-merge heuristic: the jump penalty -- one
interpretable parameter -- subsumes the old margin / min-region / merge-gap knobs.
A switch happens only when the accumulated emission gain beats the jump cost, so
near-identical references do not trigger spurious flips and thin windows cannot
drive a segmentation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .similarity import WindowSimilarity

# Defaults, both from first principles rather than fit to any dataset.
DEFAULT_JUMP_RATE = 1e-3  # prior probability of switching reference per window
_CONF = 0.9  # posterior a state must reach to bound a breakpoint confidently


@dataclass
class Segment:
    state: str
    start_window: int  # inclusive window index range, for per-segment statistics
    end_window: int
    msa_start: int
    msa_end: int
    query_start: int
    query_end: int
    n_windows: int
    mean_posterior: float
    # Query-coordinate uncertainty of this segment's LEFT breakpoint (None for the
    # first segment); the window band around the switch where neither side is sure.
    breakpoint_lo: int | None
    breakpoint_hi: int | None


def _logsumexp(a: np.ndarray, axis: int) -> np.ndarray:
    m = np.max(a, axis=axis, keepdims=True)
    m = np.where(np.isneginf(m), 0.0, m)
    return np.squeeze(m, axis=axis) + np.log(np.sum(np.exp(a - m), axis=axis))


def _estimate_identity(result: WindowSimilarity) -> float:
    finite = [b for b in result.best_sim if not math.isnan(b)]
    p = float(np.median(finite)) if finite else 0.95
    return min(max(p, 0.80), 0.999)


def segment_query(
    result: WindowSimilarity,
    *,
    identity: float | None = None,
    jump_rate: float = DEFAULT_JUMP_RATE,
    states: list[str] | None = None,
) -> tuple[list[Segment], str | None]:
    """Return ``(segments, major_state)`` from an HMM decode of ``result``."""
    labels = states or list(result.similarities)
    n_win = len(result.positions)
    if not labels or n_win == 0:
        return [], None

    p = identity if identity is not None else _estimate_identity(result)
    emission = _emission_matrix(result, labels, p)  # shape (S, W)

    s_count = len(labels)
    if s_count == 1:
        path = np.zeros(n_win, dtype=int)
        post = np.ones((n_win, 1))
    else:
        trans = _transition_matrix(s_count, jump_rate)
        path = _viterbi(emission, trans)
        post = _posterior(emission, trans)

    segments = _path_to_segments(path, post, labels, result)
    major = _major_state(path, labels)
    return segments, major


def _emission_matrix(result: WindowSimilarity, labels: list[str], p: float) -> np.ndarray:
    """Per-state per-window copying log-likelihood, comparable across states.

    Each state contributes ``k*log(p) + (n-k)*log(1-p)`` over its comparable sites,
    then is padded to the window's maximum comparable count with a neutral
    ``log(0.5)`` per missing site. Without the padding a gapped reference (n=0)
    would emit 0 -- the maximum -- and the HMM would prefer absent references; the
    padding makes a gapped state emit low, so a real reference always wins where it
    has data, while genuine identity still drives the choice between present ones.
    """
    num = np.array([result.numerators[s] for s in labels], dtype=float)  # (S, W)
    den = np.array([result.denominators[s] for s in labels], dtype=float)
    n_max = den.max(axis=0)  # most comparable sites any reference offers per window
    pp = min(max(p, 1e-6), 1.0 - 1e-6)
    return (
        num * math.log(pp)
        + (den - num) * math.log(1.0 - pp)
        + (n_max[None, :] - den) * math.log(0.5)
    )


def _transition_matrix(s_count: int, jump_rate: float) -> np.ndarray:
    jump_rate = min(max(jump_rate, 1e-9), 0.5)
    log_stay = math.log(1.0 - jump_rate)
    log_jump = math.log(jump_rate / (s_count - 1))
    trans = np.full((s_count, s_count), log_jump)
    np.fill_diagonal(trans, log_stay)
    return trans


def _viterbi(emission: np.ndarray, trans: np.ndarray) -> np.ndarray:
    s_count, n_win = emission.shape
    prior = math.log(1.0 / s_count)
    v = prior + emission[:, 0]
    back = np.zeros((n_win, s_count), dtype=int)
    for i in range(1, n_win):
        scores = v[:, None] + trans  # (prev, cur)
        back[i] = np.argmax(scores, axis=0)
        v = scores.max(axis=0) + emission[:, i]
    path = np.zeros(n_win, dtype=int)
    path[-1] = int(np.argmax(v))
    for i in range(n_win - 1, 0, -1):
        path[i - 1] = back[i, path[i]]
    return path


def _posterior(emission: np.ndarray, trans: np.ndarray) -> np.ndarray:
    s_count, n_win = emission.shape
    prior = math.log(1.0 / s_count)
    log_alpha = np.empty((n_win, s_count))
    log_alpha[0] = prior + emission[:, 0]
    for i in range(1, n_win):
        log_alpha[i] = _logsumexp(log_alpha[i - 1][:, None] + trans, axis=0) + emission[:, i]
    log_beta = np.zeros((n_win, s_count))
    for i in range(n_win - 2, -1, -1):
        log_beta[i] = _logsumexp(
            trans + (emission[:, i + 1] + log_beta[i + 1])[None, :], axis=1
        )
    log_post = log_alpha + log_beta
    log_post -= _logsumexp(log_post, axis=1)[:, None]
    return np.exp(log_post)  # (W, S)


def _path_to_segments(
    path: np.ndarray, post: np.ndarray, labels: list[str], result: WindowSimilarity
) -> list[Segment]:
    half = (result.positions[1] - result.positions[0]) if len(result.positions) > 1 else 0
    boundaries: list[int] = [0]
    for i in range(1, len(path)):
        if path[i] != path[i - 1]:
            boundaries.append(i)
    boundaries.append(len(path))

    segments: list[Segment] = []
    for j in range(len(boundaries) - 1):
        lo_i, hi_i = boundaries[j], boundaries[j + 1] - 1
        state_idx = int(path[lo_i])
        msa_start = max(0, result.positions[lo_i] - half)
        msa_end = min(result.width, result.positions[hi_i] + half)
        mean_post = float(np.mean(post[lo_i : hi_i + 1, state_idx]))
        bp_lo, bp_hi = (None, None)
        if j > 0:
            bp_lo, bp_hi = _breakpoint_interval(
                post, path, lo_i, int(path[lo_i - 1]), state_idx, result
            )
        segments.append(Segment(
            state=labels[state_idx],
            start_window=lo_i, end_window=hi_i,
            msa_start=msa_start, msa_end=msa_end,
            query_start=result.column_to_query(msa_start),
            query_end=result.column_to_query(msa_end),
            n_windows=hi_i - lo_i + 1,
            mean_posterior=round(mean_post, 4),
            breakpoint_lo=bp_lo, breakpoint_hi=bp_hi,
        ))
    return segments


def _breakpoint_interval(
    post: np.ndarray, path: np.ndarray, switch_i: int,
    prev_state: int, new_state: int, result: WindowSimilarity,
) -> tuple[int, int]:
    """Query-coord band around a switch where neither flanking state is confident."""
    lo = switch_i
    while lo > 0 and post[lo - 1, prev_state] < _CONF:
        lo -= 1
    hi = switch_i
    while hi < len(path) - 1 and post[hi, new_state] < _CONF:
        hi += 1
    return (
        result.column_to_query(result.positions[lo]),
        result.column_to_query(result.positions[hi]),
    )


def _major_state(path: np.ndarray, labels: list[str]) -> str:
    counts = np.bincount(path, minlength=len(labels))
    return labels[int(np.argmax(counts))]
