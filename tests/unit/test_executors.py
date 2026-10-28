"""In-stage parallelism: ordering, failure handling, and the sequential path."""

from __future__ import annotations

import threading
import time

import pytest

from tessera.core.executors import parallel_map


def test_results_keep_the_input_order() -> None:
    # Completion order is deliberately the reverse of submission order.
    def work(i: int) -> int:
        time.sleep(0.05 * (5 - i))
        return i * i

    assert parallel_map(work, range(5), workers=5) == [0, 1, 4, 9, 16]


def test_runs_inline_without_a_pool_for_one_worker_or_one_item() -> None:
    threads: set[int] = set()

    def work(i: int) -> int:
        threads.add(threading.get_ident())
        return i

    assert parallel_map(work, [1, 2, 3], workers=1) == [1, 2, 3]
    assert parallel_map(work, [1], workers=8) == [1]
    assert threads == {threading.get_ident()}
    assert parallel_map(work, [], workers=8) == []


def test_a_failure_drops_the_work_not_yet_started() -> None:
    """The first exception must not be preceded by the whole rest of the stage.

    The pool's context manager shuts down waiting but not cancelling, so every queued
    genome was aligned before the first failure surfaced -- on a panel of hundreds,
    hours spent on a run already lost.
    """
    started = []
    lock = threading.Lock()

    def work(i: int) -> int:
        with lock:
            started.append(i)
        if i == 0:
            raise RuntimeError("first item fails")
        time.sleep(0.2)
        return i

    with pytest.raises(RuntimeError, match="first item fails"):
        parallel_map(work, range(40), workers=2)

    # Two workers, so at most the failing task and the one already in flight, plus
    # whatever a worker picked up before the cancellation took effect. Far short of 40.
    assert len(started) < 10, started


def test_an_exception_from_submission_still_shuts_the_pool_down() -> None:
    class Hostile:
        def __iter__(self):
            yield 1
            yield 2
            raise RuntimeError("iteration fails")

    with pytest.raises(RuntimeError, match="iteration fails"):
        parallel_map(lambda i: i, Hostile(), workers=4)
