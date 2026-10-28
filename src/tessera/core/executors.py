"""Lightweight in-stage parallelism.

The parallel work here is subprocess-bound (progressiveMauve runs one process
per query), so a thread pool is the right tool: the GIL is released while each
external tool runs, and threads avoid the pickling constraints that
``multiprocessing`` imposes on closures.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import Future, ThreadPoolExecutor, as_completed


def parallel_map[T, R](
    fn: Callable[[T], R],
    items: Iterable[T],
    workers: int,
    *,
    logger: logging.Logger | None = None,
) -> list[R]:
    """Apply ``fn`` to each item, up to ``workers`` at a time, preserving order.

    Runs sequentially when ``workers <= 1`` or there is at most one item. The first
    worker exception propagates once the tasks already running have settled, and the
    tasks not yet started are dropped rather than run first. Leaving the pool's own
    context manager to shut down did the opposite: it waits without cancelling, so a
    panel whose first genome failed to align went on aligning every remaining genome
    before reporting it -- the whole cost of the stage paid for a run already lost.
    """
    items_list: Sequence[T] = list(items)
    if workers <= 1 or len(items_list) <= 1:
        return [fn(item) for item in items_list]

    results: list[R | None] = [None] * len(items_list)
    pool = ThreadPoolExecutor(max_workers=workers)
    future_to_index: dict[Future[R], int] = {}
    try:
        future_to_index = {pool.submit(fn, item): i for i, item in enumerate(items_list)}
        for future in as_completed(future_to_index):
            index = future_to_index[future]
            results[index] = future.result()
    except BaseException:
        cancelled = sum(1 for f in future_to_index if f.cancel())
        if cancelled and logger is not None:
            logger.debug("Dropped %d queued task(s) after a failure.", cancelled)
        pool.shutdown(wait=True, cancel_futures=True)
        raise
    pool.shutdown(wait=True)
    return results  # type: ignore[return-value]
