"""Lightweight in-stage parallelism.

The parallel work here is subprocess-bound (progressiveMauve runs one process
per query), so a thread pool is the right tool: the GIL is released while each
external tool runs, and threads avoid the pickling constraints that
``multiprocessing`` imposes on closures.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import TypeVar

T = TypeVar("T")
R = TypeVar("R")


def parallel_map(
    fn: Callable[[T], R],
    items: Iterable[T],
    workers: int,
    *,
    logger: logging.Logger | None = None,
) -> list[R]:
    """Apply ``fn`` to each item, up to ``workers`` at a time, preserving order.

    Runs sequentially when ``workers <= 1`` or there is at most one item. The
    first worker exception propagates after in-flight tasks settle.
    """
    items_list: Sequence[T] = list(items)
    if workers <= 1 or len(items_list) <= 1:
        return [fn(item) for item in items_list]

    results: list[R | None] = [None] * len(items_list)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        future_to_index = {pool.submit(fn, item): i for i, item in enumerate(items_list)}
        for future in as_completed(future_to_index):
            index = future_to_index[future]
            results[index] = future.result()
    return results  # type: ignore[return-value]
