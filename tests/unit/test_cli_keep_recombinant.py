"""The keep-recombinant toggle reaches FillParams from both CLI commands."""

from __future__ import annotations

from tessera.discover.iterate import FillParams


def test_default_excludes_recombinant():
    params = FillParams(query=None, collection=None, output=None)
    assert params.keep_recombinant is False


def test_for_detection_forwards_keep_recombinant():
    params = FillParams.for_detection(
        query=None, output=None, keep_recombinant=True,
    )
    assert params.keep_recombinant is True
