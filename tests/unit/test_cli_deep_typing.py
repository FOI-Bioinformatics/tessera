"""The deep-typing toggle reaches FillParams from both CLI commands."""

from __future__ import annotations

from tessera.discover.iterate import FillParams


def test_default_deep_typing_off():
    assert FillParams(query=None, collection=None, output=None).deep_typing is False


def test_for_detection_forwards_deep_typing():
    params = FillParams.for_detection(query=None, output=None, deep_typing=True)
    assert params.deep_typing is True
