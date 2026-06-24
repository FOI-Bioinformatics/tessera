"""Registry discovery and error behaviour."""

from __future__ import annotations

import pytest

from recomfi.aligners.base import Aligner
from recomfi.aligners.base import registry as aligner_registry
from recomfi.core.errors import PluginError


def test_registry_lists_the_three_aligners() -> None:
    assert {"progressivemauve", "sibeliaz", "cactus"} <= set(aligner_registry.names())


def test_registry_create_returns_aligner_instances() -> None:
    for name in ("progressivemauve", "sibeliaz", "cactus"):
        aligner = aligner_registry.create(name)
        assert isinstance(aligner, Aligner)
        assert aligner.capabilities.name == name


def test_unknown_aligner_raises_plugin_error() -> None:
    with pytest.raises(PluginError):
        aligner_registry.get("does-not-exist")
