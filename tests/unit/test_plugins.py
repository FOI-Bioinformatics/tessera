"""Registry discovery and error behaviour."""

from __future__ import annotations

import pytest

from tessera.aligners.base import Aligner
from tessera.aligners.base import registry as aligner_registry
from tessera.core.errors import PluginError

ALL_ALIGNERS = ("progressivemauve", "mafft", "minimap2", "sibeliaz", "cactus")


def test_registry_lists_all_aligners() -> None:
    assert set(ALL_ALIGNERS) <= set(aligner_registry.names())


def test_registry_create_returns_aligner_instances() -> None:
    for name in ALL_ALIGNERS:
        aligner = aligner_registry.create(name)
        assert isinstance(aligner, Aligner)
        assert aligner.capabilities.name == name


def test_unknown_aligner_raises_plugin_error() -> None:
    with pytest.raises(PluginError):
        aligner_registry.get("does-not-exist")
