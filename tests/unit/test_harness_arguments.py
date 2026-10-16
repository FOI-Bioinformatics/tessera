"""The self-contained harnesses parse their arguments (audit 4, item F1).

``run_specificity.py`` scanned ``sys.argv`` for two flags by hand: ``--help`` started a
ten-replicate run, an unknown flag was ignored and a missing value raised an IndexError.
``run_regimes.py`` and ``characterise_false_regions.py`` did the same.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_VALIDATION = Path(__file__).resolve().parents[2] / "validation"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _VALIDATION / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


HARNESSES = ["run_specificity", "run_regimes", "characterise_false_regions"]


def _no_work(monkeypatch, module) -> None:
    """Fail the test if a harness gets past argument parsing into a simulation."""
    def boom(*args, **kwargs):
        raise AssertionError("the harness started work")

    for attr in ("simulate_clonal", "simulate_recombinant", "simulate_near_identical",
                 "_load_specificity", "_scan", "scan"):
        if hasattr(module, attr):
            monkeypatch.setattr(module, attr, boom)


@pytest.mark.parametrize("name", HARNESSES)
def test_help_exits_without_running(name: str, monkeypatch, capsys) -> None:
    module = _load(name)
    _no_work(monkeypatch, module)
    with pytest.raises(SystemExit) as exc:
        module.main(["--help"])
    assert exc.value.code == 0
    assert "--reps" in capsys.readouterr().out


@pytest.mark.parametrize("name", HARNESSES)
@pytest.mark.parametrize("argv", [["--rep", "3"], ["--reps"], ["--reps", "three"],
                                  ["--reps", "0"]])
def test_bad_arguments_are_an_error(name: str, argv: list[str], monkeypatch) -> None:
    module = _load(name)
    _no_work(monkeypatch, module)
    with pytest.raises(SystemExit) as exc:
        module.main(argv)
    assert exc.value.code == 2


def test_specificity_defaults_are_unchanged() -> None:
    """Pins the defaults the hand parser had: 10 replicates, gate 2. It fails before the
    change only because parse_args is new."""
    args = _load("run_specificity").parse_args([])
    assert (args.reps, args.min_methods) == (10, 2)
