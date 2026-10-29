#!/usr/bin/env python
"""Build real donor-absent runs for issue #90's Option A measurement.

``run_hybrids.py`` has exactly one donor-absent case, ``donorabsent_rsv``. One case cannot
calibrate a statistic, and the issue says so itself. ``measure_pendant.py`` can simulate
donor-absence cheaply by dropping the donor's rows from an existing alignment, but that
proxy is biased in a way the ``rsv_a`` / ``donorabsent_rsv`` pair exposes: dropping rows
leaves a *sparser* panel, which inflates the panel's own pendant lengths -- the statistic's
denominator -- and so understates the signal. A real run instead re-recruits a panel from
the lineages that remain and re-aligns it.

This script produces those real runs. For each eligible positive case it re-uses the case's
own dataset and pinned parent pair, sets ``case_type: panel_donor_absent``, and removes

* the donor's **parent sub-clade** where the clade names are hierarchical (``A.D.1.8`` ->
  remove ``A.D.1``), which is what ``donorabsent_rsv`` does -- a sister sub-clade is as good
  a source as the donor itself, so leaving one in means the donor is not really absent;
* the donor's **own clade** where the naming is flat (``DENV4``, ``B3``, ``denovo_3``), which
  for a flat system already removes the whole lineage.

It is a **probe, not a harness case**: nothing is added to ``run_hybrids.HYBRIDS``, the
must-pass headline is untouched, and the artifacts land in their own directory so
``measure_pendant.py --data`` can read them beside the main run.

Needs an aligner on ``PATH`` and (first run only) the Nextclade datasets:

    export PATH="$PATH:/path/to/aligner-env/bin"    # append; prepending swaps the Python
    python validation/probe_donor_absent.py
    python validation/measure_pendant.py --data validation/data/donor_absent_probe
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE / "data" / "donor_absent_probe"

# Case types whose parent pair comes from ``pick_parents``, so a donor-absent variant of
# them is the same splice against a panel missing the donor's lineage. The mosaic and
# negative types have no single donor lineage to remove.
ELIGIBLE = ("single_insert", "low_div")
# Dataset options that must be carried over for the probe to read the same tree.
CARRY = ("dataset", "clade_key", "aligner", "clades", "min_divergence", "cluster_ani",
         "pair_objective", "divergence_band")


def _load_harness():
    spec = importlib.util.spec_from_file_location("run_hybrids", HERE / "run_hybrids.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["run_hybrids"] = module
    spec.loader.exec_module(module)
    return module


def resolve_parents(case: dict, hx, logger: logging.Logger) -> tuple[str, str]:
    """The ``(backbone, donor)`` clades this case would splice, without splicing it.

    Repeats the dataset/tree/typing steps ``_prepare_case`` performs before it calls
    ``pick_parents``, so the probe can name the clade to remove *before* the run that
    consumes it. Everything here is served from the Nextclade cache after the first call.
    """
    dataset = hx.resolve_dataset(Path("/dev/null"), case["dataset"], email=None, logger=logger)
    genomes = hx.build_pool(dataset, cache_dir=hx.nextclade_cache(dataset.path, dataset.tag),
                            logger=logger)
    reference = hx.read_reference(hx._download_text(dataset, "reference", logger))
    tree = json.loads(hx._download_text(dataset, "treeJson", logger))["tree"]
    tips = hx.collect_tips(tree, case.get("clade_key"))
    if case.get("cluster_ani"):
        tips = hx._denovo_type_tips(genomes, tips, case["cluster_ani"], logger)
    clade_a, clade_b, _, _ = hx.pick_parents(
        tips, reference, case.get("clades", []), logger,
        objective=case.get("pair_objective", "max"),
        floor=case.get("min_divergence", hx.MIN_DIVERGENCE))
    return clade_a, clade_b


def build_probe(case: dict, hx, logger: logging.Logger) -> dict | None:
    """The donor-absent probe case for ``case``, or ``None`` when it cannot be posed."""
    clade_a, clade_b = resolve_parents(case, hx, logger)
    # Hierarchical names: remove the donor's parent level, as donorabsent_rsv does. Flat
    # names: the donor's own clade already *is* the whole lineage.
    absent = hx._parent_clade(clade_b) or hx.base_clade(clade_b)
    if hx.clade_match(clade_a, absent):
        logger.info("[%s] skip: removing %r would take the backbone %r with it",
                    case["name"], absent, clade_a)
        return None
    probe = {k: case[k] for k in CARRY if k in case}
    probe.update({"name": f"absent_{case['name']}", "case_type": "panel_donor_absent",
                  "clades": [clade_a, clade_b], "absent_clade": absent})
    return probe


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT, help="artifact directory")
    ap.add_argument("--only", nargs="*", default=None, help="limit to these case names")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logger = logging.getLogger("probe")
    hx = _load_harness()
    args.out.mkdir(parents=True, exist_ok=True)
    hx.DATA = args.out  # artifacts land here, not beside the must-pass run

    wanted = [c for c in hx.HYBRIDS
              if c.get("tier") is None
              and c.get("case_type", "single_insert") in ELIGIBLE
              and (args.only is None or c["name"] in args.only)]
    logger.info("%d eligible positive case(s)", len(wanted))

    done, skipped, failed = [], [], []
    for case in wanted:
        name = case["name"]
        try:
            probe = build_probe(case, hx, logger)
        except hx.CaseSkipped as exc:
            logger.info("[%s] SKIP while resolving parents: %s", name, exc)
            skipped.append((name, str(exc)))
            continue
        except Exception as exc:  # a dataset problem must not lose the rest of the probe
            logger.error("[%s] ERROR while resolving parents: %s: %s",
                         name, type(exc).__name__, exc)
            failed.append((name, f"{type(exc).__name__}: {exc}"))
            continue
        if probe is None:
            skipped.append((name, "removal would take the backbone"))
            continue
        t0 = time.monotonic()
        try:
            record = hx.run_case(probe, logger)
        except hx.CaseSkipped as exc:
            logger.info("[%s] SKIP: %s", probe["name"], exc)
            skipped.append((probe["name"], str(exc)))
            continue
        except Exception as exc:
            logger.error("[%s] ERROR: %s: %s", probe["name"], type(exc).__name__, exc)
            failed.append((probe["name"], f"{type(exc).__name__}: {exc}"))
            continue
        logger.info("[%s] done in %.0fs (removed %r)",
                    probe["name"], time.monotonic() - t0, probe["absent_clade"])
        done.append((probe["name"], probe["absent_clade"], record))

    print(f"\n=== donor-absent probe: {len(done)} built, {len(skipped)} skipped, "
          f"{len(failed)} failed ===")
    for name, absent, _ in done:
        print(f"  built    {name:<28} removed {absent}")
    for name, why in skipped:
        print(f"  skipped  {name:<28} {why}")
    for name, why in failed:
        print(f"  FAILED   {name:<28} {why}")
    print(f"\nartifacts in {args.out}")
    print(f"next: python validation/measure_pendant.py --data {args.out}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
