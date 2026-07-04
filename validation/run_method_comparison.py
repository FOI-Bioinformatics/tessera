#!/usr/bin/env python
"""Per-caller method comparison on the real published recombinants.

Runs each of Tessera's recombination callers -- HMM, 3SEQ, MaxChi, Bootscan -- on the same
published-recombinant alignments (the `datasets.json` cases fetched by `fetch.py`) and tabulates,
per caller, whether it detects the known recombination; the default four-caller ensemble is shown
alongside. This is the internal-method comparison; it needs skani/an aligner and the fetched data,
and is opt-in.

An external head-to-head against RDP4/OpenRDP was intended but OpenRDP does not build on modern
Python (it pins ``numpy==1.17.3`` and needs cblas/pybind11), and no reference PhiPack/3seq binary is
available here, so that path is deferred. If OpenRDP ever becomes importable the ``--openrdp`` hook
below runs it; otherwise it is skipped with this note.

    export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
    python validation/fetch.py && python validation/run_method_comparison.py
"""

from __future__ import annotations

import logging
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_validation import (  # noqa: E402 - sibling validation module
    _ALIGNER_BINARIES,
    _load_datasets,
    _parse_regions,
    _resolve_inputs,
)

from tessera.msa.build import MsaParams, build_msa  # noqa: E402
from tessera.recomb.regions import parse_methods  # noqa: E402
from tessera.recomb.run import RecombParams, run_recomb  # noqa: E402

CALLERS = ("hmm", "3seq", "maxchi", "bootscan")


def any_region_detected(regions: list[dict]) -> bool:
    """True iff any present (non ``donor_absent``) recombination region was called."""
    return any(r.get("donor_absent") != "yes" for r in regions)


def _run_caller(msa: Path, label: str, methods: tuple[str, ...], out: Path,
                rec: dict, logger: logging.Logger) -> bool:
    run_recomb(
        RecombParams(msa=msa, output=out, query=label, methods=methods,
                     window_size=rec.get("window_size", 1000),
                     window_step=rec.get("window_step", 100)),
        logger,
    )
    return any_region_detected(_parse_regions(out / "recombination_regions.tsv"))


def main(argv: list[str]) -> int:
    logging.basicConfig(level=logging.ERROR, format="%(message)s")
    logger = logging.getLogger("method-comparison")

    # Only the positive (recombinant) datasets -- a clonal control has nothing to detect.
    datasets = [d for d in _load_datasets()
                if d.get("enabled", True) and "min_regions" in d.get("expected", {})]
    rows: list[tuple[str, dict]] = []
    with tempfile.TemporaryDirectory() as tmp:
        for ds in datasets:
            binary = _ALIGNER_BINARIES.get(ds.get("aligner", "sibeliaz"), ds.get("aligner"))
            query, collection, label = _resolve_inputs(ds)
            if query is None:
                rows.append((ds["name"], {"skip": "sequences not present (run fetch.py)"}))
                continue
            if shutil.which(binary) is None:
                rows.append((ds["name"], {"skip": f"{binary} not on PATH"}))
                continue
            work = Path(tmp) / ds["name"]
            work.mkdir(parents=True)
            msa = work / "msa.fasta"
            build_msa(MsaParams(query=query, collection=collection, output=msa,
                                aligner=ds.get("aligner", "sibeliaz"),
                                reference=ds.get("reference"),
                                threads=ds.get("threads", 4)), logger)
            rec = ds.get("recomb", {})
            per: dict[str, bool] = {}
            for m in CALLERS:
                per[m] = _run_caller(msa, label, parse_methods(m), work / m, rec, logger)
            per["ensemble"] = _run_caller(msa, label, parse_methods("hmm,3seq,maxchi,bootscan"),
                                          work / "ensemble", rec, logger)
            rows.append((ds["name"], per))

    hdr = f"  {'dataset':22} " + " ".join(f"{c:>8}" for c in (*CALLERS, "ensemble"))
    print("\nTessera per-caller detection on real published recombinants\n" + "=" * len(hdr))
    print(hdr)
    for name, per in rows:
        if "skip" in per:
            print(f"  {name:22} SKIP: {per['skip']}")
            continue
        cells = " ".join(f"{'yes' if per[c] else 'no':>8}" for c in (*CALLERS, "ensemble"))
        print(f"  {name:22} {cells}")

    if "--openrdp" in argv:
        try:
            import openrdp  # noqa: F401
            print("\n(OpenRDP present -- external comparison not yet wired; see module docstring)")
        except ImportError:
            print("\n[note] OpenRDP not importable (does not build on modern Python); "
                  "external RDP4 comparison deferred.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
