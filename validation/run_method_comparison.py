#!/usr/bin/env python
"""Method comparison on the real published recombinants: Tessera vs OpenRDP.

Two comparisons on the same published-recombinant alignments (the `datasets.json` positives,
fetched by `fetch.py`):

1. **Tessera's own callers** -- HMM, 3SEQ, MaxChi, Bootscan -- each run alone, plus the default
   four-caller ensemble, tabulating per-caller detection of the known recombination.
2. **OpenRDP** (the maintained RDP5 reimplementation, PoonLab/OpenRDP) -- its seven methods
   (GENECONV, Bootscan, MaxChi, SiScan, Chimaera, 3SEQ, RDP) run on the *same* alignment, so
   Tessera sits beside the standard methods on identical data.

OpenRDP is pure Python but pins ``numpy<2``/``h5py<3.11`` and is easiest to keep in its own
environment. The harness resolves it in this order: ``$OPENRDP_CMD`` (a shell command), then an
``openrdp`` on PATH, then ``conda run -n $OPENRDP_CONDA_ENV openrdp`` (default env ``openrdp``). If
none resolves, the OpenRDP half is skipped with a note -- Tessera's own comparison still runs.

    # once, in a dedicated env:
    #   conda create -n openrdp python=3.10 -y
    #   conda run -n openrdp pip install "numpy<2.0.0" "scipy>=1.5.0" "h5py<3.11.0" \
    #       git+https://github.com/PoonLab/OpenRDP
    export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
    python validation/fetch.py && python validation/run_method_comparison.py
    python validation/run_method_comparison.py --no-openrdp   # Tessera callers only
"""

from __future__ import annotations

import csv as _csv
import io
import logging
import os
import shlex
import shutil
import subprocess
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
# OpenRDP method names as they appear in its CSV `Method` column (lower-cased).
OPENRDP_METHODS = ("geneconv", "bootscan", "maxchi", "siscan", "chimaera", "threeseq", "rdp")
# SiScan is permutation-heavy: it timed out (>15 min) on the full-genome HIV alignment while
# detecting nothing on the datasets where it did complete, so it is opt-in (--with-siscan).
DEFAULT_OPENRDP_METHODS = ("geneconv", "bootscan", "maxchi", "chimaera", "threeseq", "rdp")
_ALLOWED_BASES = set("ATGCN-")


def any_region_detected(regions: list[dict]) -> bool:
    """True iff any present (non ``donor_absent``) recombination region was called."""
    return any(r.get("donor_absent") != "yes" for r in regions)


def sanitize_alignment(text: str) -> str:
    """Map any non-``ATGCN-`` base to ``N`` (OpenRDP rejects IUPAC ambiguity codes)."""
    out: list[str] = []
    for ln in text.splitlines():
        if ln.startswith(">"):
            out.append(ln)
        else:
            out.append("".join(c if c.upper() in _ALLOWED_BASES else "N" for c in ln.upper()))
    return "\n".join(out) + "\n"


def parse_openrdp_csv(text: str) -> list[dict]:
    """Parse OpenRDP's ``Method,Start,End,Recombinant,Parent1,Parent2,Pvalue`` CSV."""
    rows: list[dict] = []
    reader = _csv.DictReader(io.StringIO(text))
    for row in reader:
        if row.get("Method"):
            rows.append({k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()})
    return rows


def openrdp_detections(rows: list[dict], query_label: str,
                       methods: tuple[str, ...] = OPENRDP_METHODS) -> dict[str, bool]:
    """Per-method: did OpenRDP surface the known recombinant?

    A method ``detects`` when the query label appears anywhere in a triplet it reports
    (Recombinant / Parent1 / Parent2). Using the whole triplet is robust to the parent/child
    direction ambiguity that triplet methods -- and Tessera's single backbone -- share. Only the
    requested ``methods`` are reported, so a method that was not run reads as absent, not ``no``.
    """
    hits = {m: False for m in methods}
    for row in rows:
        method = (row.get("Method") or "").strip().lower()
        if method not in hits:
            continue
        triplet = {row.get("Recombinant"), row.get("Parent1"), row.get("Parent2")}
        if query_label in triplet:
            hits[method] = True
    return hits


def openrdp_command() -> list[str] | None:
    """Resolve how to invoke OpenRDP, or ``None`` if it cannot be found."""
    if os.environ.get("OPENRDP_CMD"):
        return shlex.split(os.environ["OPENRDP_CMD"])
    if shutil.which("openrdp"):
        return ["openrdp"]
    if shutil.which("conda"):
        return ["conda", "run", "-n", os.environ.get("OPENRDP_CONDA_ENV", "openrdp"), "openrdp"]
    return None


def run_openrdp(msa: Path, cmd: list[str], work: Path, query_label: str,
                methods: tuple[str, ...], logger: logging.Logger) -> dict[str, bool] | None:
    """Run OpenRDP on a sanitized copy of ``msa``; return per-method detections or ``None``."""
    clean = work / "clean.fasta"
    clean.write_text(sanitize_alignment(msa.read_text()))
    out_csv = work / "openrdp.csv"
    try:
        proc = subprocess.run([*cmd, str(clean), "-o", str(out_csv), "-m", *methods],
                              capture_output=True, text=True, timeout=1800)
    except (subprocess.TimeoutExpired, FileNotFoundError) as exc:
        logger.warning("OpenRDP run failed: %s", exc)
        return None
    if not out_csv.exists():
        logger.warning("OpenRDP produced no output: %s", proc.stderr.strip()[:200])
        return None
    return openrdp_detections(parse_openrdp_csv(out_csv.read_text()), query_label, methods)


def _run_caller(msa: Path, label: str, methods: tuple[str, ...], out: Path,
                rec: dict, logger: logging.Logger) -> bool:
    run_recomb(
        RecombParams(msa=msa, output=out, query=label, methods=methods,
                     window_size=rec.get("window_size", 1000),
                     window_step=rec.get("window_step", 100)),
        logger,
    )
    return any_region_detected(_parse_regions(out / "recombination_regions.tsv"))


def _fmt(flag: bool) -> str:
    return "yes" if flag else "no"


def main(argv: list[str]) -> int:
    logging.basicConfig(level=logging.ERROR, format="%(message)s")
    logger = logging.getLogger("method-comparison")

    orc = None if "--no-openrdp" in argv else openrdp_command()
    orc_methods = OPENRDP_METHODS if "--with-siscan" in argv else DEFAULT_OPENRDP_METHODS

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
            per: dict[str, object] = {}
            for m in CALLERS:
                per[m] = _run_caller(msa, label, parse_methods(m), work / m, rec, logger)
            per["ensemble"] = _run_caller(msa, label, parse_methods("hmm,3seq,maxchi,bootscan"),
                                          work / "ensemble", rec, logger)
            if orc is not None:
                per["openrdp"] = run_openrdp(msa, orc, work, label, orc_methods, logger)
            rows.append((ds["name"], per))

    _print_combined_table(rows, orc, orc_methods)
    if orc is None:
        print("\n[note] OpenRDP not run (--no-openrdp, or no openrdp/conda env resolved). "
              "Set $OPENRDP_CMD or install into a conda env named 'openrdp'.")
    return 0


def _print_combined_table(rows: list[tuple[str, dict]], orc: list[str] | None,
                          methods: tuple[str, ...]) -> None:
    """One table: Tessera's callers and, if run, OpenRDP's methods on the same alignments."""
    tcols = (*CALLERS, "ensemble")
    ocols = tuple(methods) if orc is not None else ()
    sep = "  |  "

    def cells(values: list[str]) -> str:
        return " ".join(f"{v:>8}" for v in values)

    tess_hdr, op_hdr = cells(list(tcols)), cells(list(ocols))
    prefix = f"  {'dataset':20} "
    group = f"  {'':20} {'Tessera'.center(len(tess_hdr))}"
    col = prefix + tess_hdr
    if ocols:
        group += f"{sep}{'OpenRDP (RDP5)'.center(len(op_hdr))}"
        col += sep + op_hdr

    print("\nMethod comparison on real published recombinants (Tessera vs OpenRDP)\n"
          + "=" * len(col))
    print(group)
    print(col)
    for name, per in rows:
        if "skip" in per:
            print(f"  {name:20} SKIP: {per['skip']}")
            continue
        line = f"  {name:20} " + cells([_fmt(bool(per[c])) for c in tcols])
        if ocols:
            det = per.get("openrdp")
            ovals = ["-" if det is None else _fmt(det[m]) for m in ocols]
            line += sep + cells(ovals)
        print(line)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
