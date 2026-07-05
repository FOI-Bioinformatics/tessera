#!/usr/bin/env python
"""GENECONV gap probe: does an OpenRDP method detect what Tessera's ensemble misses?

Tessera already ships MaxChi and Bootscan, and its HMM/similarity core covers the RDP
similarity-switch signal. The one OpenRDP method with a genuinely different signal is GENECONV
(runs of consecutive matching polymorphic sites). On the easy real recombinants Tessera detects
everything, so any benefit can only appear on the *hard* synthetic cases (short tracts, low
divergence, mosaics) or in specificity (negatives).

This harness reuses the synthetic-hybrid builder (`run_hybrids.py`), which makes those hard cases
with known ground truth, and OpenRDP (`run_method_comparison.py`). Per case it runs Tessera's
default ensemble and OpenRDP's methods on the SAME `panel.msa.fasta`, then cross-tabs against the
ground truth so we can count where GENECONV succeeds and Tessera's ensemble fails -- the value
signal that would justify writing a GENECONV caller.

Caveat: OpenRDP runs all-vs-all on the panel, so a GENECONV hit on a Tessera miss shows the signal
*exists* in the alignment, not that a query-vs-panel GENECONV caller would fire. It is a gate, not a
proof. Needs the aligner env, network (Nextclade datasets), and the `openrdp` conda env; opt-in.

    export PATH="$PATH:$HOME/miniforge3/envs/recomfi-aln/bin"
    python validation/run_method_comparison_hybrids.py                 # the hard/negative cases
    python validation/run_method_comparison_hybrids.py --adversarial   # sub-window tract sweep
    python validation/run_method_comparison_hybrids.py lowdiv_rsv      # one case by name
    python validation/run_method_comparison_hybrids.py --no-openrdp    # Tessera only
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_hybrids import (  # noqa: E402 - sibling validation module
    HYBRIDS,
    CaseSkipped,
    _build_and_score,
    _prepare_case,
    parse_regions,
)
from run_method_comparison import (  # noqa: E402 - sibling validation module
    CALLERS,
    DEFAULT_OPENRDP_METHODS,
    OPENRDP_METHODS,
    openrdp_command,
    run_openrdp,
)

from tessera.recomb.regions import DEFAULT_METHODS, parse_methods  # noqa: E402

# Case types that stress detection (positives) or specificity (negatives); the easy single-insert
# core cases are excluded unless asked for with --all, since Tessera already passes them.
HARD_CASE_TYPES = frozenset({
    "neg_pure", "neg_within", "low_div", "panel_donor_absent", "panel_equidistant",
    "mosaic", "mask_sibling",
})
NEGATIVE_CASE_TYPES = frozenset({"neg_pure", "neg_within"})


def _fmt(flag: bool) -> str:
    return "yes" if flag else "no"


def tessera_callers(regions: list[dict], true_spans: list[tuple[int, int, str]],
                    is_negative: bool) -> set[str]:
    """Which Tessera callers fired on a relevant region.

    For a positive, "relevant" means a present region overlapping a true donor span; for a negative,
    any present (donor-present) region is a false call, so all of them count.
    """
    present = [r for r in regions if r.get("donor_absent") != "yes"]
    if is_negative:
        relevant = present
    else:
        relevant = [r for r in present
                    if any(int(r["query_start"]) <= hi and int(r["query_end"]) >= lo
                           for lo, hi, _ in true_spans)]
    fired: set[str] = set()
    for r in relevant:
        fired.update(m.strip() for m in (r.get("methods") or "").split(",") if m.strip())
    return fired & set(CALLERS)


def gap_summary(rows: list[dict], methods: tuple[str, ...] = DEFAULT_OPENRDP_METHODS) -> dict:
    """Count, per OpenRDP method, rescues (detects a positive Tessera's ensemble missed) and
    false calls (fires on a negative). ``rows`` are records with keys ``is_negative``,
    ``tessera_pass`` and ``openrdp`` (a per-method dict, or ``None`` when OpenRDP did not run)."""
    pos = [r for r in rows if not r["is_negative"] and r.get("openrdp") is not None]
    neg = [r for r in rows if r["is_negative"] and r.get("openrdp") is not None]
    per_method = {}
    for m in methods:
        per_method[m] = {
            "rescues": sum(1 for r in pos if r["openrdp"].get(m) and not r["tessera_pass"]),
            "false_calls": sum(1 for r in neg if r["openrdp"].get(m)),
        }
    return {
        "n_positive": len(pos),
        "n_negative": len(neg),
        "tessera_missed": sum(1 for r in pos if not r["tessera_pass"]),
        "tessera_false_calls": sum(1 for r in neg if not r["tessera_pass"]),
        "per_method": per_method,
    }


def _select(names: list[str], include_core: bool, adversarial: bool) -> list[dict]:
    """Cases to probe. Default: the hard must-pass case types. ``--adversarial`` selects the
    adversarial tier instead (sub-window tracts expected to break Tessera's windowed callers)."""
    if adversarial:
        want = [c for c in HYBRIDS if c.get("tier") == "adversarial"]
        return [c for c in want if not names or c["name"] in names]
    out = []
    for c in HYBRIDS:
        if c.get("tier", "must_pass") != "must_pass":
            continue  # frontier / adversarial probes are scored separately
        ct = c.get("case_type", "single_insert")
        if names:
            if c["name"] in names:
                out.append(c)
        elif ct in HARD_CASE_TYPES or (include_core and ct == "single_insert"):
            out.append(c)
    return out


def _probe_case(case: dict, orc, orc_methods: tuple[str, ...],
                logger: logging.Logger) -> dict:
    setup = _prepare_case(case, logger)
    result = _build_and_score(setup, "tip", parse_methods(",".join(DEFAULT_METHODS)),
                              setup.out, logger)
    is_neg = setup.case_type in NEGATIVE_CASE_TYPES
    spans = setup.true_spans or [(setup.q_start, setup.q_end, setup.clade_b)]
    regions = parse_regions(setup.out / "recombination_regions.tsv")
    det = None
    if orc is not None:
        work = setup.out / "openrdp"
        work.mkdir(parents=True, exist_ok=True)
        det = run_openrdp(setup.out / "panel.msa.fasta", orc, work,
                          setup.query_label, orc_methods, logger)
    return {
        "name": setup.name, "case_type": setup.case_type, "is_negative": is_neg,
        "tessera_pass": bool(result["pass"]), "tessera_detected": bool(result.get("detected")),
        "tessera_callers": tessera_callers(regions, spans, is_neg), "openrdp": det,
    }


def main(argv: list[str]) -> int:
    logging.basicConfig(level=logging.WARNING, format="[%(levelname)s] %(message)s")
    logger = logging.getLogger("gap-probe")

    orc = None if "--no-openrdp" in argv else openrdp_command()
    orc_methods = OPENRDP_METHODS if "--with-siscan" in argv else DEFAULT_OPENRDP_METHODS
    names = [a for a in argv if not a.startswith("-")]
    adversarial = "--adversarial" in argv
    cases = _select(names, include_core="--all" in argv, adversarial=adversarial)

    kind = "ADVERSARIAL" if adversarial else "hard"
    print(f"\nGENECONV gap probe ({kind}) -- {len(cases)} case(s); "
          f"OpenRDP {'on' if orc else 'off'}\n" + "=" * 72)
    rows: list[dict] = []
    for case in cases:
        try:
            row = _probe_case(case, orc, orc_methods, logger)
        except CaseSkipped as exc:
            print(f"  [SKIP ] {case['name']:16} {exc}")
            continue
        except Exception as exc:  # noqa: BLE001 - report, never abort the batch
            logger.exception("[%s] ERROR", case["name"])
            print(f"  [ERROR] {case['name']:16} {type(exc).__name__}: {exc}")
            continue
        rows.append(row)
        _print_row(row, orc_methods if orc else ())

    if orc:
        _print_summary(gap_summary(rows, orc_methods))
    else:
        print("\n[note] OpenRDP off; ran Tessera only (no gap comparison).")
    return 0


def _print_row(row: dict, methods: tuple[str, ...]) -> None:
    tag = "neg" if row["is_negative"] else "pos"
    verdict = "PASS" if row["tessera_pass"] else "FAIL"
    callers = ",".join(c for c in CALLERS if c in row["tessera_callers"]) or "-"
    line = (f"  {row['name']:16} {row['case_type']:14} [{tag}] Tessera {verdict:4} "
            f"({callers})")
    if methods:
        det = row["openrdp"]
        if det is None:
            line += "  | OpenRDP did not run"
        else:
            hits = ",".join(m for m in methods if det.get(m)) or "-"
            line += f"  | OpenRDP: {hits}"
    print(line)


def _print_summary(summ: dict) -> None:
    print("\n" + "=" * 72)
    print(f"positives {summ['n_positive']} (Tessera missed {summ['tessera_missed']}), "
          f"negatives {summ['n_negative']} (Tessera false-called {summ['tessera_false_calls']})")
    print("  per OpenRDP method -- rescues (detect a positive Tessera missed) / "
          "false-calls on negatives:")
    for m, s in summ["per_method"].items():
        mark = "  <-- fills a gap" if s["rescues"] else ""
        print(f"    {m:10} rescues {s['rescues']}  false-calls {s['false_calls']}{mark}")
    gc = summ["per_method"].get("geneconv", {})
    print(f"\n  GENECONV verdict: {gc.get('rescues', 0)} rescue(s), "
          f"{gc.get('false_calls', 0)} false-call(s) on negatives.")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
