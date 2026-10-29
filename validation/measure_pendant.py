#!/usr/bin/env python
"""Measure the Option A pendant-branch statistic proposed for issue #90.

Issue #90: with the donor's whole sub-clade removed from an RSV-A panel, the scan names a
donor from a *neighbouring* sub-clade with four-caller agreement and no caveat. Four
candidate discriminators have already been measured on the issue and eliminated, all of
them functions of pairwise identity. Option A in
``docs/superpowers/specs/2026-10-29-issue-90-donor-absent-design.md`` proposes evidence of a
different kind: not how well the best reference matches, but *how the query's tract attaches
to the panel's structure*. A tract drawn from a sampled lineage should attach by a short
pendant branch; a tract from an unsampled lineage by a long one, because the sites
separating its true parent from every panel member sit on that branch.

This script measures that statistic. It does **not** implement it in the package: the design
note is explicit that the measurement comes first and decides whether anything is built.

For one region of one case the statistic is

1. restrict the alignment to the region's columns,
2. build a p-distance matrix over the panel plus the query,
3. neighbour-join it (Saitou-Nei, a few lines of numpy -- no tree library, no binary),
4. read the pendant branch by which each tip attaches, and
5. express the query tract's pendant length against the panel's own pendant lengths over the
   same columns, so the panel calibrates itself and the statistic is comparable across
   pathogens and divergences.

The harness itself has exactly one real donor-absent case, which is too few to calibrate
anything on -- the issue says so itself. ``probe_donor_absent.py`` builds the rest for real
(panel re-recruited from the surviving lineages, re-aligned), and they are the headline test
here. Each region is additionally measured under two *simulated* removals, which are cheap
but biased, and are reported separately for that reason:

``full``        the panel as run: the donor's lineage present for a positive case, really
                absent for a ``donorabsent_*``/``absent_*`` run;
``loo_exact``   the donor's own lineage's rows dropped, its sibling sub-clades retained;
``loo_parent``  the donor's whole parent sub-clade's rows dropped, which is the lineage set
                the real ``donorabsent_rsv`` case removes (``run_hybrids._parent_clade`` +
                ``clade_match``, imported here so the removal is identical).

Dropping rows from a reference-anchored alignment does not move its columns, so the simulated
panels need no realignment -- but they are biased, and measurably so. A real run re-recruits a
panel from the lineages that remain; dropping rows just leaves a *sparser* panel, whose own
pendant lengths are longer, and those are the statistic's denominator. On ``rsv_a``, whose
donor-absent counterpart exists for real, the simulation reads 1.02 where the real run reads
3.22. The script prints that contrast first, so the proxy is never read as the measurement.

The falsification criterion is the design note's: if the donor-absent values fall inside the
donor-present distribution, as the four identity discriminators did, Option A is dead and
joins them.

Needs **no aligner, no network and no downloaded data**: it reads the per-case artifacts that
``run_hybrids.py`` leaves in ``validation/data/hybrids/<case>/`` (``panel.msa.fasta``,
``lineages.tsv``, ``recombination_regions.tsv``). Run ``run_hybrids.py`` first if that
directory is empty.

    python validation/probe_donor_absent.py          # the real donor-absent runs, once
    python validation/measure_pendant.py --out validation/data/pendant.tsv
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from tessera.recomb.similarity import _canonical_mask, _read_alignment

HERE = Path(__file__).resolve().parent
HYBRIDS = HERE / "data" / "hybrids"
PROBE = HERE / "data" / "donor_absent_probe"  # written by probe_donor_absent.py
QUERY_LABEL = "hybrid"
MIN_PANEL = 4  # fewer tips than this and a pendant length is not worth reading

COLUMNS = (
    "case", "klass", "region", "variant", "n_panel", "n_dropped", "dropped_clades",
    "columns", "comparable", "true_backbone", "true_donor",
    "called_donor", "called_donor_clade", "attribution", "donor_lineage_present",
    "methods", "n_methods", "length_bp", "mean_sim_minor",
    "pendant_query", "pendant_panel_median", "pendant_panel_max",
    "panel_med_pairwise", "panel_spread_share",
    "quantile", "ratio_max", "ratio_median", "robust_z",
)
STATS = ("quantile", "ratio_max", "ratio_median", "robust_z")
# Discriminator 1 from the issue -- an absolute threshold on the named donor's identity over
# the tract -- was measured there and eliminated. It is re-run here on exactly the same case
# set as a control: a pendant statistic is only worth anything if it separates where this
# does not. Negated, because donor-absence is expected to show as *lower* identity.
CONTROL = "neg_mean_sim_minor"


def _load_harness():
    """Import ``run_hybrids`` for its clade helpers, so the simulated removal of a donor
    lineage is the same operation the real ``panel_donor_absent`` case performs."""
    spec = importlib.util.spec_from_file_location("run_hybrids", HERE / "run_hybrids.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["run_hybrids"] = module
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------- neighbour joining


def nj_pendant_lengths(d: np.ndarray) -> np.ndarray:
    """Neighbour-join ``d`` and return each original tip's pendant branch length.

    The pendant length is the branch joining a tip to its immediate parent node, fixed at
    the moment the tip is first joined (Saitou-Nei), or by the three-point formula for
    whichever tips survive into the final triplet. Negative lengths, which NJ can produce
    on a non-additive matrix, are returned as-is; the caller clips them.
    """
    n = int(d.shape[0])
    if n < 2:
        raise ValueError("neighbour joining needs at least two tips")
    big = np.full((2 * n, 2 * n), np.nan, dtype=float)
    big[:n, :n] = d
    active = list(range(n))
    pendant = np.full(n, np.nan, dtype=float)
    nxt = n

    while len(active) > 3:
        m = len(active)
        idx = np.asarray(active)
        sub = big[np.ix_(idx, idx)]
        r = sub.sum(axis=1)
        # Saitou-Nei Q: the pair minimising it is the one to join.
        q = (m - 2) * sub - r[:, None] - r[None, :]
        np.fill_diagonal(q, np.inf)
        a, b = np.unravel_index(int(np.argmin(q)), q.shape)
        i, j = active[a], active[b]
        dij = big[i, j]
        li = 0.5 * dij + (r[a] - r[b]) / (2 * (m - 2))
        if i < n:
            pendant[i] = li
        if j < n:
            pendant[j] = dij - li
        u = nxt
        nxt += 1
        for k in active:
            if k in (i, j):
                continue
            big[u, k] = big[k, u] = 0.5 * (big[i, k] + big[j, k] - dij)
        big[u, u] = 0.0
        active = [k for k in active if k not in (i, j)] + [u]

    if len(active) == 3:
        i, j, k = active
        three = {
            i: 0.5 * (big[i, j] + big[i, k] - big[j, k]),
            j: 0.5 * (big[i, j] + big[j, k] - big[i, k]),
            k: 0.5 * (big[i, k] + big[j, k] - big[i, j]),
        }
        for tip, length in three.items():
            if tip < n:
                pendant[tip] = length
    else:  # exactly two nodes left: split the remaining distance evenly
        i, j = active
        for tip in (i, j):
            if tip < n:
                pendant[tip] = 0.5 * big[i, j]
    return pendant


def pdist_matrix(
    rows: dict[str, np.ndarray], labels: list[str], start: int, end: int
) -> tuple[np.ndarray, int]:
    """Pairwise p-distance over the canonical columns of ``[start, end)``.

    Returns the matrix and the smallest number of columns any pair had in common, so a span
    too gappy to carry a tree can be reported rather than silently trusted.
    """
    stack = np.stack([rows[label][start:end] for label in labels])
    canon = _canonical_mask(stack)
    n = len(labels)
    d = np.zeros((n, n), dtype=float)
    least = end - start
    for i in range(n):
        for j in range(i + 1, n):
            both = canon[i] & canon[j]
            comparable = int(np.count_nonzero(both))
            least = min(least, comparable)
            if comparable == 0:
                d[i, j] = d[j, i] = np.nan
            else:
                diff = int(np.count_nonzero(stack[i][both] != stack[j][both]))
                d[i, j] = d[j, i] = diff / comparable
    return d, least


def pendant_stats(
    rows: dict[str, np.ndarray], panel: list[str], start: int, end: int
) -> dict[str, float] | str:
    """The statistic for one tract against one panel, or a string saying why not.

    Option A calibrates the query's pendant against *the panel's own* pendant lengths, so a
    panel with only two or three members has no distribution to calibrate against and the
    statistic is simply undefined -- which is a finding about where it can be used, not an
    inconvenience to be worked around.
    """
    if len(panel) < MIN_PANEL:
        return f"panel has {len(panel)} member(s); the statistic needs >= {MIN_PANEL}"
    labels = [*panel, QUERY_LABEL]
    d, comparable = pdist_matrix(rows, labels, start, end)
    if not np.isfinite(d).all():
        return "a pair of records shares no comparable column over the span"
    # NJ can return small negative branch lengths on a non-additive matrix; clip to zero,
    # the usual convention, which cannot manufacture a long pendant.
    pendant = np.clip(nj_pendant_lengths(d), 0.0, None)
    lq = float(pendant[-1])
    pan = pendant[:-1]
    median, top = float(np.median(pan)), float(np.max(pan))
    mad = float(np.median(np.abs(pan - median)))
    # How tightly the panel is packed: the panel's own pendant scale against its overall
    # spread. A dense panel (many close relatives) has short pendants, so a SMALL share --
    # and that is the panel on which an unsampled donor's long pendant can stand out.
    off = d[: len(panel), : len(panel)][np.triu_indices(len(panel), k=1)]
    pairwise = float(np.median(off)) if off.size else float("nan")
    return {
        "comparable": float(comparable),
        "panel_med_pairwise": pairwise,
        "panel_spread_share": median / pairwise if pairwise > 0 else float("nan"),
        "pendant_query": lq,
        "pendant_panel_median": median,
        "pendant_panel_max": top,
        "quantile": float(np.count_nonzero(pan < lq) / pan.size),
        "ratio_max": lq / top if top > 0 else float("nan"),
        "ratio_median": lq / median if median > 0 else float("nan"),
        "robust_z": (lq - median) / (1.4826 * mad) if mad > 0 else float("nan"),
    }


# --------------------------------------------------------------------------- case artifacts


@dataclass
class Case:
    name: str
    klass: str  # positive | negative | donor_absent
    msa: Path
    lineages: dict[str, str]  # panel label -> clade
    true_backbone: str
    true_donor: str
    regions: list[dict[str, str]]


def _classify(name: str) -> str:
    # "absent_*" is probe_donor_absent.py's naming for a real donor-absent run; the
    # harness's own single case predates it and is called donorabsent_rsv.
    if name == "donorabsent_rsv" or name.startswith("absent_"):
        return "donor_absent"
    if name.startswith(("neg_", "negwithin_")):
        return "negative"
    return "positive"


def load_case(case_dir: Path) -> Case | None:
    """Read one case's artifacts, or ``None`` when the case did not run (a SKIP)."""
    msa, lin, reg = (case_dir / "panel.msa.fasta", case_dir / "lineages.tsv",
                     case_dir / "recombination_regions.tsv")
    if not (msa.exists() and lin.exists() and reg.exists()):
        return None

    lineages: dict[str, str] = {}
    query_genotype = ""
    with open(lin, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["source"] == "query":
                query_genotype = row["genotype"]
            else:
                lineages[row["label"]] = row["genotype"]

    # The harness writes the query's row as "hybrid(<backbone>+<donor>)"; a negative case
    # carries an empty donor.
    inner = query_genotype[query_genotype.find("(") + 1 : query_genotype.rfind(")")]
    backbone, _, donor = inner.partition("+")

    with open(reg, encoding="utf-8", newline="") as handle:
        regions = list(csv.DictReader(handle, delimiter="\t"))

    return Case(name=case_dir.name, klass=_classify(case_dir.name), msa=msa,
                lineages=lineages, true_backbone=backbone.strip(),
                true_donor=donor.strip(), regions=regions)


# --------------------------------------------------------------------------- measurement


def panel_variants(case: Case, panel: list[str], hx) -> list[tuple[str, list[str], list[str]]]:
    """``(variant, kept labels, dropped clades)`` for the panels this case supports.

    ``loo_*`` is offered only where it is meaningful and safe: the case must be a positive
    with a donor clade distinct from its backbone, and the removal must not take the
    backbone's representatives with it.
    """
    out = [("full", panel, [])]
    if case.klass != "positive" or not case.true_donor:
        return out
    if hx.clade_match(case.true_donor, case.true_backbone):
        return out  # a within-clade splice has no separate donor lineage to remove

    for variant, absent in (("loo_exact", case.true_donor),
                            ("loo_parent", hx._parent_clade(case.true_donor))):
        if not absent or hx.clade_match(case.true_backbone, absent):
            continue  # no parent level to remove, or removing it would take the backbone
        dropped = sorted({c for label, c in case.lineages.items()
                          if label in panel and hx.clade_match(c, absent)})
        kept = [label for label in panel
                if not hx.clade_match(case.lineages.get(label, "?"), absent)]
        if not dropped or len(kept) < MIN_PANEL:
            continue
        out.append((variant, kept, dropped))
    return out


def measure_case(
    case: Case, hx, undefined: list[tuple[str, int, str, int, str]]
) -> list[dict[str, object]]:
    rows = _read_alignment(str(case.msa))
    if QUERY_LABEL not in rows:
        raise SystemExit(f"{case.name}: no {QUERY_LABEL!r} record in {case.msa}")
    panel = sorted(label for label in rows if label != QUERY_LABEL)
    variants = panel_variants(case, panel, hx)

    out: list[dict[str, object]] = []
    for n, region in enumerate(case.regions, start=1):
        start, end = int(region["msa_start"]), int(region["msa_end"])
        called = region["minor_parent"]
        called_clade = case.lineages.get(called, "?")
        for variant, kept, dropped in variants:
            stats = pendant_stats(rows, kept, start, end)
            if isinstance(stats, str):
                undefined.append((case.name, n, variant, len(kept), stats))
                continue
            out.append({
                "case": case.name, "klass": case.klass, "region": n, "variant": variant,
                "n_panel": len(kept), "n_dropped": len(panel) - len(kept),
                "dropped_clades": ",".join(dropped) or "-",
                "columns": end - start, "comparable": int(stats.pop("comparable")),
                "true_backbone": case.true_backbone or "NA",
                "true_donor": case.true_donor or "NA",
                "called_donor": called, "called_donor_clade": called_clade,
                "attribution": (
                    "NA" if not case.true_donor
                    else "correct" if hx.clade_match(called_clade, case.true_donor)
                    else "wrong"),
                # What the statistic is asked to detect: for a loo_* panel the donor's
                # lineage is gone by construction, whatever the run originally named.
                "donor_lineage_present": (
                    "no" if variant != "full" or case.klass == "donor_absent"
                    else "yes" if any(hx.clade_match(c, case.true_donor)
                                      for c in case.lineages.values())
                    else "no"),
                "methods": region["methods"],
                "n_methods": len([m for m in region["methods"].split(",") if m]),
                "length_bp": int(region["length_bp"]),
                "mean_sim_minor": float(region["mean_sim_minor"]),
                **stats,
            })
    return out


def primary(rows: list[dict[str, object]]) -> dict[str, object]:
    """The case's primary region: most callers agreeing, longest breaking the tie.

    Chosen on caller agreement and length alone -- never on the statistic under test -- so
    no case can be handed a favourable region by the choice.
    """
    return max(rows, key=lambda r: (int(r["n_methods"]), int(r["length_bp"])))


# --------------------------------------------------------------------------- reporting


@dataclass
class Separation:
    """How well one scaling separates donor-absent from donor-present values."""

    auc: float           # chance a random donor-absent value exceeds a random donor-present one
    best_cut: float      # the cut minimising total misclassifications
    best_errors: int     # how many it still makes
    cost_of_all: int     # donor-present cases falsely flagged by the cut that catches every absent
    sens_at_clean: float # share of donor-absent caught by the highest cut that flags no present


def _separation(present: np.ndarray, absent: np.ndarray) -> Separation:
    """Separation of ``absent`` (expected high) from ``present`` (expected low).

    Reports the Mann-Whitney AUC, the error-minimising cut, and the two operating points a
    decision actually turns on: what a caveat on every donor-absent case would cost in false
    caveats, and how many donor-absent cases a caveat that never fires wrongly would catch.
    """
    gt = (absent[:, None] > present[None, :]).sum()
    eq = (absent[:, None] == present[None, :]).sum()
    auc = float((gt + 0.5 * eq) / (absent.size * present.size))
    best = (10**9, float("nan"))
    for cut in np.unique(np.concatenate([present, absent])):
        errors = int((present >= cut).sum() + (absent < cut).sum())
        if errors < best[0]:
            best = (errors, float(cut))
    # Catch every donor-absent case: the cut must sit at or below the lowest of them.
    cost_of_all = int((present >= absent.min()).sum())
    # Flag no donor-present case: the cut must sit above the highest of them.
    sens_at_clean = float((absent > present.max()).sum() / absent.size)
    return Separation(auc, best[1], best[0], cost_of_all, sens_at_clean)


def _report(
    by_case: dict[str, list[dict[str, object]]],
    undefined: list[tuple[str, int, str, int, str]],
) -> None:
    # One row per (case, variant): the case's primary region, chosen within each variant.
    prim: dict[tuple[str, str], dict[str, object]] = {}
    for name, rows in by_case.items():
        for variant in {str(r["variant"]) for r in rows}:
            subset = [r for r in rows if r["variant"] == variant]
            prim[(name, variant)] = primary(subset)

    if undefined:
        print("\n=== regions where the statistic is undefined ===")
        print("Option A calibrates against the panel's own pendant lengths, so these carry no")
        print("statistic at all -- a limit on where it can be used, not a measurement.")
        seen: set[str] = set()
        for case, region, variant, n_panel, why in undefined:
            print(f"  {case:<20} region {region} [{variant:<10}] n_panel={n_panel}: {why}")
            seen.add(case)
        print(f"  -> {len(seen)} case(s) affected")

    print("\n=== the paired RSV contrast (same query, same true donor A.D.1.8) ===")
    print("case / variant                 n_panel  pendant_q  panel_med  ratio_median  robust_z"
          "  called")
    for key in (("rsv_a", "full"), ("rsv_a", "loo_exact"), ("rsv_a", "loo_parent"),
                ("donorabsent_rsv", "full")):
        r = prim.get(key)
        if r is None:
            continue
        tag = f"{key[0]} [{key[1]}]"
        real = " (real run, re-aligned)" if key[0] == "donorabsent_rsv" else ""
        print(f"{tag:<30} {r['n_panel']:>7}  {r['pendant_query']:>9.5f}  "
              f"{r['pendant_panel_median']:>9.5f}  {r['ratio_median']:>12.4f}  "
              f"{r['robust_z']:>8.3f}  {r['called_donor_clade']}{real}")

    print("\n=== primary region per case and panel ===")
    print("case\tvariant\tdonor_lineage\tattribution\tn_panel\tn_meth\tpendant_q\t"
          "panel_med\tquantile\tratio_median\trobust_z")
    for name, variant in sorted(prim):
        r = prim[(name, variant)]
        print(f"{name}\t{variant}\t{r['donor_lineage_present']}\t{r['attribution']}\t"
              f"{r['n_panel']}\t{r['n_methods']}\t{r['pendant_query']:.5f}\t"
              f"{r['pendant_panel_median']:.5f}\t{r['quantile']:.4f}\t"
              f"{r['ratio_median']:.4f}\t{r['robust_z']:.3f}")

    present = [r for r in prim.values()
               if r["variant"] == "full" and r["klass"] == "positive"
               and r["donor_lineage_present"] == "yes"]
    real_absent = [r for r in prim.values()
                   if r["klass"] == "donor_absent" and r["variant"] == "full"]

    print("\n=== why the statistic works on some panels and not others ===")
    print("ratio_median against the panel's spread share (its own median pendant over its")
    print("median pairwise distance). A DENSE panel has short pendants -- a small share -- so an")
    print("unsampled donor's long pendant stands out; a sparse panel's own pendants are already")
    print("as long as the query's, and the signal has nothing to stand out against.")
    loo = [r for r in prim.values() if r["variant"] != "full"]
    if len(loo) >= 3:
        dens = np.array([float(r["panel_spread_share"]) for r in loo])
        rat = np.array([float(r["ratio_median"]) for r in loo])
        ok = np.isfinite(dens) & np.isfinite(rat)
        if ok.sum() >= 3:
            # Spearman: rank correlation, so the heavy tail in ratio_median cannot drive it.
            rd = np.argsort(np.argsort(dens[ok])).astype(float)
            rr = np.argsort(np.argsort(rat[ok])).astype(float)
            rho = float(np.corrcoef(rd, rr)[0, 1])
            print(f"  simulated donor-absent panels: n={int(ok.sum())}, "
                  f"Spearman rho(panel_spread_share, ratio_median) = {rho:+.3f}")
        print("  case                 variant     n_panel  spread_share  ratio_median")
        for r in sorted(loo, key=lambda r: float(r["panel_spread_share"])):
            print(f"  {r['case']:<20} {r['variant']:<10} {r['n_panel']:>8}  "
                  f"{float(r['panel_spread_share']):>12.4f}  {float(r['ratio_median']):>12.4f}")

    blocks: list[tuple[str, list[dict[str, object]]]] = []
    if real_absent:
        blocks.append(("REAL donor-absent runs (panel re-recruited and re-aligned)",
                       real_absent))
    for variant, label in (("loo_parent", "simulated: parent sub-clade's rows dropped"),
                           ("loo_exact", "simulated: donor's own rows dropped, siblings kept")):
        rows = [r for r in prim.values() if r["variant"] == variant]
        if rows:
            blocks.append((label, rows))

    for label, absent in blocks:
        print(f"\n=== Option A falsification test -- {label} ===")
        print(f"donor-present n={len(present)}   donor-absent n={len(absent)}")
        for stat in (CONTROL, *STATS):
            def value(row: dict[str, object], stat: str = stat) -> float:
                if stat == CONTROL:
                    return -float(row["mean_sim_minor"])
                return float(row[stat])
            pv = np.array([value(r) for r in present if np.isfinite(value(r))])
            av = np.array([value(r) for r in absent if np.isfinite(value(r))])
            if pv.size == 0 or av.size == 0:
                continue
            sep = _separation(pv, av)
            overlap = int(((av >= pv.min()) & (av <= pv.max())).sum())
            tag = f"{stat} (control)" if stat == CONTROL else stat
            print(f"{tag:>24}: present [{pv.min():>9.4f}, {pv.max():>8.4f}] med "
                  f"{np.median(pv):>8.4f} | absent [{av.min():>8.4f}, {av.max():>9.4f}] med "
                  f"{np.median(av):>8.4f} | AUC {sep.auc:.3f}")
            print(f"{'':>24}  {overlap}/{av.size} donor-absent inside the donor-present range; "
                  f"best cut {sep.best_cut:.4f} -> {sep.best_errors}/{pv.size + av.size} wrong")
            print(f"{'':>24}  to flag every donor-absent case: {sep.cost_of_all}/{pv.size} "
                  f"donor-present cases falsely caveated; "
                  f"to flag none wrongly: {sep.sens_at_clean:.0%} of donor-absent caught")


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, nargs="+", default=[HYBRIDS, PROBE],
                    help="per-case artifact directories to read; the donor-present baseline "
                         "and the donor-absent runs normally live in different ones "
                         "(default: validation/data/hybrids and .../donor_absent_probe)")
    ap.add_argument("--out", type=Path, default=None, help="write the per-region TSV here")
    args = ap.parse_args(argv)

    roots = [d for d in args.data if d.is_dir()]
    if not roots:
        listed = ", ".join(str(d) for d in args.data)
        print(f"No case artifacts under {listed}. Run validation/run_hybrids.py first.",
              file=sys.stderr)
        return 2
    for missing in (d for d in args.data if not d.is_dir()):
        print(f"  {missing}: not present; skipped", file=sys.stderr)

    hx = _load_harness()
    undefined: list[tuple[str, int, str, int, str]] = []
    measured: list[dict[str, object]] = []
    by_case: dict[str, list[dict[str, object]]] = {}
    for case_dir in sorted((p for root in roots for p in root.iterdir() if p.is_dir()),
                           key=lambda p: p.name):
        case = load_case(case_dir)
        if case is None:
            print(f"  {case_dir.name}: no completed run; skipped", file=sys.stderr)
            continue
        rows = measure_case(case, hx, undefined)
        if rows:
            measured.extend(rows)
            by_case[case.name] = rows

    if not measured:
        print("No regions measured.", file=sys.stderr)
        return 2

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with open(args.out, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerows(measured)
        print(f"wrote {args.out} ({len(measured)} rows)")

    _report(by_case, undefined)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
