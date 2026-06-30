"""The self-contained HTML report: the verdict, cards, query mosaic, region and
coverage tables, the parent-free signal section, and ``write_html_report``.

Static styling and content (CSS, glossary, references) come from
``report_assets``; the colour palette and the embedded interactive figure come
from ``report_plots``.
"""

from __future__ import annotations

import html
import logging
from pathlib import Path

from .analyze import AnalysisResult, stats_sort_key
from .coverage import CoverageGap
from .diagnostics import RecombinationSignal
from .regions import Region
from .report_assets import _CSS, _GLOSSARY, _REFERENCES
from .report_context import ReportContext
from .report_plots import GREY, _color_map, build_interactive_figure
from .similarity import WindowSimilarity
from .typing import LineageMap, typed


def _fmt_int(n) -> str:
    return f"{int(round(float(n))):,}"


def _fmt_kb(n) -> str:
    n = float(n)
    return f"{n / 1000:.1f} kb" if abs(n) >= 1000 else f"{int(round(n))} bp"


def _swatch(color: str) -> str:
    return f'<span class="sw" style="background:{color}"></span>'


def _summary(
    result: WindowSimilarity, regions: list[Region], datasets: list[str]
) -> dict:
    """Headline numbers: query length, backbone, donor list, recombinant fraction."""
    qcum = result.query_cumulative
    query_len = int(qcum[-1]) if len(qcum) else 0
    if regions:
        major = regions[0].major_parent
    elif datasets:
        major = datasets[0]
    else:
        major = "n/a"
    present = [r for r in regions if not r.donor_absent]
    absent = [r for r in regions if r.donor_absent]
    recomb_bp = sum(max(0, r.query_end - r.query_start) for r in present)
    minors: list[str] = []
    for r in present:
        if r.minor_parent not in minors:
            minors.append(r.minor_parent)
    return {
        "query_len": query_len, "major": major,
        "n_regions": len(present), "n_absent": len(absent),
        "recomb_bp": recomb_bp,
        "pct": (100.0 * recomb_bp / query_len) if query_len else 0.0,
        "minors": minors,
        "confidence": _confidence(present, absent),
    }


def _confidence(present: list[Region], absent: list[Region]) -> str | None:
    """A plain-language confidence label for the call, from FDR / support / coverage.

    None when there is nothing to qualify (no recombination). A donor-present region
    with a tiny q-value, a clear discordant-site majority, and a well-covered donor is
    *high*; a significant but marginal or under-covered call is *moderate*; a
    donor-absent or barely-significant call is *low*.
    """
    if not present:
        return "low" if absent else None
    qs = [r.qvalue for r in present if r.qvalue is not None]
    sups = [r.support for r in present if r.support is not None]
    min_q = min(qs) if qs else 1.0
    max_sup = max(sups) if sups else 0.0
    undercovered = any(r.donor_undercovered for r in present)
    # Agreement is a first-class confidence signal: a region independently called by
    # more than one method is high-confidence even if no single method is overwhelming.
    agreement = any(len(r.methods) >= 2 for r in present)
    if not undercovered and (agreement or (min_q <= 1e-5 and max_sup >= 0.70)):
        return "high"
    if min_q <= 0.05 and not undercovered:
        return "moderate"
    return "low"


def _absent_clause(s: dict) -> str:
    if not s["n_absent"]:
        return ""
    word = "region" if s["n_absent"] == 1 else "regions"
    return (f' {s["n_absent"]} further {word} diverge from every reference &mdash; a '
            f'<strong>donor may be missing</strong> (see Reference coverage).')


def _query_lineage_clause(query: str, query_lineage: str | None) -> str:
    """A clause naming the query's own typed lineage (omitted if it equals the label)."""
    if not query_lineage or query_lineage == query:
        return ""
    return f' The query is typed as <strong>{html.escape(query_lineage)}</strong>.'


def _verdict_html(
    s: dict, query: str, colors: dict[str, str], lineage_map: LineageMap | None = None,
    query_lineage: str | None = None,
) -> str:
    major = html.escape(typed(s["major"], lineage_map))
    q_clause = _query_lineage_clause(query, query_lineage)
    if s["n_regions"] == 0:
        lead = (f'No recombination among the present references &mdash; the query is most '
                f'similar to <strong>{major}</strong> throughout.')
        return f'<p class="verdict">{lead}{q_clause}{_absent_clause(s)}</p>'
    donors = ", ".join(
        f'{_swatch(colors.get(m, GREY))}<strong>{html.escape(typed(m, lineage_map))}</strong>'
        for m in s["minors"]
    )
    word = "region" if s["n_regions"] == 1 else "regions"
    conf = s.get("confidence")
    conf_clause = (
        f' <strong>{html.escape(conf)}</strong> confidence.' if conf else ""
    )
    return (
        f'<p class="verdict">The query is a <strong>recombinant</strong>: a '
        f'<strong>{major}</strong> backbone carrying {s["n_regions"]} donor {word} from '
        f'{donors}, covering <span class="mono">{_fmt_kb(s["recomb_bp"])}</span> '
        f'(<span class="mono">{s["pct"]:.1f}%</span>) of the query.{conf_clause}'
        f'{q_clause}{_absent_clause(s)}</p>'
    )


def _cards_html(s: dict, colors: dict[str, str], lineage_map: LineageMap | None = None) -> str:
    donors = "".join(
        f'{_swatch(colors.get(m, GREY))}{html.escape(typed(m, lineage_map))} ' for m in s["minors"]
    ) or "&mdash;"
    regions_val = f'<span class="big">{s["n_regions"]}</span>'
    if s["n_absent"]:
        regions_val += f'<span class="sub">&nbsp;+{s["n_absent"]} donor-absent</span>'
    cards = [
        ("Recombinant regions", regions_val),
        ("Backbone &middot; major parent",
         f'{_swatch(colors.get(s["major"], GREY))}{html.escape(typed(s["major"], lineage_map))}'),
        ("Donor &middot; minor parent", donors),
        ("Query recombinant",
         f'<span class="big">{s["pct"]:.1f}%</span>'
         f'<span class="sub">&nbsp;{_fmt_kb(s["recomb_bp"])} of {_fmt_kb(s["query_len"])}</span>'),
    ]
    items = "".join(
        f'<div class="card"><div class="k">{k}</div><div class="v">{v}</div></div>'
        for k, v in cards
    )
    return f'<div class="cards">{items}</div>'


def _mosaic_html(
    regions: list[Region], colors: dict[str, str], s: dict, gaps: list[CoverageGap],
    lineage_map: LineageMap | None = None,
) -> str:
    query_len = s["query_len"] or 1
    backbone = colors.get(s["major"], GREY)
    segs = "".join(
        f'<div class="seg" style="left:{100.0 * r.query_start / query_len:.3f}%;'
        f'width:{100.0 * max(0, r.query_end - r.query_start) / query_len:.3f}%;'
        f'background:{colors.get(r.minor_parent, GREY)}" '
        f'title="{html.escape(r.minor_parent)}: '
        f'{_fmt_int(r.query_start)}-{_fmt_int(r.query_end)} bp"></div>'
        for r in regions if not r.donor_absent  # absent-donor regions show as gaps
    )
    segs += "".join(
        f'<div class="gap" style="left:{100.0 * g.query_start / query_len:.3f}%;'
        f'width:{100.0 * max(0, g.query_end - g.query_start) / query_len:.3f}%" '
        f'title="low coverage ({g.kind}): closest reference {html.escape(g.best_label)} '
        f'~{g.mean_best:.2f}"></div>'
        for g in gaps
    )
    axis = "".join(
        f"<span>{_fmt_kb(frac * query_len)}</span>" for frac in (0, 0.25, 0.5, 0.75, 1.0)
    )
    legend = (f'<span class="leg">{_swatch(backbone)}'
              f'{html.escape(typed(s["major"], lineage_map))} (backbone)</span>')
    legend += "".join(
        f'<span class="leg">{_swatch(colors.get(m, GREY))}'
        f'{html.escape(typed(m, lineage_map))} (donor)</span>'
        for m in s["minors"]
    )
    if gaps:
        legend += '<span class="leg"><span class="hatch"></span>low coverage</span>'
    return (
        f'<div class="mosaic"><div class="track" style="--bb:{backbone}">{segs}</div>'
        f'<div class="axis">{axis}</div><div class="legend">{legend}</div>'
        f'<p class="cap">The query genome painted by its closest reference per window: solid '
        f'backbone is the major parent; coloured segments are donor regions; hatched bands are '
        f'low-coverage stretches where even the closest reference is a poor match.</p></div>'
    )


def _regions_html(
    regions: list[Region], colors: dict[str, str], query_len: int,
    lineage_map: LineageMap | None = None,
) -> str:
    if not regions:
        return '<p class="empty">No recombinant regions were called.</p>'
    head = (
        "<tr><th>Donor (minor)</th><th>Backbone (major)</th><th>Query span (bp)</th>"
        "<th>Length</th><th>% query</th><th>Sim donor</th><th>Sim backbone</th>"
        "<th>Support</th><th>q-value</th><th>Breakpoint</th><th>Method(s)</th></tr>"
    )
    rows = ""
    for r in regions:
        qlen = max(0, r.query_end - r.query_start)
        pct = (100.0 * qlen / query_len) if query_len else 0.0
        if r.donor_absent:
            swatch = '<span class="sw hatch"></span>'
            donor = '<em>missing</em>'
            flag = ('<span class="flag" title="query is far from every reference here; '
                    'the true donor is likely absent">donor absent</span>')
        else:
            swatch = _swatch(colors.get(r.minor_parent, GREY))
            donor = html.escape(typed(r.minor_parent, lineage_map))
            flag = ('<span class="flag" title="donor is itself a poor match">low conf</span>'
                    if r.donor_undercovered else "")
        support = "&ndash;" if r.support is None else f"{r.support:.2f}"
        qval = "&ndash;" if r.qvalue is None else f"{r.qvalue:.2g}"
        if r.breakpoint_lo is None:
            bp = "&ndash;"
        elif r.breakpoint_lo == r.breakpoint_hi:
            bp = _fmt_int(r.breakpoint_lo)
        else:
            bp = f"{_fmt_int(r.breakpoint_lo)}&ndash;{_fmt_int(r.breakpoint_hi)}"
        if not r.methods:
            methods_cell = "&ndash;"
        else:
            label = ", ".join(r.methods)
            agree = (' <span class="flag" title="called by more than one method">agree</span>'
                     if len(r.methods) >= 2 else "")
            phi = (' <span class="flag" title="overlaps a parent-free Rmin interval">'
                   'PHI</span>' if r.parent_free_support else "")
            methods_cell = f'{html.escape(label)}{agree}{phi}'
        rows += (
            "<tr>"
            f'<td class="lbl">{swatch}{donor}{flag}</td>'
            f'<td class="lbl">{html.escape(typed(r.major_parent, lineage_map))}</td>'
            f'<td class="num">{_fmt_int(r.query_start)}&ndash;{_fmt_int(r.query_end)}</td>'
            f'<td class="num">{_fmt_kb(qlen)}</td>'
            f'<td class="num">{pct:.1f}%</td>'
            f'<td class="num">{r.mean_sim_minor:.3f}</td>'
            f'<td class="num">{r.mean_sim_major:.3f}</td>'
            f'<td class="num strong">{support}</td>'
            f'<td class="num">{qval}</td>'
            f'<td class="num">{bp}</td>'
            f'<td class="lbl">{methods_cell}</td>'
            "</tr>"
        )
    return f'<div class="scroll"><table class="table">{head}{rows}</table></div>'


def _method_comparison_html(
    breakdown: list[dict], methods_run: tuple[str, ...], per_major: dict[str, str],
    lineage_map: LineageMap | None = None,
) -> str:
    """A compact region x method agreement matrix; omitted for a single-method run."""
    if len(methods_run) < 2:
        return ""
    majors = ", ".join(
        f'<span class="mono">{html.escape(m)}</span> &rarr; '
        f'{html.escape(typed(per_major.get(m, "n/a"), lineage_map))}'
        for m in methods_run
    )
    intro = (
        f'<p class="cap">Each caller ran independently on the same alignment; a region '
        f'found by more than one is more trustworthy (and raises the confidence above). '
        f'Backbone per method &mdash; {majors}.</p>'
    )
    if not breakdown:
        return intro + '<p class="empty">No regions were called by any method.</p>'
    head = (
        "<tr><th>Donor (minor)</th><th>Query span (bp)</th>"
        + "".join(f"<th>{html.escape(m)}</th>" for m in methods_run)
        + "<th>PHI/Rmin</th></tr>"
    )
    rows = ""
    for b in breakdown:
        cells = "".join(
            f'<td class="num">{"&check;" if m in b["per_method_support"] else "&middot;"}</td>'
            for m in methods_run
        )
        phi = "&check;" if b["parent_free_support"] else "&middot;"
        rows += (
            "<tr>"
            f'<td class="lbl">{html.escape(typed(b["minor_parent"], lineage_map))}</td>'
            f'<td class="num">{_fmt_int(b["query_start"])}&ndash;{_fmt_int(b["query_end"])}</td>'
            f'{cells}<td class="num">{phi}</td>'
            "</tr>"
        )
    return f'{intro}<div class="scroll"><table class="table">{head}{rows}</table></div>'


def _method_section(
    method_breakdown: list[dict] | None, methods_run: tuple[str, ...],
    per_major: dict[str, str] | None, lineage_map: LineageMap | None,
) -> str:
    """Wrap the method-comparison table in a report section (empty for one method)."""
    body = _method_comparison_html(
        method_breakdown or [], methods_run, per_major or {}, lineage_map
    )
    if not body:
        return ""
    return (
        '<section class="section"><div class="eyebrow">Method comparison</div>'
        f"{body}</section>"
    )


def _winners_html(
    analysis: AnalysisResult, colors: dict[str, str], lineage_map: LineageMap | None = None
) -> str:
    items = sorted(analysis.winners_with_ties.items(), key=lambda x: x[1], reverse=True)
    if not items:
        return '<p class="empty">No window winners recorded.</p>'
    top = items[0][1] or 1
    bars = "".join(
        f'<div class="barrow"><span class="blabel">{_swatch(colors.get(label, GREY))}'
        f'{html.escape(typed(label, lineage_map))}</span><span class="btrack">'
        f'<span class="bfill" style="width:{100.0 * count / top:.1f}%;'
        f'background:{colors.get(label, GREY)}"></span></span>'
        f'<span class="bnum">{_fmt_int(count)}</span></div>'
        for label, count in items
    )
    return f'<div class="bars">{bars}</div>'


def _stats_html(
    analysis: AnalysisResult, major: str, lineage_map: LineageMap | None = None
) -> str:
    head = "".join(f"<th>{html.escape(h)}</th>" for h in analysis.stats_header)
    rows = ""
    for dataset, values in sorted(analysis.stats.items(), key=stats_sort_key, reverse=True):
        cls = ' class="hl"' if dataset == major else ""
        cells = f'<td class="lbl">{html.escape(typed(dataset, lineage_map))}</td>' + "".join(
            f'<td class="num">{html.escape(str(v))}</td>' for v in values
        )
        rows += f"<tr{cls}>{cells}</tr>"
    return f'<div class="scroll"><table class="table"><tr>{head}</tr>{rows}</table></div>'


def _methods_html(provenance: dict[str, str]) -> str:
    glossary = "".join(
        f"<dt>{html.escape(t)}</dt><dd>{html.escape(d)}</dd>" for t, d in _GLOSSARY
    )
    references = "".join(
        f"<li><strong>{html.escape(t)}.</strong> {html.escape(c)}</li>"
        for t, c in _REFERENCES
    )
    params = "".join(
        f'<tr><th>{html.escape(str(k))}</th><td class="mono">{html.escape(str(v))}</td></tr>'
        for k, v in provenance.items()
    )
    return (
        '<details class="methods"><summary>Methods &amp; glossary</summary>'
        '<p class="cap">Tessera segments the query against the reference panel with an HMM '
        '(jpHMM-style) and reports a region only when its donor beats the major parent on the '
        'sites that distinguish them (a sign test on discordant sites, immune to window '
        'overlap; Benjamini-Hochberg FDR across segments), with a posterior breakpoint '
        'interval. By default it runs an ensemble of callers (HMM and the 3SEQ triplet '
        'test) and merges their regions into one consensus, so a region found by more than '
        'one method is flagged as agreeing and treated as higher confidence (see the '
        'caller line under Run parameters). It remains an indicative screen, not a full '
        'phylogenetic test (e.g. GARD) -- confirm strong candidates.</p>'
        f'<dl class="glossary">{glossary}</dl>'
        f'<h3>References</h3><ul class="refs">{references}</ul>'
        f'<h3>Run parameters</h3><table class="kv">{params}</table></details>'
    )


def _footer_html(provenance: dict[str, str]) -> str:
    files = [
        "recombination_regions.tsv", "recombination_methods.tsv", "coverage_gaps.tsv",
        "recombination_profile.tsv", "window_winners.tsv", "similarity_stats.tsv",
        "similarity_windows.tsv", "similarity_top*.pdf", "similarity_pair.pdf",
    ]
    flist = ", ".join(f"<code>{f}</code>" for f in files)
    ver = html.escape(provenance.get("tessera version", ""))
    date = html.escape(provenance.get("date (UTC)", ""))
    return (
        f'<footer><div>Generated by Tessera <span class="mono">{ver}</span> &middot; '
        f'<span class="mono">{date}</span> UTC</div>'
        f'<div>Companion files in this folder: {flist}.</div></footer>'
    )


def _caveat_html(gaps: list[CoverageGap], threshold: float) -> str:
    if not gaps:
        return ""
    n = len(gaps)
    word = "region" if n == 1 else "regions"
    verb = "is" if n == 1 else "are"
    total = sum(max(0, g.query_end - g.query_start) for g in gaps)
    worst = min(gaps, key=lambda g: g.mean_best)
    return (
        f'<div class="caveat"><span class="ic">&#9888;</span><div>'
        f'<strong>Possible missing reference.</strong> {n} {word} '
        f'(<span class="mono">{_fmt_kb(total)}</span> of the query) {verb} poorly covered: '
        f'even the closest reference stays below <span class="mono">{threshold:.2f}</span>. '
        f'The weakest is query <span class="mono">{_fmt_int(worst.query_start)}&ndash;'
        f'{_fmt_int(worst.query_end)}</span>, where <strong>{html.escape(worst.best_label)}'
        f'</strong> reaches only <span class="mono">{worst.mean_best:.2f}</span> &mdash; its '
        f'true source may not be in the collection. See <em>Reference coverage</em> below, or '
        f'run <span class="mono">tessera find-references</span> to search NCBI.</div></div>'
    )


def _coverage_html(gaps: list[CoverageGap], threshold: float) -> str:
    intro = (
        f'<p class="cap">Stretches where even the closest reference is below the '
        f'<span class="mono">{threshold:.3f}</span> best-similarity threshold. '
        f'<strong>divergent</strong> = the query is genuinely far from every reference '
        f'(a likely missing reference); <strong>low information</strong> = too few comparable '
        f'bases to judge.</p>'
    )
    if not gaps:
        return (
            '<p class="cap">Every region has a close reference &mdash; no coverage gaps at the '
            f'<span class="mono">{threshold:.3f}</span> threshold.</p>'
        )
    head = (
        "<tr><th>Closest reference</th><th>Query span (bp)</th><th>Length</th>"
        "<th>Windows</th><th>Mean best sim</th><th>Kind</th></tr>"
    )
    rows = "".join(
        "<tr>"
        f'<td class="lbl">{html.escape(g.best_label)}</td>'
        f'<td class="num">{_fmt_int(g.query_start)}&ndash;{_fmt_int(g.query_end)}</td>'
        f'<td class="num">{_fmt_kb(g.length_bp)}</td>'
        f'<td class="num">{_fmt_int(g.n_windows)}</td>'
        f'<td class="num strong">{g.mean_best:.3f}</td>'
        f'<td class="lbl">{html.escape(g.kind.replace("_", " "))}</td>'
        "</tr>"
        for g in gaps
    )
    return f'{intro}<div class="scroll"><table class="table">{head}{rows}</table></div>'


def _signal_profile_div(signal: RecombinationSignal) -> str:
    """A compact line of the PHI profile (mean local incompatibility) vs MSA position;
    a dip marks where recombination signal concentrates."""
    if not signal.profile:
        return ""
    import plotly.graph_objects as go

    xs = [msa_pos for msa_pos, _, _ in signal.profile]
    ys = [value for _, _, value in signal.profile]
    fig = go.Figure(go.Scatter(x=xs, y=ys, mode="lines", line={"color": "#6b46c1", "width": 1}))
    fig.update_layout(
        template="plotly_white", height=220, margin={"l": 50, "r": 20, "t": 10, "b": 40},
        xaxis_title="MSA position", yaxis_title="local incompatibility",
    )
    return fig.to_html(full_html=False, include_plotlyjs=False)


def _signal_html(signal: RecombinationSignal | None, alpha: float = 0.05) -> str:
    """The parent-free recombination-signal section: PHI p-value, Rmin, profile."""
    if signal is None:
        return (
            '<p class="cap">Too few informative sites in the alignment for a parent-free '
            'recombination test.</p>'
        )
    significant = signal.phi_p < alpha
    verdict = (
        '<strong>significant recombination signal</strong>' if significant
        else 'no significant recombination signal'
    )
    intervals = ", ".join(
        f"{_fmt_int(lo)}&ndash;{_fmt_int(hi)}" for lo, hi in signal.rmin_intervals[:8]
    )
    if len(signal.rmin_intervals) > 8:
        intervals += ", &hellip;"
    intro = (
        '<p class="cap">A parent-free test that asks only whether the alignment carries '
        'recombination at all (no candidate parents needed), strongest where divergence is '
        'low or the true donor is absent. The PHI test (Bruen et al. 2006) is significant '
        'when nearby sites are more compatible than a random reordering; Rmin (Hudson &amp; '
        'Kaplan 1985) is the minimum number of recombination events the four-gamete test '
        'forces, with the intervals (query coordinates) as breakpoint candidates.</p>'
    )
    rows = (
        f'<tr><td class="lbl">PHI test</td><td class="num strong">p = {signal.phi_p:.4g}</td>'
        f'<td class="lbl">{verdict} (alpha {alpha:g}; {signal.n_informative} informative '
        f'sites)</td></tr>'
        f'<tr><td class="lbl">Min recombination events (Rmin)</td>'
        f'<td class="num strong">{signal.rmin}</td>'
        f'<td class="lbl">{"intervals " + intervals if intervals else "none"}</td></tr>'
    )
    return (
        f'{intro}<div class="scroll"><table class="table">{rows}</table></div>'
        f'{_signal_profile_div(signal)}'
    )


def write_html_report(
    result: WindowSimilarity,
    analysis: AnalysisResult,
    regions: list[Region],
    datasets: list[str],
    provenance: dict[str, str],
    output_dir: Path,
    logger: logging.Logger,
    ctx: ReportContext,
) -> Path:
    """Write a single self-contained ``report.html``."""
    gaps = ctx.gaps
    lineage_map = ctx.lineage_map
    threshold = ctx.coverage_threshold
    fig = build_interactive_figure(result, datasets, regions, gaps)
    plot_div = fig.to_html(full_html=False, include_plotlyjs="inline")

    colors = _color_map(datasets)
    s = _summary(result, regions, datasets)
    organism_html = (
        f'<div class="organism">{html.escape(ctx.organism)}</div>' if ctx.organism else ""
    )
    extras = "".join(
        f'<section class="section"><div class="eyebrow">{html.escape(title)}</div>{body}</section>'
        for title, body in (ctx.extra_sections or [])
    )

    doc = (
        '<!DOCTYPE html>\n<html lang="en"><head><meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '<link rel="icon" href="data:,">\n'
        f"<title>Tessera report: {html.escape(result.query)}</title>\n"
        f"<style>{_CSS}</style></head><body><div class=\"wrap\">"
        '<header><div class="eyebrow">Tessera &middot; recombination report</div>'
        f'<h1 class="mono">{html.escape(result.query)}</h1>'
        f"{organism_html}"
        f"{_verdict_html(s, result.query, colors, lineage_map, ctx.query_lineage)}"
        f"{_caveat_html(gaps, threshold)}</header>"
        f"{_cards_html(s, colors, lineage_map)}"
        '<section class="section"><div class="eyebrow">Query mosaic</div>'
        f"{_mosaic_html(regions, colors, s, gaps, lineage_map)}</section>"
        '<section class="section"><div class="eyebrow">Recombinant regions</div>'
        f'{_regions_html(regions, colors, s["query_len"], lineage_map)}</section>'
        f"{_method_section(ctx.method_breakdown, ctx.methods_run, ctx.per_major, lineage_map)}"
        '<section class="section"><div class="eyebrow">Reference coverage</div>'
        f"{_coverage_html(gaps, threshold)}</section>"
        f"{extras}"
        '<section class="section"><div class="eyebrow">Similarity across the alignment</div>'
        '<p class="cap">Each line is one reference\'s similarity to the query along the '
        'alignment; coloured bands are called donor regions, hatched bands are low-coverage '
        'stretches. Drag to zoom, hover for values.</p>'
        f"{plot_div}</section>"
        '<section class="section"><div class="eyebrow">Recombination signal (parent-free)</div>'
        f"{_signal_html(ctx.signal)}</section>"
        '<section class="section"><div class="eyebrow">Window winners</div>'
        '<p class="cap">Windows in which each reference is the query\'s closest match '
        '(ties included).</p>'
        f"{_winners_html(analysis, colors, lineage_map)}</section>"
        '<section class="section"><div class="eyebrow">Per-dataset similarity statistics</div>'
        f'{_stats_html(analysis, s["major"], lineage_map)}</section>'
        f'<section class="section">{_methods_html(provenance)}</section>'
        f"{_footer_html(provenance)}"
        "</div></body></html>"
    )

    out = output_dir / "report.html"
    out.write_text(doc)
    logger.info("Report written: %s", out)
    return out
