"""Static presentation constants for the HTML report: stylesheet, glossary, references.

These are pure content (no logic), kept apart from the HTML builders in
``report_html`` so the rendering code stays readable and the styling can be
edited in one place.
"""

from __future__ import annotations

_CSS = """
:root{
  --ink:#11161f;--muted:#5b6573;--faint:#8b94a3;--line:#e6e9ee;
  --panel:#f6f7f9;--paper:#fff;--bg:#fbfcfd;--accent:#2b6cb0;--hl:#eef4fb;
  --mono:ui-monospace,"SF Mono","JetBrains Mono",Menlo,Consolas,monospace;
  --sans:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font-family:var(--sans);
  font-size:15px;line-height:1.55;-webkit-font-smoothing:antialiased}
.wrap{max-width:1120px;margin:0 auto;padding:48px 28px 72px}
.mono{font-family:var(--mono);font-variant-numeric:tabular-nums}
.eyebrow{font-size:11px;font-weight:600;letter-spacing:.16em;text-transform:uppercase;color:var(--faint)}
.sw{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:6px;flex:none}
.sw.hatch{background:repeating-linear-gradient(45deg,#64748b 0 3px,#cbd5e1 3px 6px)}
header{border-bottom:1px solid var(--line);padding-bottom:26px;margin-bottom:30px}
header h1{font-size:30px;font-weight:600;margin:.35em 0 .15em;letter-spacing:-.01em;word-break:break-word}
.organism{font-size:16px;font-weight:600;color:var(--muted);margin:0 0 .55em}
.verdict{font-size:19px;line-height:1.5;margin:0;max-width:76ch}
.verdict strong{font-weight:650}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:14px;margin-bottom:40px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:15px 18px}
.card .k{font-size:11px;font-weight:600;letter-spacing:.07em;text-transform:uppercase;color:var(--faint);margin-bottom:9px}
.card .v{font-size:17px;font-weight:550;display:flex;align-items:center;flex-wrap:wrap;gap:3px}
.big{font-size:26px;font-weight:650;font-family:var(--mono)}
.card .sub{font-size:12px;color:var(--muted);font-weight:400;font-family:var(--mono)}
.section{margin:42px 0}
.section>.eyebrow{margin-bottom:14px}
.cap{color:var(--muted);font-size:13px;margin:.2em 0 1em;max-width:82ch}
.mosaic .track{position:relative;height:50px;border-radius:8px;overflow:hidden;border:1px solid var(--line);
  background:repeating-linear-gradient(90deg,transparent 0 3px,rgba(0,0,0,.015) 3px 6px),var(--bb);
  box-shadow:inset 0 1px 2px rgba(0,0,0,.05)}
.mosaic .seg{position:absolute;top:0;bottom:0;min-width:1px;box-shadow:0 0 0 1px rgba(255,255,255,.25) inset}
.mosaic .axis{display:flex;justify-content:space-between;margin-top:6px;font-family:var(--mono);font-size:11px;color:var(--faint)}
.mosaic .gap{position:absolute;top:0;bottom:0;min-width:1px;border-left:1px solid #64748b;border-right:1px solid #64748b;
  background:repeating-linear-gradient(45deg,rgba(100,116,139,.42) 0 4px,rgba(100,116,139,.14) 4px 8px)}
.mosaic .legend{display:flex;flex-wrap:wrap;gap:16px;margin-top:13px;font-size:13px}
.mosaic .leg{display:inline-flex;align-items:center;color:var(--muted)}
.mosaic .leg .hatch{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:6px;
  background:repeating-linear-gradient(45deg,#64748b 0 3px,#cbd5e1 3px 6px)}
.caveat{display:flex;gap:11px;align-items:flex-start;background:#fff7ed;border:1px solid #fed7aa;
  border-left:3px solid #dd6b20;border-radius:10px;padding:13px 16px;margin:20px 0 0;font-size:14px;color:#7c4a13;max-width:76ch}
.caveat .ic{font-weight:700;color:#dd6b20;font-size:16px;line-height:1.3}
.flag{display:inline-block;margin-left:7px;font-size:10.5px;color:#92600a;border:1px solid #e6c98a;
  background:#fdf6e7;border-radius:4px;padding:0 5px;font-weight:600;letter-spacing:.02em;vertical-align:middle}
.scroll{overflow-x:auto}
table.table{border-collapse:collapse;width:100%;font-size:14px}
table.table th{font-size:11px;font-weight:600;letter-spacing:.05em;text-transform:uppercase;color:var(--faint);
  text-align:right;padding:8px 14px;border-bottom:2px solid var(--line);white-space:nowrap}
table.table td{padding:9px 14px;border-bottom:1px solid var(--line);text-align:right}
table.table th:first-child,table.table td:first-child{text-align:left}
table.table .num{font-family:var(--mono);font-variant-numeric:tabular-nums;white-space:nowrap}
table.table .lbl{text-align:left;font-weight:500;white-space:nowrap}
table.table .strong{font-weight:650;color:var(--ink)}
table.table tr:hover td{background:var(--panel)}
table.table tr.hl td{background:var(--hl);font-weight:600}
.empty{color:var(--muted);font-style:italic}
.bars{display:flex;flex-direction:column;gap:7px;max-width:780px}
.barrow{display:grid;grid-template-columns:190px 1fr 60px;align-items:center;gap:12px}
.blabel{display:inline-flex;align-items:center;font-size:13px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.btrack{height:14px;background:var(--panel);border:1px solid var(--line);border-radius:7px;overflow:hidden}
.bfill{display:block;height:100%}
.bnum{font-family:var(--mono);font-size:13px;color:var(--muted);text-align:right}
details.methods{border:1px solid var(--line);border-radius:12px;padding:2px 18px;background:var(--paper)}
details.methods summary{cursor:pointer;font-weight:600;padding:13px 0;list-style:none;font-size:14px}
details.methods summary::-webkit-details-marker{display:none}
details.methods summary::before{content:"+";display:inline-block;width:1.3em;color:var(--faint);font-family:var(--mono)}
details.methods[open] summary::before{content:"-"}
details.methods[open]{padding-bottom:18px}
.glossary{display:grid;grid-template-columns:max-content 1fr;gap:7px 20px;margin:4px 0 20px}
.glossary dt{font-weight:600;font-size:13px}
.glossary dd{margin:0;color:var(--muted);font-size:13px}
.refs{margin:4px 0 20px;padding-left:18px}
.refs li{color:var(--muted);font-size:13px;margin:4px 0}
details.methods h3{font-size:11px;text-transform:uppercase;letter-spacing:.08em;color:var(--faint);margin:10px 0 6px}
table.kv{border-collapse:collapse}
table.kv th{text-align:left;font-weight:500;color:var(--muted);padding:3px 18px 3px 0;vertical-align:top;font-size:13px;white-space:nowrap}
table.kv td{padding:3px 0;font-size:13px}
code{font-family:var(--mono);font-size:12px;background:var(--panel);padding:1px 5px;border-radius:4px}
footer{margin-top:50px;padding-top:20px;border-top:1px solid var(--line);color:var(--faint);font-size:12.5px;
  display:flex;flex-direction:column;gap:6px}
footer code{font-size:11px}
summary:focus-visible{outline:2px solid var(--accent);outline-offset:3px}
@media (max-width:640px){.wrap{padding:32px 16px 56px}header h1{font-size:23px}.verdict{font-size:16px}
  .barrow{grid-template-columns:120px 1fr 50px}}
"""


_GLOSSARY = [
    ("Major parent (backbone)",
     "The reference the query matches in the most windows overall."),
    ("Minor parent (donor)",
     "A reference the query matches better than the backbone over a stretch of the "
     "alignment -- a candidate recombination donor."),
    ("Similarity",
     "Per-window fraction of identical canonical bases (1.0 = identical). Windows with "
     "no comparable position are ignored."),
    ("Window / step",
     "The scan slides a fixed-width window along the alignment in fixed steps; each "
     "window is scored independently."),
    ("Support",
     "The share of distinguishing (discordant) sites -- where the query matches one "
     "candidate parent but not the other -- that favour the donor. 0.5 = no "
     "preference, 1.0 = every distinguishing site favours the donor."),
    ("q-value",
     "The sign-test p-value after Benjamini-Hochberg correction across all candidate "
     "segments (false-discovery-rate control). A region is reported when q <= alpha."),
    ("Breakpoint",
     "The query position where the source switches, with a posterior-derived "
     "uncertainty interval from the HMM."),
    ("Region calling (HMM)",
     "An HMM segments the query against the reference panel (a jump rate penalises "
     "switching reference); a segment is reported as recombinant only when its donor "
     "beats the major parent on the discordant sites by a sign test at level alpha."),
    ("Method(s) / ensemble",
     "By default several callers run and their regions are merged into one consensus. "
     "A region called by more than one method (agree) is more trustworthy and raises "
     "the confidence; PHI marks regions corroborated by the parent-free Rmin signal."),
]


# Sources for the reimplemented callers and statistics (dependency-free numpy).
_REFERENCES = [
    ("HMM segmentation (jpHMM-style)",
     "Schultz A-K, Zhang M, Leitner T, et al. (2006). A jumping profile hidden Markov "
     "model and applications to recombination sites in HIV and HCV genomes. BMC "
     "Bioinformatics 7:265."),
    ("3SEQ triplet test",
     "Boni MF, Posada D, Feldman MW (2007). An exact nonparametric method for inferring "
     "mosaic structure in sequence triplets. Genetics 176(2):1035-1047."),
    ("PHI test",
     "Bruen TC, Philippe H, Bryant D (2006). A simple and robust statistical test for "
     "detecting the presence of recombination. Genetics 172(4):2665-2681."),
    ("Four-gamete test and Rmin",
     "Hudson RR, Kaplan NL (1985). Statistical properties of the number of recombination "
     "events in the history of a sample of DNA sequences. Genetics 111(1):147-164."),
    ("Benjamini-Hochberg FDR",
     "Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical "
     "and powerful approach to multiple testing. J R Stat Soc B 57(1):289-300."),
]
