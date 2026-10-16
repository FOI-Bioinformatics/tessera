"""The report headline names the run's reconciled major parent.

Audit 4, item B9 (suspected, confirmed while planning): with no region called, the headline
named the dataset that wins the most windows while run_provenance.json named the HMM's
major -- under lineage pooling a cluster's best-covering member, a different genome.
"""

from __future__ import annotations

import json
import random
import re
from pathlib import Path

from tessera.recomb.report_html import _summary
from tessera.recomb.run import RecombParams, run_recomb


def _panel(tmp_path: Path) -> Path:
    """A clonal query; its identical reference has a 3 kb gapped stretch, so a 99.6 %
    relative covers more of it and represents their pooled lineage in the HMM."""
    rng = random.Random(5)
    anc = [rng.choice("ACGT") for _ in range(12000)]

    def mut(seq: list[str], rate: float) -> list[str]:
        return [rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c
                for c in seq]

    a = mut(anc, 0.05)
    a_close = list(a)
    a_close[3000:6000] = ["-"] * 3000
    refs = {"A_close": a_close, "A_cov": mut(a, 0.004), "A_x": mut(a, 0.008),
            "A_y": mut(a, 0.008), "B": mut(anc, 0.05), "C": mut(anc, 0.05)}
    msa = tmp_path / "b9.fasta"
    msa.write_text(">query\n" + "".join(a) + "\n"
                   + "".join(f">{k}\n{''.join(v)}\n" for k, v in refs.items()))
    return msa


def test_headline_backbone_is_the_provenance_major_parent(tmp_path: Path, logger) -> None:
    out = tmp_path / "out"
    run_recomb(RecombParams(msa=_panel(tmp_path), output=out, query="query",
                            plot_format="png"), logger)
    major = json.loads((out / "run_provenance.json").read_text())["run"]["major parent"]
    verdict = re.search(r'class="verdict">(.*?)</p>', (out / "report.html").read_text())
    assert verdict is not None
    named = re.findall(r"<strong>(.*?)</strong>", verdict.group(1))
    assert named and named[0] == major, (named, major)


def test_summary_prefers_the_reconciled_major() -> None:
    class _Result:
        query_cumulative = [1000]

    assert _summary(_Result(), [], ["top_window_winner"], "reconciled")["major"] == "reconciled"
    assert _summary(_Result(), [], ["top_window_winner"])["major"] == "top_window_winner"
