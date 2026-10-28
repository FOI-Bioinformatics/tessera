"""Text IO is UTF-8 whatever the locale says.

The outputs declare UTF-8 -- report.html says so in a meta tag -- but text IO without
an explicit encoding uses the locale's. Under LC_ALL=C that made the report impossible
to write (its inlined plotly bundle is not ASCII) and a non-ASCII sequence label
impossible to read or write, which is an ordinary setting on a compute node.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

LABEL = "Ebola/Guéckédou"

# The label is spelled with escapes so the script itself is ASCII: under LC_ALL=C the
# child decodes its command line with surrogateescape, so a label passed as an argument
# would arrive carrying surrogates and fail to encode for reasons of its own.
_SCRIPT = r'''
import logging, sys
from pathlib import Path
from tessera.recomb.run import RecombParams, run_recomb

label = "Ebola/Gu\u00e9ck\u00e9dou"
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
msa = out / "panel.fasta"
src = Path(sys.argv[2]).read_text(encoding="utf-8").replace(">parent_B", ">" + label)
msa.write_text(src, encoding="utf-8")
log = logging.getLogger("t")
log.addHandler(logging.NullHandler())
run_recomb(RecombParams(msa=msa, output=out / "o", query="query", window_size=300,
                        window_step=30, plot_format="png"), log)
written = out / "o"
assert (written / "report.html").stat().st_size > 0
assert label in (written / "similarity_windows.tsv").read_text(encoding="utf-8")
print("OK")
'''


def test_a_full_run_survives_a_non_utf8_locale(tmp_path: Path, example_data: Path) -> None:
    """An alignment carrying an accented label, read and written under LC_ALL=C.

    Run in a child interpreter because the locale's encoding is fixed at startup.
    PYTHONUTF8=0 and PYTHONCOERCECLOCALE=0 defeat the two mechanisms that would
    otherwise quietly give the child UTF-8 anyway.
    """
    env = {
        "LC_ALL": "C", "LANG": "C", "PYTHONUTF8": "0", "PYTHONCOERCECLOCALE": "0",
        "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
        "MPLCONFIGDIR": str(tmp_path / "mpl"),
        "TESSERA_CACHE": str(tmp_path / "cache"),
    }
    # Carried through deliberately: without it the child imports whatever Tessera is
    # installed, which is not necessarily the tree under test.
    if os.environ.get("PYTHONPATH"):
        env["PYTHONPATH"] = os.environ["PYTHONPATH"]
    result = subprocess.run(
        [sys.executable, "-c", _SCRIPT, str(tmp_path / "run"),
         str(example_data / "divergent.msa.fasta")],
        capture_output=True, text=True, env=env, timeout=600,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "OK" in result.stdout


def test_a_non_ascii_label_round_trips_through_the_region_writer(tmp_path: Path) -> None:
    """The bytes on disk must be the UTF-8 encoding of the label, not the locale's."""
    import logging

    from tessera.recomb.regions import Region
    from tessera.recomb.report_text import write_regions_tsv

    region = Region(
        minor_parent=LABEL, major_parent="parent_A", msa_start=0, msa_end=300,
        query_start=0, query_end=300, n_windows=3, mean_sim_minor=0.99,
        mean_sim_major=0.90, margin=0.09, methods=("hmm",),
    )
    write_regions_tsv([region], tmp_path, logging.getLogger("test"))
    path = tmp_path / "recombination_regions.tsv"

    assert LABEL.encode("utf-8") in path.read_bytes()
    assert LABEL in path.read_text(encoding="utf-8")
