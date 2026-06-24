# RecomFi
RecomFi (Recombination Finder) identifies recombination events in a query sequence, contigs or genome, against a collection of reference sequences.

# Description
RecomFi is developed to identify recombination in relatively similar datasets, such as between (sub)species of a genus or family. It generates a reference-anchored "pseudo-MSA (multiple sequence alignment)" by using one sequence as a backbone. This makes RecomFi fast but limits the resolution. With the pseudo-MSA strategy the query may be a fragmented genome, for example a set of contigs, which RecomFi organizes relative to the backbone.

Recombination events are detected by sliding a window over the MSA and computing, in each window, the similarity between the query and each reference. A recombination event is indicated where the query is near reference A across most of its length but is near another reference B over a region.

RecomFi is organized as a small, modular Python package:
- Alignment is delegated to a pluggable **aligner backend** (`progressivemauve`, `sibeliaz`, or `cactus`), discovered through entry points so additional backends can be added without changing the core.
- The recombination scan uses a built-in sliding-window distance engine (no external alignment-analysis dependency).

# Installation
RecomFi needs Python (>= 3.11) and at least one aligner backend. The aligner binaries are most easily installed with conda.

```
# create an environment with an aligner backend, then install RecomFi
conda create -n recomfi -c conda-forge -c bioconda python">=3.11" mauve "boost-cpp=1.74.0"
conda activate recomfi
pip install .
```

Optional aligner backends can be added to the same environment:
```
conda install -c bioconda sibeliaz   # or: cactus
```

An `environment.yml` is provided that installs Python, all three backends, and RecomFi in one step:
```
conda env create -f environment.yml
conda activate recomfi
```

# Aligner backends
`progressivemauve` is the default and recommended backend: it produces a reference-anchored alignment, which is what the recombination scan assumes. `sibeliaz` and `cactus` are available as alternatives for cases where they fit the data better, but their output should be reviewed for the recombination use case (their coordinate anchoring differs).

# Example dataset
Find an example dataset of orthopoxvirus in `example_data/`. The query is a short-read assembly (8 contigs) of a synthetic cowpox sample with a variola segment. The collection are reference-labelled orthopoxvirus sequences from `BV-BRC.org`.

Example folder structure (query is `cowpox_with_variolaInsert.fasta.gz`):
```
.
├── collection
│   ├── camelpox.fasta.gz
│   ├── cowpox.fasta.gz
│   ├── cowpox_KC813504.fasta.gz
│   ├── monkeypox.fasta.gz
│   ├── taterapox.fasta.gz
│   ├── vaccinia.fasta.gz
│   └── variola.fasta.gz
└── cowpox_with_variolaInsert.fasta.gz
```

# Usage
Generate a multiple sequence alignment:
```
recomfi msa --query cowpox_with_variolaInsert.fasta.gz --collection collection/ --output msa.fasta

# Choose an aligner backend (default: progressivemauve) and pass tuning options:
#   recomfi msa ... --aligner sibeliaz --aligner-arg kmer=15
#   recomfi msa ... --aligner progressivemauve --aligner-arg seed_weight=11
#
# If you have a single-contig query you can use it as the backbone instead of a
# reference from the collection:
#   recomfi msa ... --query-as-backbone
```

Identify recombination events (state the query label as it appears in the MSA, i.e. the query file name without extension):
```
recomfi recomb --msa msa.fasta --query cowpox_with_variolaInsert --output recomfi_out

# The window, step, metric, number of plotted datasets and plot format are
# configurable:
#   recomfi recomb ... --window-size 1000 --window-step 100 --top-n 5 --plot-format png
#
# Region calling can be tuned (defaults derive from the window size):
#   recomfi recomb ... --min-region 1000 --margin 0.0 --merge-gap 1000
```

Run `recomfi --help`, `recomfi msa --help` or `recomfi recomb --help` for the full set of options.

# Output
RecomFi computes, in sliding windows across the MSA, the similarity of the query to each reference (1 = identical, 0 = no similarity). For each window the closest reference (or references, on a tie) is the "winner". The reference winning the most windows is the **major parent** (the backbone donor); a stretch where the query is instead closest to another reference (a **minor parent**) is reported as a recombinant region.

The terminal output (also written to the run log) summarises the window winners, the per-dataset similarity statistics, and the called regions, for example:
```
Recombination regions (major parent: cowpox_KC813504):
  Minor parent  Major parent     MSA start  MSA end  Query start  Query end  Length(bp)  Windows  Sim minor  Sim major
  ---------------------------------------------------------------------------------------------------------------------
  variola       cowpox_KC813504  60500      141000   58800        138900     80500       790      0.992      0.951
```

The **region calling is a transparent heuristic screen, not a statistical significance test** (such as 3SEQ or RDP). Treat the regions as candidates to inspect, not confirmed events.

Output files in the chosen directory:

| File | Contents |
|---|---|
| `recombination_regions.tsv` | Called regions: minor/major parent, start/end in **both MSA columns and query bases**, length, support, mean similarities |
| `similarity_windows.tsv` | Full per-window matrix: `msa_position`, `query_position`, `winner`, and one similarity column per dataset |
| `similarity_stats.tsv` | Per-dataset similarity statistics (median, windows above identity thresholds) |
| `window_winners.tsv` | Per-dataset count of windows won (ties included) |
| `similarity_top{N}.{fmt}` | Static plot of the nearest `--top-n` datasets, called regions shaded |
| `similarity_pair.{fmt}` | Static plot of the major vs leading minor parent, region shaded |
| `report.html` | Self-contained report: run provenance, the region table, the per-dataset stats, and an embedded interactive plot |

Similarity is computed only over columns where both sequences carry a canonical base (A/C/G/T); gaps, `N` and IUPAC ambiguity codes are ignored, so an `N` never counts as a match. A window with no comparable position — for example an inter-contig gap in a fragmented query — is uninformative and reported as `NA` in `similarity_windows.tsv` (and excluded from the winners, statistics and region calling).

Query coordinates are reported alongside MSA coordinates, so regions no longer need to be mapped back to the query by hand.

A plot is generated for the nearest datasets (`--top-n`, default 5), showing similarity across the MSA with called regions shaded.
![image](wiki/plot_x5.png) \
**Similarity in each window of the nearest five sequences to the query. Values towards 1 indicate high similarity. Here the query is most similar to a Cowpox sequence but has a region in the middle similar to a Variola sequence — a putative recombination event, called automatically and reported in `recombination_regions.tsv` in both MSA and query coordinates.**

The pairwise plot shows the major parent against the leading minor parent:
![image](wiki/plot_x2.png) \
**The two sequences most likely involved in the recombination, with the called region shaded.**

# Development
```
pip install -e ".[dev]"
ruff check src tests
pytest                       # add -m "not requires_binary" to skip aligner-dependent tests
```
