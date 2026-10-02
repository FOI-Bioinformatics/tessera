# Aligner backends

Tessera delegates alignment to a pluggable backend, discovered through entry points
so additional backends can be added without changing the core. All backends produce
a reference-anchored alignment, which is what the recombination scan assumes. Choose
with `--aligner` and tune with repeatable `--aligner-arg key=value`.

| Backend | Best for | Notes |
|---|---|---|
| `sibeliaz` (default) | Moderately divergent genomes, including rearrangements | Installs cleanly via conda; `kmer`, `abundance`, `bubble`, `filtermemory` |
| `mafft` | Similar, largely collinear genomes | True base-level alignment, the canonical input for the window method; adds a fragmented query with `--addfragments` and reorients reverse-strand genomes or contigs (`--adjustdirection`). `maxiterate`, `retree`, `op`, `ep`, `sixmerpair` |
| `minimap2` | Speed and assembly/contig queries | Fast assembly-to-reference projection; `preset` (default `asm20`, e.g. `asm10` for closer genomes) |
| `progressivemauve` | Genomes with large rearrangements/inversions | Tolerant but slow, heavy, and not available as a conda build on all platforms; `seed_weight`, `single` |
| `cactus` | Same-species pangenomes | Resource heavy (Toil/containers) |

`sibeliaz` is the default: it installs cleanly across platforms and, on the example
data, reproduces `progressivemauve`'s recombination coordinates. For very similar,
collinear genomes `mafft` gives the most faithful base-level signal and `minimap2`
the fastest run (and the best fit for a fragmented query); `progressivemauve` remains
an option for genomes with large rearrangements. Every backend reports the alignment in
backbone coordinates, so material inserted relative to the backbone is dropped -- `mafft`
included, which runs with `--keeplength`.

Memory with `mafft`: by default each genome is added with a pairwise local alignment
against the backbone, whose memory grows with the product of the two lengths (1.6-6 GB per
genome was measured on 8.5-30 kb pairs). When the backbone or any genome to add is longer
than 20 kb Tessera therefore passes `--6merpair` (k-mer distances), which keeps memory
small. The two modes were compared on one genome pair each: SARS-CoV-2 (30 kb, the only
one above the threshold), measles and Ebola virus gave the same alignment, and a divergent
HIV-1 pair differed in 0.3% of columns; no comparison exists at mpox size, where the
pairwise mode did not complete. `--aligner-arg sixmerpair=true` or `sixmerpair=false`
overrides the length rule. Below 20 kb, genomes are added concurrently only as far as an
estimate of the memory needed fits half of the physical memory, and `--threads` is shared
between the concurrent jobs rather than given to each. The estimate uses the host's
physical memory; under a container or scheduler memory limit lower `--threads` or pass
`sixmerpair=true`.

A multi-contig backbone is laid out as the concatenation of its contigs in the order of
its FASTA file (`sibeliaz`, `mafft`, `minimap2`, `progressivemauve`); a contig that nothing
aligns to keeps its columns, as gaps. The `cactus` backend lays contigs out in sorted-name
order. A genome the aligner cannot place against the backbone at all is kept as an all-gap
row and named in a warning, so the alignment always has one row per input genome.

Genomes are staged before alignment. A FASTA with whitespace inside its sequence lines
(trailing spaces, tabs, Windows line endings) is aligned from a cleaned copy and reported
in a warning; the input file itself is not modified.

Examples:

```
tessera msa ... --aligner minimap2 --aligner-arg preset=asm10
tessera msa ... --aligner mafft --aligner-arg maxiterate=1000
tessera msa ... --aligner sibeliaz --aligner-arg kmer=15
tessera msa ... --aligner progressivemauve --aligner-arg seed_weight=11
```

Add an optional backend to a conda environment with
`conda install -c bioconda sibeliaz` (or `cactus`). The provided `environment.yml`
installs Python, the backends, and Tessera in one step.
