# Example data

Two small, pre-aligned synthetic alignments that run directly with `tessera recomb`
-- no aligner needed -- and contrast the detection methods. Each is four sequences
(`query`, `parent_A`, `parent_B`, `outgroup`); the query is a `parent_A` backbone
carrying a `parent_B`-derived insert. Regenerate them with `python make_example.py`
(deterministic, standard library only).

## `divergent.msa.fasta` -- the easy case

Parents ~11 % apart with a large (1 kb) insert in a 3 kb genome. The default HMM caller
localizes the mosaic confidently:

```
tessera recomb --msa divergent.msa.fasta --query query --output out_divergent \
    --window-size 300 --window-step 30
```

It calls `parent_B` over the insert (q-value ~1e-29) with a sharp breakpoint; the
similarity plot shows an obvious crossover. `--method 3seq` agrees.

## `cryptic_insert.msa.fasta` -- why 3SEQ exists

Parents ~1 % apart with a short (800 bp) insert in a 10 kb genome -- only ~10
discriminating sites. A fixed base-pair window averages the tract away, so the HMM
segmentation finds **nothing**:

```
tessera recomb --msa cryptic_insert.msa.fasta --query query --output out_cryptic \
    --window-size 1000 --window-step 100                 # 0 regions
```

The 3SEQ caller pools the discriminating sites into an exact triplet test and recovers
the event (q-value ~1e-12):

```
tessera recomb --msa cryptic_insert.msa.fasta --query query --output out_cryptic_3seq \
    --window-size 1000 --window-step 100 --method 3seq    # finds parent_B insert
```

Both runs also report the parent-free PHI / Rmin signal in `recombination_profile.tsv`.
