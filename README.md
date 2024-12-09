# RecomFi
RecomFi (Recombination Finder) identify recombination events in a query sequence, contigs or genome, against a collection of reference sequences.

# Installation
Create RecomFi environment and install dependencies
```
conda create -n RecomFi progressivemauve biopython numpy pandas seaborn

# download and install custom Recan fork (required to identify recombination events)
git clone https://github.com/jaclew/recan.git
pip install recan/

# download and install repgenr (required to make multiple sequence alignment)
git clone https://github.com/FOI-Bioinformatics/repgenr.git
pip install repgenr/
```

# Description
RecomFi is developed to identify recombination in relatively similar datasets, such as between (sub)species of a genus or family. It generates a "pseudo-MSA (multiple sequence alignment)"  based on SNPs by using one sequence in the collection as a back-bone. This makes RecomFi fast but limits the resolution. With the "pseudo-MSA" strategy, the query may contain a fragmented genome, for example in the form of contigs, and RecomFi organize the contigs relative to the back-bone.

Recombination events are detected by sliding a window over the MSA, computing the distance between the query sequence and each of the collection of reference sequences. A recombination event is identified where a majority of the query is near reference sequence A and has a region where it is near another reference sequence B.

# Example dataset
Find example dataset of orthopoxvirus in `example_data/`. The query is a short-read assembly (x8 contigs) of a synthetic cowpox sample with a variola segment. The collection are reference-labelled orthopoxvirus sequences from `BV-BRC.org`

# Usage
Example folder structure. Query is `cowpox_with_variolaInsert.fasta.gz`
```
.
├── collection
│   ├── camelpox.fasta.gz
│   ├── cowpox.fasta.gz
│   ├── cowpox_KC813504.fasta.gz
│   ├── monkeypox.fasta.gz
│   ├── taterapox.fasta.gz
│   ├── vaccinia.fasta.gz
│   └── variola.fasta.gz
└── cowpox_with_variolaInsert.fasta.gz
```

Generate multiple sequence alignment
```
recomfi msa --query cowpox_with_variolaInsert.fasta.gz --collection collection/ --output msa.fasta

# Short notice on MSA:
#   If you have a single-contig query then you can apply flag --query_as_backbone to
#   use the query sequence as backbone in MSA instead of a reference from the collection.
```



Identify recombination events
```
# state the query name without extension
recomfi recomb --msa msa.fasta --query cowpox_with_variolaInsert --output recomfi_out
```

# Output
The terminal output show how the software ranks similarity (distance) of the query sequence to the collection. Similarity is computed in sliding windows across the MSA. For each window, the closest sequence in the collection (or sequences, if multiple ties) is selected as the "winner". The collection sequences stats are summarized:
```
Dataset distance winners in each window (no ties):
  Dataset          Windows
  --------------------------
  cowpox_KC813504  1077
  variola          793
  camelpox         23
  taterapox        2
  cowpox           1

Dataset distance winners in each window (ties allowed):
  Dataset          Windows
  --------------------------
  cowpox_KC813504  1236
  variola          949
  camelpox         176
  taterapox        160
  cowpox           155
  monkeypox        155
  vaccinia         152
```

The median similarity (distance; 1=identical, 0=no similarity) across all windows
```
Dataset distance stats:
  Dataset          Tot windows  Median distance  Distances >0  Distances >99%  Distances >95%  Distances >90%  Distances >80%  Distances >70%
  ---------------------------------------------------------------------------------------------------------------------------------------------
  cowpox_KC813504  2058         1.0              1884          1098            1862            1882            1884            1884
  variola          2058         0.98             1826          800             1697            1803            1818            1826
  taterapox        2058         0.98             1854          318             1784            1836            1844            1854
  camelpox         2058         0.98             1906          127             1812            1886            1896            1906
  vaccinia         2058         0.97             1744          122             1615            1729            1744            1744
  cowpox           2058         0.97             1867          25              1652            1858            1867            1867
  monkeypox        2058         0.97             1877          13              1532            1829            1862            1873
```

The above tables are written to output files `distances.tsv` and `distances_winners.tsv`.

Located in the output folder are two plots (and an equivalent interactive HTML-plot):

A plot is  generated for the x5 nearest datasets, showing similarity across the MSA.
![image](wiki/plot_x5.png) \
**The image show similarity in each window of the nearest five sequences to the query. Values towards 1 indicate high similarity. In the image, the query is most similar to a Cowpox sequence (light-blue line) but has a region at the middle that is similar to a Variola sequence (brown line). Concluding from the image is a putative recombination event at approximately 60-140 kbp in the MSA. Please keep in mind that these coordinates need to be mapped to the query sequence.**


The software tries to determine the one or two collection sequences most likely to recombine in the query:
![image](wiki/plot_x2.png) \
**The image show the similarity in each window of the nearest two sequences to the query. It may give a clearer view than the previous plot with x5 datasets.**
