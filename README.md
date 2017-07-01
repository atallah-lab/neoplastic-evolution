# sim-develop
Developments for simulating transposable element activity in genome evolution.

#### gcbasedins.r
* Insert L1 element copies in a given chromosome with probability based on local G/C content.
* Dependencies: R (>= 2.8.0), Bioconductor packages - Biostrings, BSgenome (for default hg19)
```
./gcbasedins.r <optional: input genome (fasta, default - hg19)> <output file name (fasta)>
```

