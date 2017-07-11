# sim-develop
Developments for simulating transposable element activity in genome evolution.

#### gcbasedins.r
* Insert L1 element copies in a given chromosome with probability based on local G/C content.
* Dependencies: R(>= 2.8.0, Packages - Biostrings, BSgenome (for default hg19))
```
./gcbasedins.r <optional: input genome (fasta, default - hg19)> <output file name (fasta)>
```
#### filtgff3.bash
* Takes in a .gff3 file and outputs a tab-delimited file with 3 columns: chromosome name, exon start loc, exon end loc

#### exchk.r
* Example of function for checking whether a given breakpoint lies inside an exon
* Dependencies: R(>= 3.0.0, Packages - data.table), filtgff3.bash
```
./exchk.r <input gff3> <chromosome name (no quotes needed)> <breakpoint location>
```


