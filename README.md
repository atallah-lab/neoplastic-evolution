# sim-develop
Developments for simulating transposable element activity in genome evolution.

#### sim-main.r
* Main function for simulating transposition of L1 elements based on Snap-Velcro model (Monot et al. 2013), with truncation modeled after Lee et al. 2012.
* If the file gmap.rda does not exist, this function will call genmap.r to create a map file and require arguments as detailed for genmap.r
* Dependencies: R(>= 2.8.0, Packages - Biostrings)
* *Monot et al., 2013, The Specificity and Flexibility of L1 Reverse Transcription Priming at Imperfect T-Tracts*
* *Lee et al., 2012, Landscape of Somatic Retrotransposition in Human Cancers*

#### genmap.r
* Maps potential insertion sites based on Snap-Velcro model, and generates an insertion probability distribution with respect to the site categories (closed-tight, closed-loose, open-tight, open-loose)
* Dependencies: R(>= 2.8.0, Packages - Biostrings, BSgenome (for default hg38))
```
./genmap.r <fraction of endonuclease-independent (ENi) insertions> <chromosome name (e.g. chrX)>
```

#### exchk.r
* Example of function for checking whether a given breakpoint lies inside an exon
* Dependencies: R(>= 3.0.0, Packages - data.table), filtgff3.bash
```
./exchk.r <input gff3> <chromosome name (e.g. X)> <breakpoint location>
```

#### gcbasedins.r
* Insert L1 element copies in a given chromosome with probability based on local G/C content.
* Dependencies: R(>= 2.8.0, Packages - Biostrings, BSgenome (for default hg38))
```
./gcbasedins.r <optional: input genome (fasta, default - hg19)> <output file name (fasta)>
```
#### filtgff3.bash
* Takes in a .gff3 file and outputs a tab-delimited file with 3 columns: chromosome name, exon start loc, exon end loc

#### map_L1_consensus_sequences_hg38.R
* Map an arbitrary set of L1 consensus sequences against hg38. 

#### trpd.dat
* Contains tab-delimited values approximating the distribution of truncations as shown in Lee et al. 2012 Fig 3 (A)

#### ./data/GenBank_L19092.fa 
* consensus sequence retrieved from GenBank, https://www.ncbi.nlm.nih.gov/nuccore/L19092.1?report=GenBank

#### ./data/Brouha_consensus_L1.fa 
* Single L1 consensus sequence taken from Brouha, et. al. paper. 

#### ./data/gb_AC002980.fa 
* Set of 5 L1 consensus sequences retrieved from NCBI BLAST search against hg38 for Accession No. AC002980. (07/17/2017). This was the top hit noted in Table 4 of the Brouha paper. 

#### ./data/hsflil1_8438.bed
* BED file retrieved from L1Base, http://l1base.charite.de/BED/hsflil1_8438.bed