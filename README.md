# sim-develop
An R project developed to simulate transposable element activity in genome evolution.
### Background 
---
##### What is the 'Snap-Velcro' Model?
Prior studies have demonstrated that most L1 insertions occur in sequences related to the L1 EN consensus sequence (degenerate 5′-TTTT/A- 3′ sites). The *Snap-Velcro* model describes four possible types of insertion site guided by two parameters. The *snap* parameter represents the last four nucleotides of the primer and is considered a perfect terminal match and 'closed' if it ends with four T’s. The snap is considered 'open' if it contains a mismatch in the last four T’s. The *velcro* region of the primer (6 upstream basepairs) can be classified as being either *tightly-fastened* or *loosely-fastened* based on its position-weighted T-density score. (Monot et al., 2013)

##### What is a position weighted T-density score?
The position-weighted T-density score is the sum of thymine (T) nucleotides where each 'T' is weighted in terms of the inverse proportionality of its distance from the 3’ end of the primer. The *maximum* velcro score corresponds to six consecutive T’s and is:

(1/5) + (1/6) + ... + (1/10) = **0.84563492** 

Each score for a velcro region is divided by this maximum score as a normalization. If the position-weighted T-density score of the velcro region is ≥0.5 then that region is considered tightly-fastened, otherwise it is loosely-fastened.

##### Relevance of 'Snap-Velcro' Model
Monot et al. studed the relative enrichment of L1 insertion frequencies amongst the four snap-velcro categories in humans. They recorded the number of insertions observed at each type of EN site, and divided this by the number of site-types in the reference genome (hg19 in their study). Amongst two datasets of human L1 insertions (Solyom, 2012; Lee, 2012) they found that the average relative enrichment for the closed/tight, closed/loose, open/tight, and open/loose categories is, respectively 11.55, 7.25, 1.95, 1.00. 

![Snap Velcro Model](https://github.com/atallah-lab/sim-develop/blob/comments/data/images/snap-velcro.png)


### Code
---
#### gen-sim.r
* Main function for simulating retrotransposition of L1 elements in hg38
* Insertion site probabilities based on Snap-Velcro model (Monot et al., 2013)
* Truncation modeled after Lee et al. (2012)
* Transduction modeled after Goodier et al. (2000)
* Dependencies: R(>= 2.8.0, Packages - Biostrings, BSgenome (for default hg38), GenomicRanges)
* *Monot et al., 2013, The Specificity and Flexibility of L1 Reverse Transcription Priming at Imperfect T-Tracts*
* *Lee et al., 2012, Landscape of Somatic Retrotransposition in Human Cancers*
* *Goodier et al., 2000, Transduction of 3′-flanking sequences is common in L1 retrotransposition*
```
./gen-sim.r <copy number>
```

#### chrom-sim.r
* Generates insertion sites in a single chromosome based on Snap-Velcro model
* Dependencies: R(>= 2.8.0, Packages - Biostrings, BSgenome (for default hg38))
```
./chrom-sim.r <Copy Number> <Chromosome Name (e.g. chr1)>
```

#### mapchrom_parallel.r
* Experiment in parallelization applied to mapchrom.r 
* Dependencies: R(>= 2.8.0, Packages - Biostrings, BSgenome (for default hg38), parallel, foreach, doParallel)
```
./mapchrom_parallel.r <chromosome start number (e.g. 1-23), <chromosome stop number (e.g. 2-24)>
```
#### mapgen.r
* Maps potential insertion sites. Creates an .rda file for each chromosome in a provided range.
```
./mapgen.r <start chromosome (e.g. 1)> <end chromosome (e.g. Y)>
```

#### get_sv_dist.r
* Stores the counts of each EN site category for each chromosome in a 24 x 4 array.

#### hgextract.r
* Extracts segments from a genome in ranges provided in a tab-delimited text file, and outputs the results to a Fasta file.

#### makeL1RankTable.r
* Creates a .rda file containing the loci of retrotransposition-competent L1's in hg38 and their rankings.

#### transduct_truncate.r
* Simulates truncation and transduction of L1 elements selected for retrotransposition based on ranking.

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
#### gff3_exonfilt.bash
* Takes in a .gff3 file and outputs a tab-delimited file with 3 columns: chromosome name, exon start loc, exon end loc

#### map_L1_consensus_sequences_hg38.R
* Map an arbitrary set of L1 consensus sequences against hg38. 

#### ./data/GenBank_L19092.fa 
* consensus sequence retrieved from GenBank, https://www.ncbi.nlm.nih.gov/nuccore/L19092.1?report=GenBank

#### ./data/Brouha_consensus_L1.fa 
* Single L1 consensus sequence taken from Brouha et al. paper. 

#### ./data/gb_AC002980.fa 
* Set of 5 L1 consensus sequences retrieved from NCBI BLAST search against hg38 for Accession No. AC002980. (07/17/2017). This was the top hit noted in Table 4 of the Brouha paper. 

#### ./data/hsflil1_8438.bed
* BED file retrieved from L1Base, http://l1base.charite.de/BED/hsflil1_8438.bed

#### ./data/L1PA
* directory to hold L1PA annotation output.

#### ./data/L1PA/by_chr
* all L1PA annotation output, by chromosome.

#### ./data/L1PA/L1PA.fa
* L1PA2, L1PA3, L1PA4 sequences.

#### ./data/L1PA/L1PA2.bed, L1PA3.bed, L1PA4.bed
* L1PA2|3|4 annotation output, by sequence, over entire genome.