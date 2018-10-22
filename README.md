# sim-develop
An R project developed to simulate transposable element activity in genome evolution.
### Background 
---
#### Introduction
 As reported by Monot (Monot, C. et al. 2013) L1 retrotransposons utilize an endonuclease (EN) to initiate retrotransposition by first nicking the genomic DNA. The DNA end that results is then used to prime transcription within the L1 RNA poly(A) tail. This process is known as __target-primed reverse transcription (TPRT)__. Efficient priming is detected with only 4 matching nucleotides at the primer's 3' end. Furthermore, the L1 ribonucleoprotein particle (RNP) tolerates terminal mismatches if they are compensated for within the final 10 bases of the primer by an increase in the number of matching nucleotides.

#### What's the 'Snap-Velcro' Model?
Prior studies have demonstrated that most L1 insertions occur in sequences related to the L1 EN consensus sequence (degenerate 5′-TTTT/A- 3′ sites). The *Snap-Velcro* model describes four possible insertion sites that are modeled by two key parameters. The *snap* parameter represents the last four nucleotides of the primer. It is considered to be a perfect terminal match and 'closed' if it ends with four T’s. Alternatively, the snap is 'open' if it contains a mismatch in the last four T’s. The *velcro* region of the primer (6 upstream basepairs) is classified as either *tightly-fastened* or *loosely-fastened* based on its position-weighted T-density score (Monot et al., 2013). Thus, the distribution of new L1 insertions within the human genome is the result, in part, of the specificity of the initial EN cleavage site of L1 reverse transcription initiation.

#### What's a position weighted T-density score?
A T-density is calculated by dividing the number of Ts in the oligonucleotide by the length of the oligonucleotide. The position-weighted T-density gives more weight to Ts which are closer to the 3' end of the primer. The weight is inversely proportional to the distance from the 3' end.

LOU519:
GAGCCAGGAGGAATACTTTT
1 + (1/2) + (1/3) + (1/4) + (1/7) = **2.23**

LOU541:
TTTTTTTTTTTTTTTTTTTT
1 + (1/2) + (1/3) + ... + (1/19) + (1/20) = **3.60**

The position-weighted T-density is calculated by dividing the position-weighted T-count of a primer to the maximum position-weighted T-count. Thus the position-weighted T-density of LOU519 is equal to 0.62 (2.23/3.60) and the position-weighted T-density of LOU541 is equal to 1 (3.60/3.60). 

The *maximum* velcro score corresponds to six consecutive T’s and is:

GAGCCAGCCATTTTTTTTTT
(1/5) + (1/6) + ... + (1/10) = **0.84563492** 

Each score for a velcro region is divided by this maximum score as a normalization. If the position-weighted T-density score of the velcro region is ≥0.5 then that region is considered tightly-fastened, otherwise it is loosely-fastened.

#### Relevance of 'Snap-Velcro' Model in the Simulation of Transposable Elements
Monot et al. studed the relative enrichment of L1 insertion frequencies amongst the four snap-velcro categories in humans. They recorded the number of insertions observed at each type of EN site and divided this by the number of site-types in the reference genome (hg19 in their study). Amongst two datasets of human L1 insertions (Solyom, 2012; Lee, 2012) they found that the average relative enrichment for the closed/tight, closed/loose, open/tight, and open/loose categories was, respectively 11.55, 7.25, 1.95, 1.00. 


## Source Code

#### somlin-L1.ipynb
* Notebook for simulating evolution of a **SOM**atic cell **LIN**eage, incorporating **L1** retrotransposition events and response to mutations

#### sompop-L1-v1.0.ipynb
* A version of the simulation implemented in somlin-L1.ipynb that assumes an initially stationary **SOM**atic cell **POP**ulation, incorporates cell death, and focuses on stochastic simulation of population dynamics.

#### sompop-L1-v2.4.ipynb
* A version of sompop-L1.ipynb that uses a two-hit model of driver and passenger mutations, and tracks disrupted genes for each clone of the population. Version 2.4 allows for multiple mutations per cell per time step.

#### sompop-L1-v2.3.ipynb
* A version of sompop-L1.ipynb that uses a two-hit model of driver and passenger mutations, and tracks disrupted genes for each clone of the population.

#### gen-sim.r
* Script for simulating retrotransposition of L1 elements in a single genome
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
#### gen-sim-notebook.ipynb
* This is a highly commented, notebook version of the site-generation part of gen-sim.r meant for learning and
developing the genome-level retrotransposition model.
* Dependencies: R(>= 2.8.0, Packages - Biostrings, BSgenome (for default hg38), GenomicRanges), Jupyter, IRKernel

#### sim_exon_annotation.ipynb
* A notebook for creating a table of exon loci in hg38 as well as their parent genes. This table is stored in the data file ./data/exann.rda

#### mapchrom_parallel.r
* Experiment in parallelization applied to mapchrom.r 
* Dependencies: R(>= 2.8.0, Packages - Biostrings, BSgenome (for default hg38), parallel, foreach, doParallel)
```
./mapchrom_parallel.r <chromosome start number (e.g. 1-23), <chromosome stop number (e.g. 2-24)>
```
#### mapgenome.r
* Maps potential L1 insertion sites based on the Snap-Velcro model (Monot et al.). Creates an .rda file for each chromosome in a provided range.
```
./mapgen.r <start chromosome # (e.g. 1)> <end chromosome # (e.g. 24)>
```
#### mapsequence.r
* Function to map potential L1 insertion sites within a provided DNAString object

#### get_chrom_sv_dists.r
* Stores the counts of each EN site category for each chromosome in a 24 x 4 matrix, as well as their probabilities.

#### save_chrom_sv_dists.ipynb
* The notebook used to generate the file ./data/chromSitesPd.rda, which contains the probability of insertion distributions for the various snap-velcro categories for each chromosome of hg38. This is similar to the above script, except it adds endonuclease-independent insertions as a class and accounts for their probability (assumed 0.1 in the script).

#### hgextract.r
* Extracts segments from a genome in ranges provided in a tab-delimited text file, and outputs the results to a Fasta file.

#### makeL1RankTable.r
* Creates a .rda file containing the loci of retrotransposition-competent L1's in hg38 and their rankings.

#### process_L1s.r
* Function to simulate truncation and transduction of L1 elements selected for retrotransposition based on activity ranking.

#### gcbasedins.r
* Insert L1 element copies in a given chromosome with probability based on local G/C content.
* Dependencies: R(>= 2.8.0, Packages - Biostrings, BSgenome (for default hg38))
```
./gcbasedins.r <optional: input genome (fasta, default - hg19)> <output file name (fasta)>
```
#### map_L1_consensus_sequences_hg38.R
* Map an arbitrary set of L1 consensus sequences against hg38. 

## Data

#### ./data/cancer_type_pd_cgc.rda
* Contains R data objects which are length-3 arrays storing probabilities of driver, passenger and null-effect L1 insertions in tumor suppressor genes (TSGs) for three cancer types (lung, colon, and brain), using TSGs reported in Cancer Gene Census.

#### ./data/human_tsgs_cgc.txt
* List of human tumor suppressor genes reported by Cancer Gene Census

### ./data/human_tsgs_tsgene.txt
* List of human tumor suppressor genes reported by TSGene 2.0 database.

#### ./data/root_maps
* To run somlin-L1.ipynb, Make this file a soft link to a directory for storing the chromosome map files (~ 3GB) generated by mapgenome.r

#### ./data/exann.rda
* R data file containing the loci of exons in hg38 as well as the symbols of their parent genes (generated using exon_sim_annotation.ipynb)

#### ./data/exonicvsnon_counts.rda
* R data file containing the frequencies of L1 target sites of each snap-velcro category within exonic regions and outside exonic regions of the human genome. Exons were taken from Ensembl v86 human genome annotation.

#### ./data/genevsnon_counts.rda
* R data file containing the frequencies of L1 target sites of each snap-velcro category within gene boundaries and outside genes of the human genome. Gene boundaries were taken from Ensembl v86 human genome annotation.

#### ./data/gene_pdt.rda
* R data file containing the probability of exonic L1 insertion in each gene of the human genome based on weighted L1 target site count.

#### ./data/gene_pdt_allsites.rda
* The same as gene_pdt.rda described above, though containing probabilities of insertion anywhere within the gene boundaries.

#### ./data/L1RankTable.rda
* R data file containing the loci of intact L1 elements of hg38, as well as their computed activity scores (using the script makeL1RankTable.rda)

#### ./data/L1transdpd.csv
* Probability distribution for the transduction lengths of L1 elements (Goodier et al. (2000))

#### ./data/L1truncpd.csv
* Probability distribution for the truncation fractions of L1 elements (Lee et al. (2012))

#### ./data/chrmpd.rda
* R data file containing two objects: chrmcnt - number of sites for each of the 4 snap-velcro categories for each chromosome, chrmpd - probability of insertion for each of the 4 categories for each chromosome

#### ./data/chromSitePd.rda
* R data file containing a 24x5 matrix with probabilities of selection for each insertion site class, including 0.1 probability of endonuclease-independent insertions

#### ./data/GenBank_L19092.fa 
* consensus sequence retrieved from GenBank, https://www.ncbi.nlm.nih.gov/nuccore/L19092.1?report=GenBank

#### ./data/Brouha_consensus_L1.fa 
* Hot L1 consensus sequence taken from Brouha et al. paper. 

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
