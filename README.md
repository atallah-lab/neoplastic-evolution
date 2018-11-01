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

#### src_sompop_v0_2_4.r
* An R file containing function definitions. Functions define a stochastic simulation of neoplastic evolution under selective effects of L1 retrotransposition. A diploid system is modeled, differentiating heterozygous and homozygous mutations. Time leaps are scaled to generation time of the population.

#### src_sompop_v0_2_3.r
* The preceding version of the code described directly above. Time steps are uniform in this version.

#### src_sompop_v0_1_1.r
* An earlier version of the simulation described for src_sompop_v0_2_4.r. Driver and passenger mutations have fixed effects. Heterozygous and homozygous insertions are not differentiated.

#### run_v2_4.ipynb
* Jupyter notebook for defining parameters and running sompop_v0_2_4. Single test runs and simple visualization of results, or large batches can be performed.

#### run_v2_3.ipynb
* See description for run_v2_4.ipynb.
* Applies to sompop_v0_2_3

#### run_v1_1.ipynb
* See description for run_v2_4.ipynb.
* Applies to sompop_v0_1_1


## Data

#### ./data/tumor_type_pd_cgc.rda
* Contains R data objects which are length-3 arrays storing probabilities of driver, passenger and null-effect L1 insertions in tumor suppressor genes (TSGs) for three cancer types (lung, colon, and brain), using TSGs reported in Cancer Gene Census.

#### ./data/tumor_type_driver_lists_cgc.rda
* R data file with tables of tumor suppressor genes associated with three tissue types (lung, colon, brain) by symbol and Entrez ID.

#### ./data/exonicvsnon_counts.rda
* R data file containing the frequencies of L1 target sites of each snap-velcro category within exonic regions and outside exonic regions of the human genome. Exons were taken from Ensembl v94 human genome annotation.

#### ./data/gene_pd_exon.rda
* R data file containing the probability of exonic L1 insertion in each protein-coding gene of the human genome based on weighted L1 target site count.

#### ./data/gene_pdt_allsites.rda
* The same as gene_pdt.rda described above, though containing probabilities of insertion anywhere within the gene boundaries.
