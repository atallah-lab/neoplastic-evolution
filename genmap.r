#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

#### Initialize parameters
chrnm <- args[1]
cat("\nChromosome: ",chrnm,"\n")

ptm <- proc.time()
#### Read reference genome
cat("\nReading reference genome (GRCh38)...\n")
genome<-Hsapiens
chr<-genome[[chrnm]]

#### Create EN target site annotation
cat("\nCreating EN target site annotation...\n")
# Define target sequence
tar <- DNAString("TTTT")
# Find locations of target in chromosome, one mismatch is allowed
mtchView <- matchPattern(tar,chr,max.mismatch=1)
# Get start points of all targets
insites <- start(mtchView)

#### Calculate position weighted T-density scores
cat("\nCalculating position weighted T-density scores...\n")
# Create a list of start and end points of the 6 bp upstream of each target
primrngs <- IRanges(start=insites-9,width=6)
targets <- DNAStringSet(chr,start(mtchView),end(mtchView))
prmrs <- DNAStringSet(chr,start(primrngs),end(primrngs))
tmp <- vmatchPattern("T",prmrs)
tmp <- startIndex(tmp)
primrnks <- lapply(tmp,function(x) sum(1/x)/2.45) # 2.45 is the maximum primer score, corresponding to a TTTTTT primer
remove(tmp)

#### Calculate distribution of target categories
cat("\nCalculating target category distribution...\n")

# Store indices of sites of each category
ict<-which(targets==tar & primrnks >= 0.5)
icl<-which(targets==tar & primrnks < 0.5)
iot<-which(targets!=tar & primrnks >= 0.5)
iol<-which(targets!=tar & primrnks < 0.5)

remove(genome,tar,mtchView,targets,primrngs,primrnks,prmrs,chrnm)

cat("\nSaving map file...\n")
save.image("gmap.rda")
cat("\nRunning time:\n")
proc.time() - ptm
