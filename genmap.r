#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

cat("\nGenerating insertion map...\n")

#### Initialize parameters
ENifrc <- as.numeric(args[1])
chrnm <- args[2]
cat("\nChromosome: ",chrnm,"\n")
cat("ENi insertion fraction: ",ENifrc,"\n")

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
ttight <- targets[primrnks>=0.5]
tloose <- targets[primrnks<0.5]
lenl <- length(tloose)
lenti <- length(ttight)
lenta <- length(targets)
ccl <-sum(vcountPattern("TTTT",tloose,max.mismatch=0))
cct <-sum(vcountPattern("TTTT",ttight,max.mismatch=0))
col <- lenl-ccl
cot <- lenti-cct
pd <- c(11.55*cct,7.25*ccl,1.95*cot,1*col)
pd <- (pd/sum(pd))*(1-ENifrc)

#### Store indices of sites of each category
ict<-which(targets==tar & primrnks >= 0.5)
icl<-which(targets==tar & primrnks < 0.5)
iot<-which(targets!=tar & primrnks >= 0.5)
iol<-which(targets!=tar & primrnks < 0.5)

remove(genome,chr,tar,mtchView,targets,primrngs,primrnks,ttight,tloose,lenl,lenti,lenta,ccl,cct,col,cot)

cat("\nSaving map file...")
save.image("gmap.rda")