#!/usr/bin/env Rscript

  # Computes target distributions for n chromosomes, in parallel.
  #
  # Command line args:
  #   args[i..n] list of chromosomes.  
  #
  # Returns:
  #   Separate map.rda files per chromosome.

#source("https://bioconductor.org/biocLite.R")
#require(BSgenome.Hsapiens.UCSC.hg38)

library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

# parallel packages
library(parallel)
library(foreach)
library(doParallel)

#devtools::install_github("collectivemedia/tictoc")
#devtools::install_github("eddelbuettel/rbenchmark")
#library(tictoc)
#library(rbenchmark)

args = commandArgs(trailingOnly=TRUE)

#### Read reference genome
cat("Reading reference genome (GRCh38)...\n")

cores <- detectCores() - 1 # Calculate number of cores, less one. 
cl <- makeCluster(cores, type="FORK") # Initiate cluster

registerDoParallel(cl)

### GLOBAL VARS
target <- DNAString("TTTT")
stop <- 22
chromosomes <- c('chr1', 
		'chr2', 
		'chr3', 
		'chr4', 
		'chr5', 
		'chr6', 
		'chr7', 
		'chr8',
		'chr9',  
		'chr10', 
		'chr11', 
		'chr12', 
		'chr13',
		'chr14',  
		'chr15', 
		'chr16',
		'chr17', 
		'chr18', 
		'chr19', 
		'chr20',
		'chr21',  
		'chr22', 
		'chrX',
		'chrY'
				)

mtchViews <- list() # FIXME
primrnks <- list() # FIXME

#loop in parallel
foreach(i=1:(length(chromosomes)-stop), 
		.export = c("unmasked", "matchPattern"),
		.packages = c("BSgenome", "Biostrings")) %dopar% {
		 mtchViews[[i]]  <- matchPattern(target, unmasked(Hsapiens[[i]]), max.mismatch=1)
		 primrngs        <- IRanges(start=start(mtchViews[[i]])-9, width=6) # primer range
		 targets         <- DNAStringSet(unmasked(Hsapiens[[i]]),start(mtchViews[[i]]),end(mtchViews[[i]]))
		 prmrs           <- DNAStringSet(unmasked(Hsapiens[[i]]),start(primrngs),end(primrngs))
		 tmp             <- vmatchPattern("T",prmrs)
		 tmp             <- startIndex(tmp)

	primrnks[[i]] <- parLapply(tmp,function(x) sum(1/(x+4))/0.84563492) # corresponds to a TTTTTT primer
	remove(tmp)
}

stopCluster(cl)

