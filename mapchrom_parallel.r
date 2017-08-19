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

if (length(args) < 2) {
	message("Usage: ./mapchrom_parallel.r <starting chromosome> <ending chromosome>")
	stop("Please provide start and stop integers as command line args.")
}

if (as.numeric(args[1]) < 1 | as.numeric(args[2]) > 23) {
	message("Usage: ./mapchrom_parallel.r <starting chromosome> <ending chromosome>")
	stop("Chromosome input number[s] out of range. Please try again.")
}

if (as.numeric(args[1]) > as.numeric(args[2])) {
	stop("Chromosome start number is greater than end number, please try again.")
}

cores <- detectCores() - 1 # Calculate number of cores, less one. 
cl <- makeCluster(cores, type="FORK", outfile="./output/multicorelog.txt") # Initiate cluster

#message(class(cl)) #SOCKcluster

# insert serial backend, otherwise error in repetetive tasks
registerDoSEQ()

### GLOBAL VARS
target <- DNAString("TTTT")

start <- as.numeric(args[1])
stop  <- as.numeric(args[2])

mtchViews <- list() # FIXME
primrnks  <- list() # FIXME

#loop in parallel
foreach(i=start:(stop),
	.export   = c("unmasked", "matchPattern"),
	.packages = c("BSgenome", "Biostrings")) %dopar% {
		 
	mtchViews[[i]] <- matchPattern(target, unmasked(Hsapiens[[i]]), max.mismatch=1)
	primrngs       <- IRanges(start=start(mtchViews[[i]])-9, width=6) # primer range
	targets        <- DNAStringSet(unmasked(Hsapiens[[i]]),start(mtchViews[[i]]),end(mtchViews[[i]]))
	prmrs          <- DNAStringSet(unmasked(Hsapiens[[i]]),start(primrngs),end(primrngs))
	indices        <- startIndex(vmatchPattern("T", prmrs))

	# 'parLapply' substitutes for lapply 
	primrnks[[i]]  <- parLapply(cl=cl, indices, function(x) sum(1/(x+4))/0.84563492) # corresponds to a TTTTTT primer
	}

stopCluster(cl)

