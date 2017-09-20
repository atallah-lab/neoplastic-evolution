#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# mapgen.r #NAMING

# Maps potential insertion sites based on the 'Snap-Velcro' model of Monot et al. 
# described here:

# Clement Monot, et al. (2013) "The Specificity and Flexibility of L1 Reverse 
# Transcription Priming at Imperfect T-Tracts." PLOS Genetics, 9:5.
#
# This script creates a separate R data file containing the locations of each 
# type of EN site for both sense and anti-sense strands of each chromosome.
# Locations and scores of the putative EN sites are initially stored in memory. 
#
# Usage: ./mapgen.r <start chromosome (e.g. 1)> <end chromosome (e.g. Y)>
# Output: <chromosome name>map.rda to data directory file specified on line 123.
#
# Dependencies: R(>= 2.8.0, Packages - Biostrings, BSgenome (for default hg38), 
# GenomicRanges)
# -----------------------------------------------------------------------------

#--- Load libraries
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

#--- arg[1] is starting chromosome, arg[2] is ending
args<-commandArgs(trailingOnly=TRUE)
start<-1
end<-2
options(warn=1)

#--- hg38 reference genome
cat("\nReading reference genome (GRCh38)...\n")
genome<-Hsapiens

#--- For each chromosome in hg38 range of <start> to <end>...
for (chrnm in names(genome)[start:end]) {
	start.time<-Sys.time()

	cat("\nChromosome: ",chrnm,"\n")	
	cat("Creating EN site annotation: sense...\n")

	#--- 5' to 3' target sequence 'TTTT' (derived from Monot, et. al)
	tar <- DNAString("TTTT")
	#--- Finds target starting location in chromosome, one mismatch is allowed
	mtchView <- matchPattern(tar,genome[[chrnm]],max.mismatch=1)
	#--- Calculates position weighted T-density scores
	#cat("Calculating position-weighted T-density scores...\n")
    
	#--- Start and end points of a 6 bp upstream 'T' sequence (targets)
	primrngs <- IRanges(start=start(mtchView)-6,width=6)
	targets <- DNAStringSet(genome[[chrnm]],start(mtchView),end(mtchView))
	prmrs <- DNAStringSet(genome[[chrnm]],start(primrngs),end(primrngs))
	tmp <- vmatchPattern("T",prmrs)
	tmp <- startIndex(tmp)
	#--- 0.8456 is the maximum velcro score, corresponding to 'TTTTTT' (6).
	primrnks <- lapply(tmp,function(x) sum(1/(11-x))/0.84563492) 
	rm(tmp) # cleans up

        #--- Calculates distribution of each target category. Stores indices of
        #--- sites for each (from Monot, et. al).
        #--- Example: 'ict_s' is site for 'closed' and 'tight' insertion.
        #cat("Calculating target category distribution...\n")
        ict_s<-which(targets==tar & primrnks >= 0.5)
        icl_s<-which(targets==tar & primrnks < 0.5)
        iot_s<-which(targets!=tar & primrnks >= 0.5)
        iol_s<-which(targets!=tar & primrnks < 0.5)
        insites_s<-end(mtchView)

	cat("Creating EN site annotation: anti-sense...\n")

	#--- 3' to 5' anti-sense target sequence 'AAAA' (derived from Monot, et. al)
	tar <- DNAString("AAAA")
	mtchView <- matchPattern(tar,genome[[chrnm]],max.mismatch=1)
	#--- Calculates position weighted T-density scores (See README for explanation)
	#cat("Calculating position-weighted T-density scores...\n")
	primrngs <- IRanges(start=end(mtchView)+1,width=6)
	targets <- DNAStringSet(genome[[chrnm]],start(mtchView),end(mtchView))
        prmrs <- DNAStringSet(genome[[chrnm]],start(primrngs),end(primrngs))
        tmp <- vmatchPattern("A",prmrs)
        tmp <- startIndex(tmp)
        primrnks <- lapply(tmp,function(x) sum(1/(x+4))/0.84563492) # 0.8456 is the maximum velcro score, corresponding to AAAAAA
        rm(tmp)

	#--- Calculates the distribution of each target category. Stores indices of 
	#--- sites for each category.
	#cat("Calculating target category distribution...\n")

	ict_ns<-which(targets==tar & primrnks >= 0.5)
	icl_ns<-which(targets==tar & primrnks < 0.5)
	iot_ns<-which(targets!=tar & primrnks >= 0.5)
	iol_ns<-which(targets!=tar & primrnks < 0.5)

	insites_ns<-start(mtchView)

	n <- max(length(ict_s), length(ict_ns)) #EXPLAIN
	length(ict_s) <- n                      
	length(ict_ns) <- n
        n <- max(length(icl_s), length(icl_ns))
	length(icl_s) <- n
        length(icl_ns) <- n
        n <- max(length(iot_s), length(iot_ns))
        length(iot_s) <- n
        length(iot_ns) <- n
        n <- max(length(iol_s), length(iol_ns))
        length(iol_s) <- n
        length(iol_ns) <- n	
        n <- max(length(insites_s), length(insites_ns))
        length(insites_s) <- n
        length(insites_ns) <- n

	assign(paste0(chrnm,"ict"),cbind(ict_s,ict_ns)) # 'closed' and 'tight'.
	assign(paste0(chrnm,"icl"),cbind(icl_s,icl_ns)) # 'closed' and 'loose'.
	assign(paste0(chrnm,"iot"),cbind(iot_s,iot_ns)) # 'open' and 'tight'.
	assign(paste0(chrnm,"iol"),cbind(iol_s,iol_ns)) # 'open' and 'loose'.

	assign(paste0(chrnm,"insites"),cbind(insites_s,insites_ns)) #EXPLAIN

	cat("Saving map file...\n")
	wd <- getwd()
	save(list=c(paste0(chrnm,"ict"),paste0(chrnm,"icl"),paste0(chrnm,"iot"),paste0(chrnm,"iol"),paste0(chrnm,"insites")),file=paste(wd, "/data/",chrnm,"map.rda",sep=""))
	rm(list=c(paste0(chrnm,"ict"),paste0(chrnm,"icl"),paste0(chrnm,"iot"),paste0(chrnm,"iol"),paste0(chrnm,"insites")))
	end.time <- Sys.time()
        print(end.time-start.time)
}
