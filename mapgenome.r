#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# mapgen.r #NAMING

# Maps potential insertion sites based on the 'Snap-Velcro' model described 
# in Clement Monot, et al. (2013) "The Specificity and Flexibility of L1 Reverse 
# Transcription Priming at Imperfect T-Tracts." PLOS Genetics, 9:5.
#
# This script creates a separate R data file for each chromosome containing the 
# locations of each type of EN site for both sense and anti-sense strands.
# Locations and scores of the putative EN sites are initially stored in memory. 
#
# Usage: ./mapgen.r <start chromosome # (e.g. 1)> <end chromosome # (e.g. 24)>
# Output: <chromosome name>map.rda to data directory file specified on line 123.
#
# Dependencies: R(>= 2.8.0, Packages - Biostrings, BSgenome (for default hg38), 
# GenomicRanges)
# -----------------------------------------------------------------------------

#--- Load libraries
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
map<-c()

#--- arg[1] is starting chromosome, arg[2] is ending
args<-commandArgs(trailingOnly=TRUE)
start<-args[1]
end<-args[2]
options(warn=1)

#--- hg38 reference genome
cat("\nReading reference genome (GRCh38)...\n")
genome<-Hsapiens

#--- For each chromosome in hg38 range of <start> to <end>...
for (chrnm in names(genome)[start:end]) {
	start.time<-Sys.time()

	### Sense strand analysis
	###------------------------------------------------------------------------###
	cat("\nChromosome: ",chrnm,"\n")	
	cat("Creating EN site annotation: sense...\n")

	#--- 5' to 3' target sequence 'TTTT' (derived from Monot, et. al)
	tar <- DNAString("TTTT")
	#--- Finds target starting location in chromosome, one mismatch is allowed
	mtchView <- matchPattern(tar,genome[[chrnm]],max.mismatch=1)

	#--- Calculates position weighted T-density scores
	primrngs <- IRanges(start=start(mtchView)-6,width=6) # start and end points of 6 bp upstream of targets
	targets <- DNAStringSet(genome[[chrnm]],start(mtchView),end(mtchView)) # collect targets from genome
	prmrs <- DNAStringSet(genome[[chrnm]],start(primrngs),end(primrngs)) # collect 'primers' - 6 bp upstream of each target
	#--- Create a list of lists containing locations of Ts in each primer
	tmp <- vmatchPattern("T",prmrs)
	tmp <- startIndex(tmp)
	#--- 0.8456 is the maximum velcro score, corresponding to 'TTTTTT' (6).
	primrnks <- lapply(tmp,function(x) sum(1/(11-x))/0.84563492) 
	rm(tmp) # cleans up

        #--- Stores indices of sites for each target category (from Monot, et. al).
        #--- Example: 'ict_s' is site for 'closed' and 'tight' insertion.
        ict_s<-which(targets==tar & primrnks >= 0.5)
        icl_s<-which(targets==tar & primrnks < 0.5)
        iot_s<-which(targets!=tar & primrnks >= 0.5)
        iol_s<-which(targets!=tar & primrnks < 0.5)
        insites_s<-end(mtchView)

	### Anti-sense strand analysis
	###------------------------------------------------------------------------###
	cat("Creating EN site annotation: anti-sense...\n")

	#--- 3' to 5' anti-sense target sequence 'AAAA' (derived from Monot, et. al)
	tar <- DNAString("AAAA")
	mtchView <- matchPattern(tar,genome[[chrnm]],max.mismatch=1)
	primrngs <- IRanges(start=end(mtchView)+1,width=6)
	targets <- DNAStringSet(genome[[chrnm]],start(mtchView),end(mtchView))
        prmrs <- DNAStringSet(genome[[chrnm]],start(primrngs),end(primrngs))
        tmp <- vmatchPattern("A",prmrs)
        tmp <- startIndex(tmp)
        primrnks <- lapply(tmp,function(x) sum(1/(x+4))/0.84563492)
        rm(tmp)

	ict_ns<-which(targets==tar & primrnks >= 0.5)
	icl_ns<-which(targets==tar & primrnks < 0.5)
	iot_ns<-which(targets!=tar & primrnks >= 0.5)
	iol_ns<-which(targets!=tar & primrnks < 0.5)
	insites_ns<-start(mtchView)

	### Concatenate data from both strands
	###-----------------------------------------------------------------------##

	#--- Make the arrays for each set of indices equal length, padding the shorter with NAs
	n <- max(length(ict_s), length(ict_ns))
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
	#--- Bind the arrays for both strands as two columns
	ict<-cbind(ict_s,ict_ns)
	icl<-cbind(icl_s,icl_ns)
	iot<-cbind(iot_s,iot_ns)
	iol<-cbind(iol_s,iol_ns)
	insites<-cbind(insites_s,insites_ns)
	#--- Store data in single S4 object
	map$insites <-insites
	map$ict<-ict
	map$icl<-icl
	map$iot<-iot
	map$iol<-iol

	#--- Save map under name of current chromosome
	assign(paste0(chrnm,"Map"),map) # 'closed' and 'tight'.

	cat("Saving map file...\n")
	wd <- getwd()
	save(list=(paste0(chrnm,"Map")),file=paste(wd, "/data/root_maps/",chrnm,".rda",sep=""))
	rm(list=c("insites","ict","icl","iot","iol"))
	end.time <- Sys.time()
        print(end.time-start.time)
}
