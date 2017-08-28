#!/usr/bin/env Rscript

library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
args<-commandArgs(trailingOnly=TRUE)
options(warn=1)

#### Read reference genome
cat("\nReading reference genome (GRCh38)...\n")
genome<-Hsapiens

for (chrnm in names(genome)[args[1]:args[2]]) {
	start.time<-Sys.time()
	cat("\nChromosome: ",chrnm,"\n")
	#### Create EN target site annotation
	
	cat("Creating EN site annotation: sense...\n")

	# Define target sequence
	tar <- DNAString("TTTT")
	# Find locations of target in chromosome, one mismatch is allowed
	mtchView <- matchPattern(tar,genome[[chrnm]],max.mismatch=1)
	#### Calculate position weighted T-density scores
	#cat("Calculating position-weighted T-density scores...\n")
	# Create a list of start and end points of the 6 bp upstream of each target
	primrngs <- IRanges(start=start(mtchView)-6,width=6)
	targets <- DNAStringSet(genome[[chrnm]],start(mtchView),end(mtchView))
	prmrs <- DNAStringSet(genome[[chrnm]],start(primrngs),end(primrngs))
	tmp <- vmatchPattern("T",prmrs)
	tmp <- startIndex(tmp)
	primrnks <- lapply(tmp,function(x) sum(1/(11-x))/0.84563492) # 0.8456 is the maximum velcro score, corresponding to TTTTTT
	rm(tmp)

        #### Calculate distribution of target categories
        #cat("Calculating target category distribution...\n")
        # Store indices of sites of each category
        ict_s<-which(targets==tar & primrnks >= 0.5)
        icl_s<-which(targets==tar & primrnks < 0.5)
        iot_s<-which(targets!=tar & primrnks >= 0.5)
        iol_s<-which(targets!=tar & primrnks < 0.5)
        insites_s<-end(mtchView)

	cat("Creating EN site annotation: anti-sense...\n")

	tar <- DNAString("AAAA")
	mtchView <- matchPattern(tar,genome[[chrnm]],max.mismatch=1)
	#### Calculate position weighted T-density scores
	#cat("Calculating position-weighted T-density scores...\n")
	primrngs <- IRanges(start=end(mtchView)+1,width=6)
	targets <- DNAStringSet(genome[[chrnm]],start(mtchView),end(mtchView))
        prmrs <- DNAStringSet(genome[[chrnm]],start(primrngs),end(primrngs))
        tmp <- vmatchPattern("A",prmrs)
        tmp <- startIndex(tmp)
        primrnks <- lapply(tmp,function(x) sum(1/(x+4))/0.84563492) # 0.8456 is the maximum velcro score, corresponding to AAAAAA
        rm(tmp)

	#### Calculate distribution of target categories
	#cat("Calculating target category distribution...\n")
	# Store indices of sites of each category
	ict_ns<-which(targets==tar & primrnks >= 0.5)
	icl_ns<-which(targets==tar & primrnks < 0.5)
	iot_ns<-which(targets!=tar & primrnks >= 0.5)
	iol_ns<-which(targets!=tar & primrnks < 0.5)
	insites_ns<-start(mtchView)

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

	assign(paste0(chrnm,"ict"),cbind(ict_s,ict_ns))
	assign(paste0(chrnm,"icl"),cbind(icl_s,icl_ns))
	assign(paste0(chrnm,"iot"),cbind(iot_s,iot_ns))
	assign(paste0(chrnm,"iol"),cbind(iol_s,iol_ns))
	assign(paste0(chrnm,"insites"),cbind(insites_s,insites_ns))

	cat("Saving map file...\n")
	save(list=c(paste0(chrnm,"ict"),paste0(chrnm,"icl"),paste0(chrnm,"iot"),paste0(chrnm,"iol"),paste0(chrnm,"insites")),file=paste("../Data/",chrnm,"map.rda",sep=""))
	rm(list=c(paste0(chrnm,"ict"),paste0(chrnm,"icl"),paste0(chrnm,"iot"),paste0(chrnm,"iol"),paste0(chrnm,"insites")))
	end.time <- Sys.time()
        print(end.time-start.time)
}
