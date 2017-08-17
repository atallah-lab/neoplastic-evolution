#!/usr/bin/env Rscript

library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

#### Read reference genome
cat("\nReading reference genome (GRCh38)...\n")
genome<-Hsapiens

rmlist<-c("args","rmlist","chrnm","genome","tar","mtchView","targets","primrngs","primrnks","prmrs","insites")

for (chrnm in names(genome)[11:24]) {
	start.time<-Sys.time()
	cat("\nChromosome: ",chrnm,"\n")
	#### Create EN target site annotation
	cat("Creating EN target site annotation...\n")
	# Define target sequence
	tar <- DNAString("TTTT")
	# Find locations of target in chromosome, one mismatch is allowed
	mtchView <- matchPattern(tar,genome[[chrnm]],max.mismatch=1)
	# Get start points of all targets
	insites <- start(mtchView)
	#### Calculate position weighted T-density scores
	cat("Calculating position weighted T-density scores...\n")
	# Create a list of start and end points of the 6 bp upstream of each target
	primrngs <- IRanges(start=insites-9,width=6)
	targets <- DNAStringSet(genome[[chrnm]],start(mtchView),end(mtchView))
	prmrs <- DNAStringSet(genome[[chrnm]],start(primrngs),end(primrngs))
	tmp <- vmatchPattern("T",prmrs)
	tmp <- startIndex(tmp)
	primrnks <- lapply(tmp,function(x) sum(1/(x+4))/0.84563492) # 0.8456 is the maximum primer score, corresponding to a TTTTTT primer
	remove(tmp)

	#### Calculate distribution of target categories
	cat("Calculating target category distribution...\n")
	# Store indices of sites of each category
	assign(paste0("ict",chrnm),which(targets==tar & primrnks >= 0.5))
	assign(paste0("icl",chrnm),which(targets==tar & primrnks < 0.5))
	assign(paste0("iot",chrnm),which(targets!=tar & primrnks >= 0.5))
	assign(paste0("iol",chrnm),which(targets!=tar & primrnks < 0.5))

	assign(paste0("insites",chrnm),start(mtchView))
	cat("Saving map file...\n")
	save(list=ls(1)[which(!ls(1) %in% rmlist)],file=paste("../Data/",chrnm,"map.rda",sep=""))
	rm(list=c(paste0("ict",chrnm),paste0("icl",chrnm),paste0("iol",chrnm),paste0("iot",chrnm),paste0("insites",chrnm)))
	end.time <- Sys.time()
        print(end.time-start.time)
}
