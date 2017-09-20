#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# get_sv_dist.r #NAMING
#
# This script loads separate data files that have been generated by 'mapgen.r' 
# and stores pd and counts of each EN-site category as a 24 x 4 array in file
# 'chrmpd.rda'. 
#
# Currently this probability distribution vector is generated from the count of 
# closed-snap/tight-velcro sites (type with the highest enrichment) in each 
# chromosome.
#
# Usage: ./get_sv_dist.r 
# Output: file named 'chrmpd.rda' to data directory.
#
# Dependencies: R(>= 2.8.0, Packages - Biostrings, BSgenome (for default hg38), 
# GenomicRanges)
# -----------------------------------------------------------------------------

#### Load libraries
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

#--- hg38 reference genome
genome <- Hsapiens

#--- Initialize arrays, one for chr's probability distribution, one for count.
chrmpd<-array(0,dim=c(24,4))
chrmcnt<-array(0,dim=c(24,4))

idx=1
for (chrnm in names(genome)[1:24]){

	load(paste0("./data/",chrnm,"map.rda")) # reads data file as '<chr1>map.rda'.

        ict<-get(paste0(chrnm,"ict")) # closed tight
        icl<-get(paste0(chrnm,"icl")) # closed loose
        iot<-get(paste0(chrnm,"iot")) # open tight
        iol<-get(paste0(chrnm,"iol")) # open loose
        insites<-get(paste0(chrnm,"insites"))

	chrmpd[idx,]<-c(11.55*length(which(!is.na(ict))),7.25*length(which(!is.na(icl))),1.95*length(which(!is.na(iot))),1*length(which(!is.na(iol))))
	chrmcnt[idx,]<-c(length(which(!is.na(ict))),length(which(!is.na(icl))),length(which(!is.na(iot))),length(which(!is.na(iol))))
	chrmpd[idx,]<-chrmpd[idx,]/sum(chrmpd[idx,])
	idx=idx+1
}
#--- Saves probability distribution file for all seleted chromosomes.
save(chrmpd,chrmcnt,file="./data/chrmpd.rda") 
