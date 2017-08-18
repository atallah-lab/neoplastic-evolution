#!/usr/bin/env Rscript

#### Load libraries
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

genome <- Hsapiens
chrmpd<-array(0,dim=c(24,4))
chrmcnt<-array(0,dim=c(24,4))

j=1
for (i in names(genome)[1:24]){
	load(paste("../Data/",i,"map.rda",sep=""))
	chrmpd[j,]<-c(11.55*length(ict),7.25*length(icl),1.95*length(iot),1*length(iol))
	chrmcnt[j,]<-c(length(ict),length(icl),length(iot),length(iol))
	chrmpd[j,]<-chrmpd[j,]/sum(chrmpd[j,])
	j=j+1
}
save(chrmpd,chrmcnt,file="../Data/chrmpd.rda")
rm(list=ls())
