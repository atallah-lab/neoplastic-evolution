#!/usr/bin/env Rscript

#### Load libraries
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

genome <- Hsapiens
chrmpd<-array(0,dim=c(24,4))
chrmcnt<-array(0,dim=c(24,4))

j=1
for (i in names(genome)[1:24]){

	load(paste0("../Data/",i,"map.rda"))

        ict<-get(paste0(i,"ict"))
        icl<-get(paste0(i,"icl"))
        iot<-get(paste0(i,"iot"))
        iol<-get(paste0(i,"iol"))
        insites<-get(paste0(i,"insites"))

	chrmpd[j,]<-c(11.55*length(which(!is.na(ict))),7.25*length(which(!is.na(icl))),1.95*length(which(!is.na(iot))),1*length(which(!is.na(iol))))
	chrmcnt[j,]<-c(length(which(!is.na(ict))),length(which(!is.na(icl))),length(which(!is.na(iot))),length(which(!is.na(iol))))
	chrmpd[j,]<-chrmpd[j,]/sum(chrmpd[j,])
	j=j+1
}
save(chrmpd,chrmcnt,file="../Data/chrmpd.rda")
