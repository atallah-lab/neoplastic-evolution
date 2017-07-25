#!/usr/bin/env Rscript
#### Load libraries
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

#### Initialize parameters
copyNum <- 20
ENifrc <- 0.1
genome <- Hsapiens
cat("\nCopy number: ",copyNum,"\n")
cat("ENi insertion fraction: ",ENifrc,"\n")

chrmlist<-sample(x=names(genome)[1:24],copyNum,replace=TRUE)
chrmlist<-table(chrmlist)
cat("\nChromosomes: ",names(chrmlist))

#### Load map file
load_map <- function(chrmnm) {
	mapnm<-paste(chrmnm,"gmap.rda",sep="")
	if (!file.exists(mapnm)) {
		stop("Map file does not exist")
	} else {
		cat("\nLoading map file...\n")
		load(mapnm)
	}
	remove(chrnm,ptm)
}
lapply(names(chrmlist),load_map)

ptm <- proc.time()
for (chrnm in names(chrmlist)) {

	chcopyNum<-chrmlist[[chrnm]]

	pd <- c(11.55*length(ict),7.25*length(icl),1.95*length(iot),1*length(iol))
	pd <- (pd/sum(pd))*(1-ENifrc)
	pd <- append(pd,ENifrc)
	cat("\nSite class distribution:\n",pd,"\n")

	#### Generate insertion sites
	cat("\nGenerating insertion sites...\n")
	classes <- sample(x = c(1:5),chcopyNum,replace=TRUE,prob=pd)
	sites = rep(0,chcopyNum)
#	cct<-0
#	ccl<-0
#	cot<-0
#	col<-0
#	cr <-0
	for (i in 1:chcopyNum) {
		if (classes[i]==1) {
			sites[i] <- insites[ict[runif(1,1,length(ict))]]
#			cct<-cct+1
		} else if (classes[i]==2) {
			sites[i] <- insites[icl[runif(1,1,length(icl))]]
#			ccl<-ccl+1
		} else if (classes[i]==3) {
			sites[i] <- insites[iot[runif(1,1,length(iot))]]
#			cot<-cot+1
		} else if (classes[i]==4) {
			sites[i] <- insites[iol[runif(1,1,length(iol))]]
#			col<-col+1
		} else if (classes[i]==5) {
			sites[i]<-runif(1,1,length(genome$chrnm))
#			cr<-cr+1
	}

	cat("\nInsertion sites:\n")
	cat(sites,"\n")
	#cat("Site targets:\n")
	#for (i in 1:copyNum) {
	#	print(chr[(sites[i]-3):sites[i]])
	#}
	gsites<-append(gsites,sites)
}

#### Create sequences for insertion
l1s <- readDNAStringSet("../hgL1.fa")
l1s <- l1s[floor(runif(copyNum,1,length(l1s)))]
trpd <- read.table("../L1trunc.dat",sep=",")
trfrcv = sample(x = trpd[[1]], copyNum, replace = TRUE, prob = trpd[[2]])
cat("\nTruncated fractions: ",1-trfrcv,"\n")
for (i in 1:copyNum) {
	len <- length(l1s[[i]])
	trlen <- round(len*trfrcv[i])
	l1s[[i]] <- l1s[[i]][len-trlen:len]
}
cat("\nRunning time:\n")
proc.time() - ptm
cat("\nSaving image...\n")
save(l1s,gsites,chrmlist,"sim-out.rda")

