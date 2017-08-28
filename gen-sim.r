#!/usr/bin/env Rscript
#### Load libraries
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
args<-commandArgs(trailingOnly=TRUE)

#### Initialize parameters
copyNum <- args[1] # Number of genome-wide insertions to simulate
ENifrc <- 0.1 # Fraction of Endonuclease-independent insertions
genome <- Hsapiens # Genome data structure
cat("\nCopy number: ",copyNum,"\n")
cat("ENi insertion fraction: ",ENifrc,"\n")

sites_loci<-c() # Empty arrays for insertion site chromosomes and loci
sites_chrm<-c()
sites_strand<-c()
strdict<-c("+","-")
names(strdict)<-c(1,2)

#### Sample chromosomes based on probability ranking
load("../Data/chrmpd.rda")
chrmlist<-sample(x=names(genome)[1:24],copyNum,replace=TRUE,prob=chrmcnt[,1]) # Here, simpyl the number of TTTT patterns provides the ranking 
chrmlist<-table(chrmlist)
cat("\nChromosomes: ",names(chrmlist),"\n")

#### Load map file for chosen chromosomes
for (i in names(chrmlist)) {
	cat("\nLoading map file...")
	load(paste0("../Data/",i,"map.rda"))
}
cat("\n")

ptm <- proc.time() # Begin timer

for (chrnm in names(chrmlist)) { # Loop through chromosomes in chrmlist
	
	cat("\nChromosome: ",chrnm)

	ict<-get(paste0(chrnm,"ict"))
	icl<-get(paste0(chrnm,"icl"))
	iot<-get(paste0(chrnm,"iot"))
	iol<-get(paste0(chrnm,"iol"))
	insites<-get(paste0(chrnm,"insites"))

	chrcopyNum<-chrmlist[[chrnm]]

	pd <- c(11.55*length(which(!is.na(ict))),7.25*length(which(!is.na(icl))),1.95*length(which(!is.na(iot))),1*length(which(!is.na(iol))))
	pd <- (pd/sum(pd))*(1-ENifrc)
	pd <- append(pd,ENifrc)
	cat("\nSite class distribution:\n",pd)

	#### Generate insertion sites
	classes <- sample(x = c(1:5),chrcopyNum,replace=TRUE,prob=pd)
	sites <- rep(0,chrcopyNum)
	strand <-rep(0,chrcopyNum)

	for (i in 1:chrcopyNum) {
		if (classes[i]==1) {
			tmp<-sample(c(1,2),1)
			sites[i] <- insites[ict[sample(c(1:length(which(!is.na(ict[,tmp])))),1),tmp]]
			strand[i] <- strdict[[tmp]]
		} else if (classes[i]==2) {
			tmp<-sample(c(1,2),1)
			sites[i] <- insites[icl[sample(c(1:length(which(!is.na(icl[,tmp])))),1),tmp]]
			strand[i] <- strdict[[tmp]]
		} else if (classes[i]==3) {
			tmp<-sample(c(1,2),1)
			sites[i] <- insites[iot[sample(c(1:length(which(!is.na(iot[,tmp])))),1),tmp]]
			strand[i] <- strdict[[tmp]]
		} else if (classes[i]==4) {
			tmp<-sample(c(1,2),1)
			sites[i] <- insites[iol[sample(c(1:length(which(!is.na(iol[,tmp])))),1),tmp]]
			strand[i] <- strdict[[tmp]]
		} else if (classes[i]==5) {
			sites[i]<-runif(1,1,length(genome[[chrnm]]))
			strand[i] <- strdict[[sample(c(1,2),1)]]
		}
	}

	cat("\nInsertion sites:\n")
	cat(sites,"\n")
	#cat("Site targets:\n")
	#for (i in 1:copyNum) {
	#	print(chr[(sites[i]-3):sites[i]])
	#}
	append(sites_loci,sites)
	append(sites_chrm,chrnm)
	append(sites_strand,strand)
}

rm(ict,icl,iot,iol,insites)

#### Create sequences for insertion
load("../Data/L1RankTable.rda")
L1RankTable$score[1:40] <- L1RankTable$score[1:40]/sum(L1RankTable$score[1:40])
l1indcs <- sample(x=c(1:40),copyNum,replace=TRUE,prob=L1RankTable$score[1:40])
trpd <- read.table("../Data/L1truncpd.csv",sep=",")
tdpd <- read.table("../Data/L1transdpd.txt",sep="\t")
trfrcv <- sample(x = trpd[[1]], copyNum, replace = TRUE, prob = trpd[[2]])
tdlenv <- sample(x = tdpd[[2]], copyNum, replace = TRUE, prob = tdpd[[2]])
trlenv<-rep(0,copyNum)
for (i in 1:copyNum) {
	len <- L1RankTable[[3]][l1indcs[i]]-L1RankTable[[2]][l1indcs[i]]
	trlenv[i] <- round(len*trfrcv[i])
}
gr <- GRanges(L1RankTable[[1]][l1indcs],IRanges(L1RankTable[[2]][l1indcs]+trlenv,L1RankTable[[3]][l1indcs]+tdlenv),strand=L1RankTable[[5]][l1indcs])
l1s <- getSeq(genome,gr)
cat("\nRunning time:\n")
proc.time() - ptm
cat("\nSaving image...\n")
save(tdlenv,trlenv,l1s,l1indcs,sites_loci,sites_chrm,sites_strand,file="../Data/gen-sim-out.rda")

