#!/usr/bin/env Rscript
#### Load libraries
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

#### Initialize parameters
copyNum <- 3 # Number of genome-wide insertions to simulate
ENifrc <- 0.1 # Fraction of Endonuclease-independent insertions
genome <- Hsapiens # Genome data structure
cat("\nCopy number: ",copyNum,"\n")
cat("ENi insertion fraction: ",ENifrc,"\n")

sites_loci<-c() # Empty arrays for insertion site chromosomes and loci
sites_chrm<-c()

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
	
	ict<-get(paste0(chrnm,"ict")) # Assign temporary variables to 
	icl<-get(paste0(chrnm,"icl"))
	iot<-get(paste0(chrnm,"iot"))
	iol<-get(paste0(chrnm,"iol"))
	insites<-get(paste0(chrnm,"insites"))

	chcopyNum<-chrmlist[[chrnm]]

	pd <- c(11.55*length(ict),7.25*length(icl),1.95*length(iot),1*length(iol))
	pd <- (pd/sum(pd))*(1-ENifrc)
	pd <- append(pd,ENifrc)
	cat("\nSite class distribution:\n",pd,"\n")

	#### Generate insertion sites
	cat("\nGenerating insertion sites...\n")
	classes <- sample(x = c(1:5),chcopyNum,replace=TRUE,prob=pd)
	sites = rep(0,chcopyNum)

	for (i in 1:chcopyNum) {
		if (classes[i]==1) {
			sites[i] <- insites[ict[runif(1,1,length(ict))]]
		} else if (classes[i]==2) {
			sites[i] <- insites[icl[runif(1,1,length(icl))]]
		} else if (classes[i]==3) {
			sites[i] <- insites[iot[runif(1,1,length(iot))]]
		} else if (classes[i]==4) {
			sites[i] <- insites[iol[runif(1,1,length(iol))]]
		} else if (classes[i]==5) {
			sites[i]<-runif(1,1,length(genome[[chrnm]]))
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
}

rm(ict,icl,iot,iol,insites)

#### Create sequences for insertion
load("../Data/L1RankTable.rda")
L1RankTable$score[1:40] <- L1RankTable$score[1:40]/sum(L1RankTable$score[1:40])
l1indcs <- sample(x=c(1:40),copyNum,replace=TRUE,prob=L1RankTable$score[1:40])
trpd <- read.table("../Data/L1trunc.dat",sep=",")
tdpd <- read.table("../Data/Transduction_PD.dat",sep="\t")
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
save(tdlenv,trlenv,l1s,l1indcs,sites_loci,sites_chrm,file="../Data/gen-sim-out.rda")

