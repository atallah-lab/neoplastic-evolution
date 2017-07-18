#!/usr/bin/env Rscript

#### Load libraries
library(Biostrings)
args = commandArgs(trailingOnly=TRUE)

#### Initialize parameters
copyNum <- 3
cat("\nCopy number: ",copyNum,"\n")

#### Load (existing) or generate (non-existant) map file
if (!file.exists("gmap.rda")) {
	if (length(args) < 2) {
		stop("gmap.rda does not exist, please provide ENi insertion fraction and chromosome name as arguments")
	} else {
		system(paste("./genmap.r",args[1],args[2]))
	}
} else {
	cat("\nLoading map file...\n")
	load("gmap.rda")
}

#### Generate insertion sites
cat("\nGenerating insertion sites...\n")
sites <- rep(0,copyNum)
for (i in 1:copyNum) {
	n <- runif(1,0,1)
	if (n > 0 & n <= pd[1]) {
		j <-floor(runif(1,1,length(ict)))
		sites[i]<-insites[ict[j]] }
	if (n > pd[1] & n <= sum(pd[1:2])) {
		j <- floor(runif(1,1,length(icl)))
		sites[i]<-insites[icl[j]] }
	if (n > sum(pd[1:2]) & n <= sum(pd[1:3])) {
		j <- floor(runif(1,1,length(iot)))
		sites[i]<-insites[iot[j]] }
	if (n > sum(pd[1:3]) & n <= 1-ENifrc) {
		j <- floor(runif(1,1,length(iol)))
		sites[i]<-insites[iol[j]] }
	if (n > (1-ENifrc) & n <= 1) {
		sites[1]<-floor(runif(1,1,length(chr))) }
}

cat("\nInsertion sites:\n")
cat(sites,"\n")

#### Create sequence for insertion
l1s <- readDNAStringSet("./hgL1.fa")
l1 <- l1s[floor(runif(1,1,length(l1s)))]
l1 <- l1[[1]]
l1
trpd <- read.table("./trpd.dat",sep=",")
trfrc <- trpd[[2]][round(runif(1,1,500))]
cat("\nTruncation fraction: ",trfrc,"\n\n")
trlen <- round(length(l1)*trfrc)
l1 <- l1[length(l1)-trlen:length(l1)]
l1

#### Write output fasta file
#cat("\nWriting output file (Fasta)...\n\n")
#writeXStringSet(chrm,file=outname,format="fasta")

