#!/usr/bin/env Rscript

#### Load libraries
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
args = commandArgs(trailingOnly=TRUE)

#### Initialize parameters
copyNum <- 20
ENifrc <- 0.1
chrnm <- "chr1"
genome <- Hsapiens
cat("\nCopy number: ",copyNum,"\n")
cat("\nENi insertion fraction: ",ENifrc,"\n")
mapnm <- paste(chrnm,"gmap.rda",sep="")
#### Load (existing) or generate (non-existant) map file
if (!file.exists(mapnm)) {
	if (length(args) < 2) {
		stop("gmap.rda does not exist, please provide ENi insertion fraction and chromosome name as arguments 1 and 2")
	} else {
		system(paste("./genmap.r",args[1],args[2]))
	}
} else {
	cat("\nLoading map file...\n")
	load(mapnm)
}
remove(.Random.seed)

pd <- c(11.55*length(ict),7.25*length(icl),1.95*length(iot),1*length(iol))
pd <- (pd/sum(pd))*(1-ENifrc)
pd <- append(pd,ENifrc)
cat("\nSite class distribution:\n",pd,"\n")

ptm <- proc.time()
#### Generate insertion sites
cat("\nGenerating insertion sites...\n")
classes <- sample(x = c(1:5),copyNum,replace=TRUE,prob=pd)
sites = rep(0,copyNum)
cct<-0
ccl<-0
cot<-0
col<-0
for (i in 1:copyNum) {
	if (classes[i]==1) {
		sites[i] <- insites[ict[runif(1,1,length(ict))]]
		cct<-cct+1
	} else if (classes[i]==2) {
		sites[i] <- insites[icl[runif(1,1,length(icl))]]
		ccl<-ccl+1
	} else if (classes[i]==3) {
                sites[i] <- insites[iot[runif(1,1,length(iot))]]
		cot<-cot+1
        } else if (classes[i]==4) {
                sites[i] <- insites[iol[runif(1,1,length(iol))]]
		col<-col+1
        } else if (classes[i]==5) {
		sites[i]<-runif(1,1,length(genome$chrnm))
}

cat("\nInsertion sites:\n")
cat(sites,"\n")
#cat("Site targets:\n")
#for (i in 1:copyNum) {
#	print(chr[(sites[i]-3):sites[i]])
#}

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
remove(ict,icl,iot,iol,insites,chr,len,trlen,classes,i,args)
cat("\nSaving image...\n")
save.image("sim-out.rda")
