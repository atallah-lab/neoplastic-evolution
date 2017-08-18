#!/usr/bin/env Rscript

#### Load libraries
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
args = commandArgs(trailingOnly=TRUE)

#### Initialize parameters
copyNum <- args[2]
ENifrc <- 0.1
chrnm <- args[1]
genome <- Hsapiens
cat("\nChromosome: ",chrnm,"\n")
cat("\nCopy number: ",copyNum,"\n")
cat("\nENi insertion fraction: ",ENifrc,"\n")
mapnm <- paste("../Data/",chrnm,"map.rda",sep="")
#### Load (existing) or generate (non-existant) map file
if (!file.exists(mapnm)) {
	stop("Map file does not exist")
} else {
	cat("\nLoading map file...\n")
	load(mapnm)
}
rm(.Random.seed)

ict<-get(paste0(chrnm,"ict")) # Assign temporary variables
icl<-get(paste0(chrnm,"icl"))
iot<-get(paste0(chrnm,"iot"))
iol<-get(paste0(chrnm,"iol"))
insites<-get(paste0(chrnm,"insites"))

ptm <- proc.time()
pd <- c(11.55*length(ict),7.25*length(icl),1.95*length(iot),1*length(iol))
pd <- (pd/sum(pd))*(1-ENifrc)
pd <- append(pd,ENifrc)
cat("\nSite class distribution:\n",pd,"\n")

#### Generate insertion sites
cat("\nGenerating insertion loci...\n\n")
classes <- sample(x=c(1:5),copyNum,replace=TRUE,prob=pd)
sites = rep(0,copyNum)

for (i in 1:copyNum) {
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

proc.time() - ptm
cat("\nSaving image...\n")
save(sites,chrnm,ENifrc,pd,file=paste("../Data/",chrnm,"-sim-out.rda",sep=""))
