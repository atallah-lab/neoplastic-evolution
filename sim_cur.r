#!/usr/bin/env Rscript

library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

copyNum <- 3
chrnm <- "chr1"
ENifrc <- 0.1
cat("\nCopy number: ",copyNum,"\n")
cat("Chromosome: ",chrnm,"\n")
cat("ENi insertion fraction: ",ENifrc,"\n")

#### Read reference genome
cat("\nReading reference genome (GRCh38)...\n")
genome<-Hsapiens
chr<-genome[[chrnm]]

#### Create EN target site annotation
cat("\nCreating EN target site annotation...\n")
# Define target sequence
tar <- DNAString("TTTT")
# Find locations of target in chromosome, one mismatch is allowed
mtchView <- matchPattern(tar,chr,max.mismatch=1)
# Get start points of all targets
insites <- start(mtchView)
# Create a list of start and end points of the 6 bp upstream of each target
ranges <- IRanges(start=insites-9,width=6)

#### Calculate position weighted T-density scores
cat("\nCalculating position weighted T-density scores...\n")
targets <- DNAStringSet(chr,start(mtchView),end(mtchView))
prmrs <- DNAStringSet(chr,start(ranges),end(ranges))
tmp <- vmatchPattern("T",prmrs)
tmp <- startIndex(tmp)
prmrnks <- lapply(tmp,function(x) sum(1/x)/2.45) # 2.45 is the maximum primer score, corresponding to a TTTTTT primer
remove(tmp)

#### Calculate distribution of target categories
cat("\nGenerating target category distribution...\n")
ttight <- targets[prmrnks>=0.5]
tloose <- targets[prmrnks<0.5]
lenl <- length(tloose)
lenti <- length(ttight)
lenta <- length(targets)
ccl <-sum(vcountPattern("TTTT",tloose,max.mismatch=0))
cct <-sum(vcountPattern("TTTT",ttight,max.mismatch=0))
col <- lenl-ccl
cot <- lenti-cct
pd <- c(11.55*cct,7.25*ccl,1.95*cot,1*col)
pd <- (pd/sum(pd))*(1-ENifrc)

ict<-which(targets==tar & prmrnks >= 0.5)
icl<-which(targets==tar & prmrnks < 0.5)
iot<-which(targets!=tar & prmrnks >= 0.5)
iol<-which(targets!=tar & prmrnks < 0.5)

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
cat(sites,"\n\n")

#### Create sequence for insertion
#l1s <- readDNAStringSet("./hgL1.fa")
#l1 <- l1s[floor(runif(1,1,length(l1s)))]
#pd <- read.table("./truncpd.dat")
#i <- round(runif(1,1,500))
#trfrc <- pd[[2]][i]
#trlen <- round(length(l1)*trfrc)

#### Write output fasta file
#cat("\nWriting output file (Fasta)...\n\n")
#writeXStringSet(chrm,file=outname,format="fasta")

