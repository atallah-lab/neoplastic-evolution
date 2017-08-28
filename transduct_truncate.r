#!/usr/bin/env Rscript
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

args <- commandArgs(trailingOnly=TRUE)
copyNum<-args[1]
genome<-Hsapiens

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
cat("\nTruncation fractions: ",trfrcv,"\n")
cat("\nTransduction lengths: ",trlenv,"\n")
gr <- GRanges(L1RankTable[[1]][l1indcs],IRanges(L1RankTable[[2]][l1indcs]+trlenv,L1RankTable[[3]][l1indcs]+tdlenv),strand=L1RankTable[[5]][l1indcs])
l1s <- getSeq(genome,gr)
cat("\n")
l1s
cat("\nSaving image...\n")
save(tdlenv,trlenv,l1s,l1indcs,file="trans_trunc-out.rda")
