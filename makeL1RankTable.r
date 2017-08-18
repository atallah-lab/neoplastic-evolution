#!/usr/bin/env Rscript

# Arguments: score file, l1base file, image name

library(Biostrings)
library(data.table)
args=commandArgs(trailingOnly=TRUE)

scoreTab <- read.table(args[1])
l1baseTab <- read.table(args[2])
L1RankTable<-data.table(
variable1=l1baseTab[,1][scoreTab[,1]],
variable2=l1baseTab[,2][scoreTab[,1]]+scoreTab[,2]-1,
variable3=l1baseTab[,2][scoreTab[,1]]+scoreTab[,3]-1
)
names(L1RankTable)<-c("chr","start","end")
L1RankTable$score <- scoreTab[,4]
L1RankTable$strand <- scoreTab[,5]
save.image(args[3])
