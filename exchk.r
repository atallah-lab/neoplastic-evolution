#!/usr/bin/env Rscript

##### This function will check if a generated breakpoint lies within an exon

##### Inputs :  1) input gff3
#		2) chromosome name
#		3) breakpoint location 

library(data.table)
args = commandArgs(trailingOnly=TRUE)
cat("\nInput gff3:",args[1],"\n")
cat("\nGenerating exon annotation... \n")
system(paste("./filtgff3.bash",args[1]))
cat("Done.\n")

genann = read.table(paste(args[1],".filt", sep=""))
chrmann = genann[genann$V1 == args[2],]
if (nrow(chrmann)==0){
	stop("Flawed chromosome name, 1-22, X, Y, or MT")
}

bkpt <- as.numeric(args[3])
excheck <- function(ann, bpt) {
	tmp = between(bpt,ann[,2],ann[,3])
	length(which(tmp==TRUE))
}

cat("\nBreakpoint inside exon (1/0): ",excheck(chrmann,bkpt),"\n")


