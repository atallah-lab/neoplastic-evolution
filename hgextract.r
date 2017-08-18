#!/usr/bin/env Rscript

# This script will extract segments from hg38 provided in a tab-delimited
# text file of the following format for each row: chromosomeName	startPt	endPt
# and write the extractions to a Fasta file, with the sequences simply numbered.

library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

args = commandArgs(trailingOnly=TRUE)
if (length(args) <2) {
	print("Please provide input and output file names")
} else if (length(args)>2) {
	print("Only first two arguments used")
}
tmp<-read.table(args[1])
gr <- GRanges(tmp[[1]],IRanges(tmp[[2]],tmp[[3]]))
seqs <- getSeq(Hsapiens,gr)
names(seqs)<-c(1:length(seqs))
cat("\nWriting output file (Fasta)...\n\n")
writeXStringSet(seqs,file=args[2],format="fasta")

