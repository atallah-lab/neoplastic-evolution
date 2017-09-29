#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# hgextract.r

# This script extracts segments of the hg38 reference genome that are provided 
# to it in tab-delimited form, assuming the following format for each row: 
#
#	chromosomeName		startPt		endPt
#
# It writes the extractions to a Fasta file, with each sequence numbered.
#
# Usage: ./hgextract.r <...> <...>
# Output: <...>
#
# Dependencies: R(>= 2.8.0, Packages - Biostrings, BSgenome (for default hg38), 
# GenomicRanges)
# -----------------------------------------------------------------------------

#--- Load libraries
library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

args = commandArgs(trailingOnly=TRUE)

if (length(args) <2) {
message("Usage: ./hgextract.r <> <>")
stop("Please provide input and output file names.")
} else if (length(args)>2) {
message("Usage: ./hgextract.r <> <>")
message("Only first two arguments used.")
}
tmp<-read.table(args[1])
gr <- GRanges(tmp[[1]],IRanges(tmp[[2]],tmp[[3]]))
seqs <- getSeq(Hsapiens,gr)
names(seqs)<-c(1:length(seqs))

cat("\nWriting output file (Fasta)...\n\n")
writeXStringSet(seqs,file=args[2],format="fasta")

