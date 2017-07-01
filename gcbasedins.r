#!/usr/bin/env Rscript

library(Biostrings)

####### Parse cmd line args, read input genome
args = commandArgs(trailingOnly=TRUE)
if (length(args)==1) {
	library(BSgenome.Hsapiens.UCSC.hg19)
	chrm="chr17"
	genome = DNAStringSet(Hsapiens[[chrm]])
	names(genome) = c(chrm)
	outname=args[1]
	cat("\n\nUsing hg19 chromosome 17...\n")
} else if (length(args)==2) {
	genome <- readDNAStringSet(args[1])
	chrm=names(genome)[1]
	cat("\n\nFasta file was read...\n")
	outname=args[2]
} else if (length(args)<1) {
	print("Please provide at least one argument: output file name")
}

genlen = length(genome[[chrm]])

####### gcContent function
gcContent <- function(x) {
	alf <- alphabetFrequency(x, as.prob=TRUE)
	sum(alf[c("G","C")])
}

copyNum=3
cat("\nCopy number: ", copyNum, "\n")
locGCWin=6
cat("\nLocal G/C content window:", locGCWin, "bp\n")

####### Generate random breakpoints and filter based on 
####### local G/C content (+/- locGCWin bp)
bkpts=c()
cat("\nGenerating breakpoints\n")
while (length(bkpts) < copyNum){
	i <- floor(runif(1,(locGCWin/2)+1,genlen-(locGCWin/2)))
        if (runif(1,0,1) <=  gcContent(genome[[chrm]][i-(locGCWin):i+(locGCWin/2)])) {
                bkpts <- append(bkpts, i)
		cat("\tBreakpoint accepted: ",i,"\n")
        }
#	else {
#		cat("\tBreakpoint rejected")
#	}
}

####### Insert sequences at accepted breakpoints (copy and paste type)
cat("\nInserting sequences...\n")
L1 = readDNAStringSet("~/jackgl/Data/humangenome/hgL1.fa")
L1 = L1[["L1HS\tL1\tHomo sapiens"]]
#L1 = DNAString(c("ATGC"))
for (i in bkpts){
	tmp <- append(genome[[chrm]][1:i],L1)
	genome[[chrm]] <- append(tmp,genome[[chrm]][(i+1):genlen])
}

####### Write output fasta file
cat("\nWriting output file (Fasta)...\n\n")
writeXStringSet(genome,file=outname,format="fasta")

