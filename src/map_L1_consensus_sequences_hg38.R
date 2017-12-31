library(Biostrings)
library(GenomicRanges)
library(rtracklayer)

# May need to...
# source("https://bioconductor.org/biocLite.R")
# setwd("/path/to/sim-develop/repo")

require(BSgenome.Hsapiens.UCSC.hg38)

# http://l1base.charite.de/BED/hsflil1_8438.bed
# List of known chromosomes containing L1's. 
bed_reads <- GRanges(import.bed(con="./data/hsflil1_8438.bed"))

#iranges object for each chr
iranges_list <- list()

# Tracks listed in 'hsflil1_8438.bed'. Look in ./data directory for this file. 
chromosomes <- c('chr1', 
				'chr2', 
				'chr3', 
				'chr4', 
				'chr5', 
				'chr6', 
				'chr7', 
				'chr8', 
				'chr10', 
				'chr11', 
				'chr12', 
				'chr16', 
				'chr17')

for (i in 1:length(chromosomes)) {
  iranges_list[[i]] <- ranges(bed_reads[seqnames(bed_reads) == chromosomes[[i]] ])
}

# 'match' is a list containing objects per chromosome. 
# 'match_list' is a list containing match objects. 
# The subject sequences derive from 'hsflil1_8438.bed' file (in ./data directory).
# . 
check_for_hits_against_consensus_L1s <- function(consensus, relax) {
	match <- list()
	match_list <- list()
		for (i in 1:length(chromosomes)) {
		  # get unmasked sequence for this chromosome
		  unmasked <- unmasked(Hsapiens[[i]])

		  for (j in 1:length(consensus)) { # for number of sequences in this accession number.
		  	# is there a match between accession sequence and area of this chromosome? 
			match[[i]] <- matchPattern(consensus[[j]], unmasked[start(iranges_list[[i]]):end(iranges_list[[i]])], max.mismatch=relax)
			# add views objects for this accession # to list 
			match_list[[j]] <- match
		}
	}
	return (match_list)
}

# Intact L1's taken from Brouha paper, 'Hot L1s account for the 
# bulk of retrotransposition in the human population' (Table 4). 
# The top 9 accessions were used to perform BLAST searches and 
# resulting sequences were applied over a list of known ranges in hg38. 

# BLAST searches were performed for the following accession numbers:
# gb|AC002980 (complete sequence), emb|AL512428, gb|AC021017, gb|AC004200, emb|AL356438
# emb|AL137845, emb|AL121825, gb|AC091612, emb|AL354951.

# In addition, we use two 'canonical' consensus sequences, the s0 sequence is from NCBI 
# and s1 is from the Brouha paper. 

# s0 is a consensus sequence found at: 
# 	https://www.ncbi.nlm.nih.gov/nuccore/L19092.1?report=GenBank
s0 = readDNAStringSet("./data/GenBank_L19092.fa")

# Brouha consensus sequence: 
s1 = readDNAStringSet("./data/Brouha_consensus_L1.fa")

# Sequences derived from BLAST search
s2 = readDNAStringSet("./data/gb_AC002980.fa")
#s3 = readDNAStringSet("./data/gb_AC004200.fa")
#s4 = readDNAStringSet("./data/gb_AC021017.fa")
#s5 = readDNAStringSet("./data/gb_AC091612.fa")
#s6 = readDNAStringSet("./data/emb_AL121825.fa")
#s7 = readDNAStringSet("./data/emb_AL137845.fa")
#s8 = readDNAStringSet("./data/emb_AL354951.fa")
#s9 = readDNAStringSet("./data/emb_AL356438.fa")
#s10 = readDNAStringSet("./data/emb_AL512428.fa")

# Sample run. Substitute any of the above sequences to test. 
# The second argument is the max number of acceptable mismatches. 

# EXAMPLE:
# for s0, no hits
hits <- check_for_hits_against_consensus_L1s(s0, 5)
# for s1, no hits
hits <- check_for_hits_against_consensus_L1s(s1, 5)
# for s2, hits for all 5 sequences
hits <- check_for_hits_against_consensus_L1s(s2, 5)


# ------------------- PRINT TO FILE -----------------------------------------
writeHitsToFile <- function(seqname, matches, strand, filename, append=FALSE) {
  
  hits <- data.frame(seqname=rep.int(seqname, length(matches)),
                     start=start(matches),
                     end=end(matches),
                     strand=rep.int(strand, length(matches)),
                     patternID=names(matches),
                     check.names=FALSE)
  
  write.table(hits, file=filename, append=append, 
              quote=FALSE, sep="\t",row.names=FALSE, col.names=!append)
}

# ------------------- SEARCH WHOLE GENOME AND PRINT  ------------------------
match_hg38_L1_sites <- function(patterns, x, y, append=FALSE) {
  
  library(BSgenome.Hsapiens.UCSC.hg38)
  genome <- BSgenome.Hsapiens.UCSC.hg38
  seqnames <- seqnames(genome)
  seqnames_as_string <- paste(seqnames[1:x], collapse=", ")
  
  cat("Searching", providerVersion(genome), "for L1 sites on",
      "chromosome[s]", seqnames_as_string, "...\n")
  append <- TRUE
  
  number_of_sequences <- length(seqnames[1:x])
  
  for (i in 1:number_of_sequences) {
    target_sequence <- unmasked(genome[[i]])
    chr <- paste(c("chr", i), collapse = "")
    
    cat(">>> Finding L1s for", chr, "...\n")
    
    for (i in seq_len(length(patterns))) {
      
      # setup outfile to print results
      datetime_stamp <- format(Sys.time(), "%m%d%Y")
      f <- paste(chr, datetime_stamp, sep="_")
      
      pattern <- patterns[[i]]
      matches <- matchPattern(pattern, target_sequence, max.mismatch=y)
      names(matches) <- rep.int(i, length(matches))
      
      cat(">>> Writing hits to file for", paste(toString(pattern[1:20], "...\n")), "on chromosome", chr, "...\n")
      outfilename <- paste(f, pattern[1:10], sep="_")
      outfilename <- paste(outfilename, ".txt", sep="_")
      outfile <- file(outfilename, "w")
      
      writeHitsToFile(chr, matches, "+", file=outfile, append=append)
      cat(">>> DONE.\n")
    }
    cat(">>> DONE.\n")
    close(outfile)
  }
}

# ---------------------------------------------------------------------------
# Alternative helper function that's faster than stock model.
compare.myDNA <- function(x,y){
    as.integer(x) == as.integer(y)
}

# ---------------------------------------------------------------------------
# 'subjects' is a GRanges object containing all sequences 
# from the hsflil1_8438.bed file.
subjects <- GRanges()
for (i in 1:length(chromosomes)) {
  unmasked <- unmasked(Hsapiens[[i]])
  gr <- GRanges(chromosomes[[i]], IRanges(start=start(iranges_list[[i]]),end=end(iranges_list[[i]])), strand="+" )
  subjects <- append(subjects, gr)
}