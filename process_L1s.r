##--------------------------------------------------
## This function simulates the sampling of L1s based on
## estimated activity ranking, and the processes of
## transduction and truncation before insertion
##--------------------------------------------------

process_L1s <- function(genome, L1RankTable, trpd, tdpd, copyNum) {

	tmp <- L1RankTable$score # Get score (bit score) from L1RankTable
	tmp <- tmp[1:41] # Take top 40 as active (Brouha et al.)

	# Sample copyNum L1s from the list (with replacement), based on activity ranking
	l1indcs <- sample(x=c(1:41), copyNum, replace=TRUE, prob=tmp)

	# Sample copyNum truncation fractions and transduction lengths from their respective 
	# probability densities trpd, and tdpd
	trfrc <- sample(x=trpd[[1]], copyNum, replace=TRUE, prob=trpd[[2]])
	tdlen <- sample(x=tdpd[[1]], copyNum, replace=TRUE, prob=tdpd[[2]])

	# Convert fraction of sequence left after truncation to length of truncated portion
	trlen <- rep(0,copyNum)
	for (i in 1:copyNum){
		len <- L1RankTable[[3]][l1indcs[i]]-L1RankTable[[2]][l1indcs[i]]
		trlen[i] <- len-round(len*trfrc[i])
	}

	# Obtain the sequences which will be inserted
	gr <- GRanges(L1RankTable[[1]][l1indcs],IRanges(L1RankTable[[2]][l1indcs]+trlen,L1RankTable[[3]][l1indcs]+tdlen),strand=L1RankTable[[5]][l1indcs])
	l1s <-getSeq(genome,gr)

	# Return
	# l1s - DNAStringSet object containing inserted sequences
	# l1indcs - Indices of L1s chosen from the rank table
	# tdlen - Array of sampled transduction lengths
	# trlen - Array of sampled lengths of truncated portions
	return(list(l1s,l1indcs,tdlen,trlen))
}
