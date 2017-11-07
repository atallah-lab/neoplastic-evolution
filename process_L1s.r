process_L1s <- function(genome, L1RankTable, trpd, tdpd, copyNum) {

	tmp <- L1RankTable$score # Sort table by bit score
	tmp <- tmp[order(-tmp)]
	tmp <- tmp[1:41] # Take top 40 as active (Brouha et al.)
	tmp <- (tmp-min(tmp))/(max(tmp)-min(tmp)) # Normalize bit score (0-1)
	tmp <- tmp^5.2 # Scale so that ~84% of transposition probability lies in top 6

	l1indcs <- sample(x=c(1:41), copyNum, replace=TRUE, prob=tmp)

	trfrc <- sample(x=trpd[[1]], copyNum, replace=TRUE, prob=trpd[[2]])
	tdlen <- sample(x=tdpd[[1]], copyNum, replace=TRUE, prob=tdpd[[2]])

	trlen <- rep(0,copyNum)
	for (i in 1:copyNum){
		len <- L1RankTable[[3]][l1indcs[i]]-L1RankTable[[2]][l1indcs[i]]
		trlen[i] <- len-round(len*trfrc[i])
	}

	gr <- GRanges(L1RankTable[[1]][l1indcs],IRanges(L1RankTable[[2]][l1indcs]+trlen,L1RankTable[[3]][l1indcs]+tdlen),strand=L1RankTable[[5]][l1indcs])
	l1s <-getSeq(genome,gr)
	return(list(l1s,l1indcs,tdlen,trlen))
}
