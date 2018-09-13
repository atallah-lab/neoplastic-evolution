#!/usr/bin/env Rscript

mapsequence <- function(seq) {

	start.time<-Sys.time()
	
	map<-c()
	seq <- DNAString(paste0("GGGGGGGGGG",seq,"GGGGGGGGGG")) # add non-target buffer
	
#---Create EN site annotation for sense strand

	tar <- DNAString("TTTT") # Define target pattern
	mtchView <- matchPattern(tar,seq,max.mismatch=1) # Find locations of target in chromosome, one mismatch is allowed

	#--- Calculate position weighted T-density scores
	velrngs <- IRanges(start=start(mtchView)-6,width=6) # Create a list of start and end points of the 6 bp upstream of each target
	snaps <- DNAStringSet(seq,start(mtchView),end(mtchView))
	velcros <- DNAStringSet(seq,start(velrngs),end(velrngs))
	tmp <- vmatchPattern("T",velcros)
	tmp <- startIndex(tmp)
	velrnks <- lapply(tmp,function(x) sum(1/(11-x))/0.84563492) # 0.8456 is the maximum velcro score, corresponding to TTTTTT
	rm(tmp)

        # Store indices of sites of each category
        ict_s<-which(snaps==tar & velrnks >= 0.5)
        icl_s<-which(snaps==tar & velrnks < 0.5)
        iot_s<-which(snaps!=tar & velrnks >= 0.5)
        iol_s<-which(snaps!=tar & velrnks < 0.5)
        insites_s<-end(mtchView)-10

#--- Create EN site annotation for anti-sense strand

	tar <- DNAString("AAAA")
	mtchView <- matchPattern(tar,seq,max.mismatch=1)

	#--- Calculate position weighted T-density scores
	velrngs <- IRanges(start=end(mtchView)+1,width=6)
	snaps <- DNAStringSet(seq,start(mtchView),end(mtchView))
        velcros <- DNAStringSet(seq,start(velrngs),end(velrngs))
        tmp <- vmatchPattern("A",velcros)
        tmp <- startIndex(tmp)
        velrnks <- lapply(tmp,function(x) sum(1/(x+4))/0.84563492) # 0.8456 is the maximum velcro score, corresponding to AAAAAA
        rm(tmp)

	# Store indices of sites of each category
	ict_ns<-which(snaps==tar & velrnks >= 0.5)
	icl_ns<-which(snaps==tar & velrnks < 0.5)
	iot_ns<-which(snaps!=tar & velrnks >= 0.5)
	iol_ns<-which(snaps!=tar & velrnks < 0.5)
	insites_ns<-start(mtchView)-10

#--- Store the loci for the sense and anti-sense strands in a 2-column array

	# Make loci arrays equal length (the length of which ever is longer)	
	# The last elements of the column of the strand with fewer sites are filled with NA
	n <- max(length(ict_s), length(ict_ns))
	length(ict_s) <- n                      
	length(ict_ns) <- n
    n <- max(length(icl_s), length(icl_ns))
	length(icl_s) <- n
    length(icl_ns) <- n
    n <- max(length(iot_s), length(iot_ns))
    length(iot_s) <- n
    length(iot_ns) <- n
    n <- max(length(iol_s), length(iol_ns))
    length(iol_s) <- n
    length(iol_ns) <- n	
    n <- max(length(insites_s), length(insites_ns))
    length(insites_s) <- n
    length(insites_ns) <- n

	# Bind the arrays as columns
	ict<-cbind(ict_s,ict_ns)
	icl<-cbind(icl_s,icl_ns)
	iot<-cbind(iot_s,iot_ns)
	iol<-cbind(iol_s,iol_ns)
	insites<-cbind(insites_s,insites_ns)
	map$insites<-insites
	map$ict<-ict
	map$icl<-icl
	map$iot<-iot
	map$iol<-iol

	return(map)
}
