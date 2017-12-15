#--- Load libraries and necessary data files, define variables
library(data.tree)
library(data.table)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
source("process_L1s.r")
load("./data/chrmpd.rda") # load chromosome probability distribution
load("./data/L1RankTable.rda")
load("./data/geneann.rda")
trpd <- read.table("./data/L1truncpd.csv",sep=",")
tdpd <- read.table("./data/L1transdpd.csv",sep=",")
for (i in names(Hsapiens)[1:24]){ # load all chromosome map files
	        load(paste0("./data/root_maps/",i,".rda"))
}
strdict<-c("+","-")
names(strdict)<-c(1,2)

#--- Define functions

rank_clone <- function(r, geneann, sites_chrm, sites_loci) {
	    
	gene_hits=0; # set counter to zero
	tsg_hits=0;
        for (i in 1:length(sites_chrm)) { # loop over chromosomes inserted into
		tmp=geneann[geneann$chrom==unique(sites_chrm)[i]] # reduce annotation table to entries for current chrom
            	chrmann_ntsg=tmp[tmp$istsg==0]
	        chrmann_tsg =tmp[tmp$istsg==1]
	        tmp = sites_loci[sites_chrm==unique(sites_chrm)[i]] # reduce insertion loci to entries for current chrom
		tmp_hits = between(tmp,chrmann_ntsg$start,chrmann_ntsg$end) # create logical for insertions, whether into non-tsg-gene or not
		gene_hits=gene_hits+length(which(tmp_hits==TRUE)) # count the number of non-tsg-gene insertions
		tmp_hits  = between(tmp,chrmann_tsg$start,chrmann_tsg$end) 
		tsg_hits =tsg_hits+length(which(tmp_hits==TRUE))
	}
        
        if (gene_hits > 0) {
		r=0
	} else if (tsg_hits > 0) {
	        r = r+tsg_hits*2; # TSG insertion doubles cell division rate
	}
	    
	return(r)
}

maybeTranspose <- function(node,tnum) {

	if (sample(x=c(0,1),1,prob=c(1-node$tp, node$tp))) {
		ncln=ncln+1
        	simout <- gen_sim(Hsapiens, 1)#round(runif(1,1,3)))
	        r_tmp <- rank_clone(node$r, geneann, simout$sites_chrm, simout$sites_loci)
	        node$AddChild(tnum)
		Set(node$children, r=r_tmp, tp=node$tp, l1s=node$l1s, sites=node$sites, ncells=1)
		        
	}

}

update_chrom_map <- function(chrnm,insites,ict,icl,iot,iol,sites_chrm,sites_loci,l1s) {

	chrloci = sites_loci[sites_chrm==chrnm] # Get the sites where insertions occurred in the chromosome
	chrl1s = l1s[sites_chrm==chrnm] # Get the L1s elements which were inserted

	for (i in 1:length(chrloci)) { # Loop over the simulated insertion points
		insites[which(is.na(insites))]<- -1 # Replace NA with -1
		indx <- insites>chrloci[i] # Get indices of target sites which lie downstream of the point
		insites[indx] <- insites[indx] + width(chrl1s[i]) # Shift the target sites by the length of the L1
		l1_map <- mapSeq_SV(chrl1s[i]) # Map target sites in the L1
		l1_map$insites <- l1_map$insites + chrloci[i] # Convert L1 loci to chromosome loci
		insites <- rbind(insites,l1_map$insites) # Add target sites within L1 to chrom map
		ict <- rbind(ict,l1_map$ict)
		icl <- rbind(icl,l1_map$icl)
		iot <- rbind(iot,l1_map$iot)
		iol <- rbind(iol,l1_map$iol)
	}

	return(list(insites,ict,icl,iot,iol))

}

#--- Define parameters

#--- Set simulation parameters
ENifrc<- .1       # Fraction of endonuclease-independent (random) insertions
inPopSize <- 10   # Initial number of cells in root clone
inDivRate <- 1    # Initial division rate
intp <- 0.5       # Initial probability of transposition

NT <- 10           # Number of time iterations

#--- Clone tree creation
CellPop <- Node$new(1)
CellPop$ncells <- inPopSize
CellPop$r <- inDivRate
CellPop$tp <- intp
CellPop$l1s <- DNAStringSet()
CellPop$sites <- c()

for (i in 2:NT) {
	    
    CellPop$Do(maybeTranspose,i)
    CellPop$Do(function(node) node$ncells <- node$ncells*(2^node$r))
                   
}

save(CellPop,file="./data/CellPop_out.rda")
