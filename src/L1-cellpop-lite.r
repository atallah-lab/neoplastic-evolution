#!/usr/bin/env Rscript

# Arguments:
#
# (1) - run type ('single', 'batch', 'endless', see below)
# (2) - number of timesteps
# (3) - output file name
# 
# 'single' runs the simulation once for the specified number of timesteps and saves the result
# 'batch' runs the simulation a number of times (see batch section) while varying parameters, and saves results separately
# 'endless' runs the simulation endlessly (until manual termination) and saves the results to the same file after each timestep 

#--- Load libraries and necessary data files
######################################################################################
library(data.tree)
library(data.table)
load('../data/mutation_type_pd.rda')

args<-commandArgs(trailingOnly=TRUE)

#--- Define functions
######################################################################################
birthrate <- function(node, type, sd, sp) {
        nd <- 0; np <- 0
        if (type==1)     { b <- ((1+sd)^(node$nd+1))/((1+sp)^(node$np  )); nd <- 1 }
        else if (type==0){ b <- ((1+sd)^(node$nd  ))/((1+sp)^(node$np+1)); np <- 1 }
        else             { b <- node$B }
        return(list(b,nd,np))
            
}

maybeTranspose <- function(node, sd, sp) {

    if (node$B==0){ # if the birth rate of the clone is zero, skip the node
        return()
    }
    
    # increase the number of cells by the existing number * the birth rate
    nc <- tail(node$ncells,n=1) + round(tail(node$ncells,n=1)*node$B)
    
    # sample from binomial distribution for number of transpositions
    if (nc < 4.2e9) {ntrans <- rbinom(1,nc,mu)} # rbinom() fails for large n
    else {ntrans <- nc*mu} # If n is too large, use the expected number of events (mean of distribution)
    if (ntrans > 0) {
        
        nc <- nc-ntrans
        types <- sample(sample(1:-1,ntrans,replace=TRUE,prob=pdfinal)) # sample mutation types from pdfinal
        for (i in 1:ntrans) {
            l<<-l+1
            tmp <- birthrate(node,types[i],sd,sp)
            node$AddChild(l, ncells=1, B=tmp[[1]], nd=node$nd+tmp[[2]], np=node$np+tmp[[3]])
        }
        
    }   
    node$ncells <- append(node$ncells,nc)

}



#--- Set simulation parameters
######################################################################################

N0 <- 1           # Initial number of cells in root clone
B0 <- 1
mu <- 0.2         # Probability of transposition / timestep of a single cell

NT <- args[2]     # Number of time steps

#--- Generate clone tree
######################################################################################

if (args[1]=='batch') {
    sv <- c(0.000, 0.001, 0.010, 0.100, 1.000) # array of selection strengths over 4 orders of magnitude
	nrun <- 0
	for (sdi in 1:5) {
		for (spi in 1:5) {
			
		    l<-1 # Clone counter
		    N <- N0 # total population size

		    CellPop <- Node$new(1)
            CellPop$ncells <- c(N0) # Set initial number of cells of clone
            CellPop$B <- B0 # Set initial birth rate of clone
            CellPop$np <- 0 # Set initial number of drivier mutations
            CellPop$nd <- 0 # Set initial number of passenger mutations

            ptm <- proc.time()
            for (i in 2:NT) {

                t <- Traverse(CellPop,traversal='pre-order',filterFun=function(x) tail(x$ncells,n=1) > 0)    
                lapply(t,maybeTranspose,sv[sdi],sv[spi])

                N <- append(N,sum(vapply(CellPop$Get('ncells'),tail,n=1L,FUN.VALUE = numeric(1))))

                
            }
            (proc.time() - ptm)[3]


		    save(CellPop,N,file=paste0(args[3],nrun,".rda"))
		    rm(CellPop,N)
		    nrun <- nrun+1

		}

	}
} else if (args[1] == 'single') {

    l<-1 # Clone counter
    N <- N0 # total population size

    CellPop <- Node$new(1) # Initialize data.tree as single node
    CellPop$ncells <- c(N0) # Set initial number of cells of clone
    CellPop$B <- B0 # Set initial birth rate of clone
    CellPop$np <- 0 # Set initial number of drivier mutations
    CellPop$nd <- 0 # Set initial number of passenger mutations


    ptm <- proc.time()
    for (i in 2:NT) {

        t <- Traverse(CellPop,traversal='pre-order',filterFun=function(x) tail(x$ncells,n=1) > 0)    
        lapply(t,maybeTranspose,0.1,0.001)

        N <- append(N,sum(vapply(CellPop$Get('ncells'),tail,n=1L,FUN.VALUE = numeric(1))))

        
    }
    (proc.time() - ptm)[3]

	save(CellPop,N, file=paste0(args[3]))

} else if (args[1] == 'endless') {
    
    l<-1 # Clone counter
    N <- N0 # total population size

    CellPop <- Node$new(1) # Initialize data.tree as single node
    CellPop$ncells <- c(N0) # Set initial number of cells of clone
    CellPop$B <- B0 # Set initial birth rate of clone
    CellPop$np <- 0 # Set initial number of drivier mutations
    CellPop$nd <- 0 # Set initial number of passenger mutations
   
    while (1) {
	
	ptm <- proc.time()

        t <- Traverse(CellPop,traversal='pre-order',filterFun=function(x) tail(x$ncells,n=1) > 0)    
        lapply(t,maybeTranspose,0.1,0.001)

        N <<- append(N,sum(vapply(CellPop$Get('ncells'),tail,n=1L,FUN.VALUE = numeric(1))))   
            
	save(CellPop,N, file=paste0(args[3]))              

	print(proc.time()-ptm)

}
}




