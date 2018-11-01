birthrate <- function(nd, np, sd, sp, tau) { # probability of birth per timestep
    return((((1+sd)^nd)/((1+sp)^np))*tau)
}

delta_ncells <- function(B, D, ncells) { # change in number of cells for a clone
    return (max(ncells + rbinom(1,ncells,min(B,1))-rbinom(1,ncells,min(D,1)),0))
}

get_mu_i <- function(B, mu) {return(mu*B)} # mutation rate of clone i: proportional to birth rate

get_nins <- function(ncells, mu_i) {return (rbinom(1,ncells,mu_i))} # number of insertions in current clone

# This function simulates the evolution of a cell population subject to L1 insertions over time.
# By default, a single cell is initialized with a driver insertion.
#
# Inputs:
#     N0 - initial population size (# cells)
#     mu - mutation rate (# L1 transpositions / birth)
#     tau - time resolution (# number of timesteps / birth of normal cell)
#     NT - number of timesteps to simulate
#     sd - driver selection strength (positive change in birth rate for each accumulated driver)
#     sp - passenger selection strength (negative change in birth rate for each accumulated passenger)
#     nclones - Buffer size of population data object; represents the max possible number of clones in the population
#     pd_mut - Discrete probability distribution of mutation types assumed to be in the order: driver, passenger, null, and to sum to 1
#
# Outputs:
#     Pop - data.table object containing a row for each clone in the final population, along
#           with columns for number of drivers and passengers, birth rate, mutation rate, number of cells
#     N - Array containing the population size over time

run_sim <- function(N0, mu, tau, NT, sd, sp, nclones, pd_mut, logpath) {

    Pop <- data.table(ncells=rep(0,nclones),nd=rep(0,nclones),np=rep(0,nclones),B=rep(0,nclones),mu_i=rep(0,nclones))
    # Populations are initialized with a driver mutation in a single cell
    Pop[1:2,c('ncells','nd','np'):=list(c(N0-1,1),c(0,1),c(0,0))]
    Pop[1:2,B := mapply(birthrate,nd,np,sd,sp,tau)]
    Pop[1:2,mu_i := mapply(get_mu_i, B, mu)]
    
    bkup_Pop <- data.table(ncells=rep(0,nclones),nd=rep(0,nclones),np=rep(0,nclones),B=rep(0,nclones),mu_i=rep(0,nclones))

    N <- rep(0,NT) # allocate array for population sizes
    e <- exp(1) # define Euler's number
    pd_mut <- pd_mut[1:2]/sum(pd_mut[1:2]) # make sure probability distribution is properly scaled
    write('Initialized...',file=logpath,append=TRUE)
    
    ptm <- proc.time()
    for (ii in 1:NT) {
        
        if(ii %in% c(round(NT/4),round(NT/4*2),round(NT/4*3),NT)) { # Print progress at 25% completed intervals
            write(paste0(toString(ii/NT*100),'% done | ',format((proc.time()-ptm)[1],nsmall=3),' (s)'),file=logpath,append=TRUE)            
        }   
        
        N[ii] <- sum(Pop$ncells) # get current number of cells
        if (N[ii]>=1e2*N0 || N[ii]<1) {break} # Simulation stops if population has grown by 3X or died
        D <- log(1 + (e-1)*N[ii]/N0)*tau        # compute death rate
        
        clog <- Pop$ncells>0 # get logical array for indices of active (# cells >0) clones
        
        nins <- sum(unlist(mapply(get_nins,Pop$ncells[clog],Pop$mu_i[clog],SIMPLIFY=FALSE))) # Get number of exonic insertions
        if (nins > 0) {
            
            rownew <- which(Pop$ncells==0)[1] # find first row of the data table with ncells==0

            types <- sample(1:2,nins,replace=TRUE,prob=pd_mut) # sample mutation types
            nmu <- length(types) # total number of passenger and drivers
            # List clone id of each cell; sample without replacement
            sampctr <- sample(rep(1:nclones,Pop$ncells),nmu,replace=FALSE,prob=rep(Pop$mu_i,Pop$ncells))
            ctab <- table(sampctr)
            cids <- as.integer(names(ctab)) # get row ids of sampled clones
            set(Pop,cids,1L,Pop[cids,1L] - as.integer(ctab)) # remove cells from sampled clones

            # Populate the new rows representing new clones
            if(rownew+nmu-1>nrow(Pop)){
                Pop<-rbind(Pop,bkup_Pop)
                write('Increased size of pop. object',file=logpath,append=TRUE)
            }
            Pop[rownew:(rownew+nmu-1), c("ncells","nd","np"):=list(
                1, 
                Pop[sampctr]$nd+((types==1)*1), 
                Pop[sampctr]$np+((types==2)*1))]
            
            Pop[rownew:(rownew+nmu-1), B := mapply(birthrate, nd, np, sd, sp, tau)]
            Pop[rownew:(rownew+nmu-1), mu_i := mapply(get_mu_i, B, mu)]
        }
        
        Pop[Pop$ncells>0, ncells:=mapply(delta_ncells, B, D, ncells)] # update number of cells for all clones
        Pop <- Pop[order(Pop$ncells,decreasing=TRUE),] # order data.table by ncells
    }
    print(proc.time() - ptm)

    return(list(Pop,N))

}