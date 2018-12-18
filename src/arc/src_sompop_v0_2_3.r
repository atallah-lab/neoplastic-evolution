#--- birthrate()
birthrate <- function(nd_het, np_het, nd_hom, np_hom, sld, slp, spd, spp, tau) { # probability of birth per timestep   
    return(((((1+sld)^nd_het)*((1+spd)^nd_hom))/(((1+slp)^np_het)*((1+spp)^np_hom)))*tau)
}

#--- delta_ncells()
delta_ncells <- function(B, D, ncells) { # change in number of cells for a clone
    return (max(ncells + rbinom(1,ncells,min(B,1))-rbinom(1,ncells,min(D,1)),0))
}

#--- get_mu_i()
get_mu_i <- function(B, mu) {return(mu*B)} # mutation rate of clone i: proportional to birth rate

#--- get_nins()
get_nins <- function(ncells, mu_i) {return (rbinom(1,ncells,mu_i))} # number of insertions in current clone

#--- update_mcount()
update_mcount <- function(nd_het,np_het,nd_hom,np_hom,typemut,typegene) {
    
    if (typemut==2 & typegene==0) {return(list(nd_het,np_het+1,nd_hom,np_hom)) # If heterozygous passenger mutation (most likely)
    } else if (typemut==1) {return(c(nd_het,np_het,nd_hom,np_hom)) # If no-effect (redundant) mutation
    } else if (typemut==2 & typegene==1) {return(list(nd_het+1,np_het,nd_hom,np_hom)) # If heterozygous driver mutation 
    } else if (typemut==3 & typegene==0) {return(list(nd_het,np_het-1,nd_hom,np_hom+1)) # If homozygous passenger
    } else if (typemut==3 & typegene==1) {return(list(nd_het-1,np_het,nd_hom+1,np_hom))} # If homozygous driver
    
}

#--- update_genes_m()
# This function will update the gene lists for a given clone, for a male genome
# type 1 - no effect
# type 2 - heterozygous effect
# type 3 - homozygous effect
update_genes_m <- function(genes_het,genes_hom) {
    if (tail(genes_het,1) %in% genes_hom) { # If the last gene in list is in list of homoz. disrupted genes, discard it
        genes_het <- head(genes_het,length(genes_het)-1)
        type<-1
    } else {
        if (tail(genes_het,1) %in% head(genes_het,length(genes_het)-1)) { # If the last gene is in the list of heteroz. disrupted genes
            if (tail(genes_het,1) %in% xy_genes) { # If it's an X or Y chrom. gene
                genes_het <- head(genes_het,length(genes_het)-1)
                type<-1
            }
            else if (sample(c(0,1),1)) { # With 50% probability, assume it's the other copy, add the gene to the homoz. list, and discard it (all instances) from heteroz. list
                genes_hom <- append(genes_hom,tail(genes_het,1))
                genes_het <- genes_het[genes_het!=tail(genes_het,1)]
                type<-3
            } else { # Else, assume it hit already disrupted copy, and discard from list
                genes_het <- head(genes_het,length(genes_het)-1)
                type<-1
            }
        } else {
            if (tail(genes_het,1) %in% xy_genes) { # If it's an X or Y chrom. gene
                type<-3
            } else {
                type<-2
            }
        }
    }
      return(list(genes_het,genes_hom,type))
    
}

#--- update_genes_f()
# This function will update the gene lists for a given clone, for a female genome
# type 1 - no effect
# type 2 - heterozygous effect
# type 3 - homozygous effect
update_genes_f <- function(genes_het,genes_hom) {
    if (tail(genes_het,1) %in% genes_hom) { # If the last gene in list is in list of paired disrupted genes, discard it
        genes_het <- head(genes_het,length(genes_het)-1)
        type<-1
    } else {
        if (tail(genes_het,1) %in% head(genes_het,length(genes_het)-1)) { # If the last gene is in the list of lone disrupted genes
            if (sample(c(0,1),1)) { # With 50% probability, assume it's the other copy, add the gene to the paired list, and discard it (all instances) from lone list
                genes_hom <- append(genes_hom,tail(genes_het,1))
                genes_het <- genes_het[genes_het!=tail(genes_het,1)]
                type<-3
            } else { # Else, assume it hit already disrupted copy, and discard from list
                genes_het <- head(genes_het,length(genes_het)-1)
                type<-1
            }
        } else {
            type<-2
        }
    }
      return(list(genes_het,genes_hom,type))
    
}

#--- run_sim()
run_sim <- function(N0, mu, tau, NT, sld, slp, spd, spp, gender, driverGene, geneList, nclones, logpath) {

    if (gender=='male') {
        gene_pd <- gene_pd_m
        update_genes <- update_genes_m
    } else if (gender=='female') {
        gene_pd <- gene_pd_f
        update_genes <- update_genes_f
    } else {stop('Argument gender must be \'male\' or \'female\'.')}
    
    gene_pd$type[gene_pd$gene_sym %in% geneList] <- 1
    write(paste0(toString(length(which(gene_pd$type==1))),' out of ',length(geneList),' driver genes found in gene_pd.'),file=logpath,append=TRUE)
    
    # Allocate population
    Pop <- data.table(ncells=rep(0,nclones),
                      B=rep(0,nclones),
                      mu_i=rep(0,nclones),
                      nd_het=rep(0,nclones),
                      np_het=rep(0,nclones),
                      nd_hom=rep(0,nclones),
                      np_hom=rep(0,nclones),
                      genes_het=rep(list(''),nclones),
                      genes_hom=rep(list(''),nclones),
                      lasttype=rep(0,nclones))
    # Initialize population with a heterozygous driver mutation
#     Pop[1:2,c('ncells','nd_het','np_het','nd_hom','np_hom','genes_het','genes_hom','lasttype'):=list(c(N0-1,1),
#                                                                                  c(0,1),
#                                                                                  c(0,0),
#                                                                                  c(0,0),
#                                                                                  c(0,0),
#                                                                                  list(c(''),c(driverGene)),
#                                                                                  list(c(''),c('')),
#                                                                                  c(0,1))]
    
    # Initialize population with no mutations
    Pop[1,c('ncells','nd_het','np_het','nd_hom','np_hom','genes_het','genes_hom','lasttype'):=list(c(N0),
                                                                                 c(0),
                                                                                 c(0),
                                                                                 c(0),
                                                                                 c(0),
                                                                                 list(c('')),
                                                                                 list(c('')),
                                                                                 c(0))]
    
#     # Assign birth and insertion rates
    Pop[1:2, B := mapply(birthrate, nd_het, np_het, nd_hom, np_hom, sld, slp, spd, spp, tau)]
    Pop[1:2, mu_i := mapply(get_mu_i, B, mu)]
    # For if the population data table needs to be enlarged
    bkup_Pop <- data.table(ncells=rep(0,nclones),
                      B=rep(0,nclones),
                      mu_i=rep(0,nclones),
                      nd_het=rep(0,nclones),
                      np_het=rep(0,nclones),
                      nd_hom=rep(0,nclones),
                      np_hom=rep(0,nclones),
                      genes_het=rep(list(''),nclones),
                      genes_hom=rep(list(''),nclones),
                      lasttype=rep(0,nclones))
    
    N <- rep(0,NT) # Allocate array for population size time series
    genes <- character(nrow(gene_pd))
    write('Initialized...',file=logpath,append=TRUE)

    ptm <- proc.time()
    for (ii in 1:NT) { # Loop over time steps
        if(ii %in% c(round(NT/4),round(NT/4*2),round(NT/4*3),NT)) { # Print progress at 25% completed intervals
            write(paste0(toString(ii/NT*100),'% done | ',format((proc.time()-ptm)[1],nsmall=3),' (s)'),file=logpath,append=TRUE)            
        }
        
        N[ii] <- sum(Pop$ncells) # Get current number of cells
        if (N[ii]>=3*N0 || N[ii]<1) {break} # Simulation stops if population has grown by 3X or died
        D <- N[ii]*tau/N0        # Compute death rate
        
        clog <- Pop$ncells>0 # Get logical array for indices of active (# cells >0) clones
        
        nins <- sum(unlist(mapply(get_nins,Pop$ncells[clog],Pop$mu_i[clog],SIMPLIFY=FALSE))) # Get number of exonic insertions
        if (nins > 0) {
            
            rownew <- which(Pop$ncells==0)[1] # Find first row of the data table with ncells==0
            
            gene_ids <- sample(1:nrow(gene_pd),nins,replace=TRUE,prob=gene_pd$p)
            gene_list <- gene_pd$gene_sym[gene_ids]
            genes<-append(genes,gene_list)
            genetypes <- gene_pd$type[gene_ids]

            clonesWIns <- sample(rep(1:nclones,Pop$ncells), nins, replace=FALSE, prob=rep(Pop$mu_i,Pop$ncells)) # List clone id (row number) of each cell; sample without replacement
            ctab <- table(clonesWIns)
            cids <- as.integer(names(ctab)) # Get row ids of sampled clones
            set(Pop,cids,1L,Pop[cids,1L] - as.integer(ctab)) # Remove cells from sampled clones

            # Populate the new rows representing new clones
            if(rownew+nins-1>nrow(Pop)){ # Enlarge the population data object if needed
                Pop<-rbind(Pop,bkup_Pop)
                write('Increased size of pop. object',file=logpath,append=TRUE)
            }
            new_inds <- rownew:(rownew+nins-1) # Row indices of new rows            
            Pop[new_inds, c("ncells","nd_het","np_het","nd_hom","np_hom","genes_het","genes_hom","lasttype"):=list(1, 
                                                                                                Pop$nd_het[clonesWIns], 
                                                                                                Pop$np_het[clonesWIns], 
                                                                                                Pop$nd_hom[clonesWIns], 
                                                                                                Pop$np_hom[clonesWIns], 
                                                                                                mapply(append,Pop$genes_het[clonesWIns],gene_list,SIMPLIFY=FALSE),
                                                                                                Pop$genes_hom[clonesWIns],
                                                                                                genetypes)]
            # Update gene lists
            tmp1 <- t(mapply(update_genes,Pop$genes_het[new_inds],Pop$genes_hom[new_inds],SIMPLIFY=TRUE))
            Pop$genes_het[new_inds] <- tmp1[,1]
            Pop$genes_hom[new_inds] <- tmp1[,2]
            # Update insertion counts
            tmp2<-t(mapply(update_mcount,Pop$nd_het[new_inds],Pop$np_het[new_inds],Pop$nd_hom[new_inds],Pop$np_hom[new_inds],tmp1[,3],Pop$lasttype[new_inds]))
            Pop[new_inds,c("nd_het","np_het","nd_hom","np_hom"):=list(unlist(tmp2[,1]),unlist(tmp2[,2]),unlist(tmp2[,3]),unlist(tmp2[,4]))]
            # Update birth and insertion rates
            Pop[new_inds, B := mapply(birthrate, nd_het, np_het, nd_hom, np_hom, sld, slp, spd, spp, tau)]
            Pop[new_inds, mu_i := mapply(get_mu_i, B, mu)]
        }
        
        Pop[Pop$ncells>0, ncells:=mapply(delta_ncells, B, D, ncells)] # Update number of cells for all clones
        Pop <- Pop[order(Pop$ncells,decreasing=TRUE),] # Order data.table by ncells

    }
    print(proc.time() - ptm)
    genes <- genes[!is.na(genes)]
    Pop <- Pop[,1:9]
    return(list(Pop,N,genes))

}