birthrate <- function(nd_het, np_het, nd_hom, np_hom, sld, slp, spd, spp) { # cell cycles per time step per cell 
    return(((((1+sld)^nd_het)*((1+spd)^nd_hom))/(((1+slp)^np_het)*((1+spp)^np_hom))))
}

delta_ncells <- function(B, D, ncells, tau, genTime) { # change in number of cells
    return (max(ncells + rpois(1,ncells*B*tau*genTime)-rpois(1,ncells*D),0))
}

get_mu_i <- function(B, mu, tau, genTime) {return(mu*B*tau*genTime)} # insertions per time step per cell

get_nins <- function(ncells, mu_i) {return (rpois(1,ncells*mu_i))} # number of insertions in clone

##### update_genes_f() #####
# This function will update the gene lists for a given clone, for a female genome
# types 1 - no effect
# types 2 - heterozygous effect
# types 3 - homozygous effect
update_genes_f <- function(new_genes,genes_het,genes_hom) {
    
    types <- rep(0,length(new_genes))
    jj <- 1
    for (ii in new_genes) {
        if (ii %in% genes_hom) { # If the new gene is in list of homozygous genes, ignore it
            types[jj]<-1
        } else {
            if (ii %in% genes_het) { # If the new gene is in the list of heterozygous genes
                # With 50% probability, assume it's the other copy, add the gene to the homoz. list, and 
                # discard it from heteroz. list
                if (runif(1,c(0,1))>=0.5) {
                    genes_hom <- append(genes_hom,ii)
                    genes_het <- genes_het[genes_het!=ii]
                    types[jj]<-3
                } else { # Else, assume it hit already disrupted copy, and discard from list
                    types[jj]<-1
                }
            } else {
                genes_het <- append(genes_het,ii)
                types[jj]<-2
            }
        }
        jj <- jj+1
    }
    return(list(genes_het,genes_hom,list(''),types))
    
}

##### update_genes_m() #####
# This function will update the gene lists for a given clone, for a male genome
# type 1 - no effect
# type 2 - heterozygous effect
# type 3 - homozygous effect
# type 4 - X/Y gene effect
update_genes_m <- function(new_genes,genes_het,genes_hom) {
    
    types <- rep(0,length(new_genes))
    jj <- 1
    for (ii in new_genes) {
        if (ii %in% genes_hom) { # If the new gene is in the list of homoz. genes, ignore it
            types[jj] <-1
        } else if (ii %in% xy_genes) { # If not in homoz. gene list, but in X/Y gene list
            genes_hom <- append(genes_hom,ii)
            types[jj] <- 4 # FIX ME
        } else {
            if (ii %in% genes_het) { # If the new gene is in the list of heteroz. genes
                # With 50% probability, assume it's the other allele, add the gene to the homoz. list, and 
                # discard it from heteroz. list
                if (runif(1,c(0,1))>=0.5) {
                    genes_hom <- append(genes_hom,ii)
                    genes_het <- genes_het[genes_het!=ii]
                    types[jj]<-3
                } else { # Else, assume it hit already disrupted copy, and discard from list
                    types[jj]<-1
                }
                
            } else {
                genes_het <- append(genes_het,ii)
                types[jj]<-2
            }
        }
        jj <- jj+1
    }
    return(list(genes_het,genes_hom,list(''),types))   
    
}    

##### update_mcount_f() #####
update_mcount_f <- function(nd_het,np_het,nd_hom,np_hom,typemut,typegene) {
    
    for (ii in 1:length(typemut)) {
        if (typemut[ii]==2 & typegene[ii]==0) { # If heterozygous passenger mutation (most likely)
            np_het <- np_het+1
        } else if (typemut[ii]==1) { # If no-effect (redundant) mutation
        } else if (typemut[ii]==2 & typegene[ii]==1) { # If heterozygous driver mutation 
            nd_het <- nd_het+1
        } else if (typemut[ii]==3 & typegene[ii]==0) { # If homozygous passenger
            np_het <- np_het-1
            np_hom <- np_hom+1
        } else if (typemut[ii]==3 & typegene[ii]==1) { # If homozygous driver
            nd_het <- nd_het-1
            nd_hom <- nd_hom+1
        }
    }
    return(list(nd_het,np_het,nd_hom,np_hom))
}

##### update_mcount_m() #####
update_mcount_m <- function(nd_het,np_het,nd_hom,np_hom,typemut,typegene) {
    
    for (ii in 1:length(typemut)) {
        if (typemut[ii]==2 & typegene[ii]==0) { # If heterozygous passenger mutation (most likely)
            np_het <- np_het+1
        } else if (typemut[ii]==1) { # If no-effect (redundant) mutation
        } else if (typemut[ii]==2 & typegene[ii]==1) { # If heterozygous driver mutation 
            nd_het <- nd_het+1
        } else if (typemut[ii]==3 & typegene[ii]==0) { # If homozygous passenger
            np_het <- np_het-1
            np_hom <- np_hom+1
        } else if (typemut[ii]==3 & typegene[ii]==1) { # If homozygous driver
            nd_het <- nd_het-1
            nd_hom <- nd_hom+1
        } else if (typemut[ii]==4 & typegene[ii]==0) { # If X/Y passenger
            np_hom <- np_hom+1
        } else if (typemut[ii]==4 & typegene[ii]==1) { # If X/Y driver
            nd_hom <- nd_hom+1
        }
    }
    return(list(nd_het,np_het,nd_hom,np_hom))
}

##### sompop() #####
sompop <- function(N0, mui, tau, NT, sld, slp, spd, spp, gender, inGene, geneListFile, nclones, logpath, mode, parallel) {

    if (parallel) {
        mapply = mcmapply
    }
    if (gender=='male') {
        mu <- mui*pd_exvsnon_m[1] # scale insertion rate by probability of exonic insertion
        gene_pd <- gene_pd_m
        update_genes <- update_genes_m
        update_mcount <- update_mcount_m
    } else if (gender=='female') {
        mu <- mui*pd_exvsnon_f[1]
        gene_pd <- gene_pd_f
        update_genes <- update_genes_f
        update_mcount <- update_mcount_f
    } else {stop('Argument gender must be \'male\' or \'female\'.')}
        
    # Initialize driver genes (assign their type as 1)
    geneList <- read.csv(geneListFile,header=F)$V1
    gene_pd$type[gene_pd$gene_id %in% geneList] <- 1
    
    # Write to log file
    write(paste0('\n',Sys.time(),'\n',toString(length(which(gene_pd$type==1))),' out of ',length(geneList),
                 ' driver genes found in annotation'),
          file=logpath,
          append=TRUE)
    write(paste0('Using parallel: ',parallel),file=logpath,append=TRUE)
    write(paste0('Initial size: ',N0),file=logpath,append=TRUE)
    write(paste0('Insertion rate: ',mui),file=logpath,append=TRUE)
    write(paste0('Number time steps: ',NT),file=logpath,append=TRUE)
    write(paste0('Gender: ',gender),file=logpath,append=TRUE)
    write(paste0('Selection coefficients (sD, sP, sd, sp): ',spd,', ',spp,', ',sld,', ',slp),file=logpath,append=TRUE)
    write(paste0('Driver gene list file: ',geneListFile),file=logpath,append=TRUE)
    write(paste0('Initial mutations: ',mode),file=logpath,append=TRUE)
    if (mode!='none') {
        # Check if the given initially mutated driver gene is in the gene_pd annotation
        if (length(which(gene_pd$gene_id==inGene))==0) {
            write(paste0('ERROR: Initial mutated gene: ',inGene,' not found. Quitting.'),file=logpath,append=TRUE)
            stop(paste0('ERROR: Initial mutated gene: ',inGene,' not found. Quitting.'))
        }
        write(paste0('Initial mutated gene: ',inGene,
                                    ' (', gene_pd$gene_sym[gene_pd$gene_id==inGene],')'),
                                    file=logpath,append=TRUE)
    }
    
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
                      genes_new=rep(list(''),nclones),
                      new_types=rep(list(),nclones))
    # Initialize population
    if (mode=='none') {
        # Initialize population with no mutations
        Pop[1,c('ncells','nd_het','np_het','nd_hom','np_hom','genes_het','genes_hom','genes_new','new_types'):=
             list(c(N0),
             c(0),
             c(0),
             c(0),
             c(0),
             list(c('')),
             list(c('')),
             list(c('')),
             list(c()))]
    } else if (mode=='single') {
        # Initialize population with a 1-cell homozygous driver mutation
        Pop[1:2,c('ncells','nd_het','np_het','nd_hom','np_hom','genes_het','genes_hom','genes_new','new_types'):=list(c(N0-1,1),
                                                                                     c(0,0),
                                                                                     c(0,0),
                                                                                     c(0,1),
                                                                                     c(0,0),
                                                                                     list(c(''),c('')),
                                                                                     list(c(''),c(inGene)),
                                                                                     list(c(''),c('')),
                                                                                     list(c(),c()))]
        
    } else if (mode=='inherited') {
        # Initialize all cells with heterozygous driver mutation
        Pop[1,c('ncells','nd_het','np_het','nd_hom','np_hom','genes_het','genes_hom','genes_new','new_types'):=
             list(c(N0),
             c(1),
             c(0),
             c(0),
             c(0),
             list(c(inGene)),
             list(c('')),
             list(c('')),
             list(c()))]
    } else {
        stop("Error: incorrect mode argument. Options are: 'none', 'single', 'inherited'")
    }
            
    # Assign birth and insertion rates
    Pop[1:2, B := mapply(birthrate, nd_het, np_het, nd_hom, np_hom, sld, slp, spd, spp)]
    Pop[1:2, mu_i := mapply(get_mu_i, B, mu, tau, 1)]
    
    N <- rep(0,NT) # Allocate array for population size time series
    genTime <- rep(0,NT) # Allocate array for generation time factor
    genes <- character(nrow(gene_pd))
    e <- exp(1) # define Euler's number
    write('Initialized...',file=logpath,append=TRUE)

    ptm <- proc.time()
    for (ii in 1:NT) { # Loop over time steps
        if(ii %in% c(round(NT/4),round(NT/4*2),round(NT/4*3),NT)) { # Print progress at 25% completed intervals
            write(paste0(toString(ii/NT*100),'% done | ',format((proc.time()-ptm)[3],nsmall=3),' (s)'),file=logpath,append=TRUE)            
        }
        
        clog <- Pop$ncells>0 # Get logical array for indices of active (# cells >0) clones
        
        N[ii] <- sum(Pop$ncells) # Get current number of cells
        if (N[ii]>=3*N0 || N[ii]<1) {break} # Simulation stops if population has grown by 3X or died
        genTime[ii] <- 1/(sum(Pop$B[clog]*Pop$ncells[clog])/sum(Pop$ncells[clog])) # Get generation length
        D <- N[ii]*tau*genTime[ii]/N0 # Linear death rate function (cell deaths per time-step per cell)
#         D <- log(1 + (e-1)*N[ii]/N0)*tau*genTime[ii] # Log death rate function

        nins <- sum(unlist(mapply(get_nins,Pop$ncells[clog],Pop$mu_i[clog],SIMPLIFY=FALSE))) # Get number of exonic insertions
        if (nins > 0) {
            
            rownew <- which(Pop$ncells==0)[1] # Find first row of the data table with ncells==0
            
            # Sample cells for mutations, with replacement, with probability determined by mu
            cellsWIns <- table(sample(1:sum(Pop$ncells[clog]), nins, replace=TRUE, prob=rep(Pop$mu_i[clog], Pop$ncells[clog])))
            # Get clone ID for each cell with mutation(s)
            clonesWIns <- rep(1:length(which(clog)), Pop$ncells[clog])[as.integer(names(cellsWIns))]
            ctab <- table(clonesWIns)
            cids <- as.integer(names(ctab)) # Get row ids of sampled clones
            set(Pop,cids,1L,Pop[cids,1L] - as.integer(ctab)) # Remove cells from sampled clones
            
            # Sample genes for mutations, with replacement, with probability determined by target site distribution
            gene_ids <- sample(1:nrow(gene_pd),nins,replace=TRUE,prob=gene_pd$p)
            gene_list <- gene_pd$gene_id[gene_ids]
            genes <- append(genes,gene_list)
            gene_types <- gene_pd$type[gene_ids]
            
            # Get list of genes for each cell
            tmp1 <- rep(list(),length(cellsWIns))
            tmp2 <- rep(list(),length(cellsWIns))
            for (jj in 1:length(cellsWIns)) {
                tmp1[[jj]] <- head(gene_list,as.integer(cellsWIns)[jj])
                tmp2[[jj]] <- head(gene_types,as.integer(cellsWIns)[jj])
                gene_list <- tail(gene_list,length(gene_list)-as.integer(cellsWIns)[jj])
                gene_types <- tail(gene_types,length(gene_types)-as.integer(cellsWIns)[jj])
            }
            gene_list <- tmp1
            gene_types <- tmp2
            
            new_inds <- rownew:(rownew+length(clonesWIns)-1)
            Pop[new_inds, 
                c("ncells","nd_het","np_het","nd_hom","np_hom","genes_het","genes_hom","genes_new","new_types"):=list(1, 
                                                                                        Pop$nd_het[clonesWIns], 
                                                                                        Pop$np_het[clonesWIns], 
                                                                                        Pop$nd_hom[clonesWIns], 
                                                                                        Pop$np_hom[clonesWIns], 
                                                                                        Pop$genes_het[clonesWIns],
                                                                                        Pop$genes_hom[clonesWIns],
                                                                                        gene_list,
                                                                                        gene_types)]
            
            # Update gene lists
            tmp1 <- t(mapply(update_genes,Pop$genes_new[new_inds],Pop$genes_het[new_inds],Pop$genes_hom[new_inds],SIMPLIFY=TRUE))
            Pop$genes_het[new_inds] <- tmp1[,1]
            Pop$genes_hom[new_inds] <- tmp1[,2]
            Pop$genes_new[new_inds] <- tmp1[,3]
            
            # Update insertion counts
            tmp2<-t(mapply(update_mcount,
                           Pop$nd_het[new_inds],
                           Pop$np_het[new_inds],
                           Pop$nd_hom[new_inds],
                           Pop$np_hom[new_inds],
                           tmp1[,4],
                           Pop$new_types[new_inds]))
            Pop[new_inds,c("nd_het","np_het","nd_hom","np_hom"):=
                list(unlist(tmp2[,1]),
                     unlist(tmp2[,2]),
                     unlist(tmp2[,3]),
                     unlist(tmp2[,4]))]
            
            # Update birth and insertion rates
            Pop[new_inds, B := mapply(birthrate, nd_het, np_het, nd_hom, np_hom, sld, slp, spd, spp)]
            Pop[new_inds, mu_i := mapply(get_mu_i, B, mu, tau, genTime[ii])]
        }
        
        Pop[Pop$ncells>0, ncells:=mapply(delta_ncells, B, D, ncells, tau, genTime[ii])] # Update number of cells for all clones
        Pop <- Pop[order(Pop$ncells,decreasing=TRUE),] # Order data.table by ncells
    }
    print(proc.time() - ptm)
    genes <- genes[!is.na(genes)]
    Pop <- Pop[,1:9]
    return(list(Pop,N,genes,genTime))

}
