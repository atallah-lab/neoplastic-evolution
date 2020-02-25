birthrate <- function(nd_het, np_het, nd_hom, np_hom, sld, slp, spd, spp) { # cell cycles per time step per cell 
    return(((((1+sld)^nd_het)*((1+spd)^nd_hom))/(((1+slp)^np_het)*((1+spp)^np_hom))))
}

get_mu_i <- function(B, mu, tau, genTime) {return(mu*B*tau*genTime)} # insertions per time step per cell

get_nmut <- function(ncells, mu_i) {return (rpois(1,ncells*mu_i))} # number of insertions in clone

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
update_mcount_f <- function(nd_het,np_het,nd_hom,np_hom,typemut,typegene,typemut2,targettypemut2) {
    
    for (ii in 1:length(typemut)) {
        if (typemut2[ii]!=targettypemut2) {
            next
        }
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
update_mcount_m <- function(nd_het,np_het,nd_hom,np_hom,typemut,typegene,typemut2,targettypemut2) {
    
    for (ii in 1:length(typemut)) {
        if (typemut2[ii]!=targettypemut2) {
            next
        }
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
sompop <- function(N0, mu1, mu2, tau, NT, sld, slp, spd, spp, gender, inGene, geneListFile, nclones, logpath) {
    
    if (gender=='male') {
        mu1 <- mu1*pd_exvsnon_m_l1[1] # scale insertion rate by probability of exonic insertion
        gene_pd_l1 <- gene_pd_m_l1
        gene_pd_bg <- gene_pd_m_bg
        update_genes <- update_genes_m
        update_mcount <- update_mcount_m
    } else if (gender=='female') {
        mu1 <- mu1*pd_exvsnon_f_l1[1]
        gene_pd_l1 <- gene_pd_f_l1
        gene_pd_bg <- gene_pd_f_bg
        update_genes <- update_genes_f
        update_mcount <- update_mcount_f
    } else {stop('Argument gender must be \'male\' or \'female\'.')}
        
    # Initialize driver genes (assign their type as 1)
    geneList <- read.csv(geneListFile,header=F)$V1
    gene_pd_l1$type[gene_pd_l1$gene_id %in% geneList] <- 1
    gene_pd_bg$type[gene_pd_bg$gene_id %in% geneList] <- 1
    
    
    # Write to log file
    write(paste0('\n',Sys.time(),'\n',toString(length(which(gene_pd_l1$type==1))),' out of ',length(geneList),
                 ' driver genes found in annotation'),
          file=logpath,
          append=TRUE)
    write(paste0('Initial size: ',N0),file=logpath,append=TRUE)
    write(paste0('Insertion rate: ',mu1),file=logpath,append=TRUE)
    write(paste0('Background rate: ',mu2),file=logpath,append=TRUE)
    write(paste0('Number time steps: ',NT),file=logpath,append=TRUE)
    write(paste0('Gender: ',gender),file=logpath,append=TRUE)
    write(paste0('Selection coefficients (sD, sP, sd, sp): ',spd,', ',spp,', ',sld,', ',slp),file=logpath,append=TRUE)
    write(paste0('Driver gene list file: ',geneListFile),file=logpath,append=TRUE)
    
    # Allocate population
    Pop <- data.table(ncells=rep(0,nclones),
                      B=rep(0,nclones),
                      mu1_i=rep(0,nclones),
                      mu2_i=rep(0,nclones),
                      nd_het=rep(0,nclones),
                      np_het=rep(0,nclones),
                      nd_hom=rep(0,nclones),
                      np_hom=rep(0,nclones),
                      nd_het_bg=rep(0,nclones),
                      np_het_bg=rep(0,nclones),
                      nd_hom_bg=rep(0,nclones),
                      np_hom_bg=rep(0,nclones),
                      genes_het=rep(list(''),nclones),
                      genes_hom=rep(list(''),nclones),
                      genes_new=rep(list(''),nclones),
                      new_types=rep(list(),nclones),
                      mut_types=rep(list(),nclones))
    
    # Initialize population with no mutations
    Pop[1,c('ncells','nd_het','np_het','nd_hom','np_hom','nd_het_bg','np_het_bg','nd_hom_bg','np_hom_bg','genes_het','genes_hom','genes_new','new_types','mut_types'):=
         list(c(N0),
         c(0),
         c(0),
         c(0),
         c(0),
         c(0),
         c(0),
         c(0),
         c(0),
         list(c('')),
         list(c('')),
         list(c('')),
         list(c()),
         list(c()))]
    
            
    # Assign initial birth and insertion rates
    Pop[1:2, B := mapply(birthrate, nd_het, np_het, nd_hom, np_hom, sld, slp, spd, spp)]
    Pop[1:2, mu1_i := mapply(get_mu_i, B, mu1, tau, 1)]
    Pop[1:2, mu2_i := mapply(get_mu_i, B, mu2, tau, 1)]
        
    N <- rep(0,NT) # Allocate array for population size time series
    genTime <- rep(0,NT) # Allocate array for generation time factor
    genes <- character(nrow(gene_pd_l1))
    #e <- exp(1) # define Euler's number (for log death rate)
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
        #D <- log(1 + (e-1)*N[ii]/N0)*tau*genTime[ii] # Log death rate function
        
        if (ii==NT) {
            break
        }
        
        #### Simulate mutations
        nins <- sum(unlist(mapply(get_nmut,Pop$ncells[clog],Pop$mu1_i[clog],SIMPLIFY=FALSE))) # Get number of exonic insertions
        nbg <- sum(unlist(mapply(get_nmut,Pop$ncells[clog],Pop$mu2_i[clog],SIMPLIFY=FALSE))) # Get number of null background mutations
        nmut <- nins + nbg
        
        if (nmut > 0) {
            
            rownew <- which(Pop$ncells==0)[1] # Find first row of the data table with ncells==0
            
            # Sample cells for mutations, with replacement, with probability determined by mu
            cellst <- table(sample(1:sum(Pop$ncells[clog]), nmut, replace=TRUE, prob=rep(Pop$mu_i[clog], Pop$ncells[clog])))
            # Get clone ID for each cell with mutations
            clones <- rep(which(clog), Pop$ncells[clog])[as.integer(names(cellst))]
            # Create a table of the number of mutated cells per clone
            ctab <- table(clones)
            # Get the row IDs of the clones
            cids <- as.integer(names(ctab))
            # Remove cells from sampled clones
            set(Pop,cids,1L,Pop[cids,1L] - as.integer(ctab))
            
            # Sample genes for mutations, with replacement
            gene_ids <- sample(1:nrow(gene_pd_l1),nins,replace=TRUE,prob=gene_pd_l1$p)
            gene_ids_bg <- sample(1:nrow(gene_pd_bg),nbg,replace=TRUE,prob=gene_pd_bg$p)
            gene_ids <- append(gene_ids, gene_ids_bg)
            mutation_types <- c(rep(1,nins), rep(2,nbg))
            # Shuffle insertion and background mutations
            randinds <- sample(1:length(gene_ids),length(gene_ids),replace=FALSE)
            gene_ids <- gene_ids[randinds]
            mutation_types <- mutation_types[randinds]
            gene_list <- gene_pd_l1$gene_id[gene_ids]
            genes <- append(genes,gene_list)
            gene_types <- gene_pd_l1$type[gene_ids]
                                    
            # Get list of genes for each cell
            tmp1 <- rep(list(),length(cellst))
            tmp2 <- rep(list(),length(cellst))
            tmp3 <- rep(list(),length(cellst))
            for (jj in 1:length(cellst)) {
                tmp1[[jj]] <- head(gene_list,as.integer(cellst)[jj])
                tmp2[[jj]] <- head(gene_types,as.integer(cellst)[jj])
                tmp3[[jj]] <- head(mutation_types,as.integer(cellst)[jj])
                gene_list <- tail(gene_list,length(gene_list)-as.integer(cellst)[jj])
                gene_types <- tail(gene_types,length(gene_types)-as.integer(cellst)[jj])
                mutation_types <- tail(mutation_types,length(mutation_types)-as.integer(cellst)[jj])
            }
            gene_list <- tmp1
            gene_types <- tmp2
            mutation_types <- tmp3
                        
            new_inds <- rownew:(rownew+length(clones)-1)
            Pop[new_inds, 
                c("ncells",
                  "nd_het","np_het","nd_hom","np_hom",
                  "nd_het_bg","np_het_bg","nd_hom_bg","np_hom_bg",
                  "genes_het","genes_hom","genes_new","new_types","mut_types"):=list(1, 
                        Pop$nd_het[clones], 
                        Pop$np_het[clones], 
                        Pop$nd_hom[clones], 
                        Pop$np_hom[clones], 
                        Pop$nd_het_bg[clones], 
                        Pop$np_het_bg[clones], 
                        Pop$nd_hom_bg[clones], 
                        Pop$np_hom_bg[clones], 
                        Pop$genes_het[clones],
                        Pop$genes_hom[clones],
                        gene_list,
                        gene_types,
                        mutation_types)]
                        
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
                           Pop$new_types[new_inds],
                           Pop$mut_types[new_inds],
                           1)) # 1 indicates exonic insertion
            Pop[new_inds,c("nd_het","np_het","nd_hom","np_hom"):=
                list(unlist(tmp2[,1]),
                     unlist(tmp2[,2]),
                     unlist(tmp2[,3]),
                     unlist(tmp2[,4]))]
            
            # Update background null mutation counts
            tmp2<-t(mapply(update_mcount,
                           Pop$nd_het_bg[new_inds],
                           Pop$np_het_bg[new_inds],
                           Pop$nd_hom_bg[new_inds],
                           Pop$np_hom_bg[new_inds],
                           tmp1[,4],
                           Pop$new_types[new_inds],
                           Pop$mut_types[new_inds],
                           2)) # 2 indicates background mutation
            Pop[new_inds,c("nd_het_bg","np_het_bg","nd_hom_bg","np_hom_bg"):=
                list(unlist(tmp2[,1]),
                     unlist(tmp2[,2]),
                     unlist(tmp2[,3]),
                     unlist(tmp2[,4]))]
            
            # Update birth and insertion rates
            Pop[new_inds, B := mapply(birthrate, nd_het+nd_het_bg, np_het+np_het_bg, nd_hom+nd_hom_bg, np_hom+np_hom_bg, sld, slp, spd, spp)]
            Pop[new_inds, mu1_i := mapply(get_mu_i, B, mu1, tau, genTime[ii])]
            Pop[new_inds, mu2_i := mapply(get_mu_i, B, mu2, tau, genTime[ii])]
            
        }
        
        clog <- Pop$ncells>0 # Get logical array for indices of still active (# cells >0) clones
        
        #### Simulate births
        # Sample total number of births in population
        nbirths <- rpois(1, sum(Pop$ncells*Pop$B)*tau*genTime[ii])  
        # Sample cells with replacement for having a birth, with probability determined by their birth rate
        cells <- sample(1:sum(Pop$ncells[clog]), nbirths, replace=TRUE, prob=rep(Pop$B[clog], Pop$ncells[clog]))
        # Create a table of the number of births per clone
        ctab <- table(rep(which(clog), Pop$ncells[clog])[cells])
        cids <- as.integer(names(ctab)) # Get row ids of the clones
        set(Pop,cids,1L,Pop[cids,1L] + as.integer(ctab)) # Add cells to sampled clones
        
        #### Simulate deaths
        ndeaths <- min(rpois(1, N[ii]*D),sum(Pop$ncells)) # Sample total number of deaths in population
        cells <- sample(1:sum(Pop$ncells[clog]), ndeaths, replace=FALSE) # Sample cells w/o replacement for death
        ctab <- table(rep(which(clog), Pop$ncells[clog])[cells]) # Create a table of the number of deaths per clone
        cids <- as.integer(names(ctab))
        set(Pop,cids,1L,Pop[cids,1L] - as.integer(ctab)) # Subtract cells from sampled clones
        
        # Order population
        Pop <- Pop[order(Pop$ncells,decreasing=TRUE),] # Order data.table by ncells
                
    }
    print(proc.time() - ptm)
    genes <- genes[!is.na(genes)]
    Pop <- Pop[,1:14]
    return(list(Pop,N,genes,genTime))

}
