#!/usr/bin/env Rscript

library(data.table)

args=commandArgs(trailingOnly=TRUE)

system(paste0("grep ID=CDS ../Data/humangenome/Homo_sapiens.GRCh38.89.gff3 > tmp"))
cdsann <- fread("tmp") # extract all lines from .gff3 describing genes into a table
system("rm tmp")
cat("Coding region list from hg38 .gff3 file")
names(cdsann) <- c('chrom','source','type','start','end','score','strand','phase','attributes')
head(cdsann)

tsg_table <- fread("https://bioinfo.uth.edu/TSGene/All_down_exp_TSGs_pan-cancer.txt")

j=1;
tsg_inds<-rep(0,nrow(tsg_table)) # allocate memory for array of TSG indices in the list of genes
for (i in tsg_table$GeneName){

        tmp <- grep(paste0("Name=",i,";"),cdsann$attributes) # search for "Name=<gene name>" in the attributes column of the .gff3
        if (length(tmp)==0){ # if it was not found, place a NA for index of the current TSG
                tsg_inds[j]=NA;
	} else { # if found, store its index in the list
	        tsg_inds[j] <- tmp
        }
        j<-j+1 # increment tsg_inds array counter

        # print(grep(paste0("Name=",i,";"),cdsann$V9)) # print to screen the index of the current TSG (not necessary)
}

cdsann$istsg=0;
cdsann$istsg[tsg_inds]=1;

save(cdsann, file=args[1])
