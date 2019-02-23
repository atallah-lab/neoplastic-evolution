#!/usr/bin/env Rscript
library(optparse)
library(data.table)
library(parallel)

option_list = list(
        make_option(c("--drivers"), type="character", default="../data/gene_lists/drivers_lusc_filt.txt",
                help="List of driver genes: text file with single column of driver gene ENSEMBL IDs\n\t\t[default = %default (Lung Squamous Cell Carcinoma)]"),
        
        make_option(c("-g", "--gender"), type="character", default="female",
                help="Gender of host: 'male' or 'female' [default = %default]"),
        
        make_option(c("-n", "--initial_size"), type="integer", default=1e3,
                help="Number of cells in initial population [default = %default]"),
        
        make_option(c("-m", "--insertion_rate"), type="double", default=1,
                help="Average number of L1 insertions per cell cycle [default = %default]"),
        
        make_option(c("-D", "--driver_strength_hom"), type="double", default=0.1,
                help="Selection coefficient of homozygous driver mutations [default = %default]"),
        
        make_option(c("-P", "--passenger_strength_hom"), type="double", default=0.005,
                help="Selection coefficient of homozygous passenger mutations [default = %default]"),
        
        make_option(c("-d", "--driver_strength_het"), type="double", default=NULL,
                help="Selection coefficient of heterozygous driver mutations [default = sD/10]"),

        make_option(c("-p", "--passenger_strength_het"), type="double", default=NULL,
                help="Selection coefficient of heterozygous passenger mutations [default = sP/10]"),

        make_option(c("-t", "--number_timesteps"), type="integer", default=1e3,
                help="Number of time steps to simulate\n\t\t[default = %default]\n\t\t(Note: by default 1 time step = 1 cell generation)"),
        
        make_option(c("--tau"), type="double", default=1,
                help="Time resolution: number of generations simulated per timestep (can be < 1) [default = %default]"),

        make_option(c("-l", "--log_path"), type="character", default="./log.txt",
                help="Path to log file [default = %default]"),
    
        make_option(c("--parallel"), type="logical", default=FALSE,
                help="Whether to use parallel processing [default = %default]"),

        make_option(c("--initial_mut"), type="character", default="none",
                help="'none' - no mutations in initial population\n\t\t'single' - 1 gene is homozygously disrupted in a single cell of population (flag --initial_gene specifies this gene)\n\t\t'inherited' - 1 gene is heterozygously disrupted in all cells of initial population (flag --initial_gene specifies this gene)\n\t\t[default = %default]"),

        make_option(c("--initial_gene"), type="character", default="ENSG00000141510",
                help="Gene disrupted at initialization if using --im 'single' or 'inherited' [default = %default (TP53 gene)]"),

        make_option(c("-f", "--data_folder"), type="character", default="../data/genomes/hg38/",
                help="Path to directory storing data files for simulated genome [default = %default]"),

        make_option(c("-o", "--output_file"), type="character", default="../data/output/test_out.rda",
                help="Path to output data file [default = %default]")

)
opt = parse_args(OptionParser(option_list=option_list))

load(paste0(opt$data_folder,'/exonicvsnon_counts.rda')) # Data objects holding probability of exonic vs non-exonic insertion
load(paste0(opt$data_folder,'/gene_pd_exon.rda')) # Data objects holding probability of insertion in each gene
source('./src_sompop_v0_2_5.r')
xy_genes <- gene_pd_m$gene_id[gene_pd_m$chrom %in% c('X','Y')] # Getting list of genes on chroms. X and Y by symbol
maxNClones <- opt$initial_size*4;
if (is.null(opt$driver_strength_het)) {opt$driver_strength_het=opt$driver_strength_hom/10}
if (is.null(opt$passenger_strength_het)) {opt$passenger_strength_het=opt$passenger_strength_hom/10}

out <- sompop(
			opt$initial_size, # Initial population size
			opt$insertion_rate, # Mutation rate
			opt$tau, # Time resolution
			opt$number_timesteps, # Number of time steps
			opt$driver_strength_het, # Het. driver coefficient
			opt$passenger_strength_het, # Het. passenger coefficient
			opt$driver_strength_hom, # Hom. driver coefficient
			opt$passenger_strength_hom, # Hom. passenger coefficient
			opt$gender, # Gender
			opt$initial_gene, # Initially disrupted gene
			opt$drivers, # File path: List of driver genes
			maxNClones, # Maximum number of clones simulatable
			opt$log_path, # Log file path
            opt$initial_mut, # How to initialize mutations in population
            opt$parallel
)

if (!is.null(opt$output_file)) {
	Pop <- out[[1]]
	N <- out[[2]]
	mut_genes <- out[[3]]
	gen_time <- out[[4]]
	mu <- opt$insertion_rate
	N0 <- opt$initial_size
	sD <- opt$driver_strength_hom
	sP <- opt$passenger_strength_hom
	sd <- opt$driver_strength_het
	sp <- opt$passenger_strength_het
	save(Pop,N,mut_genes,gen_time,N0,mu,sD,sP,sd,sp,file=opt$output_file)
}














