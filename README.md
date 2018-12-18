# sim-develop
An R project developed to simulate transposable element activity in population dynamics.

## Source Code

#### ./src/run_v2_5.r
* An R script for running sompop_v0_2_5 simulation.

  
```{r}
Usage: ./run_v2_5.r [options]


Options:
	-d DRIVERS, --drivers=DRIVERS
		List of driver genes: text file with single column of driver gene ENSEMBL IDs
		[default = ../data/gene_lists/drivers_lusc_filt.txt (Lung Squamous Cell Carcinoma)]

	-g GENDER, --gender=GENDER
		Gender of host: 'male' or 'female' [default = female]

	-n INITIAL_SIZE, --initial_size=INITIAL_SIZE
		Number of cells in initial population [default = 1000]

	-mu INSERTION_RATE, --insertion_rate=INSERTION_RATE
		Average number of L1 insertions per cell cycle [default = 1]

	-sD DRIVER_STRENGTH_HOM, --driver_strength_hom=DRIVER_STRENGTH_HOM
		Selection coefficient of homozygous driver mutations [default = 0.1]

	-sP PASSENGER_STRENGTH_HOM, --passenger_strength_hom=PASSENGER_STRENGTH_HOM
		Selection coefficient of homozygous passenger mutations [default = 0.005]

	-sd DRIVER_STRENGTH_HET, --driver_strength_het=DRIVER_STRENGTH_HET
		Selection coefficient of heterozygous driver mutations [default = sD/10]

	-sp PASSENGER_STRENGTH_HET, --passenger_strength_het=PASSENGER_STRENGTH_HET
		Selection coefficient of heterozygous passenger mutations [default = sP/10]

	-t NUMBER_TIMESTEPS, --number_timesteps=NUMBER_TIMESTEPS
		Number of time steps to simulate
		[default = 1000]
		(Note: by default 1 time step = 1 cell generation)

	--tau=TAU
		Time resolution: number of time steps each generation is divided into [default = 1]

	-l LOG_PATH, --log_path=LOG_PATH
		Path to log file [default = ./log.txt]

	--initial_mut=INITIAL_MUT
		'none' - no mutations in initial population
		'single' - 1 gene is heterozygously disrupted in a single cell of population (flag --initial_gene specifies this gene)
		'inherited' - 1 gene is heterozygously disrupted in all cells of initial population (flag --initial_gene specifies this gene)
		[default = none]

	--initial_gene=INITIAL_GENE
		Gene disrupted at initialization if using --im 'single' or 'inherited' [default = ENSG00000141510 (TP53 gene)]

	-df DATA_FOLDER, --data_folder=DATA_FOLDER
		Path to directory storing data files for simulated genome [default = ../data/genomes/hg38/]

	-o OUTPUT_FILE, --output_file=OUTPUT_FILE
		Path to output data file [default = ../data/output/test_out.rda]

	-h, --help
		Show this help message and exit
```

#### ./src/run_v2_5.ipynb
* Jupyter notebook for defining parameters and running sompop_v0_2_5. Single test runs and simple visualization of results, or batch runs can be performed.

#### ./src/src_sompop_v0_2_5.r
* An R file containing function definitions for simulation

## Data

#### ./data/genomes
* Directory for storing sub-folders for each host genome that can be simulated. Contains only hg38 by default.
* Sub-folders hold data files required for simulation

#### ./data/output
* Default folder for storing output data

#### ./data/gene_lists
* Folder storing text files with lists of driver genes by ENSEMBL ID for simulating different tumor types

#### ./data/genomes/hg38/drivers_XXXX.txt
* Text file containing a column of selected tumor suppressor genes for the tumor type XXXX (Lung squamous cell - LUSC, Glioblastoma multiforme - GBM, Colorectal adenocarcinoma - COREAD). Present driver gene lists were selected from IntOGen database: https://www.intogen.org/search

#### ./data/genomes/hg38/drivers_XXXX_filt.txt
* Same as above except some filters have been applied to the lists

#### ./data/genomes/hg38/exonicvsnon_counts.rda
* R data file containing the frequencies of L1 target sites of each snap-velcro category within exonic regions and outside exonic regions of the human genome. Exons were taken from Ensembl v94 human genome annotation.

#### ./data/genomes/hg38/gene_pd_exon.rda
* R data file containing the probability of exonic L1 insertion in each protein-coding gene of the human genome based on weighted L1 target site count.

#### ./data/genomes/hg38/gene_pdt_allsites.rda
* The same as gene_pdt.rda described above, though containing probabilities of insertion anywhere within the gene boundaries.
