# ZIKV evolutionary dynamics
## Data
### Sequence data
`thai_fp_brazil_sequences_align.fasta`:
This file contains the aligned ZIKV sequences used for analysis in this study. Sequence names are formatted as the original sequence accession ID, followed by an underscore (_) and the sequence collection date.

`ZIKV_sequence_metadata.xlsx`: **Supplementary Data S1**

This metadata file contains the original accession ID, updated sequence ID used in the FASTA file, sequence collection date, country of collection, institutions generating the sequences, and factors related to sequence generation process. 

### Time-resolved phylogenetic trees

`time_tree.treefile`: The time-resolved phylogenetic tree for sequences from Thailand, French Polynesia and Brazil.

`thai_time_tree.treefile`: The subtree used for lineage detection analysis, containing only Thai sequences and excluding six sequences as described in the manuscript draft.


## Step by step instructions running the analysis
- Software requirements for all analysis
  - Step 1. Install `R`.
  - Step 2. Install `RStudio` or another R-compatible development environment.

### Lineage detection
- Software requirement
  - Install the required R packages: `phytools`, `data.table`.

- File requirement  
  The following files are required and are provided in this repository `lineage detection` folder.

  - Phylowave codebase files (_Lefrancq, Noémie, et al., Nature_): `2_1_Index_computation_20251129.R`, `2_2_Lineage_detection_20260127.R`
  - Input tree file: `thai_time_tree.treefile`
  - Analysis script: `phylowave_zikv.R`

- To run the lineage detection algorithm:
  - Step 1. Ensure that the required R packages are installed.
  - Step 2. Open `phylowave_zikv.R` in RStudio and run the R script.
  - Step 3. The analysis will generate the following output files in the working directory: `split_phylowave.rds`, `dataset_with_nodes_phylowave.rds`.
   
  Copies of these output files are also provided for reference and comparison.

### Viral fitness estimation
- Software requirement
  - Install the required R packages: `phytools`, `data.table`, `cmdstanr`, `stringr`

- File requirement  
  The following files are required and are provided in this repository `fitness estimation` folder.
  
  - Input tree file: `thai_time_tree.treefile`
  - Input Phylowave results file (output from running the lineage detection): `dataset_with_nodes_phylowave.rds`
  - Bayesian model stan file: `Model_lineage_fitness_tstart.stan`
  - Analysis script: `fitness_estimation.R`

- To run the fitness estimation:
  - Step 1. Ensure that the required R packages are installed. 
  - Step 2. Open `fitness_estimation.R` in RStudio and run the R script.
  - Step 3. The analysis will generate the following output file in the working directory:   
    The stan model output: `res_fitness.rds`  
    The posterior estimates of fitness: `figuredata_fitness_est.rds`, `figuredata_fitness_data.rds`, 
`relative_fitness_refgroup1.rds`, `relative_fitness_refparent.rds`, `relative_fitness_refgroup2.rds`, `time_varying_fitness.rds`

  Copies of these output files are also provided for reference and comparison.

### Figure generation

- Software requirement
  - Install the required R packages: `ape`, `readxl`, `data.table`, `treeio`, `ggtree`, `ggplot2`, `ggpubr`, `colorspace`, `ggh4x`, `stringr`, `cowplot`

- File requirement  
  The following files are required and are provided in this repository `figures` folder.
  
  - Input tree file: `time_tree.treefile`, `thai_time_tree.treefile`
  - Sequence meta data: `ZIKV_sequence_metadata.xlsx`
  - Lineage detection output: `dataset_with_nodes_phylowave.rds`
  - Fitness estimation output: `figuredata_fitness_est.rds`, `figuredata_fitness_data.rds`, 
`relative_fitness_refparent.rds`, `relative_fitness_refgroup2.rds`, `time_varying_fitness.rds`
  - R files recording lineage-defining substitutions: `aa_mutation_table.rds`, `utr_mutation_table.rds`, `aa_table_columnname.rds`, `utr_table_columnname.rds`
  - Phylowave codebase files (_Lefrancq, Noémie, et al., Nature_): `2_1_Index_computation_20251129.R`
  - Analysis script: `figures.R`

- To generate the figures in main manuscript:
  - Step 1. Ensure that the required R packages are installed. 
  - Step 2. Open `figures.R` in RStudio and run the R script.
  - Step 3. The analysis will generate the 3 figures in main manuscript.










