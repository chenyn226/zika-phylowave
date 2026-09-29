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
  - Step 2. Install `RStudio` or another `R`-compatible development environment.

### Lineage detection
- Software requirement
  - Install the required `R` packages: `phytools`, `data.table`.

- File requirement  
  The following files are required and are provided in this repository `lineage detection` folder.

  - Phylowave codebase files (_Lefrancq, Noémie, et al., Nature_): `2_1_Index_computation_20251129.R`, `2_2_Lineage_detection_20260127.R`
  - Input tree file: `thai_time_tree.treefile`
  - Analysis script: `phylowave_zikv.R`

- To run the lineage detection algorithm:
  - Step 1. Ensure that the required `R` packages are installed.
  - Step 2. Open `phylowave_zikv.R` in `RStudi`o and run the `R` script.
  - Step 3. The analysis will generate the following output files in the working directory: `split_phylowave.rds`, `dataset_with_nodes_phylowave.rds`.  
  Copies of these output files are also provided for reference and comparison.





