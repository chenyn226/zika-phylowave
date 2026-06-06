# ZIKV evolutionary dynamics
## Data
### ZIKV sequence data
`thai_fp_brazil_sequences_align.fasta`:
This file contains the aligned ZIKV sequences used for analysis in this study. Sequence names are formatted as the original sequence accession ID, followed by an underscore (_) and the sequence collection date.

`sequence_collection_info.rds`:
This metadata file contains the original accession ID, updated sequence ID used in the FASTA file, sequence collection date, and country of collection for each sequence. 

### Time-resolved phylogenetic trees
`time_tree.treefile`:
The time-resolved phylogenetic tree for sequences from Thailand, French Polynesia and Brazil.

`thai_time_tree.treefile`:
The subtree used for lineage detection analysis, containing only Thai sequences and excluding six sequences as described in the manuscript draft.

### Lineage detection results
The output files generated from the lineage detection analysis, using the implementation provided in the accompanying code file. 

`dataset_with_nodes_phylowave.rds`

`sequence_collection_info.rds`

## Code

### Files loaded from phylowave codebase
_Lefrancq, Noémie, et al. "Learning the fitness dynamics of pathogens from phylogenies." Nature 637.8046 (2025): 683-690._

`2_1_Index_computation_20251129.R`

`2_2_Lineage_detection_20260127.R`

### Lineage detection algorithm implementation
`phylowave_zikv.R`

### Viral fitness estimation 
`fitness_estimation.R`

`Model_lineage_fitness_tstart.stan`


