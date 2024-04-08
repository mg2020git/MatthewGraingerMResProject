# Matthew Shaun Grainger's CMEE MRes Research Project Repository
This project investigates the reproducibility of microbial community assembly using network analysis. Raw data was sourced from Replaying the tape of ecology to domesticate wild microbiota https://doi.org/10.1101/2023.07.07.548163 (Pascual-Garcia et al., 2023), which is currently undergoing review. This data consists of communities sampled from 275 tree holes. Each community has 1 starting sample and 4 final samples from after a period of growth within a controlled environment. These samples were also grouped into community classes via their pairwise beta-diversity similarity. Here, we construct 4 co-occurrence networks - 1 across starting samples, 1 across starting samples whilst taking into account community classes, 1 across final samples, and 1 across final samples whilst taking into account community classes.

## Languages
- Python 3.10.12
- R 4.1.2
- Julia 1.10.2
- Perl 5.34.0

## Dependencies
### Dependencies within Python
- pandas 2.1.3

### Dependencies within R
- ggplot2 3.4.4

### FlashWeave
PLACEHOLDER

### functionInk
PLACEHOLDER

## Installation
PLACEHOLDER

## Usage
PLACEHOLDER

## Repository contents

### Data
- seqtable_readyforanalysis.csv : ASV table from https://doi.org/10.1101/2023.07.07.548163
- metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv : the associated metadata to the ASV table.
- taxa_wsp_readyforanalysis.csv : taxonomic information.
- starting_asv_table.tsv : ASV table for only the starting communities.
- starting_metadata.tsv : metadata for only the starting communities.
- final_asv_table.tsv : ASV table for only the final communities.
- final_metadata.tsv : metadata for only the final communities.
For each of the 4 networks (starting communities ignoring classes, starting communities accounting for classes, final communities ignoring classes, and final communities accounting for classes) the following data files can be found. Here, '#' represents either 'starting', 'starting_classes', 'final', or 'final_classes', respectively.
- #_network_output.edgelist : output for applying FlashWeave to the respecitve ASV table.
- #_network_data.tsv : FlashWeave output converted into functionInk input format.

### Code
- GraingerResearchProjectPipeline.ipynb : current iteration of my project's associated Jupyter notebook. Includes the steps of the project's pipeline.
- /functionInk : clone of the functionInk repository, available from: https://github.com/apascualgarcia/functionInk

### Results
- Clusters-NL_Average_StopStep-#_network_data.tsv : functionInk output describing the contents of each cluster.
- #_network_cytoscape_viz_1.png : Cytoscape visualisations of networks, created by following the functionInk vignette.
- #_Plot_PartitionDensityVsStep.pdf : plots of the partition density at each clustering step, during the running of functionInk.
- HistCompact-NL_Average_NoStop_#_network_data.tsv : simplified description of the clustering process. Output from functionInk.
- HistCompact-NL_Average_StopStep-#_network_data.tsv simplified description of the clustering process, stopping at the step at which the maximum partition density is reached. Output from functionInk.
- HistExtend-NL_Average_NoStop_#_network_data.tsv detailed description of the clustering process. Output from functionInk.
- HistExtend-NL_Average_StopStep_#_network_data.tsv detailed description of the clustering process, stopping at the step at which the maximum partition density is reached. Output from functionInk.
- Nodes-Similarities_#_network_data.tsv : functionInk output describing similarities between nodes.
- Partition-NL_Average-StopStep-#_network_data.tsv : functionInk output describing which cluster each ASV belongs to.

## Contact
Matthew Shaun Grainger
matthew.grainger20@imperial.ac.uk
