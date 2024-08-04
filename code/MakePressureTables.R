

# Load packages
library(dplyr)
# Source function
source("clean_ASV_table.R")

# Import ASV table
ASV.table=read.table(file = "../data/seqtable_readyforanalysis.csv", sep="\t")
# Import interaction table
Interactions_table = read.table(file = "../data/network_data_no_network.tsv", sep="\t")
# Import metadata
sample_md <-read.table(file = "../data/metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv", sep="\t", header=TRUE)

# Clean interactions_table
Interactions_table <- Interactions_table[,c(1,2,3)] # Get rid of type column
Interactions_table <- Interactions_table %>% # Get rid of rows including partition variable
  filter(startsWith(V1, "A") & startsWith(V2, "A"))
colnames(Interactions_table) <- c("ASV_A", "ASV_B", "Interaction") # Changing column names

# Clean ASV table
nreads = 10000 # minimum number of reads to consider a sample
exclude_exp = c("4M") # A vector of characters with the experiments that should be excluded
match_exp = TRUE # Set to true if only starting communities that were resurrected should be included
output.label = "Time0D_7D_matched"
clean.data.list = clean_ASV_table(ASV.table ,sample_md , match_exp = TRUE, exclude_exp = c("4M"),
                                  nreads = 10000) # Outputs list containing cleaned asv table and meta data
ASV.table = clean.data.list$ASV.table # Cleaned asv table
sample_md = clean.data.list$sample_md # cleaned meta data
id.t7 = grep("7D", sample_md$Experiment) # Indices of rows which are final communities
samples.t7 = as.character(sample_md$sampleid[id.t7]) # Sampleids of final communities
ASV.table = ASV.table[samples.t7, ] # ASV table containing only final communities

# Data frame for pressure values
pressure_table <- matrix(0, nrow = nrow(ASV.table), ncol = ncol(ASV.table)) # Same dimensions as ASV table
pressure_table <- data.frame(pressure_table) # Making into data frame
colnames(pressure_table) <- colnames(ASV.table) # Same columns as asv table (each ASV is a column)
rownames(pressure_table) <- rownames(ASV.table) # Same rows as asv table (each sample is a row)

# Filling in pressure table
for (i in 1:nrow(Interactions_table)){ # For every interaction
  asv1 <- Interactions_table[i,1]
  asv2 <- Interactions_table[i,2]
  interaction <- Interactions_table[i,3]
  for (sample in rownames(ASV.table)) { # For every sample
    if (asv1 %in% colnames(ASV.table) & asv2 %in% colnames(ASV.table)){ # If both of the interacting ASVs are in the ASV table
      pressure_table[sample, asv1] <- pressure_table[sample, asv1] + (interaction*ASV.table[sample,asv2])
      pressure_table[sample, asv2] <- pressure_table[sample, asv2] + (interaction*ASV.table[sample,asv1])
    }
  }
}

# Writing as .csv
write.csv(pressure_table, file = "../data/pressure_table.csv", row.names = TRUE)















