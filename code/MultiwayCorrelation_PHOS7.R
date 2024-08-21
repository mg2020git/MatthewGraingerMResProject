

# Row: ASV
# Column: node purity, max depth, mean rel abun, mean pressure, betweenness, degree, clustering coeff, neighbourhood connectivity

# Load packages
library(dplyr)
library(psych)
library(GGally)
# Source function
source("clean_ASV_table.R")

# Importing data
ASV.table=read.table(file = "../data/seqtable_readyforanalysis.csv", sep="\t") # ASV table
sample_md <-read.table(file = "../data/metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv", sep="\t", header=TRUE) # Meta data to clean asv table
pressure_table=read.table(file = "../data/pressure_table.csv", sep=",", header = TRUE) # pressure table
network_metrics=read.table(file = "../data/combined_network_metrics.csv", sep = ",", header = TRUE) # network metrics table
# Importing RF importance data
rf_importance = read.table(file = "../results/Optimised_PHOS7_ASV/RFEVarImp_Class-log.mP7.norm_ASV_AllRepTrain.tsv")

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
# Get average relative abundance of every ASV
ASV.table <- t(ASV.table)
ASV.table <- as.data.frame(ASV.table)
ASV.table$mean_abundance <- rowMeans(ASV.table)
ASV.table$ASV <- rownames(ASV.table)  # Create a new column from row names
ASV.table <- ASV.table[, c("ASV", colnames(ASV.table))]  # Move the new column to the first position
ASV.table <- ASV.table %>% select(ASV, mean_abundance) # only keep asvs and their mean abundance


# Get average pressure of every ASV
pressure_table <- t(pressure_table) # transpose
pressure_table <- as.data.frame(pressure_table) # make back into data frame from matrix
colnames(pressure_table) <- pressure_table[1,] # making ifrst row into column names
pressure_table <- pressure_table[-1,] # getting rid of first row
asv_names_pressure <- row.names(pressure_table)
pressure_table <- as.data.frame(apply(pressure_table, 2, function(x) as.numeric(x))) # making data frame numeric
pressure_table$mean_pressure <- rowMeans(pressure_table) # calculating mean pressure of each asv
row.names(pressure_table) <- asv_names_pressure
pressure_table$ASV <- rownames(pressure_table) # making asvs into column
pressure_table <- pressure_table %>% select(ASV, mean_pressure) # only keeping asvs and their mean pressures

# Clean network metrics table
network_metrics <- network_metrics[, c(9, 1:8, 10:20)] # Reorder
network_metrics <- network_metrics[, c("name", "BetweennessCentrality", "ClusteringCoefficient", "Degree", "NeighborhoodConnectivity")]
colnames(network_metrics) <- c("ASV", "BetweennessCentrality", "ClusteringCoefficient", "Degree", "NeighborhoodConnectivity")

# Clean rf importance table
colnames(rf_importance) <- rf_importance[1,] # making ifrst row into column names
rf_importance <- rf_importance[-1,] # getting rid of first row
rf_importance <- rf_importance[, c("variable", "mse_increase","node_purity_increase","mean_min_depth", "times_a_root")]
colnames(rf_importance) <- c("ASV", "mse_increase","node_purity_increase","mean_min_depth", "times_a_root")
rf_importance <- subset(rf_importance, mse_increase >= 10^-15) # Getting rid of ASVs with low importance

# Merging data
correlation_df <- merge(rf_importance, ASV.table, by = "ASV")
correlation_df <- merge(correlation_df, pressure_table, by = "ASV")
correlation_df <- merge(correlation_df, network_metrics, by = "ASV")
correlation_df <- correlation_df[,-1]
correlation_df <- as.data.frame(apply(correlation_df, 2, function(x) as.numeric(x))) # making data frame numeric
# Removing: time_a_root, NeighborhoodConnectivity
correlation_df <- correlation_df[, -c(4, 10)]

# Creating correlation plot using Spearmans
#ggpairs(correlation_df, 
#        upper = list(continuous = wrap("cor", method = "spearman")), # Display Spearman correlation
#        lower = list(continuous = "points"), # Show scatter plots
#        diag = list(continuous = "densityDiag")) # Show density plots on the diagonal

# Correlation plot using Spearmans
pairs.panels(correlation_df, method = "spearman", ellipses = FALSE, stars=TRUE)
# Creating correlation plot using Pearsons no log
pairs.panels(correlation_df, ellipses = FALSE, stars=TRUE)

# Log normalising all but pressure
pressure <- correlation_df$mean_pressure
correlation_df <- correlation_df[,-5]
correlation_df <- log(correlation_df + 0.0000000001)
correlation_df$mean_pressure <- pressure
# Creating another correlation plot using Pearsons
pairs.panels(correlation_df, ellipses = FALSE, stars=TRUE)

# log nromalising all
#correlation_df <- correlation_df[,-8]
#pressure <- log(pressure + 0.0000000001)
#correlation_df$mean_pressure <- pressure
#pairs.panels(correlation_df, ellipses = FALSE, stars=TRUE)

# Correlation test of degrees
cor_test <- cor.test(correlation_df$Degree, correlation_df$mse_increase, method = "pearson")










