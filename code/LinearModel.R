


# LOADING LIBRARIES 
library(dplyr)
library(stringr)
library(ggplot2)
library(randomForest)
library(randomForestExplainer)
library(ggrepel)
library(reshape2)  # Or use tidyr's pivot_longer

#ASV TABLE
ASV.table=read.table(file = "../data/seqtable_readyforanalysis.csv", sep="\t")
# Metadata
sample_md <-read.table(file = "../data/metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv", sep="\t", header=TRUE)
# Function data
data.func <- read.csv(file = "../data/20151016_Functions_remainder.csv")
# Cluster data
cluster_data <- read.csv(file = "../data/max_tot_ext_network_table.tsv", sep = "\t")


# SOURCING FUNCTIONS
source("clean_ASV_table.R")

# CLEANING DATA - Removing experiment 4M, keeping only final samples, and only keeping samples present in the metadata
clean.data.list = clean_ASV_table(ASV.table ,sample_md , match_exp = TRUE, exclude_exp = c("4M"),
                                  nreads = 10000) # Outputs list containing cleaned asv table and meta data
ASV.table = clean.data.list$ASV.table
sample_md = clean.data.list$sample_md




##### BUILDING REGRESSION DATA #####
## Modifying ASV table
# Keeping only final communities
id.t7 = grep("7D", sample_md$Experiment) # Indices of rows whcih are final communities
samples.t7 = as.character(sample_md$sampleid[id.t7]) # Sampleids of final communities
xIn.tmp = ASV.table[samples.t7, ] # ASV table containing only final communities
# Transposing asv table for merging
xIn.tmp <- t(xIn.tmp) # Transposing ASV table so that there is one column with all ASVs
xIn.tmp <- as.data.frame(xIn.tmp) # Making back into data frame type object
xIn.tmp$ASV = rownames(xIn.tmp) # Making ASV column from row names
rownames(xIn.tmp) = seq(1,dim(xIn.tmp)[1]) # Row names are now numbers
xIn.tmp = xIn.tmp[,c(dim(xIn.tmp)[2],1:(dim(xIn.tmp)[2]-1))] # reorder so that the first column is ASV
# Preparing cluster data for merging
cluster_data = cluster_data[,1:2] # Getting rid of 'set' column
cluster_counts <- table(cluster_data$functionInk) # Sizes of each cluster

# ONLY INCLUDING CLUSTERS 36
cluster_data <- cluster_data[cluster_data$functionInk ==36, ]

# Merging
xIn.tmp = merge(xIn.tmp, cluster_data, by = "ASV", all.x = FALSE) # Merging wiht cluster data by ASV
xIn.tmp <- xIn.tmp[, -1] # Get rid of ASV column
xIn.tmp <- xIn.tmp %>%
  group_by(functionInk) %>%
  summarize(across(everything(), \(x) sum(x, na.rm = TRUE)), .groups = 'drop')
# Transposing back
xIn.tmp <- t(xIn.tmp) 
xIn.tmp <- as.data.frame(xIn.tmp) # Making back into data frame type object
xIn.tmp$sampleid = rownames(xIn.tmp) # create a new column for sampleid
xIn.tmp = xIn.tmp[,c(dim(xIn.tmp)[2],1:(dim(xIn.tmp)[2]-1))] # reorder so that the first column is sampleid
rownames(xIn.tmp) = seq(1,dim(xIn.tmp)[1]) # Row names are now numbers
col_names <- as.character(xIn.tmp[xIn.tmp$sampleid == 'functionInk', -1]) # Making the clusters into column names
names(xIn.tmp)[-1] <- col_names # Making the clusters into column names
xIn.tmp <- xIn.tmp[xIn.tmp$sampleid != 'functionInk', ] # Removing row that contains cluster names

## Modifying function data
# Changing units
data.func$mgCO2.7 <- data.func$mgCO2.7 * 1000 # Converting CO2 from miligr. to microgr.
names(data.func)[names(data.func) == "mgCO2.7"] <- "mCO2.7" # Changing column name accordingly
data.func$ATP7 <- data.func$ATP7 / 1000 # Changing nanomolar to micromolar
data.func$ATP14 <- data.func$ATP14 / 1000
# Normalising by number of cells
data.func$ATP7.norm <- data.func$ATP7 / data.func$CPM7
data.func$mG7.norm <- data.func$mG7 / data.func$CPM7
data.func$mN7.norm <- data.func$mN7 / data.func$CPM7
data.func$mX7.norm <- data.func$mX7 / data.func$CPM7
data.func$mP7.norm <- data.func$mP7 / data.func$CPM7
data.func$mCO2.7.norm <- data.func$mCO2.7 / data.func$CPM7
# Taking log of normalised functions and cell count
data.func$log.ATP7.norm <- log(data.func$ATP7.norm + 0.00001)
data.func$log.mG7.norm <- log(data.func$mG7.norm + 0.00001)
data.func$log.mN7.norm <- log(data.func$mN7.norm + 0.00001)
data.func$log.mX7.norm <- log(data.func$mX7.norm + 0.00001)
data.func$log.mP7.norm <- log(data.func$mP7.norm + 0.00001)
data.func$log.mCO2.7.norm <- log(data.func$mCO2.7.norm + 0.00001)
# data.func.x has only the chosen function and condition e.g. ATP7
data.func.x <- data.func
data.func.x$sampleid <- paste0(data.func$Community, ".", data.func$Replicate) # adding sampleid column to data.func.x

# Merging
DATA = merge(xIn.tmp, data.func.x, by = "sampleid", all.x = FALSE) # Merging wiht cluster data by ASV
DATA <- DATA[, c("36", "log.ATP7.norm", "log.mG7.norm", "log.mN7.norm", "log.mX7.norm", "log.mP7.norm", "log.mCO2.7.norm")]
names(DATA)[1] <- "Abundance"



################################### MAKE PLOT

# lienar models
atplm <- lm(log.ATP7.norm ~ Abundance, data = DATA)
xylolm <- lm(DATA$log.mX7.norm ~ DATA$Abundance, data = DATA)
gluclm <- lm(DATA$log.mG7.norm ~ DATA$Abundance, data = DATA)
chitlm <- lm(DATA$log.mN7.norm ~ DATA$Abundance, data = DATA)
phoslm <- lm(DATA$log.mP7.norm ~ DATA$Abundance, data = DATA)
colm <- lm(DATA$log.mCO2.7.norm ~ DATA$Abundance, data = DATA)

DATA2 <- DATA[, 2:7]

# Base plot with the first model's data
plot(DATA$Abundance, DATA$log.ATP7.norm, pch = 16, col = "black", xlab = "Abundance", ylab = "Function", ylim = range(DATA2))

# Add regression lines for each model
abline(atplm, col = "red", lwd = 2)
abline(xylolm, col = "blue", lwd = 2)
abline(gluclm, col = "green", lwd = 2)
abline(chitlm, col = "purple", lwd = 2)
abline(phoslm, col = "orange", lwd = 2)
abline(colm, col = "brown", lwd = 2)

# Optionally, add a legend to identify the models
legend("topleft", legend = c("ATP", "Xylosidase", "Glucosidase", "Chitinase", "Phosphatase", "Carbon dioxide"),
       col = c("red", "blue", "green", "purple", "orange", "brown"), lty = 1, lwd = 2)





# Function to create a color with transparency
color_with_transparency <- function(color, alpha) {
  rgb(t(col2rgb(color) / 255), alpha = alpha, maxColorValue = 1)
}

# Plot with transparent points
plot(DATA$Abundance, DATA$log.ATP7.norm, pch = 16, col = color_with_transparency("black", 0.25), 
     xlab = "Abundance", ylab = "Function", ylim = range(DATA2))

# Add regression lines for each model
abline(atplm, col = "red", lwd = 2)
abline(xylolm, col = "blue", lwd = 2)
abline(gluclm, col = "green", lwd = 2)
abline(chitlm, col = "purple", lwd = 2)
abline(phoslm, col = "orange", lwd = 2)
abline(colm, col = "brown", lwd = 2)

# Optionally, add a legend to identify the models
legend("topleft", legend = c("ATP", "Xylosidase", "Glucosidase", "Chitinase", "Phosphatase", "Carbon dioxide"),
       col = c("red", "blue", "green", "purple", "orange", "brown"), lty = 1, lwd = 2)




