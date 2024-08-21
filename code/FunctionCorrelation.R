


##### First find out which #####




# LOADING LIBRARIES 
library(dplyr)
library(stringr)
library(ggplot2)
library(randomForest)
library(randomForestExplainer)
library(psych)
# Sourcing functions
source("clean_ASV_table.R")

# Importing function data
data.func <- read.csv(file = "../data/20151016_Functions_remainder.csv")
# ASV table
ASV.table=read.table(file = "../data/seqtable_readyforanalysis.csv", sep="\t")
# Metadata
sample_md <-read.table(file = "../data/metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv", sep="\t", header=TRUE)

# --- Reshaping data parameters
nreads = 10000 # minimum number of reads to consider a sample
exclude_exp = c("4M") # A vector of characters with the experiments that should be excluded
match_exp = TRUE # Set to true if only starting communities that were resurrected should be included
output.label = "Time0D_7D_matched" 

# CLEANING DATA - Removing experiment 4M, keeping only final samples, and only keeping samples present in the metadata
clean.data.list = clean_ASV_table(ASV.table ,sample_md , match_exp = TRUE, exclude_exp = c("4M"),
                                  nreads = 10000) # Outputs list containing cleaned asv table and meta data
ASV.table = clean.data.list$ASV.table
sample_md = clean.data.list$sample_md

##### BUILDING REGRESSION DATA #####
## Modifying ASV table
id.t7 = grep("7D", sample_md$Experiment) # Indices of rows whcih are final communities
samples.t7 = as.character(sample_md$sampleid[id.t7]) # Sampleids of final communities
xIn.tmp = ASV.table[samples.t7, ] # ASV table containing only final communities
# Adding column with the logarithm of the number of reads in each sample
xIn.tmp$counts = rowSums(xIn.tmp) # add the sum of the counts for each sample as an additional predictor
xIn.tmp$counts = log(xIn.tmp$counts) # log it to avoid having too high numbers
xIn.tmp$sampleid = rownames(xIn.tmp) # create a new column for sampleid
xIn = xIn.tmp # ASV table with only final samples, and a column for the log of the number of reads for each sample
rownames(xIn) = seq(1,dim(xIn)[1]) # Row names are now numbers
xIn = xIn[,c(dim(xIn)[2],1:(dim(xIn)[2]-1))] # reorder so that the first column is sampleid


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

# log but not normalised 
data.func$ATP7.LOG <- log(data.func$ATP7 + 0.00001)
data.func$mG7.LOG <- log(data.func$mG7 + 0.00001)
data.func$mN7.LOG <- log(data.func$mN7 + 0.00001)
data.func$mX7.LOG <- log(data.func$mX7 + 0.00001)
data.func$mP7.LOG <- log(data.func$mP7 + 0.00001)
data.func$mCO2.7.LOG <- log(data.func$mCO2.7 + 0.00001)
data.func$CPM7.LOG <- log(data.func$CPM7 + 0.00001)

data.func.x <- data.func

data.func.x$sampleid <- paste0(data.func$Community, ".", data.func$Replicate) # adding sampleid column to data.func.x
sample.id <- intersect(data.func.x$sampleid, xIn.tmp$sampleid)  # The sampleid of samples in the function data and the asv table
# IMPORTANT NOTE: Only 1035/1402 samples in xIn.tmp are present in data.func
data.func.x <- data.func.x[data.func.x$sampleid %in% sample.id, ] # Only keeping samples in both the asv table and the function data
xIn <- xIn[xIn$sampleid %in% sample.id, ] # Only keeping sampleids that are in both the function data and the asv table
xIn <- xIn[order(xIn$sampleid), ] # Sorting so in same order as function data
data.func.x <- data.func.x[order(data.func.x$sampleid), ] # sorting so in same order as asv table
xIn <- data.frame(xIn[, -which(names(xIn) == "sampleid")]) # get rid of sampleid column from asv table
data.func.x <- data.func.x[, -which(names(data.func.x) == "sampleid")] # get rid of sampleid column from function data
#yIn = as.numeric(data.func.x) # response vector

# Splitting into unaltered, nromalised, logged, and log normalised
functions_unaltered <- data.func.x[,c("ATP7","mG7", "mX7", "mN7", "mP7", "mCO2.7", "CPM7")]
functions_normalised <- data.func.x[,c("ATP7.norm","mG7.norm", "mX7.norm", "mN7.norm", "mP7.norm", "mCO2.7.norm", "CPM7")]
functions_log <- data.func.x[,c("ATP7.LOG","mG7.LOG", "mX7.LOG", "mN7.LOG", "mP7.LOG", "mCO2.7.LOG", "CPM7.LOG")]
functions_log_norm <- data.func.x[,c("log.ATP7.norm","log.mG7.norm", "log.mX7.norm", "log.mN7.norm", "log.mP7.norm", "log.mCO2.7.norm", "CPM7.LOG")]

################### Plots of distribution and correlation of functions ###################
# Unaltered functions spearman
pairs.panels(functions_unaltered, method = "spearman", ellipses = FALSE, stars=TRUE) # NOT NORMAL DISTR
# Normalised function spearman
pairs.panels(functions_normalised, method = "spearman", ellipses = FALSE, stars=TRUE) # NOT NORMAL DISTR

# Log functions pearson
pairs.panels(functions_log, method = "pearson", ellipses = FALSE, stars=TRUE)
# log normalised pearsons
pairs.panels(functions_log_norm, method = "pearson", ellipses = FALSE, stars=TRUE)



################ COMBINING FUNCTIONS ################
functions_log_norm_nocount <- data.func.x[,c("log.ATP7.norm","log.mG7.norm", "log.mX7.norm", "log.mN7.norm", "log.mP7.norm", "log.mCO2.7.norm")]
# Assuming `log_normalized_functions` is your data frame or matrix of log-normalized functions
pca_result <- prcomp(functions_log_norm_nocount, scale. = TRUE)
combined_function_metric <- pca_result$x[, 1]  # Extracting PC1 scores























