###################################
#     random_forest_cfl.R  
##################################
# Author: Alberto Pascual-García
# Copyright (c)  Alberto Pascual-García,  2024
# Web:  apascualgarcia.github.io
# 
# Date: 2024-05-13
# Script Name: random_forest.R
# Modified by Matthew Shaun Grainger





##### GENERAL SETUP #####
rm(list = ls()) # Clear workspace

# LOADING LIBRARIES 
library(dplyr)
library(stringr)
library(randomForest)
library(foreach)
library(doParallel)

# SELECT FUNCTIONS AND CONDITIONS OF INTEREST TO OPTIMISE FOR
select.func="log.ATP" # determine the function you want to predict. Name should match those from
# the column name of the function data frame.
select.cond=c("7.norm") # determine the function conditions you want to predict. Names should 
#match those used to differentiate the different conditions within the function df.
predictor.unit = "ASV" # ASV or cluster as predictor - only for use in naming output files

# SOURCING FUNCTIONS
source("/rds/general/user/mg2020/home/researchprojecthpc/clean_ASV_table.R")

## READING INPUT FILES
file.ASV = "/rds/general/user/mg2020/home/researchprojecthpc/seqtable_readyforanalysis.csv"
file.taxa = "/rds/general/user/mg2020/home/researchprojecthpc/taxa_wsp_readyforanalysis.csv"
file.Meta = "/rds/general/user/mg2020/home/researchprojecthpc/metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv"
file.func = "/rds/general/user/mg2020/home/researchprojecthpc/20151016_Functions_remainder.csv"
# ASV Table
ASV.table=read.table(file = file.ASV, sep="\t")
# Metadata
sample_md <-read.table(file = file.Meta, sep="\t", header=TRUE)
# Function data
data.func <- read.csv(file = file.func)

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
names(data.func)[names(data.func) == "mgCO2.7"] <- "μgCO2.7" # Changing column name accordingly
data.func$ATP7 <- data.func$ATP7 / 1000 # Changing nanomolar to micromolar
data.func$ATP14 <- data.func$ATP14 / 1000
# Normalising by number of cells
data.func$ATP7.norm <- data.func$ATP7 / data.func$CPM7
data.func$mG7.norm <- data.func$mG7 / data.func$CPM7
data.func$mN7.norm <- data.func$mN7 / data.func$CPM7
data.func$mX7.norm <- data.func$mX7 / data.func$CPM7
data.func$mP7.norm <- data.func$mP7 / data.func$CPM7
data.func$μgCO2.7.norm <- data.func$μgCO2.7 / data.func$CPM7
# Taking log of normalised functions and cell count
data.func$log.ATP7.norm <- log(data.func$ATP7.norm + 0.00001)
data.func$log.mG7.norm <- log(data.func$mG7.norm + 0.00001)
data.func$log.mN7.norm <- log(data.func$mN7.norm + 0.00001)
data.func$log.mX7.norm <- log(data.func$mX7.norm + 0.00001)
data.func$log.mP7.norm <- log(data.func$mP7.norm + 0.00001)
data.func$log.μgCO2.7.norm <- log(data.func$μgCO2.7.norm + 0.00001)
# Selecting only the function and condition of interest e.g. ATP7
functions.list <- colnames(data.func) # Columns in function data
data.func.x <- data.frame(data.func[,grep(select.func, functions.list), drop = FALSE]) # Dataframe with only selected function
functions.list <- colnames(data.func.x) # Functions of interest with different conditions e.g. ATP7 and ATP14
data.func.x <- data.func.x[,grep(select.cond, functions.list), drop = FALSE] # Selects only the chosen condition e.g. ATP7
# data.func.x has only the chosen function and condition e.g. ATP7
data.func.x$sampleid <- paste0(data.func$Community, ".", data.func$Replicate) # adding sampleid column to data.func.x
sample.id <- intersect(data.func.x$sampleid, xIn.tmp$sampleid)  # The sampleid of samples in the function data and the asv table
# IMPORTANT NOTE: Only 1035/1402 samples in xIn.tmp are present in data.func
data.func.x <- data.func.x[data.func.x$sampleid %in% sample.id, ] # Only keeping samples in both the asv table and the function data
xIn <- xIn[xIn$sampleid %in% sample.id, ] # Only keeping sampleids that are in both the function data and the asv table
xIn <- xIn[order(xIn$sampleid), ] # Sorting so in same order as function data
data.func.x <- data.func.x[order(data.func.x$sampleid), ] # sorting so in same order as asv table
xIn <- data.frame(xIn[, -which(names(xIn) == "sampleid")]) # get rid of sampleid column from asv table
data.func.x <- data.func.x[, -which(names(data.func.x) == "sampleid")] # get rid of sampleid column from function data
yIn = as.numeric(data.func.x) # response vector






##### OPTIMISING NUMBER OF VARIABLES #####
set.seed(3062024) # today

# Parameters for RF
ntree.fixed <- 10000 # THIS WAS FOUND TO BE SUFFICIENTLY OPTIMISED
Nrand=25 # number of test random forests for each variable number
# Numbers of variables to test (10 different numbers in total). There are 1459 total variables.
nvar.test <- c(10,20,30,50,100,200,500,750,1000,1400)

# Function that takes a number of trees, and performs RF using the data set with fixed number of variables
random_forest_varnum <- function(varnum_input){ # yIn, xIn, and ntree are already defined so they do not need to be inputs for the function
  RF.tmp <- randomForest(y=yIn,x=xIn,
                         importance=T, proximity = T, 
                         ntree=ntree.fixed, mtry = varnum_input)
  mean_OOB_mse_tmp <- mean(RF.tmp$mse)
  return(mean_OOB_mse_tmp)
}
  
# Setting up how different jobs allocated
iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX")) # Gets the PBS environment variable
nvar.tmp <- nvar.test[iter] # The number of variables for the current job
  
# Vector that stores results
mse_results <- numeric(Nrand)

# Parallel backend
num_cores <- 25 # CHANGE DEPENDING ON NUMBER OF CPUS ALLOCATED IN PBS SHELL SCRIPT
cluster_of_cores <- makeCluster(num_cores) # Creates a cluster of these cores
registerDoParallel(cluster_of_cores) # Registers parallel backend with foreach package

# Perform the tests for given number of trees in parallel
mse_results <- foreach(i=1:Nrand, .combine='c', .packages='randomForest') %dopar% {
  random_forest_varnum(nvar.tmp)
}

# Stop parallel cluster
stopCluster(cluster_of_cores)

# Save results in file
mse_results <- c(nvar.tmp, mse_results)
saveRDS(mse_results, file=paste0(select.func, select.cond, "_", predictor.unit, "_test_var_mses_", nvar.tmp, "_vars.rds"))
# This file contains the results vector.
# The first item is the number of variables
# Each of the other items is an MSE value for one of the 25 test RFs of that variable number


  
  
