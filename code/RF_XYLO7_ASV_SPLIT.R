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

# generate plots for the 10 most important variables. This is controlled by
# variable Nsel in the section plots below

# SELECT FUNCTIONS AND CONDITIONS OF INTEREST
select.func="log.mX" # determine the function you want to predict. Name should match those from
# the column name of the function data frame.
select.cond=c("7.norm") # determine the function conditions you want to predict. Names should 
#match those used to differentiate the different conditions within the function df.
predictor.unit = "ASV" # ASV or cluster as predictor - only for use in naming output files
replicate.train = "Split" # All replicates trained using OOB or used one replicate to predict other replicates. Only used for naming.
just_func <- gsub("log.", "", select.func) # For plot titles
just_cond <- gsub(".norm", "", select.cond) # for plot titles

# LOADING LIBRARIES 
library(dplyr)
library(stringr)
library(ggplot2)
library(randomForest)
library(randomForestExplainer)

# SOURCING FUNCTIONS
source("/rds/general/user/mg2020/home/researchprojecthpc/clean_ASV_table.R")
#source("../code/clean_ASV_table.R")

## READING INPUT FILES
file.ASV = "/rds/general/user/mg2020/home/researchprojecthpc/seqtable_readyforanalysis.csv"
file.taxa = "/rds/general/user/mg2020/home/researchprojecthpc/taxa_wsp_readyforanalysis.csv"
file.Meta = "/rds/general/user/mg2020/home/researchprojecthpc/metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv"
file.func = "/rds/general/user/mg2020/home/researchprojecthpc/20151016_Functions_remainder.csv"
#file.ASV = "../data/seqtable_readyforanalysis.csv"
#file.taxa = "../data/researchprojecthpc/taxa_wsp_readyforanalysis.csv"
#file.Meta = "../data/metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv"
#file.func = "../data/20151016_Functions_remainder.csv"

select.class = select.func

# --- Reshaping data parameters
nreads = 10000 # minimum number of reads to consider a sample
exclude_exp = c("4M") # A vector of characters with the experiments that should be excluded
match_exp = TRUE # Set to true if only starting communities that were resurrected should be included
output.label = "Time0D_7D_matched" 

# READING INPUT FILES
# ASV Table
ASV.table=read.table(file = file.ASV, sep="\t")
# Metadata
sample_md <-read.table(file = file.Meta, sep="\t", header=TRUE)
# Function data
data.func <- read.csv(file = file.func)

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
# Splitting into 1df for training (rep1 samples) and 1df for testing (rep2,3,4 samples)
id.rep1 = grep("Rep1", sample_md$replicate) # Indices of rows which are rep1 communities
samples.rep1 = as.character(sample_md$sampleid[id.rep1]) # Sampleids of rep1 communities
xIn.train <- xIn.tmp[samples.rep1, ] # ASV table containing only rep1 communities
id.rep234 = sample_md$replicate %in% c("Rep2", "Rep3", "Rep4") # Indices of rep2,3,4
samples.rep234 = as.character(sample_md$sampleid[id.rep234]) # samples ids rep2,3,4
xIn.test <- xIn.tmp[samples.rep234, ] # Only rep2,3,4 communities
# create a new column for sampleid
xIn.train$sampleid = rownames(xIn.train)
xIn.test$sampleid = rownames(xIn.test)
# Formatting
rownames(xIn.train) = seq(1,dim(xIn.train)[1]) # Row names are now numbers
xIn.train = xIn.train[,c(dim(xIn.train)[2],1:(dim(xIn.train)[2]-1))] # reorder so that the first column is sampleid
rownames(xIn.test) = seq(1,dim(xIn.test)[1]) # Row names are now numbers
xIn.test = xIn.test[,c(dim(xIn.test)[2],1:(dim(xIn.test)[2]-1))] # reorder so that the first column is sampleid

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
# Splitting function data for the training (rep1) and testing (rep2,3,4) data
sample.id.train <- intersect(data.func.x$sampleid, xIn.train$sampleid)  # The sampleid of samples in the function data and the training data
sample.id.test <- intersect(data.func.x$sampleid, xIn.test$sampleid) # The sampleid of samples in the function data and the testing data
data.func.x.train <- data.func.x[data.func.x$sampleid %in% sample.id.train, ] # Only keeping samples in both the training data and the function data
data.func.x.test <- data.func.x[data.func.x$sampleid %in% sample.id.test, ] # Only keeping samples in both the testing data and the function data
xIn.train <- xIn.train[xIn.train$sampleid %in% sample.id.train, ] # Only keeping sampleids that are in both the function data and the training data
xIn.train <- xIn.train[order(xIn.train$sampleid), ] # Sorting so in same order as function data
data.func.x.train <- data.func.x.train[order(data.func.x.train$sampleid), ] # sorting so in same order as training data
xIn.test <- xIn.test[xIn.test$sampleid %in% sample.id.test, ] # Only keeping sampleids that are in both the function data and the testing data
xIn.test <- xIn.test[order(xIn.test$sampleid), ] # Sorting so in same order as function data
data.func.x.test <- data.func.x.test[order(data.func.x.test$sampleid), ] # sorting so in same order as testing data
xIn.train <- data.frame(xIn.train[, -which(names(xIn.train) == "sampleid")]) # get rid of sampleid column from training data
xIn.test <- data.frame(xIn.test[, -which(names(xIn.test) == "sampleid")]) # get rid of sampleid column from testing data
data.func.x.train <- data.func.x.train[, -which(names(data.func.x.train) == "sampleid")] # get rid of sampleid column from function data
data.func.x.test <- data.func.x.test[, -which(names(data.func.x.test) == "sampleid")] # get rid of sampleid column from function data
yIn.train = as.numeric(data.func.x.train) # response vector training
yIn.test = as.numeric(data.func.x.test) # response vector testing
## NOTE: Using replicate 1 communities for training, rep 2, 3, 4 for predicting ##








##### TRAINING RF WITH OPTIMAL PARAMETERS #####
set.seed(3062024)
optimal_ntree <- 10000 # Optimal number of trees (found using another script)
optimal_nvar <- 750 # Optimal number of variables (found using another script)
# Name of output file
fileOut=paste("OPTIMALTRAINED_Class-",select.class,select.cond,
              "_",predictor.unit,"_",replicate.train,
              ".RDS",sep="")
# Performing random forest
RF.out <-randomForest(y=yIn.train,x=xIn.train,
                        importance=T, proximity = T, 
                        ntree=optimal_ntree, mtry=optimal_nvar) # Running random forest model and storing it as RF.out
# Saving to a file
saveRDS(RF.out, file = fileOut)

##### Comparing predictions and actual values #####
predictions <- predict(RF.out, xIn.test)
mse_of_predictions <- mean((predictions - yIn.test)^2)
rsq_of_predictions <- cor(predictions, yIn.test)^2
cat("Mean Squared Error (MSE):", mse_of_predictions, "\n")
cat("R-squared:", rsq_of_predictions, "\n")


# Save predictions, mse and R-squared value to files
predictions_output_file <- paste("SPLIT_Predictions_", select.class, select.cond, 
                                 "_", predictor.unit,".csv", sep="")
rsq_mse_output_file <- paste("SPLIT_RSquaredMSE_", select.class, select.cond,
                             "_", predictor.unit,".txt", sep="")
predictions_df <- data.frame(sampleid = sample.id.test, Actual = yIn.test, Predicted = predictions)
write.csv(predictions_df, predictions_output_file, row.names = FALSE)
cat("R-squared:", rsq_of_predictions, "\n", "MSE:", mse_of_predictions, "\n", file = rsq_mse_output_file)






##### Basic overview of the RF model trained on only rep1 #####
## Error plot
plotOut=paste("SPLITPlot_Error_Class-",select.class,select.cond,
              "_",predictor.unit,"_",replicate.train,
              ".pdf",sep="") # Name of file
plotTitle =  paste("Predicting ",just_func,just_cond,
                   " using ",predictor.unit, " abundances"
                   ,sep="") # Plot title name
pdf(plotOut)
plot(RF.out) # Plot RF to file
dev.off()

## Random forest explanation (overview) # MOST IMPORTANT THING # REMEMBER TO MANUALLY RENAME OUTPUT FILE SO NOT OVERWRITTEN
#fileOut=paste("SPLITRFExplained-",select.class,select.cond,
#              "_",predictor.unit,"_",replicate.train,
#              ".html",sep="") # Name of file
#explain_forest(RF.out, path = fileOut)


