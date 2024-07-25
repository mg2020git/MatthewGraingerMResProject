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

# SELECT FUNCTIONS AND CONDITIONS OF INTEREST
select.func="log.mX" # determine the function you want to predict. Name should match those from
# the column name of the function data frame.
select.cond=c("7.norm") # determine the function conditions you want to predict. Names should 
#match those used to differentiate the different conditions within the function df.
predictor.unit = "Cluster" # ASV or cluster as predictor - only for use in naming output files
replicate.train = "AllRepTrain" # All replicates trained using OOB or used one replicate to predict other replicates. Only used for naming.
just_func <- gsub("log.", "", select.func) # For plot titles
just_cond <- gsub(".norm", "", select.cond) # for plot titles
select.class = select.func

# LOADING LIBRARIES 
library(dplyr)
library(stringr)
library(ggplot2)
library(randomForest)
library(randomForestExplainer)

# Load in optimal RF
RF.out <- readRDS("/rds/general/user/mg2020/home/researchprojecthpc/OPTIMALRF_Class-log.mX7.norm_Cluster_AllRepTrain.RDS")


##### Analyse interactions #####
## randomForestExplainer importance of variables
importance.df <- measure_importance(RF.out) # Measure importance of variables in MULTIPLE METRICS
rownames(importance.df) <- NULL
## Most important variables # DECIDE WHICH MEASURES OF IMPORTANCE TO USE + DECIDE HOW MANY VARIABLES TO USE (PROB BASED ON OPTIMAL NUM VAR)
ASV.top <- important_variables(importance.df, k=20, measures = c("mean_min_depth", "times_a_root", "no_of_trees", "mse_increase", "node_purity_increase")) 
fileOut=paste("RFETopVar_",select.class,select.cond,
              "_",predictor.unit,"_",replicate.train,
              ".tsv",sep="") # Name of file
write.table(ASV.top,file=fileOut,sep = "\t",quote=FALSE, row.names = FALSE) # Write to file
## Interactions between most important variables
interaction_df <- min_depth_interactions(RF.out, vars = ASV.top)
fileOut=paste("RFEInteractions_Class-",select.class,select.cond,
              "_",predictor.unit,"_",replicate.train,
              ".tsv",sep="") # Name of file
write.table(interaction_df,file=fileOut,sep = "\t",quote=FALSE, row.names = FALSE) # Write to file
## Plot of interactions between most important variables
plotOut=paste("Plot_Interactions_Class-",select.class,select.cond,
              "_",predictor.unit,"_",replicate.train,
              ".pdf",sep="") # Name of plot file
pdf(plotOut,width=10) 
plot_min_depth_interactions(interaction_df)
dev.off()



