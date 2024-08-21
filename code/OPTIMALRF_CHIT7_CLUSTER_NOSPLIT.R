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
select.func="log.mN" # determine the function you want to predict. Name should match those from
# the column name of the function data frame.
select.cond=c("7.norm") # determine the function conditions you want to predict. Names should 
#match those used to differentiate the different conditions within the function df.
predictor.unit = "Cluster" # ASV or cluster as predictor - only for use in naming output files
replicate.train = "AllRepTrain" # All replicates trained using OOB or used one replicate to predict other replicates. Only used for naming.
just_func <- gsub("log.", "", select.func) # For plot titles
just_cond <- gsub(".norm", "", select.cond) # for plot titles

# LOADING LIBRARIES 
library(dplyr)
library(stringr)
library(ggplot2)
library(randomForest)
library(randomForestExplainer)
library(ggrepel)

# SOURCING FUNCTIONS
source("/rds/general/user/mg2020/home/researchprojecthpc/clean_ASV_table.R")

## READING INPUT FILES
file.ASV = "/rds/general/user/mg2020/home/researchprojecthpc/seqtable_readyforanalysis.csv"
file.taxa = "/rds/general/user/mg2020/home/researchprojecthpc/taxa_wsp_readyforanalysis.csv"
file.Meta = "/rds/general/user/mg2020/home/researchprojecthpc/metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv"
file.func = "/rds/general/user/mg2020/home/researchprojecthpc/20151016_Functions_remainder.csv"
file.cluster = "/rds/general/user/mg2020/home/researchprojecthpc/max_tot_ext_network_table.tsv"

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
# Cluster data
cluster_data <- read.csv(file = file.cluster, sep = "\t")

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

cluster_data <- cluster_data[cluster_data$functionInk %in% names(cluster_counts[cluster_counts >= 2]), ]

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
# Adding column with the logarithm of the number of reads in each sample
xIn.tmp$Counts = rowSums(xIn.tmp[,-1]) # add the sum of the counts for each sample as an additional predictor
xIn.tmp$Counts = log(xIn.tmp$Counts) # log it to avoid having too high numbers
xIn = xIn.tmp # ASV table with only final samples, and a column for the log of the number of reads for each sample

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
# Selecting only the function and condition of interest e.g. ATP7
functions.list <- colnames(data.func) # Columns in function data
data.func.x <- data.frame(data.func[,grep(select.func, functions.list), drop = FALSE]) # Dataframe with only selected function
functions.list <- colnames(data.func.x) # Functions of interest with different conditions e.g. ATP7 and ATP14
data.func.x <- data.func.x[,grep(select.cond, functions.list), drop = FALSE] # Selects only the chosen condition e.g. ATP7
# data.func.x has only the chosen function and condition e.g. ATP7
data.func.x$sampleid <- paste0(data.func$Community, ".", data.func$Replicate) # adding sampleid column to data.func.x
sample.id <- intersect(data.func.x$sampleid, xIn$sampleid)  # The sampleid of samples in the function data and the asv table
# IMPORTANT NOTE: Only 1035/1402 samples in xIn.tmp are present in data.func
data.func.x <- data.func.x[data.func.x$sampleid %in% sample.id, ] # Only keeping samples in both the asv table and the function data
xIn <- xIn[xIn$sampleid %in% sample.id, ] # Only keeping sampleids that are in both the function data and the asv table
xIn <- xIn[order(xIn$sampleid), ] # Sorting so in same order as function data
data.func.x <- data.func.x[order(data.func.x$sampleid), ] # sorting so in same order as asv table
xIn <- xIn[, -which(names(xIn) == "sampleid")] # get rid of sampleid column from asv table
data.func.x <- data.func.x[, -which(names(data.func.x) == "sampleid")] # get rid of sampleid column from function data
yIn = as.numeric(data.func.x) # response vector








##### RUNNING RF WITH OPTIMAL PARAMETERS #####
set.seed(3062024)
optimal_ntree <- 10000 # Optimal number of trees (found using another script)
optimal_nvar <- 50 # Optimal number of variables (found using another script)
# Name of output file
fileOut=paste("OPTIMALRF_Class-",select.class,select.cond,
              "_",predictor.unit,"_",replicate.train,
              ".RDS",sep="")
# Performing random forest
RF.out <-randomForest(y=yIn,x=xIn,
                        importance=T, proximity = T, 
                        ntree=optimal_ntree, mtry=optimal_nvar) # Running random forest model and storing it as RF.out
# Saving to a file
saveRDS(RF.out, file = fileOut)





################################# ANALYSE #################################
## randomForestExplainer importance of variables
importance.df = measure_importance(RF.out) # Measure importance of variables in MULTIPLE METRICS
fileOut=paste("RFEVarImp_Class-",select.class,select.cond,
              "_",predictor.unit,"_",replicate.train,
              ".tsv",sep="") # Name of file
rownames(importance.df) <- NULL
write.table(importance.df,file=fileOut,sep = "\t",quote=FALSE, row.names = FALSE) # Write to file

# Measure importance of variables by MSE, Node purity and mean min depth
top10asvs <- important_variables(importance.df,
                                 k = 10,
                                 measures = c("mean_min_depth", "mse_increase", "node_purity_increase"),
                                 ties_action = "draw")
importance.df$colour_group <- ifelse(importance.df$variable %in% top10asvs, "Top 10 ASVs", "Other ASVs")
top_10_data <- importance.df[importance.df$variable %in% top10asvs, ]

## Plot of multiway importance
p <- ggplot(importance.df, aes(x = mean_min_depth, y = node_purity_increase)) +
  geom_point(aes(color = colour_group, size = mse_increase), shape = 16) +  # Map color and size
  scale_color_manual(
    name = "Predictor importance",  # Update legend title
    values = c("Top 10 ASVs" = "blue", "Other ASVs" = "black"),  # Define colors
    labels = c("Other Predictors", "Top 10 Predictors"),  # Match the labels to the colors
    guide = guide_legend(
      order = 1,  # Ensure this guide is first
      override.aes = list(shape = 16)  # Ensure shape consistency
    )
  ) +
  scale_size_continuous(
    name = "MSE Increase",  # Update legend title
    range = c(2, 10),       # Size range for points
    limits = c(0, 1)        # Set limits for the size scale
  ) +
  geom_label_repel(
    data = top_10_data,
    aes(x = mean_min_depth, y = node_purity_increase, label = variable, color = colour_group),
    size = 5, 
    box.padding = 0.5, 
    point.padding = 0.5, 
    force = 50,  
    segment.color = "grey50",
    segment.size = 0.5,
    fill = "white",
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  xlim(0, 30) +  # Set x-axis limits
  ylim(0, 350) + # Set y-axis limits
  xlab("Mean minimum depth") +  # Update x-axis label
  ylab("Node purity increase") +  # Update y-axis label
  ggtitle("Explaining variation in the chitinase concentration per cell between microbial community samples using the abundances of subgroups") +  # Set the plot title dynamically
  theme_minimal() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),  # Title size and style
    axis.title.x = element_text(size = 18),  # X-axis title size
    axis.title.y = element_text(size = 18),  # Y-axis title size
    axis.text.x = element_text(size = 14),   # X-axis text size
    axis.text.y = element_text(size = 14),    # Y-axis text size
    legend.title = element_text(size = 16),  
    legend.text = element_text(size = 14),   
    legend.key.size = unit(1.5, "lines") 
  )

pdf("Plot_RFEMultiwayImp_CHIT_CLUSTER.pdf", height =6)
print(p)
dev.off() # Plots multiway importance

## Error plot
plotOut=paste("Plot_Error_Class-",select.class,select.cond,
              "_",predictor.unit,"_",replicate.train,
              ".pdf",sep="") # Name of file
plotTitle =  paste("Predicting ",just_func,just_cond,
                   " using ",predictor.unit, " abundances"
                   ,sep="") # Plot title name
pdf(plotOut)
plot(RF.out) # Plot RF to file
dev.off()

## Variables in order of mse_increase
sorted_importance_df_mse <- importance.df[order(-importance.df$mse_increase), ]
sorted_importance_df_mse <- sorted_importance_df_mse[,c("variable","mse_increase")]
fileOut=paste("varMSEINCR_Class-",select.class,select.cond,
              "_",predictor.unit,"_",replicate.train,
              ".tsv",sep="") # Name of file
write.table(sorted_importance_df_mse,file=fileOut,sep = "\t",quote=FALSE, row.names = FALSE) # Write to file

## Variables in order of node_purity_increase
sorted_importance_df_nodepurity <- importance.df[order(-importance.df$node_purity_increase), ]
sorted_importance_df_nodepurity <- sorted_importance_df_nodepurity[,c("variable","node_purity_increase")]
fileOut=paste("varNODEPURINCR_Class-",select.class,select.cond,
              "_",predictor.unit,"_",replicate.train,
              ".tsv",sep="") # Name of file
write.table(sorted_importance_df_nodepurity,file=fileOut,sep = "\t",quote=FALSE, row.names = FALSE) # Write to file

## Random forest explanation (overview) # MOST IMPORTANT THING # REMEMBER TO MANUALLY RENAME OUTPUT FILE SO NOT OVERWRITTEN
#fileOut=paste("RFExplained-",select.class,select.cond,
#              "_",predictor.unit,"_",replicate.train,
#              ".html",sep="") # Name of file
#explain_forest(RF.out, path = fileOut)

