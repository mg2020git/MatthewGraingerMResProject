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
select.func="log.mP" # determine the function you want to predict. Name should match those from
# the column name of the function data frame.
select.cond=c("7.norm") # determine the function conditions you want to predict. Names should 
#match those used to differentiate the different conditions within the function df.
predictor.unit = "ASV" # ASV or cluster as predictor - only for use in naming output files
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
#
#source("clean_ASV_table.R")

## READING INPUT FILES
file.ASV = "/rds/general/user/mg2020/home/researchprojecthpc/seqtable_readyforanalysis.csv"
file.taxa = "/rds/general/user/mg2020/home/researchprojecthpc/taxa_wsp_readyforanalysis.csv"
file.Meta = "/rds/general/user/mg2020/home/researchprojecthpc/metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv"
file.func = "/rds/general/user/mg2020/home/researchprojecthpc/20151016_Functions_remainder.csv"
#
#file.ASV = "../data/seqtable_readyforanalysis.csv"
#file.taxa = "../data/taxa_wsp_readyforanalysis.csv"
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
# Network data
#network_data <- read.csv(file = "../data/combined_network_metrics.csv")
#
network_data <- read.csv(file = "/rds/general/user/mg2020/home/researchprojecthpc/combined_network_metrics.csv")

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
xIn.tmp$Counts = rowSums(xIn.tmp) # add the sum of the counts for each sample as an additional predictor
xIn.tmp$Counts = log(xIn.tmp$Counts) # log it to avoid having too high numbers
xIn.tmp$sampleid = rownames(xIn.tmp) # create a new column for sampleid
xIn = xIn.tmp # ASV table with only final samples, and a column for the log of the number of reads for each sample
rownames(xIn) = seq(1,dim(xIn)[1]) # Row names are now numbers
xIn = xIn[,c(dim(xIn)[2],1:(dim(xIn)[2]-1))] # reorder so that the first column is sampleid

# Only keeping ASVs with the top 5% number of interactions (greater than degree of 21)
network_data <- network_data[,c("name", "Degree")] # only keeping asv and degree columns
colnames(network_data)[colnames(network_data) == "name"] <- "ASV" # renaming name column to asv
percentile_95 <- quantile(network_data$Degree, 0.95) # 95th percentile
network_data <- network_data[network_data$Degree > 21, ] # Only keep ASVs in top 5% of number of interactions
network_data <- network_data$ASV # Get rid of degree column
network_data <- intersect(network_data, colnames(xIn))
xIn <- xIn[, c(network_data, "sampleid")] # Only keep ASVs in xIn in top 5% degree

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








##### RUNNING RF WITH OPTIMAL PARAMETERS #####
set.seed(3062024)
optimal_ntree <- 10000 # Optimal number of trees (found using another script)
optimal_nvar <- ncol(xIn) # number of columns (asvs) in predictor data frame
# Name of output file
fileOut=paste("DEGREERF_Class-",select.class,select.cond,
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
importance.df <- measure_importance(RF.out) # Measure importance of variables in MULTIPLE METRICS
fileOut=paste("DEGREERFEVarImp_Class-",select.class,select.cond,
              "_",predictor.unit,"_",replicate.train,
              ".tsv",sep="") # Name of file
rownames(importance.df) <- NULL
write.table(importance.df,file=fileOut,sep = "\t",quote=FALSE, row.names = FALSE) # Write to file
# importance.df <- read.table("DEGREERFEVarImp_Class-log.mP7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")

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
  ggtitle("Importance of taxa in explaining differences\nin phosphatase production between microbial\ncommunity samples") +  # Set the plot title dynamically
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title.x = element_text(size = 14),  # X-axis title size
    axis.title.y = element_text(size = 14),  # Y-axis title size
    axis.text.x = element_text(size = 12),   # X-axis text size
    axis.text.y = element_text(size = 12),    # Y-axis text size
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),   
    legend.key.size = unit(1.5, "lines") 
  )

pdf("DEGREEPlot_RFEMultiwayImp_PHOS_ASV.pdf", width = 7, height = 6)
print(p)
dev.off() # Plots multiway importance

## Error plot
plotOut=paste("DEGREEPlot_Error_Class-",select.class,select.cond,
              "_",predictor.unit,"_",replicate.train,
              ".pdf",sep="") # Name of file
plotTitle =  paste("Predicting ",just_func,just_cond,
                   " using ",predictor.unit, " abundances"
                   ,sep="") # Plot title name
pdf(plotOut)
plot(RF.out) # Plot RF to file
dev.off()


## Random forest explanation (overview) # MOST IMPORTANT THING 
#fileOut=paste("RFExplained-",select.class,select.cond,
#              "_",predictor.unit,"_",replicate.train,
#              ".html",sep="") # Name of file
#explain_forest(RF.out, path = fileOut)


