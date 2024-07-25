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

#setwd("~/home/matthew/Documents/ResearchProject/ResearchProjectRepository/code") # MAKE SURE THIS IS THE WORKING DIRECTORY
# setwd("../../code") # Setting working directory
rm(list = ls()) # Clear workspace

# RANDOM FOREST PARAMETERS
optimization=0 # Should RF parameters be optimized? (=1)
run_RF=1 # Should RF be run (=1) or just read from file (=0)
type_RF="regression" # Should RF be a "classification" or a "regression"
ntree.min=10000 # number of trees, mandatory if optimization = 0 (and only used in that case)
mtry.min="default" # number of variables randomly selected to build the trees. Fix to "default" 
# if you don't have an informed guess
partial.plots=0 # should partial plots be generated (=1), fix to 0 otherwise. It will
# generate plots for the 10 most important variables. This is controlled by
# variable Nsel in the section plots below

# SELECT FUNCTIONS AND CONDITIONS OF INTEREST
select.func="log.ATP" # determine the function you want to predict. Name should match those from
# the column name of the function data frame.
select.cond=c("7.norm") # determine the function conditions you want to predict. Names should 
#match those used to differentiate the different conditions within the function df.
predictor.unit = "ASV" # ASV or cluster as predictor - only for use in naming output files
replicate.train = "AllRepTrain" # All replicates trained using OOB or used one replicate to predict other replicates. Only used for naming.
just_func <- gsub("log.", "", select.func) # For plot titles
just_cond <- gsub(".norm", "", select.cond) # for plot titles



select.class = select.func

# --- Input and output files
file.ASV = "seqtable_readyforanalysis.csv"
file.taxa = "taxa_wsp_readyforanalysis.csv"
file.Meta = "metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv"
file.func = "20151016_Functions_remainder.csv"

# Directories
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/code/")[[1]][1]
dirSrc=paste(this.dir,"/code/",sep="") # Directory where the code is
dirASV=paste(this.dir,"/data/",sep="") # Dir of ASV table
dirMeta=paste(this.dir,"/data/",sep="") # Dir of metadata
dirFunc=paste(this.dir,"/data/",sep="") #Dir of response function data frame
dirOut=paste(this.dir,"/results/fnl_random_forest",sep="") # Dir of output data

# --- Reshaping data parameters
nreads = 10000 # minimum number of reads to consider a sample
exclude_exp = c("4M") # A vector of characters with the experiments that should be excluded
match_exp = TRUE # Set to true if only starting communities that were resurrected should be included
output.label = "Time0D_7D_matched" 

# LOADING LIBRARIES 
library(tidyverse)
library(stringr)
library(ggplot2)
library(caret)
library(randomForest)
library(randomForestExplainer)

# SOURCING FUNCTIONS
source(paste0(dirSrc, "clean_ASV_table.R"))

# READING INPUT FILES
# ASV Table
setwd(dirASV)
ASV.table=read.table(file = file.ASV, sep="\t")
# Metadata
setwd(dirMeta)
sample_md <-read.table(file = file.Meta, sep="\t", header=TRUE)
# Function data
setwd(dirFunc)
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
xIn.tmp$sampleid = rownames(xIn.tmp) # create a new column for sampleid
xIn = xIn.tmp # ASV table with only final samples, and a column for the log of the number of reads for each sample
rownames(xIn) = seq(1,dim(xIn)[1]) # Row names are now numbers
xIn = xIn[,c(dim(xIn)[2],1:(dim(xIn)[2]-1))] # reorder so that the first column is sampleid
# Splitting into 2dfs: one for rep1 and one for rep2,3,4

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






##### OPTIMISING RF PARAMETERS #####
set.seed(3062024) # today

# Estimate the optimal RF parameters:
# the two steps (ntree and mtry) should possibly be iterated
dir.create(dirOut) # creating new directory for rf outputs
setwd(dirOut) # setting this as the working directory

# OPTIMISING NUMBER OF TREES
mtry.default=sqrt(dim(xIn)[2]) # default number of variables taken to build the trees
if(optimization == 1){ # If want to optimise rf parameters
  range="large"
  Nrand=25 # number of test random forests for each tree number (each value in ntree.test)
  if(range == "small"){ # Determining numbers of trees to test
    ntree.test=seq(from=100,to=500,by=20) # small range
  }else{
    ntree.test=c(100,250,500, seq(from=1000,to=11000,by=2000)) # large (e.g. 100 trees, 250 trees, 500 trees...)
    ntree.test=c(ntree.test,14000,18000,22000) # very large
    #ntree.test=c(ntree.test,12000,14000,16000,18000,20000) # very large
  }
  mtry.test=mtry.default #300 # if different than default this was fixed after one iteration with the next step below
  OOB=matrix(0,nrow=length(ntree.test),ncol=Nrand) # matrix with a cell for each of the 25 tests of each number of trees for RF
  i=0
  # Loop performing RF Nrand times for each number of trees in ntree.test, and storing it in OOB data frame+
  for(ntree.tmp in ntree.test){ # for each number of trees
    i=i+1
    for(k in 1:Nrand){ # for each test
      RF.tmp <-randomForest(y=yIn,x=xIn,
                            importance=T, proximity = T, 
                            ntree=ntree.tmp,mtry = mtry.test) # perform RF
      
      if (type_RF == "classification") {
        OOB[i,k]=mean(RF.tmp$err.rate[,1])
      } else if (type_RF == "regression") {
        OOB[i,k]=mean(RF.tmp$mse) # mean square error for regression
      } else {
        warning("type not recognised")
      }
      
    }
  }
  OOB.mean=rowMeans(OOB) # mean squared error for each number of trees
  OOB.std=apply(OOB,1, sd, na.rm = TRUE)
  OOB.df=data.frame(cbind(ntree.test,OOB.mean,OOB.std)) # Data frame with different numbers of trees, their means, and their sds
  #rownames(OOB.df)=ntree.test
  OOB.min.id=which.min(OOB.mean) # id of number of trees with minimum mse
  ntree.min=OOB.df$ntree.test[OOB.min.id] # number of trees with minimum mse
  ymin=OOB.df$OOB.mean-OOB.df$OOB.std/sqrt(Nrand) 
  ymax=OOB.df$OOB.mean+OOB.df$OOB.std/sqrt(Nrand)
  fileOut=paste("optimization_ntree_Class-",select.class,select.cond,"_mtry",trunc(mtry.test),
                "_",range,"_",predictor.unit,"_",replicate.train, ".csv",sep="")
  write.table(OOB.df,file = fileOut,sep="\t",quote = FALSE,row.names = FALSE) # Creating file with the different numbers of trees and mse etc
  plotOut=paste("Plot_optimization_ntree_Class-",select.class,select.cond,"_mtry",trunc(mtry.test),
                "_",range,"_",predictor.unit,"_",replicate.train,".pdf",sep="")
  pdf(plotOut) # Creating file with plot of mse against number of trees
  g=ggplot()+
    geom_vline(xintercept = ntree.min,linetype = 'dotted', col = 'red')+
    geom_point(data=OOB.df,aes(x=ntree.test,y=OOB.mean))+
    ylab("Mean OOB error")+xlab("Number of trees")+
    scale_y_continuous(trans='log10')+
    geom_errorbar(aes(x=ntree.test,ymin=ymin,ymax=ymax))+
    theme_bw()
  print(g)
  dev.off()
}else # if we do not optimize
  if(mtry.min == "default"){ # we need to give a value if the user choose a default value
    mtry.min = mtry.default
}

## OPTIMISING NUMBER OF VARIABLE
if(optimization == 1){ # If want to optimise rf parameters
  range="large"
  Nrand=25 # number of test random forests for each number of variables
  # Now we fix the optimal ntree and look for the optimization of mtry
  ntree.in=ntree.min
  #ntree.in = 6000 # same order of magnitude than the minimum.
  mtry.test=seq(from=mtry.default/2,to= dim(xIn)[2],by=mtry.default/2)
  #OOB=vector(mode="numeric",length=length(mtry.test))
  OOB=matrix(0,nrow=length(mtry.test),ncol=Nrand) # data frame with a number of tests for each different numbers of parameters
  i=0
  for(mtry.tmp in mtry.test){ # For every number of parameters
    i=i+1
    for(k in 1:Nrand){ # For every test
      RF.tmp <-randomForest(y=yIn,x=xIn,
                            importance=T, proximity = T, 
                            ntree=ntree.in,mtry=mtry.tmp) # Random forest
      
      if (type_RF == "classification") {
        OOB[i,k]=mean(RF.tmp$err.rate[,1])
      } else if (type_RF == "regression") {
        OOB[i,k]=mean(RF.tmp$mse) # mse if regression
      } else {
        warning("type not recognised")
      }
      
    }
  }
  OOB.mean=rowMeans(OOB)
  OOB.std=apply(OOB,1, sd, na.rm = TRUE)
  OOB.df=data.frame(cbind(mtry.test,OOB.mean,OOB.std))
  #rownames(OOB.df)=ntree.test
  OOB.min.id=which.min(OOB.mean)
  mtry.min=OOB.df$mtry.test[OOB.min.id]
  ymin=OOB.df$OOB.mean-OOB.df$OOB.std/sqrt(Nrand)
  ymax=OOB.df$OOB.mean+OOB.df$OOB.std/sqrt(Nrand)
  fileOut=paste("optimization_mtry_Class-",select.class,select.cond,"_ntree-",trunc(ntree.min),
                "_",range,"_",predictor.unit,"_",replicate.train,".csv",sep="")
  write.table(OOB.df,file = fileOut,sep="\t",quote = FALSE,row.names = FALSE) # File with different numbers of parameters and mse etc
  plotOut=paste("Plot_optimization_mtry_Class-",select.class,select.cond,
                "_ntree",ntree.in,"_",predictor.unit,"_",replicate.train,".pdf",sep="")
  pdf(plotOut) # Plot of different numbers of paramters and associated mse
  g=ggplot()+
    geom_point(data=OOB.df,aes(x=mtry.test,y=OOB.mean))+
    ylab("Mean OOB error")+xlab("Number of variables")+
    geom_errorbar(aes(x=mtry.test,ymin=ymin,ymax=ymax))
  print(g)
  dev.off()
}else{ # if we do not optimize
  if(mtry.min == "default"){ # we need to give a value if the user choose a default value
    mtry.min = mtry.default
  }
}




##### RUNNING RF WITH OPTIMAL (OR DEFAULT) PARAMETERS #####
ntree.in=ntree.min # Takes optimised number of trees or default if running without optimisation
mtry.in=mtry.min # Takes optimised number of variables or default if running without optimisation
fileOut=paste("RandForestOut_Class-",select.class,select.cond,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              "_",predictor.unit,"_",replicate.train,
              ".RDS",sep="") # File name
if(run_RF == 1){
  RF.out <-randomForest(y=yIn,x=xIn,
                        importance=T, proximity = T, 
                        ntree=ntree.in,mtry=mtry.in) # Running random forest model and storing it as RF.out
  # --- Have a look at the output
  RF.out
  
  if (type_RF == "classification") {
    mean(RF.out$err.rate[,1]) # this is the mean OOB
  } else if (type_RF == "regression") {
    mean(RF.out$mse) # this is the mean OOB
  } else {
    warning("type not recognised")
  }
  
  saveRDS(RF.out, file = fileOut)
}else{
  RF.out = readRDS(file = fileOut)
}

RF.out 






##### ANALYSE #####

setwd(dirOut)
## Extract important variables
# Measure of importance of variables (IN NO ORDER)
importance.df = measure_importance(RF.out) # Measure importance of variables in MULTIPLE METRICS
fileOut=paste("varImpExt_Class-",select.class,select.cond,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              "_",predictor.unit,"_",replicate.train,
              ".tsv",sep="") # Name of file
rownames(importance.df) <- NULL
write.table(importance.df,file=fileOut,sep = "\t",quote=FALSE, row.names = FALSE) # Write to file
# Most important variables
ASV.top = important_variables(importance.df) # Most important variables (BASED UPON COMBINED IMPORTANCE METRICS IN UNKNOWN WAY)
fileOut=paste("varTop_Class-",select.class,select.cond,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              "_",predictor.unit,"_",replicate.train,
              ".tsv",sep="") # Name of file
write.table(ASV.top,file=fileOut,sep = "\t",quote=FALSE, row.names = FALSE) # Write to file
## Top ASVs in mse_increase
sorted_importance_df_mse <- importance.df[order(-importance.df$mse_increase), ]
fileOut=paste("varTopMSEINCR_Class-",select.class,select.cond,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              "_",predictor.unit,"_",replicate.train,
              ".tsv",sep="") # Name of file
write.table(sorted_importance_df_mse,file=fileOut,sep = "\t",quote=FALSE, row.names = FALSE) # Write to file
## Top ASVs in node_purity_increase
sorted_importance_df_nodepurity <- importance.df[order(-importance.df$node_purity_increase), ]
fileOut=paste("varTopNODEPURINCR_Class-",select.class,select.cond,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              "_",predictor.unit,"_",replicate.train,
              ".tsv",sep="") # Name of file
write.table(sorted_importance_df_nodepurity,file=fileOut,sep = "\t",quote=FALSE, row.names = FALSE) # Write to file


## Variable importance again, this time using Caret function NOT SEEMINLY VERY USEFUL
fileOut=paste("varImpCaret_Class-",select.class,select.cond,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              "_",predictor.unit,"_",replicate.train,
              ".tsv",sep="") # Name of file
varImp.out=varImp(RF.out) # Taking importance of variables again (using diff package to above) NO INFO ON IMPORTANCE METRICS
varImp.out=varImp.out[which(rowSums(varImp.out) != 0),] # Keep only variables with an importance that is not 0
write.table(varImp.out,file=fileOut,sep = "\t",quote=FALSE) # Write to file



## Variable importance plots % Inc MSE
# 1 - using varImpPlot() NOT VERY USEFUL AS DON'T KNOW MEASURE OF IMPORTANCE
plotOut=paste("Plot_varImp_Class-",select.class,select.cond,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              "_",predictor.unit,"_",replicate.train, "_function1",
              ".pdf",sep="") # Name of plot file
pdf(plotOut,width=10)
varImpPlot(RF.out, type=1,main="")# ,xlab="Variable importance",title="") # THIS IS JUST A SIMPLE, ALTERNATIVE WAY OF VISUALISING THE SAME INFO AS THE BELOW GGPLOT 
dev.off()
# 2 - Instead of using varImpPlot(), get variable importance from the model fit. Variable importance plot using %IncMSE as MeanDecreaseAccuracy not available
# USING importance() which is different to measure_importance() although not sure if diff calculations or same
# Size indicates node purity
plotOut=paste("Plot_varImp_Class-",select.class,select.cond,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              "_",predictor.unit,"_",replicate.train, "_function2",
              ".pdf",sep="") # Name of plot file
pdf(plotOut,width=10)
ImpData <- as.data.frame(importance(RF.out)) # Calculating importance of variables again (3rd time?) using different function
ImpData$Var.Names <- row.names(ImpData) # Names of variables as column
names(ImpData)[names(ImpData) == "%IncMSE"] <- "PercentageIncMSE" # Renaming to avoid syntax errors
ImpData.sort.idx = sort(ImpData$PercentageIncMSE, 
                        decreasing = TRUE, index.return = T)
quantile(ImpData.sort.idx$x) # Take quantile
Imp.Data.sort = ImpData[ImpData.sort.idx$ix, ] # Reorders original df in sorted order
Imp.Data.sort = Imp.Data.sort[which(Imp.Data.sort$PercentageIncMSE > 20),] # Reducing clutter by only selecting more important ASVs

ggplot(Imp.Data.sort, aes(x=Var.Names, y=PercentageIncMSE)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=PercentageIncMSE), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) # Plotting variable importance
dev.off()

# Plot variable importance, multiway
plotOut=paste("Plot_MindepthNodepurMseinc_Class-",select.class,select.cond,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              "_",predictor.unit,"_",replicate.train,
              ".pdf",sep="") # Plot file name
plotTitle =  paste("Predicting ",just_func,just_cond,
                   " using ",predictor.unit, " abundances"
                   ,sep="") # Plot title name
pdf(plotOut, height =6)
p = plot_multi_way_importance(
  importance_frame = importance.df,
  x_measure = "mean_min_depth",
  y_measure = "node_purity_increase",
  size_measure = "mse_increase",
  min_no_of_trees = 0,
  no_of_labels = 10,
  main = plotTitle
)
p = p + xlab("Average depth") +ylab("Node purity increase") +
  labs(
    size = "MSE Increase (%)",          # Change legend title for size
  )
p
dev.off() # Plots multiway importance

## Random forest explanation (overview)                                 # MOST IMPORTANT THING
# REMEMBER TO MANUALLY RENAME OUTPUT FILE SO NOT OVERWRITTEN
explain_forest(RF.out) # gives distribution of minimal depth

## Error plot
plotOut=paste("Plot_Error_Class-",select.class,select.cond,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              "_",predictor.unit,"_",replicate.train,
              ".pdf",sep="") # Name of file
plotTitle =  paste("Predicting ",just_func,just_cond,
                   " using ",predictor.unit, " abundances"
                   ,sep="") # Plot title name
pdf(plotOut)
plot(RF.out) # Plot RF to file
dev.off()


## Plotting interactions between top 10 variables (according to mean minimal depth and number of trees in which appeared)
ASV.top.10 = important_variables(importance.df, k = 10, measures = c("mean_min_depth", "no_of_trees")) # Most important variables accroding to mean minimal depth and number of trees in which appeared
interactions_df <- min_depth_interactions(RF.out, ASV.top.10)
fileOut=paste("VariableInteractions_Class-",select.class,select.cond,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              "_",predictor.unit,"_",replicate.train,
              ".tsv",sep="") # Name of file
write.table(ASV.top,file=fileOut,sep = "\t",quote=FALSE, row.names = FALSE) # Write to file

plotOut=paste("Plot_Interactions-",select.class,select.cond,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              "_",predictor.unit,"_",replicate.train,
              ".pdf",sep="") # Name of file
pdf(plotOut)
plot_min_depth_interactions(interactions_df)
dev.off()







## Partial plots sorted by importance
# ... these plots may take a long time
if(partial.plots == 1){
  Nsel=10 # select only 10 vars
  imp <- importance(RF.out) # Importance of variables
  impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)] # Ordered with most important at top (I THINK, DIDN'T CHECK)
  
  impvar=impvar[1:Nsel] # 10 most important var
  Nlev=levels(yIn)
  
  # Interpret partial plots
  # https://stats.stackexchange.com/questions/147763/meaning-of-y-axis-in-random-forest-partial-dependence-plot
  # p=seq(from=0.005, to=0.995, by=0.005) # to understand the plots
  # plot(p,log(p/(1-p))) # plot the logit function
  for (i in seq_along(impvar)) {
    var.lab=str_replace_all(impvar[i],pattern="[[:punct:]]",replacement = "")
    var.lab=str_replace_all(var.lab,pattern = " ",replacement = "_")
    plotOut=paste("Plot_PartialDependence_",trunc(ntree.in),"_mtry",trunc(mtry.in),
                  "Class-",select.class,select.cond,"_Var-",var.lab,
                  "_",predictor.unit,"_",replicate.train,
                  ".pdf",sep="")
    pdf(file=plotOut,width=12)
    op <- par(mfrow=c(2, 3))
    for(level in levels(yIn)){
      partialPlot(RF.out, xIn, impvar[i], xlab=impvar[i],
                  which.class = level,
                  main=paste("Partial Dependence for class", level)) #ylim=c(30, 70))
    }
    par(op)
    dev.off()
  }
}

