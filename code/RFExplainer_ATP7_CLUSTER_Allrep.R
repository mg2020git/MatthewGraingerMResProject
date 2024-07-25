# Loading packages
library(randomForest)
library(randomForestExplainer)

# Reading in RF model
RF.out <- readRDS("../results/Optimised_ATP7_Cluster/OPTIMALRF_Class-log.ATP7.norm_Cluster_AllRepTrain.RDS")

# Random forest explanation (overview)
fileOut=paste("OptimalRFExplained_ATP7_Cluster_AllRep.html",sep="") # Name of file
explain_forest(RF.out, path = fileOut)