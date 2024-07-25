# Loading packages
library(randomForest)
library(randomForestExplainer)

# Reading in RF model
RF.out <- readRDS("../results/Optimised_XYLO7_ASV/OPTIMALRF_Class-log.mX7.norm_ASV_AllRepTrain.RDS")

# Random forest explanation (overview)
fileOut=paste("OptimalRFExplained_XYLO7_ASV_AllRep.html",sep="") # Name of file
explain_forest(RF.out, path = fileOut)