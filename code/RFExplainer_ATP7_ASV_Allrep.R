# Loading packages
library(randomForest)
library(randomForestExplainer)

# Reading in RF model
RF.out <- readRDS("../results/Optimised_ATP7_ASV/OPTIMALRF_Class-log.ATP7.norm_ASV_AllRepTrain.RDS")

# Random forest explanation (overview)
fileOut=paste("OptimalRFExplained_ATP7_ASV_AllRep.html") # Name of file
explain_forest(RF.out, path = fileOut)