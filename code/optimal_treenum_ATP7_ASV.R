


# Loading files containing mses for the 25 tests of each tree number
data_100_trees <- readRDS("../data/Optimise_RF_ASV_ATP7/test_mses_100_trees.rds")
data_250_trees <- readRDS("../data/Optimise_RF_ASV_ATP7/test_mses_250_trees.rds")
data_500_trees <- readRDS("../data/Optimise_RF_ASV_ATP7/test_mses_500_trees.rds")
data_1000_trees <- readRDS("../data/Optimise_RF_ASV_ATP7/test_mses_1000_trees.rds")
data_3000_trees <- readRDS("../data/Optimise_RF_ASV_ATP7/test_mses_3000_trees.rds")
data_5000_trees <- readRDS("../data/Optimise_RF_ASV_ATP7/test_mses_5000_trees.rds")
data_7000_trees <- readRDS("../data/Optimise_RF_ASV_ATP7/test_mses_7000_trees.rds")
data_9000_trees <- readRDS("../data/Optimise_RF_ASV_ATP7/test_mses_9000_trees.rds")
data_11000_trees <- readRDS("../data/Optimise_RF_ASV_ATP7/test_mses_11000_trees.rds")
data_14000_trees <- readRDS("../data/Optimise_RF_ASV_ATP7/test_mses_14000_trees.rds")
data_18000_trees <- readRDS("../data/Optimise_RF_ASV_ATP7/test_mses_18000_trees.rds")
data_22000_trees <- readRDS("../data/Optimise_RF_ASV_ATP7/test_mses_22000_trees.rds")

# Tree OOB data frame
OOB_df <- data.frame(trees_100 = data_100_trees,
                     trees_250 = data_250_trees,
                     trees_500 = data_500_trees,
                     trees_1000 = data_1000_trees,
                     trees_3000 = data_3000_trees,
                     trees_5000 = data_5000_trees,
                     trees_7000 = data_7000_trees,
                     trees_9000 = data_9000_trees,
                     trees_11000 = data_11000_trees,
                     trees_14000 = data_14000_trees,
                     trees_18000 = data_18000_trees,
                     trees_22000 = data_22000_trees)

# Getting rid of tree numbers row
OOB_df <- OOB_df[-1,]
# Transposing
OOB_df <- data.frame(t(OOB_df))
# Calculating mean mse across the 25 tests for each tree number
OOB_df$MeanMSE <- rowMeans(OOB_df)

######### TREE NUMBER WITH LOWEST MSE IS 22000 trees (the highest number of trees)
# BUT 10,000 clearly enough as not much difference after 9,000 trees




