


# Loading files containing mses for the 25 tests of each variable number
data_10_vars <- readRDS("../data/Optimise_RF_ASV_ATP7/log.ATP7.norm_ASV_test_var_mses_10_vars.rds")
data_20_vars <- readRDS("../data/Optimise_RF_ASV_ATP7/log.ATP7.norm_ASV_test_var_mses_20_vars.rds")
data_30_vars <- readRDS("../data/Optimise_RF_ASV_ATP7/log.ATP7.norm_ASV_test_var_mses_30_vars.rds")
data_50_vars <- readRDS("../data/Optimise_RF_ASV_ATP7/log.ATP7.norm_ASV_test_var_mses_50_vars.rds")
data_100_vars <- readRDS("../data/Optimise_RF_ASV_ATP7/log.ATP7.norm_ASV_test_var_mses_100_vars.rds")
data_200_vars <- readRDS("../data/Optimise_RF_ASV_ATP7/log.ATP7.norm_ASV_test_var_mses_200_vars.rds")
data_500_vars <- readRDS("../data/Optimise_RF_ASV_ATP7/log.ATP7.norm_ASV_test_var_mses_500_vars.rds")
data_750_vars <- readRDS("../data/Optimise_RF_ASV_ATP7/log.ATP7.norm_ASV_test_var_mses_750_vars.rds")
data_1000_vars <- readRDS("../data/Optimise_RF_ASV_ATP7/log.ATP7.norm_ASV_test_var_mses_1000_vars.rds")
data_1400_vars <- readRDS("../data/Optimise_RF_ASV_ATP7/log.ATP7.norm_ASV_test_var_mses_1400_vars.rds")

# Tree OOB data frame
OOB_df <- data.frame(vars_10 = data_10_vars,
                     vars_20 = data_20_vars,
                     vars_30 = data_30_vars,
                     vars_50 = data_50_vars,
                     vars_100 = data_100_vars,
                     vars_200 = data_200_vars,
                     vars_500 = data_500_vars,
                     vars_750 = data_750_vars,
                     vars_1000 = data_1000_vars,
                     vars_1400 = data_1400_vars
                     )

# Getting rid of variable numbers row
OOB_df <- OOB_df[-1,]
# Transposing
OOB_df <- data.frame(t(OOB_df))
# Calculating mean mse across the 25 tests for each variable number
OOB_df$MeanMSE <- rowMeans(OOB_df)

######### VARIABLE NUMBER WITH LOWEST MSE IS 500 variables
# Not huge difference overall.
# Strange pattern - no clear pattern of reduced mse as increase variable, seems a bit random
# 750 second lowest MSE, seems reasonable




