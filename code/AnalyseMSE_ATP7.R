
fileOut=paste("varMSEINCR_Class-",select.class,select.cond,
              "_",predictor.unit,"_",replicate.train,
              ".tsv",sep="") # Name of file


# Loading in MSE increase values for ASVs and clusters - each df has column for variable (asvs or clusters) and column for mse increase when removed
asv_mse_data <- read.table("../data/varMSEINCR_Class-log.ATP7.norm_ASV_AllRepTrain.tsv", sep="")
cluster_mse_data <- read.table("../data/varMSEINCR_Class-log.ATP7.norm_Cluster_AllRepTrain.tsv", sep="")
cluster_data = read.table("../data/max_tot_ext_network_table.tsv", sep = "\t")

## Merging ASV mse data with cluster data (cluster for each asv)
# Making it so that ASV variables have exact same name as in cluster data by adding 'ASV_' to them
for (i in 1:nrow(asv_mse_data)){
  if (is.numeric(asv_mse_data[i, 1])){
    asv_mse_data[i,1] <- paste0("ASV_", asv_mse_data[i,1])
  }
}
# Making name of variable column into 'ASV'
names(asv_mse_data)[names(asv_mse_data) == "variable"] <- "ASV"
# Getting rid of 'set' column from cluster data
cluster_data = cluster_data[,1:2]
# Merging by ASV such that ave column for cluster, column for ASV, column for mse increase
merged_mse_data <- merge(asv_mse_data, cluster_data, by = "ASV", all.x = TRUE)
# Remove ASV column
merged_mse_data <- merged_mse_data[, c("functionInk", "mse_increase")]
# Summing mses for all ASVs in each cluster
merged_mse_data <- merged_mse_data %>%
  group_by(functionInk) %>%
  summarize(across(everything(), \(x) sum(x, na.rm = TRUE)), .groups = 'drop')

# merged_mse_data has one column with each cluster, and one column with its mse_increase as a sum of the mse_increase of each of its ASVs from the ASV random forest
# cluster_mse_data has one column with each cluster, and one column with its mse_increase as was found from the cluster random forest

## Merging these two data frames
# Making column names the same or different
names(cluster_mse_data)[names(cluster_mse_data) == "variable"] <- "functionInk"
names(cluster_mse_data)[names(cluster_mse_data) == "mse_increase"] <- "mse_increase_cluster"
names(merged_mse_data)[names(merged_mse_data) == "mse_increase"] <- "mse_increase_sum_asvs"
# Merging
merged_2_mse_data <- merge(asv_mse_data, cluster_mse_data, by = "functionInk", all.y = TRUE)
# Column for difference
merged_2_mse_data$cluster_minus_asvsum_mse <- merged_2_mse_data$mse_increase_cluster - merged_2_mse_data$mse_increase_sum_asvs

# merged_2_mse_data has a column for clusters, a column for mse increase of clusters, and column for mse increase of sum of ASVs in clusters
# Also has column for difference between the two mse columns

