

# Getting top 10 vars from the following RFs:
# - ASV abundance RFs for each function
# - high degree abundance RFs for each function
# - subgroup abundance RFs for each function
# - large subgroup abundance RFs for each function


# Importing importance frames
importance.df.1 <- read.delim("../results/Optimised_ATP7_ASV/RFEVarImp_Class-log.ATP7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.2 <- read.delim("../results/Optimised_XYLO7_ASV/RFEVarImp_Class-log.mX7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.3 <- read.delim("../results/Optimised_GLUC7_ASV/RFEVarImp_Class-log.mG7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.4 <- read.delim("../results/Optimised_CHIT7_ASV/RFEVarImp_Class-log.mN7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.5 <- read.delim("../results/Optimised_PHOS7_ASV/RFEVarImp_Class-log.mP7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.6 <- read.delim("../results/Optimised_CO7_ASV/RFEVarImp_Class-log.μgCO2.7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")

importance.df.1a <- read.delim("../results/DegreeRF_ATP7/DEGREERFEVarImp_Class-log.ATP7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.2a <- read.delim("../results/DegreeRF_XYLO7/DEGREERFEVarImp_Class-log.mX7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.3a <- read.delim("../results/DegreeRF_GLUC7/DEGREERFEVarImp_Class-log.mG7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.4a <- read.delim("../results/DegreeRF_CHIT7/DEGREERFEVarImp_Class-log.mN7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.5a <- read.delim("../results/DegreeRF_PHOS7/DEGREERFEVarImp_Class-log.mP7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.6a <- read.delim("../results/DegreeRF_CO7/DEGREERFEVarImp_Class-log.mCO2.7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")

importance.df.1b <- read.delim("../results/Optimised_ATP7_Cluster/RFEVarImp_Class-log.ATP7.norm_Cluster_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.2b <- read.delim("../results/Optimised_XYLO7_Cluster/RFEVarImp_Class-log.mX7.norm_Cluster_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.3b <- read.delim("../results/Optimised_GLUC7_Cluster/RFEVarImp_Class-log.mG7.norm_Cluster_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.4b <- read.delim("../results/Optimised_CHIT7_Cluster/RFEVarImp_Class-log.mN7.norm_Cluster_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.5b <- read.delim("../results/Optimised_PHOS7_Cluster/RFEVarImp_Class-log.mP7.norm_Cluster_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.6b <- read.delim("../results/Optimised_CO7_Cluster/RFEVarImp_Class-log.μgCO2.7.norm_Cluster_AllRepTrain.tsv", header = TRUE, sep = "\t")

importance.df.1c <- read.delim("../results/cluster large/RFEVarImp_Class-log.ATP7.norm_Cluster_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.2c <- read.delim("../results/cluster large/RFEVarImp_Class-log.mX7.norm_Cluster_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.3c <- read.delim("../results/cluster large/RFEVarImp_Class-log.mG7.norm_Cluster_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.4c <- read.delim("../results/cluster large/RFEVarImp_Class-log.mN7.norm_Cluster_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.5c <- read.delim("../results/cluster large/RFEVarImp_Class-log.mP7.norm_Cluster_AllRepTrain.tsv", header = TRUE, sep = "\t")
importance.df.6c <- read.delim("../results/cluster large/RFEVarImp_Class-log.mCO2.7.norm_Cluster_AllRepTrain.tsv", header = TRUE, sep = "\t")


# CALCULATING TOP 10 MOST IMPORTANT VARIABLES
top10asvs.1 <- important_variables(importance.df.1, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.2 <- important_variables(importance.df.2, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.3 <- important_variables(importance.df.3, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.4 <- important_variables(importance.df.4, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.5 <- important_variables(importance.df.5, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.6 <- important_variables(importance.df.6, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")

top10asvs.1a <- important_variables(importance.df.1a, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.2a <- important_variables(importance.df.2a, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.3a <- important_variables(importance.df.3a, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.4a <- important_variables(importance.df.4a, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.5a <- important_variables(importance.df.5a, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.6a <- important_variables(importance.df.6a, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")

top10asvs.1b <- important_variables(importance.df.1b, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.2b <- important_variables(importance.df.2b, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.3b <- important_variables(importance.df.3b, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.4b <- important_variables(importance.df.4b, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.5b <- important_variables(importance.df.5b, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.6b <- important_variables(importance.df.6b, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")

top10asvs.1c <- important_variables(importance.df.1c, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.2c <- important_variables(importance.df.2c, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.3c <- important_variables(importance.df.3c, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.4c <- important_variables(importance.df.4c, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.5c <- important_variables(importance.df.5c, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")
top10asvs.6c <- important_variables(importance.df.6c, k = 10, measures = c("mean_min_depth", "mse_increase", "node_purity_increase"), ties_action = "draw")


# MAKING DATA FRAME WITH TOP 10 MOST IMPORTANT VARIABLES IN EACH RF in alphabetical order
AllImportantVars <- data.frame(Alltaxa_ATP = top10asvs.1,
                               Alltaxa_XYLO = top10asvs.2,
                               Alltaxa_GLUC = top10asvs.3,
                               Alltaxa_CHIT = top10asvs.4,
                               Alltaxa_PHOS = top10asvs.5,
                               Alltaxa_CO = top10asvs.6,
                               
                               Alltaxa_ATP = top10asvs.1a,
                               Alltaxa_XYLO = top10asvs.2a,
                               Alltaxa_GLUC = top10asvs.3a,
                               Alltaxa_CHIT = top10asvs.4a,
                               Alltaxa_PHOS = top10asvs.5a,
                               Alltaxa_CO = top10asvs.6a,
                               
                               Alltaxa_ATP = top10asvs.1b,
                               Alltaxa_XYLO = top10asvs.2b,
                               Alltaxa_GLUC = top10asvs.3b,
                               Alltaxa_CHIT = top10asvs.4b,
                               Alltaxa_PHOS = top10asvs.5b,
                               Alltaxa_CO = top10asvs.6b,
                               
                               Alltaxa_ATP = top10asvs.1c,
                               Alltaxa_XYLO = top10asvs.2c,
                               Alltaxa_GLUC = top10asvs.3c,
                               Alltaxa_CHIT = top10asvs.4c,
                               Alltaxa_PHOS = top10asvs.5c,
                               Alltaxa_CO = top10asvs.6c)


# Getting number of times each variable appears across all RFs
all_values <- unlist(AllImportantVars)
value_counts <- as.data.frame(table(all_values))

# Comparing only the interaction and all taxa
ImportantASVs <- data.frame(Alltaxa_ATP = top10asvs.1,
                               Alltaxa_XYLO = top10asvs.2,
                               Alltaxa_GLUC = top10asvs.3,
                               Alltaxa_CHIT = top10asvs.4,
                               Alltaxa_PHOS = top10asvs.5,
                               Alltaxa_CO = top10asvs.6,
                               
                               Alltaxa_ATP = top10asvs.1a,
                               Alltaxa_XYLO = top10asvs.2a,
                               Alltaxa_GLUC = top10asvs.3a,
                               Alltaxa_CHIT = top10asvs.4a,
                               Alltaxa_PHOS = top10asvs.5a,
                               Alltaxa_CO = top10asvs.6a)
all_values_asvs <- unlist(ImportantASVs)
asv_counts <- as.data.frame(table(all_values_asvs))

# Finding out which RFs taxa were important in
# The value to search for e.g. "ASV_29" or "36"
value_to_search <- "ASV_58" # The value to search for
logical_matrix <- AllImportantVars == value_to_search # logical matrix where TRUE indicates the presence of the value
columns_with_value <- colnames(AllImportantVars)[apply(logical_matrix, 2, any)] # column names where the value is found
print(columns_with_value)




# Converting subgroups into taxa
cluster_data <- read.csv(file = "../data/max_tot_ext_network_table.tsv", sep = "\t") # data of which asvs in subgroups
# Only keeping rows where functionInk column is in important subgroups
top10clusterasvs_1b <- subset(cluster_data, functionInk %in% top10asvs.1b)
top10clusterasvs_1b <- top10clusterasvs_1b[, "ASV"]
top10clusterasvs_2b <- subset(cluster_data, functionInk %in% top10asvs.2b)
top10clusterasvs_2b <- top10clusterasvs_2b[, "ASV"]
top10clusterasvs_3b <- subset(cluster_data, functionInk %in% top10asvs.3b)
top10clusterasvs_3b <- top10clusterasvs_3b[, "ASV"]
top10clusterasvs_4b <- subset(cluster_data, functionInk %in% top10asvs.4b)
top10clusterasvs_4b <- top10clusterasvs_4b[, "ASV"]
top10clusterasvs_5b <- subset(cluster_data, functionInk %in% top10asvs.5b)
top10clusterasvs_5b <- top10clusterasvs_5b[, "ASV"]
top10clusterasvs_6b <- subset(cluster_data, functionInk %in% top10asvs.6b)
top10clusterasvs_6b <- top10clusterasvs_6b[, "ASV"]

max_length <- max(length(top10clusterasvs_1b), length(top10clusterasvs_2b), 
                  length(top10clusterasvs_3b), length(top10clusterasvs_4b), 
                  length(top10clusterasvs_5b), length(top10clusterasvs_6b))
length(top10clusterasvs_1b) <- max_length
length(top10clusterasvs_2b) <- max_length
length(top10clusterasvs_3b) <- max_length
length(top10clusterasvs_4b) <- max_length
length(top10clusterasvs_5b) <- max_length
length(top10clusterasvs_6b) <- max_length
subgroup_asvs <- data.frame(ATPSub = top10clusterasvs_1b,
                            XYLOSub = top10clusterasvs_2b,
                            GLUCSub = top10clusterasvs_3b,
                            CHITSub = top10clusterasvs_4b,
                            PHOSSub = top10clusterasvs_5b,
                            COSub = top10clusterasvs_6b)
all_values_subgroups <- unlist(subgroup_asvs)
subgroup_counts <- as.data.frame(table(all_values_subgroups))


# Finding out which RFs subgroup taxa were important in
# The value to search for e.g. "ASV_29" or "36"
value_to_search <- "ASV_29" # The value to search for
logical_matrix <- subgroup_asvs == value_to_search # logical matrix where TRUE indicates the presence of the value
columns_with_value <- colnames(subgroup_asvs)[apply(logical_matrix, 2, any)] # column names where the value is found
print(columns_with_value)




########## TAXA INFO - INTERESTED IN TAXA THAT ARE IN BOTH INTERACTING AND ALL, OR THAT ARE IN ALL TAXA BUT NOT IN INTERACTING

### TAXA IN ALL TAXA, INTERACTING, IMPORTANT SUBGROUP
# asv29 in all except all taxa carbon dioxide # BOTH ATP,XYLO,GLUC,CHIT,PHOS - subgroup 36
# asv23 in all except all taxa carbon dioxide # BOTH ATP,XYLO,GLUC,CHIT,PHOS - subgroup 36
# asv44 in all except all taxa atp # BOTH XYLO,GLUC,CHIT,PHOS,CO - subgroup 161
# asv 24 found in all except both atp rfs # BOTH XYLO,GLUC,CHIT,PHOS,CO. SUBGROUP NOT IMPORTANT FOR GLUC - subgroup 119
# asv6 in all taxa co and in interacting atp,xylo,chit,co # BOTH CO - subgroup 148

### TAXA IN ALL TAXA, INTERACTING, NO SUBGROUP
# asv18 in all except all taxa atp # BOTH XYLO,GLUC,CHIT,PHOS,CO
# asv9 found in all except all taxa atp, all taxa gluc, interacting taxa atp # BOTH XYLO,CHIT,PHOS,CO
# asv1 in all taxa xylo, chit, co and interacting taxa xylo, chit, phos # BOTH XYLO,CHIT
# asv16 in co and interacting co # BOTH CO
# asv32 in all taxa chit, and in interacting taxa atp,xylo,gluc,chit,phos BOTH CHIT

### TAXA IN ALL TAXA, NOT INTERACTING, IMPORTANT SUBGROUP
# asv28 in all taxa atp, cylo, gluc, chit, phot, co but not in any of interacting taxa # IN ATP,XYLO,GLUC,CHIT,PHOS,CO NO INTERACTING - subgroup 161
# asv30 in all taxa atp,xylo,gluc,chit,co not in any of interacting # ATP,XYLO,GLUC,CHIT,COS NO INTERACTING - subgroup 161
# asv58 in all taxa atp and xylo # ATP, XYLO NO INTERACCTING - subgroup 177
# 104 in alltaxa phos # PHOS NO INT - subgroup 115
# 107 in all taxa phos # PHOS NO INT - subgroup 198
# 14 in all taxa co # CO NO INT - subgroup 170
# 17 in all taxa atp # ATP NO INT - subgroup 7
# 38 in all taxa gluc # GLUC NO INT - subgroup 151
# 42 in all taxa phos # PHOS NO INT - subgroup 187
# 50 in all taxa atp # ATP NO INT - subgroup 155
# 80 alltaxa phos # PHOS NO INT - subgroup 141
# 92 all taxa atp # ATP NO INT - subgroup 43

### TAXA IN ALL TAXA, NOT INTERACTING, NO SUBGROUP
# 2 in all taxa atp # ATP NO INT
# 35 in all taxa co # CO NO INT
# 43 in all taxa gluc # GLUC NO INT

### Non-important taxa
# asv7 in interacting atp, xylo, gluc, chit
# asv274 in interacting phos and co
# 45 in INTERACTING TAXA gluc
# 5 in INTERACTING taxa atp
# 54 interact atp
# 55 interact gluc
# 8 interact atp

### IMPORTANT SUBGROUPS
## RF 6 ##
# 155 for all functions
# 36 for all function
# 161 for all functions
## RF 5 ##
# 148 for all except ATP
## RF 4 ##
# 119 for XYLO, CHIT, PHOS, CO
# 177 for ATP, XYLO, CHIT, PHOS
# 43 FOR ATP, XYLO, CHIT, PHOS
## RF 3 ##
# 2 FOR XYLO, GLUC, CHIT
# 84 FOR XYLO, CHIT, PHOS
# COunts for XYLO, GLUC, CHIT
## RF 2 ##
# 187 GLUC, PHOS
# 7 ATP, CO
## RF 1 ##
# 108 GLUC
# 141 PHOS
# 151 GLUC
# 170 CO
# 179 GLUC
# 185 CO
# 199 CO
# 24 ATP
# 30 CO
# 49 ATP
# 5 ATP
# 74 ATP



### IMPORTANT SUBGROUPS, WITHOUT IMPORTANT TAXA
subgroups_with_imptaxa <- c("7", "36", "43", "115", "119", "141", "148", "151", "155", "161", "170", "177", "187", "198")
ImportantSubs <- data.frame(ATP = top10asvs.1b,
                            XYLO = top10asvs.2b,
                            GLUC = top10asvs.3b,
                            CHIT = top10asvs.4b,
                            PHOS = top10asvs.5b,
                            CO = top10asvs.6b)
ImportantSubs1 <- as.data.frame(apply(ImportantSubs, 2, function(x) ifelse(x %in% subgroups_with_imptaxa, NA, x)))
ImportantSubs2 <- as.data.frame(apply(ImportantSubs, 2, function(x) ifelse(x %in% subgroups_with_imptaxa, x, NA)))
# sub 24, 49, 5, 74 important ATP
# 84, 2, counts important xylo
# 108, counts, 2, 179 sub important gluc
# 84, 2, counts important chit
# 84 important phos
# 30, 199, 185 important co
# Total: 2, 5, 24, 30, 49, 74, 84, 108, 179, 185, 199
# 2 size 3
# 5 size 2
# 24 size 3
# 30 size 3
# 49 size 2
# 74 size 2
# 84 size 2
# 108 size 3
# 179 size 4
# 185 size 4
# 199 size 2

# 7 s2
# 36 s9
# 43 s3
# 115 s2
# 119 s2
# 141 s2
# 148 s2
# 151 s2
# 155 s4
# 161 s5
# 170 s2
# 177 s2
# 187 s2
# 198 s1





# Only comparing diff functions clusters
KeySubs <- data.frame(Alltaxa_ATP = top10asvs.1b,
                               Alltaxa_XYLO = top10asvs.2b,
                               Alltaxa_GLUC = top10asvs.3b,
                               Alltaxa_CHIT = top10asvs.4b,
                               Alltaxa_PHOS = top10asvs.5b,
                               Alltaxa_CO = top10asvs.6b)
# Getting number of times each variable appears across all RFs
key_sub_values <- unlist(KeySubs)
KeySubs_counts <- as.data.frame(table(key_sub_values))
# Finding out which RFs subgroup taxa were important in
# The value to search for e.g. "ASV_29" or "36"
value_to_search <- "36" # The value to search for
logical_matrix <- KeySubs == value_to_search # logical matrix where TRUE indicates the presence of the value
columns_with_value <- colnames(KeySubs)[apply(logical_matrix, 2, any)] # column names where the value is found
print(columns_with_value)






network_data <- read.csv(file = "../data/combined_network_metrics.csv")
cluster_sizes <- as.data.frame(table(network_data$functionInk)) # number of times each cluster appears in table
colnames(cluster_sizes) <- c("Subgroup", "Size") # renaming
cluster_sizes <- cluster_sizes[cluster_sizes$Size > 1, ]



source("../code/clean_ASV_table.R")
ASV.table=read.table(file = "../data/seqtable_readyforanalysis.csv", sep="\t")
sample_md <-read.table(file = "../data/metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv", sep="\t", header=TRUE)
# --- Reshaping data parameters
nreads = 10000 # minimum number of reads to consider a sample
exclude_exp = c("4M") # A vector of characters with the experiments that should be excluded
match_exp = TRUE # Set to true if only starting communities that were resurrected should be included
output.label = "Time0D_7D_matched" 
# CLEANING DATA - Removing experiment 4M, keeping only final samples, and only keeping samples present in the metadata
clean.data.list = clean_ASV_table(ASV.table ,sample_md , match_exp = TRUE, exclude_exp = c("4M"),
                                  nreads = 10000) # Outputs list containing cleaned asv table and meta data
ASV.table = clean.data.list$ASV.table
sample_md = clean.data.list$sample_md
id.t7 = grep("7D", sample_md$Experiment) # Indices of rows whcih are final communities
samples.t7 = as.character(sample_md$sampleid[id.t7]) # Sampleids of final communities
xIn.tmp = ASV.table[samples.t7, ] # ASV table containing only final communities
TotalAbunds = rowSums(xIn.tmp) # add the sum of the counts for each sample as an additional predictor
MeanAbunds = rowMeans(xIn.tmp)
meanabundance = mean(MeanAbunds)



















# IMPORTANT TAXA (IN ALL TAXA RFS)
"
24 important taxa
28: 20
30: 20
23: 36
24: 24
29: 64
44: 23
18: 31
9: 27
1: 61
58: 20
107: 17
14: 32
16: 36
17: 12
2: 0
32: 25
38: 17
42: 17
43: 18
50: 0
52: 19
6: 30
80: 20
92: 19
"



















