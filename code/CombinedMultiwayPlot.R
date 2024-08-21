
library(ggplot2)
library(patchwork)
library(ggrepel)


########################### ATP #########################################################################
importance.df <- read.table("../results/Optimised_ATP7_ASV/RFEVarImp_Class-log.ATP7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")

# Measure importance of variables by MSE, Node purity and mean min depth
top10asvs <- important_variables(importance.df,
                                 k = 10,
                                 measures = c("mean_min_depth", "mse_increase", "node_purity_increase"),
                                 ties_action = "draw")
importance.df$colour_group <- ifelse(importance.df$variable %in% top10asvs, "Top 10 ASVs", "Other ASVs")
top_10_data <- importance.df[importance.df$variable %in% top10asvs, ]

## Plot of multiway importance
p1 <- ggplot(importance.df, aes(x = mean_min_depth, y = node_purity_increase)) +
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
    size = 5,  # Adjusted size to match font size 12
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
  scale_y_continuous(
    limits = c(0, 350),  # Set y-axis limits
    breaks = seq(0, 350, by = 50)  # Set y-axis intervals to go up by 50
  ) +
  xlab("Mean minimum depth") +  # Update x-axis label
  ylab("Node purity increase") +  # Update y-axis label
  ggtitle("A) ATP") +  # Add the title
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title.x = element_text(size = 14),  # X-axis title size
    axis.title.y = element_text(size = 14),  # Y-axis title size
    axis.text.x = element_text(size = 12),   # X-axis text size
    axis.text.y = element_text(size = 12),   # Y-axis text size
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),   
    legend.key.size = unit(1.5, "lines"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)  # Add border around the plot
  )

########################### Xylosidase #########################################################################
importance.df <- read.table("../results/Optimised_XYLO7_ASV/RFEVarImp_Class-log.mX7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")

# Measure importance of variables by MSE, Node purity and mean min depth
top10asvs <- important_variables(importance.df,
                                 k = 10,
                                 measures = c("mean_min_depth", "mse_increase", "node_purity_increase"),
                                 ties_action = "draw")
importance.df$colour_group <- ifelse(importance.df$variable %in% top10asvs, "Top 10 ASVs", "Other ASVs")
top_10_data <- importance.df[importance.df$variable %in% top10asvs, ]

## Plot of multiway importance
p2 <- ggplot(importance.df, aes(x = mean_min_depth, y = node_purity_increase)) +
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
    size = 5,  # Adjusted size to match font size 12
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
  scale_y_continuous(
    limits = c(0, 350),  # Set y-axis limits
    breaks = seq(0, 350, by = 50)  # Set y-axis intervals to go up by 50
  ) +
  xlab("Mean minimum depth") +  # Update x-axis label
  ylab("Node purity increase") +  # Update y-axis label
  ggtitle("B) Xylosidase") +  # Add the title
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title.x = element_text(size = 14),  # X-axis title size
    axis.title.y = element_text(size = 14),  # Y-axis title size
    axis.text.x = element_text(size = 12),   # X-axis text size
    axis.text.y = element_text(size = 12),   # Y-axis text size
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),   
    legend.key.size = unit(1.5, "lines"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)  # Add border around the plot
  )

########################### Glucosidase #########################################################################
importance.df <- read.table("../results/OptimisedGLUC7_ASV/RFEVarImp_Class-log.mG7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")

# Measure importance of variables by MSE, Node purity and mean min depth
top10asvs <- important_variables(importance.df,
                                 k = 10,
                                 measures = c("mean_min_depth", "mse_increase", "node_purity_increase"),
                                 ties_action = "draw")
importance.df$colour_group <- ifelse(importance.df$variable %in% top10asvs, "Top 10 ASVs", "Other ASVs")
top_10_data <- importance.df[importance.df$variable %in% top10asvs, ]

## Plot of multiway importance
p3 <- ggplot(importance.df, aes(x = mean_min_depth, y = node_purity_increase)) +
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
    size = 5,  # Adjusted size to match font size 12
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
  scale_y_continuous(
    limits = c(0, 350),  # Set y-axis limits
    breaks = seq(0, 350, by = 50)  # Set y-axis intervals to go up by 50
  ) +
  xlab("Mean minimum depth") +  # Update x-axis label
  ylab("Node purity increase") +  # Update y-axis label
  ggtitle("C) Glucosidase") +  # Add the title
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title.x = element_text(size = 14),  # X-axis title size
    axis.title.y = element_text(size = 14),  # Y-axis title size
    axis.text.x = element_text(size = 12),   # X-axis text size
    axis.text.y = element_text(size = 12),   # Y-axis text size
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),   
    legend.key.size = unit(1.5, "lines"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)  # Add border around the plot
  )

########################### Chitinase #########################################################################
importance.df <- read.table("../results/Optimised_CHIT7_ASV/RFEVarImp_Class-log.mN7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")

# Measure importance of variables by MSE, Node purity and mean min depth
top10asvs <- important_variables(importance.df,
                                 k = 10,
                                 measures = c("mean_min_depth", "mse_increase", "node_purity_increase"),
                                 ties_action = "draw")
importance.df$colour_group <- ifelse(importance.df$variable %in% top10asvs, "Top 10 ASVs", "Other ASVs")
top_10_data <- importance.df[importance.df$variable %in% top10asvs, ]

## Plot of multiway importance
p4 <- ggplot(importance.df, aes(x = mean_min_depth, y = node_purity_increase)) +
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
    size = 5,  # Adjusted size to match font size 12
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
  scale_y_continuous(
    limits = c(0, 350),  # Set y-axis limits
    breaks = seq(0, 350, by = 50)  # Set y-axis intervals to go up by 50
  ) +
  xlab("Mean minimum depth") +  # Update x-axis label
  ylab("Node purity increase") +  # Update y-axis label
  ggtitle("D) Chitinase") +  # Add the title
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title.x = element_text(size = 14),  # X-axis title size
    axis.title.y = element_text(size = 14),  # Y-axis title size
    axis.text.x = element_text(size = 12),   # X-axis text size
    axis.text.y = element_text(size = 12),   # Y-axis text size
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),   
    legend.key.size = unit(1.5, "lines"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)  # Add border around the plot
  )

########################### Phosphatase #########################################################################
importance.df <- read.table("../results/Optimised_PHOS7_ASV/RFEVarImp_Class-log.mP7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")

# Measure importance of variables by MSE, Node purity and mean min depth
top10asvs <- important_variables(importance.df,
                                 k = 10,
                                 measures = c("mean_min_depth", "mse_increase", "node_purity_increase"),
                                 ties_action = "draw")
importance.df$colour_group <- ifelse(importance.df$variable %in% top10asvs, "Top 10 ASVs", "Other ASVs")
top_10_data <- importance.df[importance.df$variable %in% top10asvs, ]

## Plot of multiway importance
p5 <- ggplot(importance.df, aes(x = mean_min_depth, y = node_purity_increase)) +
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
    size = 5,  # Adjusted size to match font size 12
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
  scale_y_continuous(
    limits = c(0, 350),  # Set y-axis limits
    breaks = seq(0, 350, by = 50)  # Set y-axis intervals to go up by 50
  ) +
  xlab("Mean minimum depth") +  # Update x-axis label
  ylab("Node purity increase") +  # Update y-axis label
  ggtitle("E) Phosphatase") +  # Add the title
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title.x = element_text(size = 14),  # X-axis title size
    axis.title.y = element_text(size = 14),  # Y-axis title size
    axis.text.x = element_text(size = 12),   # X-axis text size
    axis.text.y = element_text(size = 12),   # Y-axis text size
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),   
    legend.key.size = unit(1.5, "lines"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)  # Add border around the plot
  )

########################### Carbon dioxide #########################################################################
importance.df <- read.table("../results/Optimised_CO7_ASV/RFEVarImp_Class-log.mCO2.7.norm_ASV_AllRepTrain.tsv", header = TRUE, sep = "\t")

# Measure importance of variables by MSE, Node purity and mean min depth
top10asvs <- important_variables(importance.df,
                                 k = 10,
                                 measures = c("mean_min_depth", "mse_increase", "node_purity_increase"),
                                 ties_action = "draw")
importance.df$colour_group <- ifelse(importance.df$variable %in% top10asvs, "Top 10 ASVs", "Other ASVs")
top_10_data <- importance.df[importance.df$variable %in% top10asvs, ]

## Plot of multiway importance
p6 <- ggplot(importance.df, aes(x = mean_min_depth, y = node_purity_increase)) +
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
    size = 5,  # Adjusted size to match font size 12
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
  scale_y_continuous(
    limits = c(0, 350),  # Set y-axis limits
    breaks = seq(0, 350, by = 50)  # Set y-axis intervals to go up by 50
  ) +
  xlab("Mean minimum depth") +  # Update x-axis label
  ylab("Node purity increase") +  # Update y-axis label
  ggtitle("F) Carbon dioxide") +  # Add the title
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Title size and style
    axis.title.x = element_text(size = 14),  # X-axis title size
    axis.title.y = element_text(size = 14),  # Y-axis title size
    axis.text.x = element_text(size = 12),   # X-axis text size
    axis.text.y = element_text(size = 12),   # Y-axis text size
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),   
    legend.key.size = unit(1.5, "lines"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)  # Add border around the plot
  )

###############################################################################################################
###############################################################################################################

combined_plot <- (p1 | p2 | p3) / (p4 | p5 | p6) +
  plot_layout(guides = 'collect') # Collects legends and applies a common legend

# Modify each plot to remove the legend
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")

# Combine all plots with one shared legend
combined_plot <- (p1 | p2 | p3) / (p4 | p5 | p6) +
  plot_layout(guides = 'collect') +
  plot_annotation(title = 'Combined Plots') # Optional: Add a title to the combined plot









