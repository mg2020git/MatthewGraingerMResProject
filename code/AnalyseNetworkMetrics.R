



# Loading packages
library(ggplot2)

# Loading data
network_data <- read.csv(file = "../data/combined_network_metrics.csv")
# Calculate the 95th percentile
percentile_95 <- quantile(network_data$Degree, 0.95)

# Save plot to PDF
pdf("DegreeDistr.pdf")
hist(network_data$Degree, 
     main = "Distribution of number of interactions",
     xlab = "Number of interactions", 
     ylab = "Frequency", 
     col = "lightblue", 
     border = "black",
     xlim = c(0, 65), 
     ylim = c(0, 500), 
     breaks = 20,
     las = 1,                # Make axis labels horizontal
     cex.axis = 1.2,         # Axis label size
     cex.lab = 1.5,          # Axis title size
     cex.main = 1.6,         # Title size
     font.main = 2)          # Bold title
axis(1, at = seq(0, 65, by = 5), cex.axis = 1.2)  # Customize x-axis intervals
# Add a vertical line at the 95th percentile
abline(v = percentile_95, col = "red", lwd = 2, lty = 2)
# Add annotation
text(x = percentile_95 + 2, y = 450, 
     labels = paste("95th percentile:", round(percentile_95, 2)), 
     col = "red", cex = 1.2, adj = 0)
dev.off()














# Calculate the 95th percentile of the original (non-log) degree values
percentile_95 <- quantile(network_data$Degree, 0.95)

# Create a histogram without plotting it to extract mids (bin centers) and counts
hist_data <- hist(network_data$Degree, plot = FALSE, breaks = 20)

# Log-transform the degree values (x-axis)
log_degrees <- log(hist_data$mids)

# Log-transform the frequency counts (y-axis)
log_counts <- log(hist_data$counts + 1)  # Add 1 to avoid log(0) issue

# Fit a linear model to the log-log data
linear_model <- lm(log_counts ~ log_degrees)

# Calculate the log of the 95th percentile of the degree
log_percentile_95 <- log(percentile_95)

# Save plot to PDF
pdf("DegreeDistr_LOG_with_LinearFit.pdf")
plot(log_degrees, log_counts, 
     type = "p",             # Plot points
     pch = 20,               # Solid circles for points
     col = "black",          # Color of the points
     xlab = "Log(Number of interactions)", 
     ylab = "Log(Frequency)",
     las = 1,                # Make axis labels horizontal
     cex.axis = 1.2,         # Axis label size
     cex.lab = 1.5,          # Axis title size
     ylim = c(0, 7),
     xlim = c(1.0, 5.0),
     main = "Log-Log Plot of Distribution of number\nof interactions with Power-Law Fit",
     cex.main = 1.6,         # Title size
     font.main = 2)          # Bold title
# Add the linear model fit line
abline(linear_model, col = "blue", lwd = 2)
# Add a vertical line at the log-transformed 95th percentile
abline(v = log_percentile_95, col = "red", lwd = 2, lty = 2)
# Add annotation for the log-transformed 95th percentile
text(x = log_percentile_95 + 0.1, y = max(log_counts) - 0.5, 
     labels = paste("95th percentile:", round(log_percentile_95, 2)), 
     col = "red", cex = 1.2, adj = 0)
# Display the slope (alpha) and intercept on the plot
coeffs <- coef(linear_model)
text(x = min(log_degrees) + 0.1, y = min(log_counts) + 0.5, 
     labels = paste("Log(Frequency) =", round(coeffs[2], 2), "* Log(Degree) +", round(coeffs[1], 2)),
     col = "blue", cex = 1.2, adj = 0)
dev.off()





























########### Cluster size distribution
cluster_sizes <- as.data.frame(table(network_data$functionInk)) # number of times each cluster appears in table
colnames(cluster_sizes) <- c("Subgroup", "Size") # renaming
cluster_sizes <- cluster_sizes[cluster_sizes$Size > 1, ]
subgroup_percentile_95 <- quantile(cluster_sizes$Size, 0.95) # 95th percentile of cluster sizes

p2 <- hist(cluster_sizes$Size, 
           main = NULL,
           xlab = "Subgroup size (number of ASVs in subgroup)", 
           ylab = "Frequency", 
           col = "lightblue", 
           border = "black",
           xlim = c(0, 20), 
           ylim = c(0, 120), 
           breaks = 20,
           las = 1,                # Make axis labels horizontal
           cex.axis = 1.2,         # Increase axis label size
           cex.lab = 1.5,          # Increase axis title size
           cex.main = 1.6,         # Increase title size
           font.main = 2,          # Make title bold
           axes = FALSE)           # Suppress default axes

# Add custom x-axis with intervals by 2
axis(1, at = seq(0, 20, by = 2), cex.axis = 1.2)  

# Add custom y-axis (if needed)
axis(2, at = seq(0, 200, by = 20), cex.axis = 1.2)  

# Add a vertical line at the 95th percentile
abline(v = 8, col = "red", lwd = 2, lty = 2)

# Add annotation
text(x = subgroup_percentile_95 + 1, y = 150, 
     labels = paste("95th percentile:", round(subgroup_percentile_95)), 
     col = "red", cex = 1.2, adj = 0)

pdf("ClusterSizeDistr.pdf")
print(p2) 
dev.off()
















