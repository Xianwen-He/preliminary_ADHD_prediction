### This script aims to perform a T-test on SC Data.
### @Yu Chen, Dec 7, 2024


library(ggplot2)
library(reshape2)
library(dplyr)
library(Matrix)
library(glmnet)
library(caret)


cleaned_data <- readRDS("./Data/connectome.dat.rds") %>% 
  filter(!is.na(DSM_Adh_Raw)) %>% 
  filter(!sapply(sc_matrix, function(matrix) is.null(matrix) || all(is.na(matrix)))) %>% 
  arrange(DSM_Adh_Raw)

row_index <- 9  
col_index <- 24  

# Safely extract the specific entry from valid matrices
entries <- sapply(cleaned_data$sc_matrix, function(mat)  mat[row_index, col_index])

# Remove NA values
entries <- na.omit(entries)

# Calculate the mean and variance of the extracted entries
entry_mean <- mean(entries)
entry_variance <- var(entries)

# Print the results
cat("Mean:", entry_mean, "\n")
cat("Variance:", entry_variance, "\n")

# Plot the distribution
hist(entries, breaks = 30, main = paste("Distribution of Matrix Entry (", row_index, ", ", col_index, ")", sep = ""),
     xlab = "Value", col = "skyblue", border = "white")

png(filename = "./Data/img/distribution_histogram.png", width = 800, height = 600)
hist(entries, breaks = 30, main = paste("Distribution of Matrix Entry (", row_index, ", ", col_index, ")", sep = ""),
     xlab = "Value", col = "skyblue", border = "white")
dev.off()

# Create output directory
output_dir <- "./Data/img"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load and clean the data
cleaned_data <- readRDS("./Data/connectome.dat.rds") %>% 
  filter(!is.na(DSM_Adh_Raw)) %>% 
  filter(!sapply(sc_matrix, function(matrix) is.null(matrix) || all(is.na(matrix)))) %>% 
  arrange(DSM_Adh_Raw)

# Define combinations of group sizes
control_ranges <- list(c(0, 100), c(100, 300), c(100, 500))
adhd_ranges <- list(c(nrow(cleaned_data) - 100, nrow(cleaned_data)), 
                    c(nrow(cleaned_data) - 100, nrow(cleaned_data)),
                    c(nrow(cleaned_data) - 100, nrow(cleaned_data)))

# Iterate over combinations of group sizes
results <- list()

for (k in seq_along(control_ranges)) {
  control_range <- control_ranges[[k]]
  adhd_range <- adhd_ranges[[k]]
  
  # Select groups
  control_group_matrices <- cleaned_data[control_range[1]:control_range[2], ]$sc_matrix
  adhd_group_matrices <- cleaned_data[adhd_range[1]:adhd_range[2], ]$sc_matrix
  
  # Compute mean matrices and hypo-connectivity
  mean_control <- Reduce("+", control_group_matrices) / length(control_group_matrices)
  mean_adhd <- Reduce("+", adhd_group_matrices) / length(adhd_group_matrices)
  hypo_connectivity <- mean_control - mean_adhd
  
  # Compute p-values
  p_values <- matrix(NA, nrow = 68, ncol = 68)
  for (i in 1:68) {
    for (j in 1:68) {
      control_values <- sapply(control_group_matrices, function(mat) mat[i, j])
      adhd_values <- sapply(adhd_group_matrices, function(mat) mat[i, j])
      p_values[i, j] <- wilcox.test(control_values, adhd_values, var.equal = FALSE)$p.value
    }
  }
  
  # Adjust p-values
  adjusted_p_values <- p.adjust(as.vector(p_values), method = "fdr")
  adjusted_p_values_matrix <- matrix(adjusted_p_values, nrow = 68, ncol = 68)
  
  # Prepare data for visualization
  adjusted_p_values_df <- melt(adjusted_p_values_matrix)
  colnames(adjusted_p_values_df) <- c("Row", "Column", "Adjusted_P_Value")
  significant_regions <- adjusted_p_values_matrix < 0.05
  significant_p_values_df <- adjusted_p_values_df %>% 
    mutate(Significant = ifelse(Adjusted_P_Value < 0.05, "Yes", "No"))
  
  # Save heatmap
  heatmap_plot <- ggplot(adjusted_p_values_df, aes(x = Column, y = Row, fill = Adjusted_P_Value)) +
    geom_tile() +
    scale_fill_gradient(low = "red", high = "white", na.value = "grey50") +
    scale_y_reverse() +
    coord_fixed() +
    theme_minimal() +
    labs(
      title = paste("Adjusted P-values Heatmap (Control:", control_range[1], "-", control_range[2], 
                    ", ADHD:", adhd_range[1], "-", adhd_range[2], ")"),
      x = "Column",
      y = "Row",
      fill = "FDR Adjusted P-Value"
    )
  ggsave(
    filename = file.path(output_dir, paste0("heatmap_DSM_Adh_Raw_", k, ".png")),
    plot = heatmap_plot, width = 8, height = 6, dpi = 200
  )
  
  # Save significant regions plot
  significant_plot <- ggplot(significant_p_values_df, aes(x = Column, y = Row, fill = Significant)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes" = "red", "No" = "grey")) +
    scale_y_reverse() +
    coord_fixed() +
    theme_minimal() +
    labs(
      title = paste("Significant Regions After FDR Adjustment (Control:", control_range[1], "-", control_range[2], 
                    ", ADHD:", adhd_range[1], "-", adhd_range[2], ")"),
      x = "Column",
      y = "Row",
      fill = "Significant"
    )
  ggsave(
    filename = file.path(output_dir, paste0("significant_DSM_Adh_Raw_", k, ".png")),
    plot = significant_plot, width = 8, height = 6, dpi = 200
  )
  
  significant_indices <- which(adjusted_p_values_matrix < 0.05, arr.ind = TRUE)
  # Save boxplots for significant indices
  for (f in seq_len(nrow(significant_indices))) {
    significant_entry <- significant_indices[f, ]
    row_idx <- significant_entry[1]
    col_idx <- significant_entry[2]
    
    # Extract values for the significant entry
    control_values <- sapply(control_group_matrices, function(mat) mat[row_idx, col_idx])
    adhd_values <- sapply(adhd_group_matrices, function(mat) mat[row_idx, col_idx])
    
    # Combine the values into a data frame for plotting
    boxplot_data <- data.frame(
      Group = c(rep("Control", length(control_values)), rep("ADHD", length(adhd_values))),
      Value = c(control_values, adhd_values)
    )
    
    # Create boxplot
    boxplot <- ggplot(boxplot_data, aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot() +
      theme_minimal() +
      labs(
        title = paste("Boxplot of Entry (", row_idx, ", ", col_idx, ")", sep = ""),
        x = "Group",
        y = "Value"
      ) +
      scale_fill_manual(values = c("Control" = "lightblue", "ADHD" = "red"))
    
    ggsave(
      filename = file.path(output_dir, paste0("boxplot_DSM_Adh_Raw_", k, "_", row_idx, "_", col_idx, ".png")),
      plot = boxplot, width = 6, height = 4, dpi = 200
    )
  }
  
  # Store results
  results[[k]] <- list(
    control_range = control_range,
    adhd_range = adhd_range,
    significant_indices = significant_indices
  )
}


# for ReadEng_AgeAdj

output_dir <- "./Data/img"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load and clean the data
cleaned_data <- readRDS("./Data/connectome.dat.rds") %>% 
  filter(!is.na(ReadEng_AgeAdj)) %>% 
  filter(!sapply(sc_matrix, function(matrix) is.null(matrix) || all(is.na(matrix)))) %>% 
  arrange(ReadEng_AgeAdj)

# Define combinations of group sizes
control_ranges <- list(c(0, 100), c(100, 300), c(100, 500))
adhd_ranges <- list(c(nrow(cleaned_data) - 100, nrow(cleaned_data)), 
                    c(nrow(cleaned_data) - 100, nrow(cleaned_data)),
                    c(nrow(cleaned_data) - 100, nrow(cleaned_data)))

# Iterate over combinations of group sizes
results <- list()

for (k in seq_along(control_ranges)) {
  control_range <- control_ranges[[k]]
  adhd_range <- adhd_ranges[[k]]
  
  # Select groups
  control_group_matrices <- cleaned_data[control_range[1]:control_range[2], ]$sc_matrix
  adhd_group_matrices <- cleaned_data[adhd_range[1]:adhd_range[2], ]$sc_matrix
  
  # Compute mean matrices and hypo-connectivity
  mean_control <- Reduce("+", control_group_matrices) / length(control_group_matrices)
  mean_adhd <- Reduce("+", adhd_group_matrices) / length(adhd_group_matrices)
  hypo_connectivity <- mean_control - mean_adhd
  
  # Compute p-values
  p_values <- matrix(NA, nrow = 68, ncol = 68)
  for (i in 1:68) {
    for (j in 1:68) {
      control_values <- sapply(control_group_matrices, function(mat) mat[i, j])
      adhd_values <- sapply(adhd_group_matrices, function(mat) mat[i, j])
      p_values[i, j] <- wilcox.test(control_values, adhd_values, var.equal = FALSE)$p.value
    }
  }
  
  # Adjust p-values
  adjusted_p_values <- p.adjust(as.vector(p_values), method = "fdr")
  adjusted_p_values_matrix <- matrix(adjusted_p_values, nrow = 68, ncol = 68)
  
  # Prepare data for visualization
  adjusted_p_values_df <- melt(adjusted_p_values_matrix)
  colnames(adjusted_p_values_df) <- c("Row", "Column", "Adjusted_P_Value")
  significant_regions <- adjusted_p_values_matrix < 0.05
  significant_p_values_df <- adjusted_p_values_df %>% 
    mutate(Significant = ifelse(Adjusted_P_Value < 0.05, "Yes", "No"))
  
  # Save heatmap
  heatmap_plot <- ggplot(adjusted_p_values_df, aes(x = Column, y = Row, fill = Adjusted_P_Value)) +
    geom_tile() +
    scale_fill_gradient(low = "red", high = "white", na.value = "grey50") +
    scale_y_reverse() +
    coord_fixed() +
    theme_minimal() +
    labs(
      title = paste("Adjusted P-values Heatmap (Low score:", control_range[1], "-", control_range[2], 
                    ", High score:", adhd_range[1], "-", adhd_range[2], ")"),
      x = "Column",
      y = "Row",
      fill = "FDR Adjusted P-Value"
    )
  ggsave(
    filename = file.path(output_dir, paste0("heatmap_ReadEng_AgeAdj_", k, ".png")),
    plot = heatmap_plot, width = 8, height = 6, dpi = 200
  )
  
  # Save significant regions plot
  significant_plot <- ggplot(significant_p_values_df, aes(x = Column, y = Row, fill = Significant)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes" = "red", "No" = "grey")) +
    scale_y_reverse() +
    coord_fixed() +
    theme_minimal() +
    labs(
      title = paste("Significant Regions After FDR Adjustment (Low score:", control_range[1], "-", control_range[2], 
                    ", High score:", adhd_range[1], "-", adhd_range[2], ")"),
      x = "Column",
      y = "Row",
      fill = "Significant"
    )
  ggsave(
    filename = file.path(output_dir, paste0("significant_ReadEng_AgeAdj_", k, ".png")),
    plot = significant_plot, width = 8, height = 6, dpi = 200
  )
  
  # Save boxplots for significant indices
  significant_indices <- which(significant_regions, arr.ind = TRUE)
  for (f in seq_len(nrow(significant_indices))) {
    significant_entry <- significant_indices[f, ]
    row_idx <- significant_entry[1]
    col_idx <- significant_entry[2]
    
    # Extract values for the significant entry
    control_values <- sapply(control_group_matrices, function(mat) mat[row_idx, col_idx])
    adhd_values <- sapply(adhd_group_matrices, function(mat) mat[row_idx, col_idx])
    
    # Combine the values into a data frame for plotting
    boxplot_data <- data.frame(
      Group = c(rep("Control", length(control_values)), rep("ADHD", length(adhd_values))),
      Value = c(control_values, adhd_values)
    )
    
    # Create boxplot
    boxplot <- ggplot(boxplot_data, aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot() +
      theme_minimal() +
      labs(
        title = paste("Boxplot of Entry (", row_idx, ", ", col_idx, ")", sep = ""),
        x = "Group",
        y = "Value"
      ) +
      scale_fill_manual(values = c("Control" = "lightblue", "ADHD" = "red"))
    
    ggsave(
      filename = file.path(output_dir, paste0("boxplot_ReadEng_AgeAdj_", k, "_", row_idx, "_", col_idx, ".png")),
      plot = boxplot, width = 6, height = 4, dpi = 200
    )
  }
  
  # Store results
  results[[k]] <- list(
    control_range = control_range,
    adhd_range = adhd_range,
    significant_indices = significant_indices
  )
}




# for ProcSpeed_AgeAdj
# Create output directory
output_dir <- "./Data/img"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load and clean the data
cleaned_data <- readRDS("./Data/connectome.dat.rds") %>% 
  filter(!is.na(ProcSpeed_AgeAdj)) %>% 
  filter(!sapply(sc_matrix, function(matrix) is.null(matrix) || all(is.na(matrix)))) %>% 
  arrange(ProcSpeed_AgeAdj)

# Define combinations of group sizes
control_ranges <- list(c(0, 100), c(100, 300), c(100, 500))
adhd_ranges <- list(c(nrow(cleaned_data) - 100, nrow(cleaned_data)), 
                    c(nrow(cleaned_data) - 100, nrow(cleaned_data)),
                    c(nrow(cleaned_data) - 100, nrow(cleaned_data)))

# Iterate over combinations of group sizes
results <- list()

for (k in seq_along(control_ranges)) {
  control_range <- control_ranges[[k]]
  adhd_range <- adhd_ranges[[k]]
  
  # Select groups
  control_group_matrices <- cleaned_data[control_range[1]:control_range[2], ]$sc_matrix
  adhd_group_matrices <- cleaned_data[adhd_range[1]:adhd_range[2], ]$sc_matrix
  
  # Compute mean matrices and hypo-connectivity
  mean_control <- Reduce("+", control_group_matrices) / length(control_group_matrices)
  mean_adhd <- Reduce("+", adhd_group_matrices) / length(adhd_group_matrices)
  hypo_connectivity <- mean_control - mean_adhd
  
  # Compute p-values
  p_values <- matrix(NA, nrow = 68, ncol = 68)
  for (i in 1:68) {
    for (j in 1:68) {
      control_values <- sapply(control_group_matrices, function(mat) mat[i, j])
      adhd_values <- sapply(adhd_group_matrices, function(mat) mat[i, j])
      p_values[i, j] <- wilcox.test(control_values, adhd_values, var.equal = FALSE)$p.value
    }
  }
  
  # Adjust p-values
  adjusted_p_values <- p.adjust(as.vector(p_values), method = "fdr")
  adjusted_p_values_matrix <- matrix(adjusted_p_values, nrow = 68, ncol = 68)
  
  # Prepare data for visualization
  adjusted_p_values_df <- melt(adjusted_p_values_matrix)
  colnames(adjusted_p_values_df) <- c("Row", "Column", "Adjusted_P_Value")
  significant_regions <- adjusted_p_values_matrix < 0.05
  significant_p_values_df <- adjusted_p_values_df %>% 
    mutate(Significant = ifelse(Adjusted_P_Value < 0.05, "Yes", "No"))
  
  # Save heatmap
  heatmap_plot <- ggplot(adjusted_p_values_df, aes(x = Column, y = Row, fill = Adjusted_P_Value)) +
    geom_tile() +
    scale_fill_gradient(low = "red", high = "white", na.value = "grey50") +
    scale_y_reverse() +
    coord_fixed() +
    theme_minimal() +
    labs(
      title = paste("Adjusted P-values Heatmap (Low score:", control_range[1], "-", control_range[2], 
                    ", High score:", adhd_range[1], "-", adhd_range[2], ")"),
      x = "Column",
      y = "Row",
      fill = "FDR Adjusted P-Value"
    )
  ggsave(
    filename = file.path(output_dir, paste0("heatmap_ProcSpeed_AgeAdj_", k, ".png")),
    plot = heatmap_plot, width = 8, height = 6, dpi = 200
  )
  
  # Save significant regions plot
  significant_plot <- ggplot(significant_p_values_df, aes(x = Column, y = Row, fill = Significant)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes" = "red", "No" = "grey")) +
    scale_y_reverse() +
    coord_fixed() +
    theme_minimal() +
    labs(
      title = paste("Significant Regions After FDR Adjustment (Low score:", control_range[1], "-", control_range[2], 
                    ", High score:", adhd_range[1], "-", adhd_range[2], ")"),
      x = "Column",
      y = "Row",
      fill = "Significant"
    )
  ggsave(
    filename = file.path(output_dir, paste0("significant_ProcSpeed_AgeAdj_", k, ".png")),
    plot = significant_plot, width = 8, height = 6, dpi = 200
  )
  
  # Save boxplots for significant indices
  significant_indices <- which(significant_regions, arr.ind = TRUE)
  for (f in seq_len(nrow(significant_indices))) {
    significant_entry <- significant_indices[f, ]
    row_idx <- significant_entry[1]
    col_idx <- significant_entry[2]
    
    # Extract values for the significant entry
    control_values <- sapply(control_group_matrices, function(mat) mat[row_idx, col_idx])
    adhd_values <- sapply(adhd_group_matrices, function(mat) mat[row_idx, col_idx])
    
    # Combine the values into a data frame for plotting
    boxplot_data <- data.frame(
      Group = c(rep("Low score", length(control_values)), rep("High score", length(adhd_values))),
      Value = c(control_values, adhd_values)
    )
    
    # Create boxplot
    boxplot <- ggplot(boxplot_data, aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot() +
      theme_minimal() +
      labs(
        title = paste("Boxplot of Entry (", row_idx, ", ", col_idx, ")", sep = ""),
        x = "Group",
        y = "Value"
      ) +
      scale_fill_manual(values = c("Low score" = "lightblue", "High score" = "red"))
    
    ggsave(
      filename = file.path(output_dir, paste0("boxplot_ProcSpeed_AgeAdj_", k, "_", row_idx, "_", col_idx, ".png")),
      plot = boxplot, width = 6, height = 4, dpi = 200
    )
  }
  
  # Store results
  results[[k]] <- list(
    control_range = control_range,
    adhd_range = adhd_range,
    significant_indices = significant_indices
  )
}


# height 
output_dir <- "./Data/img"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load and clean the data
cleaned_data <- readRDS("./Data/connectome.dat.rds") %>% 
  filter(!is.na(DSM_Adh_Raw)) %>% 
  filter(!sapply(sc_matrix, function(matrix) is.null(matrix) || all(is.na(matrix)))) %>% 
  arrange(Height)

# Define combinations of group sizes
control_ranges <- list(c(0, 100), c(100, 300), c(100, 500))
adhd_ranges <- list(c(nrow(cleaned_data) - 100, nrow(cleaned_data)), 
                    c(nrow(cleaned_data) - 100, nrow(cleaned_data)),
                    c(nrow(cleaned_data) - 100, nrow(cleaned_data)))

# Iterate over combinations of group sizes
results <- list()

for (k in seq_along(control_ranges)) {
  control_range <- control_ranges[[k]]
  adhd_range <- adhd_ranges[[k]]
  
  # Select groups
  control_group_matrices <- cleaned_data[control_range[1]:control_range[2], ]$sc_matrix
  adhd_group_matrices <- cleaned_data[adhd_range[1]:adhd_range[2], ]$sc_matrix
  
  # Compute mean matrices and hypo-connectivity
  mean_control <- Reduce("+", control_group_matrices) / length(control_group_matrices)
  mean_adhd <- Reduce("+", adhd_group_matrices) / length(adhd_group_matrices)
  hypo_connectivity <- mean_control - mean_adhd
  
  # Compute p-values
  p_values <- matrix(NA, nrow = 68, ncol = 68)
  for (i in 1:68) {
    for (j in 1:68) {
      control_values <- sapply(control_group_matrices, function(mat) mat[i, j])
      adhd_values <- sapply(adhd_group_matrices, function(mat) mat[i, j])
      p_values[i, j] <- wilcox.test(control_values, adhd_values, var.equal = FALSE)$p.value
    }
  }
  
  # Adjust p-values
  adjusted_p_values <- p.adjust(as.vector(p_values), method = "fdr")
  adjusted_p_values_matrix <- matrix(adjusted_p_values, nrow = 68, ncol = 68)
  
  # Prepare data for visualization
  adjusted_p_values_df <- melt(adjusted_p_values_matrix)
  colnames(adjusted_p_values_df) <- c("Row", "Column", "Adjusted_P_Value")
  significant_regions <- adjusted_p_values_matrix < 0.05
  significant_p_values_df <- adjusted_p_values_df %>% 
    mutate(Significant = ifelse(Adjusted_P_Value < 0.05, "Yes", "No"))
  
  # Save heatmap
  heatmap_plot <- ggplot(adjusted_p_values_df, aes(x = Column, y = Row, fill = Adjusted_P_Value)) +
    geom_tile() +
    scale_fill_gradient(low = "red", high = "white", na.value = "grey50") +
    scale_y_reverse() +
    coord_fixed() +
    theme_minimal() +
    labs(
      title = paste("Adjusted P-values Heatmap (Low score:", control_range[1], "-", control_range[2], 
                    ", High score:", adhd_range[1], "-", adhd_range[2], ")"),
      x = "Column",
      y = "Row",
      fill = "FDR Adjusted P-Value"
    )
  ggsave(
    filename = file.path(output_dir, paste0("heatmap_Height_", k, ".png")),
    plot = heatmap_plot, width = 8, height = 6, dpi = 200
  )
  
  # Save significant regions plot
  significant_plot <- ggplot(significant_p_values_df, aes(x = Column, y = Row, fill = Significant)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes" = "red", "No" = "grey")) +
    scale_y_reverse() +
    coord_fixed() +
    theme_minimal() +
    labs(
      title = paste("Significant Regions After FDR Adjustment (Low score:", control_range[1], "-", control_range[2], 
                    ", High score:", adhd_range[1], "-", adhd_range[2], ")"),
      x = "Column",
      y = "Row",
      fill = "Significant"
    )
  ggsave(
    filename = file.path(output_dir, paste0("significant_Height_", k, ".png")),
    plot = significant_plot, width = 8, height = 6, dpi = 200
  )
  
  # Save boxplots for significant indices
  significant_indices <- which(significant_regions, arr.ind = TRUE)
  for (f in seq_len(nrow(significant_indices))) {
    significant_entry <- significant_indices[f, ]
    row_idx <- significant_entry[1]
    col_idx <- significant_entry[2]
    
    # Extract values for the significant entry
    control_values <- sapply(control_group_matrices, function(mat) mat[row_idx, col_idx])
    adhd_values <- sapply(adhd_group_matrices, function(mat) mat[row_idx, col_idx])
    
    # Combine the values into a data frame for plotting
    boxplot_data <- data.frame(
      Group = c(rep("Low score", length(control_values)), rep("High score", length(adhd_values))),
      Value = c(control_values, adhd_values)
    )
    
    # Create boxplot
    boxplot <- ggplot(boxplot_data, aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot() +
      theme_minimal() +
      labs(
        title = paste("Boxplot of Entry (", row_idx, ", ", col_idx, ")", sep = ""),
        x = "Group",
        y = "Value"
      ) +
      scale_fill_manual(values = c("Low score" = "lightblue", "High score" = "red"))
    
    ggsave(
      filename = file.path(output_dir, paste0("boxplot_Height_", k, "_", row_idx, "_", col_idx, ".png")),
      plot = boxplot, width = 6, height = 4, dpi = 200
    )
  }
  
  # Store results
  results[[k]] <- list(
    control_range = control_range,
    adhd_range = adhd_range,
    significant_indices = significant_indices
  )
}



#SSAGA_Alc_12_Frq_Drk

output_dir <- "./Data/img"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load and clean the data
cleaned_data <- readRDS("./Data/connectome.dat.rds") %>% 
  filter(!is.na(DSM_Adh_Raw)) %>% 
  filter(!sapply(sc_matrix, function(matrix) is.null(matrix) || all(is.na(matrix)))) %>% 
  arrange(SSAGA_Alc_12_Frq_Drk)

# Define combinations of group sizes
control_ranges <- list(c(0, 100), c(100, 300), c(100, 500))
adhd_ranges <- list(c(nrow(cleaned_data) - 100, nrow(cleaned_data)), 
                    c(nrow(cleaned_data) - 100, nrow(cleaned_data)),
                    c(nrow(cleaned_data) - 100, nrow(cleaned_data)))

# Iterate over combinations of group sizes
results <- list()

for (k in seq_along(control_ranges)) {
  control_range <- control_ranges[[k]]
  adhd_range <- adhd_ranges[[k]]
  
  # Select groups
  control_group_matrices <- cleaned_data[control_range[1]:control_range[2], ]$sc_matrix
  adhd_group_matrices <- cleaned_data[adhd_range[1]:adhd_range[2], ]$sc_matrix
  
  # Compute mean matrices and hypo-connectivity
  mean_control <- Reduce("+", control_group_matrices) / length(control_group_matrices)
  mean_adhd <- Reduce("+", adhd_group_matrices) / length(adhd_group_matrices)
  hypo_connectivity <- mean_control - mean_adhd
  
  # Compute p-values
  p_values <- matrix(NA, nrow = 68, ncol = 68)
  for (i in 1:68) {
    for (j in 1:68) {
      control_values <- sapply(control_group_matrices, function(mat) mat[i, j])
      adhd_values <- sapply(adhd_group_matrices, function(mat) mat[i, j])
      p_values[i, j] <- wilcox.test(control_values, adhd_values, var.equal = FALSE)$p.value
    }
  }
  
  # Adjust p-values
  adjusted_p_values <- p.adjust(as.vector(p_values), method = "fdr")
  adjusted_p_values_matrix <- matrix(adjusted_p_values, nrow = 68, ncol = 68)
  
  # Prepare data for visualization
  adjusted_p_values_df <- melt(adjusted_p_values_matrix)
  colnames(adjusted_p_values_df) <- c("Row", "Column", "Adjusted_P_Value")
  significant_regions <- adjusted_p_values_matrix < 0.05
  significant_p_values_df <- adjusted_p_values_df %>% 
    mutate(Significant = ifelse(Adjusted_P_Value < 0.05, "Yes", "No"))
  
  # Save heatmap
  heatmap_plot <- ggplot(adjusted_p_values_df, aes(x = Column, y = Row, fill = Adjusted_P_Value)) +
    geom_tile() +
    scale_fill_gradient(low = "red", high = "white", na.value = "grey50") +
    scale_y_reverse() +
    coord_fixed() +
    theme_minimal() +
    labs(
      title = paste("Adjusted P-values Heatmap (Low score:", control_range[1], "-", control_range[2], 
                    ", High score:", adhd_range[1], "-", adhd_range[2], ")"),
      x = "Column",
      y = "Row",
      fill = "FDR Adjusted P-Value"
    )
  ggsave(
    filename = file.path(output_dir, paste0("heatmap_SSAGA_Alc_12_Frq_Drk_", k, ".png")),
    plot = heatmap_plot, width = 8, height = 6, dpi = 200
  )
  
  # Save significant regions plot
  significant_plot <- ggplot(significant_p_values_df, aes(x = Column, y = Row, fill = Significant)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes" = "red", "No" = "grey")) +
    scale_y_reverse() +
    coord_fixed() +
    theme_minimal() +
    labs(
      title = paste("Significant Regions After FDR Adjustment (Low score:", control_range[1], "-", control_range[2], 
                    ", High score:", adhd_range[1], "-", adhd_range[2], ")"),
      x = "Column",
      y = "Row",
      fill = "Significant"
    )
  ggsave(
    filename = file.path(output_dir, paste0("significant_SSAGA_Alc_12_Frq_Drk_", k, ".png")),
    plot = significant_plot, width = 8, height = 6, dpi = 200
  )
  
  # Save boxplots for significant indices
  significant_indices <- which(significant_regions, arr.ind = TRUE)
  for (f in seq_len(nrow(significant_indices))) {
    significant_entry <- significant_indices[f, ]
    row_idx <- significant_entry[1]
    col_idx <- significant_entry[2]
    
    # Extract values for the significant entry
    control_values <- sapply(control_group_matrices, function(mat) mat[row_idx, col_idx])
    adhd_values <- sapply(adhd_group_matrices, function(mat) mat[row_idx, col_idx])
    
    # Combine the values into a data frame for plotting
    boxplot_data <- data.frame(
      Group = c(rep("Low score", length(control_values)), rep("High score", length(adhd_values))),
      Value = c(control_values, adhd_values)
    )
    
    # Create boxplot
    boxplot <- ggplot(boxplot_data, aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot() +
      theme_minimal() +
      labs(
        title = paste("Boxplot of Entry (", row_idx, ", ", col_idx, ")", sep = ""),
        x = "Group",
        y = "Value"
      ) +
      scale_fill_manual(values = c("Low score" = "lightblue", "High score" = "red"))
    
    ggsave(
      filename = file.path(output_dir, paste0("boxplot_SSAGA_Alc_12_Frq_Drk_", k, "_", row_idx, "_", col_idx, ".png")),
      plot = boxplot, width = 6, height = 4, dpi = 200
    )
  }
  
  # Store results
  results[[k]] <- list(
    control_range = control_range,
    adhd_range = adhd_range,
    significant_indices = significant_indices
  )
}


cleaned_data <- readRDS("./Data/connectome.dat.rds") %>% 
  filter(!is.na(DSM_Adh_Raw)) %>% 
  filter(!sapply(sc_matrix, function(matrix) is.null(matrix) || all(is.na(matrix))))

# Define variable names
variables <- c("ReadEng_AgeAdj", "IWRD_TOT", "ProcSpeed_AgeAdj", "DSM_Adh_Raw","VSPLOT_TC","Height","BMI","SSAGA_Income","SSAGA_Educ","SSAGA_Alc_12_Frq_Drk","SSAGA_Times_Used_Stimulants","SSAGA_ChildhoodConduct")

# Initialize results
variance_explained <- data.frame(Variable = character(), Train = numeric(), Test = numeric())

# Loop through each variable
for (var in variables) {
  # Filter out missing or constant values for the current variable
  valid_data <- cleaned_data %>% filter(!is.na(.data[[var]])) 
  #%>% filter(sd(.data[[var]]) > 0)
  response <- valid_data[[var]]
  scmatrix <- do.call(rbind, lapply(valid_data$sc_matrix, as.vector))
  
  # Create training and testing datasets (random split for randomness control)
  set.seed(123) # Set seed for reproducibility
  train_index <- sample(seq_len(nrow(scmatrix)), size = 0.8 * nrow(scmatrix))
  train_scmatrix <- scmatrix[train_index, ]
  train_response <- response[train_index]
  test_scmatrix <- scmatrix[-train_index, ]
  test_response <- response[-train_index]
  
  # Perform LASSO regression on the training set
  lasso_model <- cv.glmnet(as.matrix(train_scmatrix), train_response, alpha = 1)
  
  # Model Summary
  best_lambda <- lasso_model$lambda.min
  predicted_train <- predict(lasso_model, as.matrix(train_scmatrix), s = best_lambda)
  predicted_test <- predict(lasso_model, as.matrix(test_scmatrix), s = best_lambda)
  
  # Variance explained
  variance_explained_train <- 1 - sum((train_response - predicted_train)^2) / sum((train_response - mean(train_response))^2)
  variance_explained_test <- 1 - sum((test_response - predicted_test)^2) / sum((test_response - mean(test_response))^2)
  
  # Store results
  variance_explained <- rbind(variance_explained, 
                              data.frame(Variable = var, Test = variance_explained_test))
}

#Train = variance_explained_train
# Plot variance explained as a barplot
variance_explained_melted <- reshape2::melt(variance_explained, id.vars = "Variable", variable.name = "Dataset", value.name = "VarianceExplained")

varexplained=ggplot(variance_explained_melted, aes(x = reorder(Variable, -VarianceExplained), y = VarianceExplained, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Variance Explained by LASSO Model for Each Variable",
    x = "Variable",
    y = "Variance Explained",
    fill = "Dataset"
  )

ggsave(filename='./Data/img/varexplained.png', plot=varexplained,
       width=8, height = 6, dpi=200)


### Analyze the point [9, 24]

# Extract the [9,24] values from the control group matrices
control_values <- sapply(control_group_matrices, function(mat) mat[9, 24])
mean_control <- mean(control_values, na.rm = TRUE)

# Extract the [9,24] values from the ADHD group matrices
adhd_values <- sapply(adhd_group_matrices, function(mat) mat[9, 24])
mean_adhd <- mean(adhd_values, na.rm = TRUE)

# Print the means
cat("Mean of [9,24] in Control Group:", mean_control, "\n")
cat("Mean of [9,24] in ADHD Group:", mean_adhd, "\n")

var(control_values, na.rm = TRUE)
var(adhd_values, na.rm = TRUE)
t.test(control_values, adhd_values, var.equal = FALSE)

# plot [9, 24]
values_df <- data.frame(
  Group = c(rep("Control", length(control_values)), rep("ADHD", length(adhd_values))),
  Value = c(control_values, adhd_values)
)
sc.sign.plot <- ggplot(values_df, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Distribution of [9,24] Values",
    x = "Group",
    y = "[9,24] Value"
  ) +
  theme_minimal()
ggsave(filename='./Data/img/sc_sign_plot.png', plot=sc.sign.plot,
       width=8, height = 6, dpi=200)
