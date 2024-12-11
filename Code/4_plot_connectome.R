### This script aims to generate plots on connectome data
### Coded by @Yu Chen, Nov, 2024

library(ggplot2)
library(tidyverse)
library(plotly)

connectome.dat <- readRDS("./Data/connectome.dat.rds")


### plot sample fc and sc 
library(dplyr)
library(ggplot2)

# Load and clean the data
cleaned_data <-  readRDS("./Data/connectome.dat.rds") %>% 
  filter(!is.na(DSM_Adh_Raw)) %>% 
  filter(!sapply(sc_matrix, function(matrix) is.null(matrix) || all(is.na(matrix))))

# Extract the first SC matrix
sc_matrix1 <- matrix(unlist(cleaned_data$sc_matrix[[2]]), nrow = 68, ncol = 68)

# Convert the matrix to a long-format data frame
sc_matrix_df <- data.frame(
  Row = rep(1:nrow(sc_matrix1), each = ncol(sc_matrix1)),  # Row indices
  Column = rep(1:ncol(sc_matrix1), times = nrow(sc_matrix1)),  # Column indices
  Value = as.vector(sc_matrix1)  # Flatten the matrix into a vector
)

# Plot and save the SC matrix heatmap
sc_plot <- ggplot(sc_matrix_df, aes(x = Row, y = Column, fill = log10(Value))) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "red", midpoint = 0, na.value = "white") +
  scale_y_reverse() +  # Reverse the y-axis to maintain a matrix-like orientation
  coord_fixed() +  # Ensure square tiles
  theme_minimal() +
  labs(
    title = "Heatmap of sc_matrix[1]",
    x = "Row",
    y = "Column",
    fill = "Log10(Value)"
  ) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5)
  )

# Save the SC matrix heatmap
ggsave("./Data/img/sample_sc.png", plot = sc_plot, width = 8, height = 6)


# Function to extract and visualize upper diagonal of FC matrices
visualize_fc_matrix <- function(fc_matrices, index = 1, save_path = NULL) {
  # Extract the innermost matrix
  fc_matrix <- fc_matrices[[index]]
  while (is.list(fc_matrix)) {
    fc_matrix <- fc_matrix[[1]]
  }
  fc_matrix <- as.matrix(fc_matrix)
  if (dim(fc_matrix)[1] != 68 || dim(fc_matrix)[2] != 68) {
    stop("The extracted matrix is not 68x68.")
  }
  
  # Extract the upper diagonal
  upper_tri_indices <- upper.tri(fc_matrix, diag = TRUE)
  
  # Create a data frame for visualization
  fc_matrix_df <- data.frame(
    Row = row(fc_matrix)[upper_tri_indices],
    Column = col(fc_matrix)[upper_tri_indices],
    Value = fc_matrix[upper_tri_indices]
  )
  
  # Plot the heatmap for the upper diagonal
  fc_plot <- ggplot(fc_matrix_df, aes(x = Column, y = Row, fill = Value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey") +
    scale_y_reverse() +  # Reverse the y-axis for matrix-like orientation
    coord_fixed() +  # Ensure square tiles
    theme_minimal() +
    labs(
      title = paste("Heatmap of Upper Diagonal of FC Matrix (Index:", index, ")"),
      x = "Column",
      y = "Row",
      fill = "Value"
    ) +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Save the plot if save_path is provided
  if (!is.null(save_path)) {
    ggsave(save_path, plot = fc_plot, width = 8, height = 6)
  }
  
  return(fc_plot)
}

# Plot and save the FC matrix heatmap
fc_plot <- visualize_fc_matrix(cleaned_data$fc_matrix, index = 1, save_path = "./Data/img/sample_fc.png")



### Plot PCA data

# function to process the data and generate plot
pca_projection_plot <- function(data, variable_name, split_value = NULL, pc_cols = c(1, 2), use_fcpcall = TRUE){
  # Step 1: Select the data set based on the `use_fcpcall` flag
  chosen_dataset <- if (use_fcpcall) "fcpcall" else "scpcall"
  
  # Step 2: Remove rows with NA in the second column of the chosen dataset
  data <- data[!is.na(data[[chosen_dataset]][, 2]), ]
  
  # Step 3: Check the type of the variable and split into groups
  variable <- data[[variable_name]]
  if (is.factor(variable) || is.character(variable) || length(unique(variable)) == 2) {
    # Binary variable case
    levels <- unique(variable)
    if (length(levels) != 2) {
      stop("Binary variables should have exactly two unique levels.")
    }
    group1 <- data %>% filter(.data[[variable_name]] == levels[1])
    group2 <- data %>% filter(.data[[variable_name]] == levels[2])
    group1_label <- paste(variable_name, "=", levels[1])
    group2_label <- paste(variable_name, "=", levels[2])
  } else {
    # Continuous variable case
    if (is.null(split_value)) {
      stop("For continuous variables, a split_value must be provided.")
    }
    group1 <- data %>% filter(.data[[variable_name]] <= split_value)
    group2 <- data %>% filter(.data[[variable_name]] > split_value)
    group1_label <- paste(variable_name, "<=", split_value)
    group2_label <- paste(variable_name, ">", split_value)
  }
  
  # Step 4: Extract PCA data
  extract_pca_data <- function(group, label, pc_cols, dataset) {
    pca_scores <- group[[dataset]] %>% bind_rows()  # Combine all PCA data.frames
    pca_scores <- pca_scores[, pc_cols]  # Subset columns using indices
    colnames(pca_scores) <- c("PCX", "PCY")  # Rename columns
    pca_scores <- pca_scores %>% mutate(Group = label)  # Add group label
    return(pca_scores)
  }
  
  group1_pca <- extract_pca_data(group1, group1_label, pc_cols, chosen_dataset)
  group2_pca <- extract_pca_data(group2, group2_label, pc_cols, chosen_dataset)
  
  # Step 5: Combine both groups into a single data frame
  pca_combined <- bind_rows(group1_pca, group2_pca)
  
  # Step 6: Plot the PCA projection
  ggplot(pca_combined, aes(x = PCX, y = PCY, color = Group)) +
    geom_point(alpha = 0.7, size = 3) +
    labs(
      title = paste("2D Projection of", chosen_dataset, "Data by", variable_name, "Groups"),
      x = paste0("Principal Component ", pc_cols[1]-1),
      y = paste0("Principal Component ", pc_cols[2]-1),
      color = "Group"
    ) +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(hjust = 0.5)
    )
}

# plot fc-pca
for (i in 3:6){
  img = pca_projection_plot(data = connectome.dat, variable_name = "DSM_Adh_Raw", split_value = 11, pc_cols = c(2, i), use_fcpcall = T)
  ggsave(filename=paste0('./Data/img/fcpca_plot_1vs', as.character(i-1), '.png'), plot=img, width=8, height = 6, dpi=200)
}

# plot for sc-pca
for (i in 3:6){
  img=pca_projection_plot(data = connectome.dat, variable_name = "DSM_Adh_Raw", split_value = 11, pc_cols = c(2, i), use_fcpcall = F)
  ggsave(filename=paste0('./Data/img/scpca_plot_2vs', as.character(i-1), '.png'), plot=img, width=8, height = 6, dpi=200)
}

# plot for sc-pca for height 
for (i in 3:6){
  img=pca_projection_plot(data = connectome.dat, variable_name = "Height", split_value = 70, pc_cols = c(2, i), use_fcpcall = F)
  ggsave(filename=paste0('./Data/img/scpca__height_plot_2_vs', as.character(i-1), '.png'), plot=img, width=8, height = 6, dpi=200)
}

# plot for sc-pca for SAGA_Alc_12_Frq_Drk
for (i in 3:6){
  img=pca_projection_plot(data = connectome.dat, variable_name = "SSAGA_Alc_12_Frq_Drk", split_value = 3, pc_cols = c(2, i), use_fcpcall = F)
  ggsave(filename=paste0('./Data/img/scpca__SAGA_Alc_12_Frq_Drk_plot_2_vs', as.character(i-1), '.png'), plot=img, width=8, height = 6, dpi=200)
}

for (i in 3:6){
  img=pca_projection_plot(data = connectome.dat, variable_name = "VSPLOT_TC", split_value = 19, pc_cols = c(2, i), use_fcpcall = F)
  ggsave(filename=paste0('./Data/img/scpca__VSPLOT_TC_plot_2_vs', as.character(i-1), '.png'), plot=img, width=8, height = 6, dpi=200)
}


### 3D Visualization

# function for plotting
pca_projection_plot_3d <- function(data, variable_name, split_value = NULL, pc_cols = c(1, 2, 3), use_fcpcall = TRUE) {
  # Step 1: Select the dataset based on the `use_fcpcall` flag
  chosen_dataset <- if (use_fcpcall) "fcpcall" else "scpcall"
  
  # Step 2: Remove rows with NA in the second column of the chosen dataset
  data <- data[!is.na(data[[chosen_dataset]][, 2]), ]
  
  # Step 3: Check the type of the variable and split into groups
  variable <- data[[variable_name]]
  if (is.factor(variable) || is.character(variable) || length(unique(variable)) == 2) {
    # Binary variable case
    levels <- unique(variable)
    if (length(levels) != 2) {
      stop("Binary variables should have exactly two unique levels.")
    }
    group1 <- data %>% filter(.data[[variable_name]] == levels[1])
    group2 <- data %>% filter(.data[[variable_name]] == levels[2])
    group1_label <- paste(variable_name, "=", levels[1])
    group2_label <- paste(variable_name, "=", levels[2])
  } else {
    # Continuous variable case
    if (is.null(split_value)) {
      stop("For continuous variables, a split_value must be provided.")
    }
    group1 <- data %>% filter(.data[[variable_name]] <= split_value)
    group2 <- data %>% filter(.data[[variable_name]] > split_value)
    group1_label <- paste(variable_name, "<=", split_value)
    group2_label <- paste(variable_name, ">", split_value)
  }
  
  # Step 4: Extract PCA data
  extract_pca_data <- function(group, label, pc_cols, dataset) {
    pca_scores <- group[[dataset]] %>% bind_rows()  # Combine all PCA data.frames
    pca_scores <- pca_scores %>%
      select(PCX = pc_cols[1], PCY = pc_cols[2], PCZ = pc_cols[3]) %>%  # Select specified PCs for 3D projection
      mutate(Group = label)  # Add group label
    return(pca_scores)
  }
  
  group1_pca <- extract_pca_data(group1, group1_label, pc_cols, chosen_dataset)
  group2_pca <- extract_pca_data(group2, group2_label, pc_cols, chosen_dataset)
  
  # Step 5: Combine both groups into a single data frame
  pca_combined <- bind_rows(group1_pca, group2_pca)
  
  # Step 6: Perform PCA on all combined PC scores (if not using precomputed PCs)
  
  # Step 7: Plot the PCA projection in 3D
  plot_ly(data = pca_combined, 
          x = ~PCX, y = ~PCY, z = ~PCZ, 
          color = ~Group, 
          colors = c("blue", "yellow"),
          type = 'scatter3d', mode = 'markers') %>%
    layout(title = paste("3D PCA Projection of", chosen_dataset, "Data by", variable_name, "Groups"),
           scene = list(
             xaxis = list(title = paste0("Principal Component ", pc_cols[1])),
             yaxis = list(title = paste0("Principal Component ", pc_cols[2])),
             zaxis = list(title = paste0("Principal Component ", pc_cols[3]))
           ))
}

# plot SC PCA in 3D
# This works in a interactive way
pca_projection_plot_3d(data = connectome.dat, variable_name = "DSM_Adh_Raw",
                       split_value = 11, pc_cols = c(1, 2, 3), use_fcpcall = F)

