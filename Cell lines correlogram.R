# Install and load packages
install.packages(c("SummarizedExperiment", "tidyr", "dplyr", "tibble", "ggplot2", 
                   "pheatmap", "ComplexHeatmap", "ggpubr", "reshape2", 
                   "umap", "class", "cluster", "viridis", "boot", "EnhancedVolcano", 
                   "corrplot", "ggpubr"))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

library(SummarizedExperiment)
library(MultiAssayExperiment)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(ggpubr)
library(reshape2)
library(PCAtools)
library(umap)
library(class)
library(cluster)
library(viridis)
library(boot)
library(EnhancedVolcano)
library(corrplot)
library(ggpubr)
library(stringr)
library(readxl)


# Read in data
normalized_ab <- read.csv("C:/Users/lcorvele/R/Correlogram/Pharmacoepigenetics/Cell_lines_sPTM_filtered.csv")

## CORRELATION PLOT

# Remove rows where the Peptide column contains "Unmod"
#normalized_ab <- normalized_ab[!grepl("Unmod", normalized_ab$PTM), ]

# Remove rows where the Peptide column contains "H1" or "H2"
#normalized_ab <- normalized_ab[!grepl("H1|H2", normalized_ab$PTM), ]

# Transpose rows and columns
normalized_ab <- t(normalized_ab)

# Set the first row as column names
colnames(normalized_ab) <- normalized_ab[1, ]

# Remove the first row (now column names) from the transposed matrix
normalized_ab <- normalized_ab[-1, ]

# Convert non-numeric columns to numeric
normalized_ab <- as.data.frame(normalized_ab)

# Convert all values in the data frame to numeric
normalized_ab_numeric <- apply(normalized_ab, 2, as.numeric)

# Change all NA values to zero
normalized_ab_numeric[is.na(normalized_ab_numeric)] <- 0

# Create a sample correlation matrix 
correlation_matrix <- cor(normalized_ab_numeric)  

# Check for any non-finite values in the correlation matrix
if (!all(is.finite(correlation_matrix))) {
  # If there are any non-finite values, replace them with 0
  correlation_matrix[!is.finite(correlation_matrix)] <- 0
}

# Perform hierarchical clustering on rows of the correlation matrix
hc_rows <- hclust(dist(correlation_matrix))

# Get the order from hierarchical clustering
order_rows <- hc_rows$order

# Reorder rows and columns based on hierarchical clustering order
correlation_matrix_reordered <- correlation_matrix[order_rows, order_rows]

# Write CSV with row names
write.csv(correlation_matrix_reordered, file = "correlation_matrix_cell_lines_filtered.csv")

# Get the number of columns in normalized_ab
n <- ncol(normalized_ab_numeric)

# Initialize a matrix to store p-values
p_mat <- matrix(NA, nrow = n, ncol = n)

# Set row and column names of p_mat to match the column names of normalized_ab
colnames(p_mat) <- rownames(p_mat) <- colnames(normalized_ab)

# Calculate p-values
for(i in 1:n){
  for(j in 1:n){
    if(i != j){
      test_result <- cor.test(normalized_ab_numeric[, i], normalized_ab_numeric[, j])
      p_mat[i, j] <- test_result$p.value
    } else {
      p_mat[i, j] <- 0  # Diagonal elements are not applicable
    }
  }
}

# Write CSV with p values
write.csv(p_mat, file = "pval_correlation_matrix_cell_lines.csv")

###########
# FOR THE SECOND CORRELOGRAM
###########

# Read the CSV file with row names
correlation_matrix_csv <- read.csv("correlation_matrix_PDX_filtered.csv", row.names = 1)

# Match the row names with the ordered names
ordered_names <- rownames(correlation_matrix_csv)

# Reorder the rows of the correlation matrix based on the matched row names
correlation_matrix_reordered <- correlation_matrix[ordered_names, ordered_names]

###############
#CREATE THE CORRELATION PLOT
###############
# Generate a custom color palette
num_colors <- 200
custom_colors <- colorRampPalette(c("blue", "white", "red"))(num_colors)

png("ptm-ptm_cor_cell_lines_v2.png", 
    height = 3000,  # height in pixels
    width = 3000,   # width in pixels
    res = 300)      # resolution in DPI
corrplot(correlation_matrix_reordered, 
         method = "square",
         type = "lower",
         tl.pos = "ld",
         tl.cex = 0.9,
         tl.col = "black",
         col = rev(COL2('RdBu', 200)),
         diag = TRUE,
         cl.cex = 1.5,
         order = "original")
dev.off()

