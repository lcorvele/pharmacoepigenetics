# Load the required libraries
install.packages("readxl")
library(readxl)
library(pals)
library(pheatmap)
library(ComplexHeatmap)
library(dplyr)
library(extrafont)

#Read the Excel file into data frames
IC50 <- read_excel("C:/Users/lcorvele/R/Heatmap/Pharmacoepigenetics/Drug_sensitivity/AUC IC50 cell lines.xlsx", sheet = 1)
AUC <- read_excel("C:/Users/lcorvele/R/Heatmap/Pharmacoepigenetics/Drug_sensitivity/AUC IC50 cell lines.xlsx", sheet = 2)
PTM <- read_excel("C:/Users/lcorvele/R/Heatmap/Pharmacoepigenetics/Drug_sensitivity/AUC IC50 cell lines.xlsx", sheet = 4)
IC50 <- as.data.frame(IC50)
AUC <- as.data.frame(AUC)
PTM <- as.data.frame(PTM)

#Use the first column (PTM) as rownames of data frame
rownames(IC50) <- IC50[,1]
rownames(AUC) <- AUC[,1]
rownames(PTM) <- PTM[,1]

#Remove the first column from the dataframe, only numeric values (values to plot) remain
IC50 <- IC50[,-1]
AUC <- AUC[,-1]
PTM <- PTM[,-1]

# Convert all values in the data frame to numeric
for (i in seq_along(PTM)) {
  PTM[,i] <- as.numeric(as.character(PTM[,i]))
}

# Calculate the Spearman correlation between each column of IC50 with all columns of PTM
spearman_cor <- sapply(seq_along(IC50), function(i) {
  sapply(PTM, function(col) cor(IC50[,i], col, method = "spearman", use = "pairwise.complete.obs"))
})
colnames_IC50 <- colnames(IC50)
colnames(spearman_cor) <- colnames_IC50

# Calculate the Spearman correlation between each column of AUC with all columns of PTM
spearman_cor_AUC <- sapply(seq_along(AUC), function(i) {
  sapply(PTM, function(col) cor(AUC[,i], col, method = "spearman", use = "pairwise.complete.obs"))
})
colnames_AUC <- colnames(AUC)
colnames(spearman_cor_AUC) <- colnames_AUC

# Save the correlations to a CSV file
write.csv(spearman_cor, file = "./Pharmacoepigenetics/Drug_sensitivity/spearman_correlations_IC50_peptidoforms_cell_lines.csv", row.names = TRUE)
write.csv(spearman_cor_AUC, file = "./Pharmacoepigenetics/Drug_sensitivity/spearman_correlations_AUC_peptidoforms_cell_lines.csv", row.names = TRUE)

# OPTIONAL: Remove rows where the Peptide column contains "Unmod"
spearman_cor_df <- as.data.frame(spearman_cor)
spearman_cor_df <- spearman_cor_df[!grepl("Unmod", spearman_cor_df$PTM), ]

# Create an empty list to store significant correlations
significant_correlations <- list()
significant_correlations_AUC <- list()

# Perform correlation tests for each column in IC50
for (i in seq_along(colnames(IC50))) {
  for (j in seq_along(colnames(PTM))) {
    # Perform Spearman correlation test
    correlation_test <- cor.test(IC50[,i], PTM[,j], method = "spearman", use = "pairwise.complete.obs")
    
    # Check if p-value is significant (e.g., less than 0.05)
    if (correlation_test$p.value < 0.05) {
      # Store significant correlations in the list
      significant_correlations[[paste(colnames(IC50)[i], colnames(PTM)[j], sep = "_")]] <- correlation_test
    }
  }
}

# Perform correlation tests for each column in AUC
for (i in seq_along(colnames(AUC))) {
  for (j in seq_along(colnames(PTM))) {
    # Perform Spearman correlation test
    correlation_test <- cor.test(AUC[,i], PTM[,j], method = "spearman", use = "pairwise.complete.obs")
    
    # Check if p-value is significant (e.g., less than 0.05)
    if (correlation_test$p.value < 0.05) {
      # Store significant correlations in the list
      significant_correlations_AUC[[paste(colnames(AUC)[i], colnames(PTM)[j], sep = "_")]] <- correlation_test
    }
  }
}

# Extract relevant information from each htest object
significant_correlations_df <- lapply(significant_correlations, function(test) {
  data.frame(
    correlation_coefficient = test$estimate,
    p_value = test$p.value
  )
})
significant_correlations_AUC_df <- lapply(significant_correlations_AUC, function(test) {
  data.frame(
    correlation_coefficient = test$estimate,
    p_value = test$p.value
  )
})

# Convert the list of dataframes to a single dataframe
significant_correlations_df <- do.call(rbind, significant_correlations_df)
significant_correlations_AUC_df <- do.call(rbind, significant_correlations_AUC_df)

# Add column names
colnames(significant_correlations_df) <- c("correlation_coefficient", "p_value")
colnames(significant_correlations_AUC_df) <- c("correlation_coefficient", "p_value")

# Save the significant correlations to a CSV file
write.csv(significant_correlations_df, file = "significant_correlations_cell_lines_peptidoforms_IC50.csv", row.names = TRUE)

###############################
#Heatmap annotation (metadata)
###############################

#Read data from Excel sheet, convert to dataframe
data_for_annotation <- as.data.frame(read_excel("C:/Users/lcorvele/R/Heatmap/Pharmacoepigenetics/Drug_sensitivity/AUC IC50 cell lines.xlsx", col_names = T, sheet = 5))
colnames(data_for_annotation) <- c('Drugs', 'Drug class')
rownames(data_for_annotation) <- colnames(spearman_cor)


#Remove the first column from the dataframe, only numeric values (values to plot) remain
data_for_annotation <- data_for_annotation[,-1, drop=FALSE]

#define colors for each levels of qualitative variables
colours_drugclass <- list('Drug class' = c('Anthracycline' = 'deeppink', 'DNMT inhibitor' = 'cyan2', 'HDAC inhibitor' = 'blueviolet'))

#create the heatmap annotation
ha <- HeatmapAnnotation(df = data_for_annotation,
                        which = 'col',
                        col = colours_drugclass,
                        annotation_width = unit(c(1, 4), 'cm'),
                        gap = unit(1, 'mm'),
                        name = "Annotation_Name",
                        annotation_legend_param = list(title_gp = gpar(fontsize = 11, fontface = "bold"), 
                                                       labels_gp = gpar(fontsize = 11)))

###############################
#Plot heatmap
###############################

#define color map
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("blue","white","red"))
col_fun(seq(-1, 0, 1))

#create heatmap object
hm <- ComplexHeatmap::Heatmap(as.matrix(spearman_cor),
                              col = col_fun,
                              cluster_columns = T,
                              cluster_rows = T,
                              column_dend_height = unit(5, "mm"), 
                              heatmap_legend_param = list(title = "Correlation factor",
                                                          title_gp = gpar(fontsize = 11, fontface = "bold"), 
                                                          labels_gp = gpar(fontsize = 11)),
                              width = unit(8, "cm"), 
                              height = unit(40, "cm"),
                              #column_km = 2,
                              top_annotation = ha
                              )

#create heatmap png
png("heatmap_peptidoform_IC50_correlation_cell_lines.png", height=1300, width=1100)
draw(hm, heatmap_legend_side = "left", annotation_legend_side = "left")
dev.off()