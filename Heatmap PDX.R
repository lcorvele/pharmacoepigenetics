#Setting up environment
#Install packages first before loading

library(readxl)
library(pals)
library(pheatmap)
library(ComplexHeatmap)
library(dplyr)
library(extrafont)

###############################
#Prepare the data for plotting
###############################

#Read data from Excel sheet 2, convert to dataframe
data_for_heatmap <- as.data.frame(read_excel("C:/Users/lcorvele/R/Heatmap/Pharmacoepigenetics/Baseline/250613_PDX_sPTM_data.xlsx", col_names = T, sheet = 2))

#Convert the column Protein from character to vector, this preserves the plotting order as it was in the Excel file
data_for_heatmap$PTM <- factor(data_for_heatmap$PTM, levels=unique(data_for_heatmap$PTM))

# Remove rows where the Peptide column contains "Fo" or "Deam" (optional)
data_for_heatmap <- data_for_heatmap[!grepl("Fo|Deam", data_for_heatmap$PTM), ]

# Remove rows where the Peptide column contains "H1" or "H2" (optional)
#data_for_heatmap <- data_for_heatmap[!grepl("H1|H2", data_for_heatmap$PTM), ]

# Remove rows where the Peptide column contains "Unmod" (optional)
#data_for_heatmap <- data_for_heatmap[!grepl("Unmod", data_for_heatmap$PTM), ]

#Convert back to dataframe
data_for_heatmap <- as.data.frame(data_for_heatmap)

# Replace NA values with 0
data_for_heatmap[data_for_heatmap == "NA"] <- NA
data_for_heatmap[is.na(data_for_heatmap)] <- 0

# Convert character data to numeric, excluding the first column
data_for_heatmap[, -1] <- sapply(data_for_heatmap[, -1], as.numeric)

#Use the first column (PTM) as rownames of data frame
rownames(data_for_heatmap) <- data_for_heatmap[,1]

#Remove the first column from the dataframe, only numeric values (values to plot) remain
data_for_heatmap <- data_for_heatmap[,-1]

## Only if you start with replicate data

# Extract condition from sample names
conditions <- sapply(strsplit(names(data_for_heatmap), "_"), function(x) x[1])

# Calculate the average for every condition
condition_averages <- sapply(unique(conditions), function(cond) {
  cols <- which(conditions == cond)
  rowMeans(data_for_heatmap[, cols, drop = FALSE])
})

# Create a dataframe from the average values
average_data <- as.data.frame(condition_averages)

# Write the dataframe to a CSV file
write.csv(average_data, "250613_PDX_sPTM_average_data_TALLonly.csv", row.names = TRUE)

## Calculate z scores

# Calculate z-scores for the data
z_scores <- t(apply(average_data, 1, scale))

# Extract the column names
column_names <- colnames(average_data)

# Assign the column names to the z-score data
colnames(z_scores) <- column_names

# Data matrix with z-scores
data_for_heatmap <- as.matrix(z_scores)


###############################
#Heatmap annotation (metadata)
###############################

#Read data from Excel sheet, convert to dataframe
data_for_annotation <- as.data.frame(read_excel("C:/Users/lcorvele/R/Heatmap/Pharmacoepigenetics/Baseline/Metadata_PDX.xlsx", col_names = T, sheet = 1))
colnames(data_for_annotation) <- c('PDX', 'Subtype', 'CIMP')
rownames(data_for_annotation) <- colnames(data_for_heatmap)


#Remove the first column from the dataframe, only numeric values (values to plot) remain
data_for_annotation <- data_for_annotation[,-1]

#define colors for each levels of qualitative variables
colours_subtype <- list('Subtype' = c('TAL1' = 'deeppink', 'TLX3' = 'cyan2', 'HOXA' = 'blueviolet', 'ETP-TALL' = 'chartreuse2', 'No data' = 'darkgrey', 'NKX2.1' = 'gold'))
colours_cimp <- list('CIMP' = c('low' = '#156082', 'high' = '#E97132'))
all_colours <- c(colours_subtype, colours_cimp)

#create the heatmap annotation
ha <- HeatmapAnnotation(df = data_for_annotation,
                        which = 'col',
                        col = all_colours,
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
col_fun = colorRamp2(c(-4, 0, 4), c("blue","white","red"))
col_fun(seq(-4, 0, 4))

#create heatmap object
hm <- ComplexHeatmap::Heatmap(as.matrix(data_for_heatmap),
                              col = col_fun,
                              cluster_columns = T,
                              cluster_rows = T,
                              column_dend_height = unit(15, "mm"), 
                              heatmap_legend_param = list(title = "PTM level",
                                title_gp = gpar(fontsize = 11, fontface = "bold"), 
                                labels_gp = gpar(fontsize = 11)),
                              width = unit(10, "cm"), 
                              height = unit(30, "cm"),
                              column_split = 4,
                              top_annotation = ha
                              )

#create heatmap png
png("250613_heatmap_PDX_sPTM_zscores_unmod_TALL_test.png", height=1100, width=1000)
draw(hm, heatmap_legend_side = "left", annotation_legend_side = "left")
dev.off()
