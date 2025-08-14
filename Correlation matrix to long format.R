# Load required library
library(igraph)
library(reshape2)
library(viridis)
library(qgraph)

#########################
# Preprocess data from correlation matrix
#########################

# Create a correlation matrix (replace this with your actual correlation matrix)
correlation_matrix <- read.csv("C:/Users/lcorvele/R/Correlogram/Pharmacoepigenetics/correlation_matrix_cell_lines_filtered.csv", row.names = 1)
# Read p-values from the CSV file
p_values <- read.csv("C:/Users/lcorvele/R/Correlogram/Pharmacoepigenetics/pval_correlation_matrix_cell_lines.csv", row.names = 1)

# Convert row names into a column
correlation_matrix$PTM <- rownames(correlation_matrix)
rownames(correlation_matrix) <- NULL
p_values$PTM <- rownames(p_values)
rownames(p_values) <- NULL

# Melt correlation matrix into long format
correlation_long <- melt(correlation_matrix, id.vars = "PTM")
p_values_long <- melt(p_values, id.vars = "PTM")

# Merge p-values with correlation_long based on source and target columns
correlation_merged <- merge(correlation_long, p_values_long, by = c("PTM", "variable"))

# Remove self-correlations and correlations with absolute values below a threshold (optional)
threshold <- 0.05
correlation_merged <- correlation_merged[correlation_merged$PTM != correlation_merged$variable & abs(correlation_merged$value.y) <= threshold, ]

# Rename columns
colnames(correlation_merged) <- c("source", "target", "correlation", "pvalue")

# Write correlation_merged to a CSV file
write.csv(correlation_merged, file = "correlation_cell_lines_merged.csv", row.names = FALSE)

#########################
# Create g object and set parameters
#########################

# Create a directed graph
g <- graph_from_data_frame(correlation_merged, directed = FALSE)

# Find histone protein information from the PTM labels by extracting everything before 'K' followed by digits
# The regular expression now accounts for periods in the histone names
histone_protein <- gsub("([H][0-9A-Z.]+)K[0-9]+.*", "\\1", V(g)$name)

# Generate the viridis color palette for the unique histone proteins
unique_histones <- unique(histone_protein)
histone_colors <- viridis(length(unique_histones), option = "D")

# Create a named vector to map proteins to the colors
names(histone_colors) <- unique_histones

# Now use the 'histone_colors' to set the color attribute for your nodes
node_colors <- histone_colors[histone_protein]

# Calculate node degrees
node_degrees <- degree(g)

# Set node color based on histone protein, Set node size based on node degrees
V(g)$color <- node_colors
V(g)$size <- node_degrees * 1  # Scale the size for better visualization

# Find the range of absolute correlation values
min_corr <- min(abs(E(g)$correlation))
max_corr <- max(abs(E(g)$correlation))

# Define the new range for edge widths
min_width <- 0.5  # This ensures that the smallest edge is still visible
max_width <- 5

# Set edge attributes (color and width) based on correlation
E(g)$color <- ifelse(E(g)$correlation < 0, adjustcolor("blue", alpha.f = 0.3), adjustcolor("red", alpha.f = 0.3))
E(g)$width <- (abs(E(g)$correlation) - min_corr) / (max_corr - min_corr) * (max_width - min_width) + min_width

# Set the Fruchterman-Reingold layout
e <- as_edgelist(g,names=FALSE)
layout <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
                                            area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1))

#########################
# PLOT DIAGRAM
#########################

# Start PNG device with desired dimensions
png("network_plot.png", width = 2000, height = 2000, res = 300)

# Plot the graph with Fruchterman-Reingold layout and straight edges
plot(g, 
     layout = layout, 
     edge.curved = 0, 
     edge.arrow.size = 0.5, 
     vertex.label.dist = 0, 
     vertex.label.cex = 0.5, 
     vertex.label.degree = 0, 
     vertex.label.rot = 45, 
     vertex.label.family = "Arial", 
     vertex.label.cex = 1.1, 
     vertex.label.color = "black",
     vertex.color = ifelse(!is.na(node_colors), node_colors, "grey"), 
     vertex.frame.color = "transparent",
     #vertex3d = TRUE
)

# Turn off the PNG device, which will save the plot to the file
dev.off()