# Pharmacoepigenetics
Custom R-scripts used in manuscript: Profiling histone post-translational modifications to identify signatures of epigenetic drug response in T-cell acute lymphoblastic leukemia

## msqrob2PTM folder:
Uses the raw peptidoform abundances exported from Progenesis QIP software as input to obtain single PTM values. 
Can also be used for differential analysis. 

Reference:
msqrob2PTM: Differential Abundance and Differential Usage Analysis of MS-Based Proteomics Data at the Posttranslational Modification and Peptidoform Level
Demeulemeester, Nina et al.
Molecular & Cellular Proteomics, Volume 23, Issue 2, 100708

## Heatmap folder:
Baseline: uses sPTM or peptidoform values for heatmap visualisation. (Fig 1, 4B)

Drug sens: calculates pearson coefficients between sensitivity metriccs (AUC or IC50) and sPTM values and creates a heatmap. (Fig 2A, 5A)

## Covariation analysis folder:
Calculate pearson correlation coefficients between PTM pairs and creates a correlogram using the corrplot package. (Fig 6A-B)
Second script is optional to parse the data to a long format that can be uploaded in Cytoscape. (Fig 6C-D)

Reference:
Wei T, Simko V (2024). R package 'corrplot': Visualization of a Correlation Matrix. (Version 0.95), https://github.com/taiyun/corrplot.

