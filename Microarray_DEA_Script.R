if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("oligo", "limma", "GEOquery", "annotate", "hgu133plus2.db"))
BiocManager::install("pd.hta.2.0")
library(oligo)
library(limma)
library(pd.hta.2.0)

# Read and normalize CEL files
#Read CEL Files
targets <- read.delim("target.txt", stringsAsFactors = FALSE)
cel_files <- targets$FileName
cel_files

# Read HTA-2.0 CEL files data using oligo package
rawData <- read.celfiles(cel_files)
exp.rawData <- exprs(rawData)


#=====================#
# RMA Normalization
#=====================#
# RMA normalization to transcript/gene level
normData  <- rma(rawData)

# Get expression matrix (normaized)
exprs_matrix <- exprs(normData)

# Save the normalized expression matrix
write.csv(exprs_matrix, "ExpSet_PostNorm_exprs_matrix.csv", quote = F)



#----------------------#
#   Box Plot
#----------------------#
#par(mfrow=c(1,2)) # run this line of code only if need 2 figures in same plot

# Get short sample names before first underscore
short_names <- sapply(sampleNames(rawData), 
                      function(x) strsplit(x, "_")[[1]][1])

# Boxplot before Normalization
#tiff(file = "Boxplot_Pre-Normalization.tiff", width = 800, height = 600)
#par(mar = c(10, 4, 4, 2))
png("Boxplot_Pre-Normalization.png")
par(mfrow = c(1, 1))
boxplot(rawData, target = "core",
        col = ifelse(short_names == "Case", "red", "blue"),
        names = short_names,
        las = 2, cex.axis = 0.9, #axis label text
        main = "Boxplot Pre-Normalization",
        ylab = "Intensities")
dev.off()


#----Boxplot after Normalization
#tiff(file="Boxplot_Post-Normalization.tiff", bg="transparent", width=400, height=500)
png("Boxplot_Post-Normalization.png")
par(mar = c(10, 4, 4, 2)) # This sets the plot margins
boxplot(exprs_matrix,
        names = short_names,
        col = ifelse(short_names == "Case", "red", "blue"),
        las = 2, #vertical x-axis label text
        cex.axis = 0.9, #axis label text
        main="Boxplot Post-Normalization", ylab="Intensities")
dev.off()



#--------------------------#
# Dimensionality reduction
#  PCA Plot
#-------------------------#
library(gplots)
library(ggplot2)
# Transpose expression matrix: rows = samples
exprs_t <- t(exprs_matrix)

# Run PCA
pca_result <- prcomp(exprs_t, scale. = TRUE)

# Calculate variance explained
pca_var <- pca_result$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)

# Prepare data for plotting
pca_df <- data.frame(
  Sample = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = short_names
)

# PCA plot with % variance and cluster ellipses (no labels)
png("PCA.png")
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.9) +
  stat_ellipse(type = "t", linetype = "dashed", level = 0.95) +
  labs(
    title = "PCA of RMA-normalized Expression Data",
    x = paste0("PC1 (", pca_var_perc[1], "%)"),
    y = paste0("PC2 (", pca_var_perc[2], "%)"),
    color = "Sample Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")
dev.off()


#===============================================#
# R Script for Gene Expression Analysis         #
# ----------------------------------------------#
# DEG Analysis - Part 2 (Downstream Analysis)
#     1. define conditions
#     2. Design matrix for model
#     3. Perform linear modelling (using limma package)
#     4. Optimize the model & perform statistical analysis for comparision
#     5. Retrieve list of expressed genes

#=======================#
#   DEG Identification
#=======================#
# Model Matrix Design
# Let's create a model matrix using the factor() function to represent the condition labels
condition <- factor(c("Case", "Case", "Case", "Case", "Case", 
                      "Fn", "Fn", "Fn", "Fn", "Fn"), 
                      levels = c("Fn", "Case")) # building the model by treating 'Case' as the reference (baseline) i.e. +veFC => gene is upregulated in Fn compared to Case
design <- model.matrix(~0 + condition)
# Now assign conditions to the columns of the model matrix
colnames(design) <- levels(condition)
design

# Fits a linear model for each gene based on the given series of arrays
# It estimates the relationship between gene expression and conditions
# This creates a model with difference between Case - Fn
fit <- lmFit(exprs_matrix, design) # exprs_matrix1[,1:ncol(exprs_matrix1)]

# Contrast Matrix Design
# Define the specific comparison between conditions to analyze
cont.matrix = makeContrasts(Fn-Case, levels=design)

# Fitting model with Contrasts(2 Groups), apply the defined contrast to the previously fitted model (fit)
fit2 <- contrasts.fit(fit, cont.matrix)

# Model optimization / Empirical Bayes Moderation
# Improves the estimation of variances for genes with low expression
# Computes moderated t-statistics and log-odds (B-stats) of differential expression by empirical Bayes shrinkage of the standard errors towards a common value
# contrast-specific information from fit2 is incorporated into the fit object, which is then passed to eBayes()
fit2 <- eBayes(fit2)  

# Get top DEGs
top_genes <- topTable(fit2, adjust="BH", sort.by="logFC", number = Inf) #number=100000);
head(top_genes)
write.csv(top_genes, "Result_Table_logFCsorted.csv", quote = F)

# after Manual Annaotation
library(dplyr)
de_annotated <- read.csv("annotated_topTable.csv", h = T, sep = ",")
deg_filtered_FN <- filter(de_annotated, P.Value < 0.05 & abs(logFC) > 0.5)
cat("Number of DEGs:", nrow(deg_filtered_FN), "\n")
write.csv(deg_filtered_FN, "DEG_filtered_logFC.5.csv", row.names = FALSE)




#================================================#
# DEG Viz (Volcano plot)
#================================================#
#install.packages("gdata")
#install.packages("gplots")
library(gdata)
library(gplots)

png(filename = "VolcanoPlot_FC_2.png")
with(top_genes, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot"))
with(subset(top_genes, P.Value < 0.05 & logFC > 0.5 ), points(logFC, -log10(P.Value), pch=20, col="red"))
with(subset(top_genes, P.Value < 0.05 & logFC < -0.5), points(logFC, -log10(P.Value), pch=20, col="green"))
dev.off()



#-----------------------------#
#  HEATMAP for filtered DEGs
#----------------------------#
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Filter required number of DEGs to make heatmap
de_annotated <- read.csv("annotated_topTable.csv", h = T, sep = ",")
deg_filtered_FN <- filter(de_annotated, P.Value < 0.05 & abs(logFC) > 3)
deg_filtered_FN$SYMBOL <- trimws(deg_filtered_FN$SYMBOL)
cat("Number of DEGs:", nrow(deg_filtered_FN), "\n")
write.csv(deg_filtered_FN, "DEG_filtered_heatmap.csv", row.names = FALSE)

# Annotated normalized counts data is in file "annotated_ExpSet_PostNorm_exprs_matrix.csv"
# and it has a column named 'SYMBOL' and samples column
annotated_norm_counts <- read.csv("annotated_ExpSet_PostNorm_exprs_matrix.csv", h = T, sep = ",")
# delete unneccesary columns X & probeset_id
annotated_norm_counts$X <- NULL
annotated_norm_counts$probeset_id <- NULL

# HEATMAP matrix preparation #-------------------------#
# Subset normalized counts (assumes annotated_norm_counts includes SYMBOL column)
heatmap_data <- annotated_norm_counts[annotated_norm_counts$SYMBOL %in% deg_filtered_FN$SYMBOL, ]
heatmap_data
# Remove rows with NA in SYMBOL
heatmap_data <- heatmap_data[!is.na(heatmap_data$SYMBOL), ]
# Set SYMBOL as rownames and remove the SYMBOL column
heatmap_data <- distinct(heatmap_data, SYMBOL, .keep_all = TRUE)
nrow(heatmap_data)
rownames(heatmap_data) <- heatmap_data$SYMBOL
heatmap_data <- heatmap_data[, -which(colnames(heatmap_data) == "SYMBOL")]
# Convert to matrix (adjust columns as per sample count - eg. Fn dataset has total 10 samples, so 1:10)
heatmap_matrix <- as.matrix(heatmap_data[, 1:10])
# SAMPLE ANNOTATION #-------------------------#
annotation_col <- data.frame(Condition = factor(condition, levels = unique(condition)))
rownames(annotation_col) <- colnames(heatmap_matrix)

ann_colors <- list(Condition = c("Case" = "darkred", "Fn" = "#1f78b4"))

# PLOT HEATMAP for filtered DEGs #-------------------------#
library(pheatmap)
library(RColorBrewer)
png(filename = "Heatmap_logFC_3.png")
pheatmap(heatmap_matrix,
         scale = "row",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         fontsize_row = 6,
         border_color = NA,
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         cluster_cols = F,
         cluster_rows = TRUE,
         color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
         main = "Heatmap of Filtered DEGs (|logFC| > 3, p < 0.05)",
         labels_col = gsub("_.+$", "", colnames(heatmap_matrix)),
         show_colnames = TRUE,
         fontsize_col = 8 # Adjust as needed
)
dev.off()

