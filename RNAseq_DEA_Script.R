if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma", force = TRUE, update = TRUE, ask = FALSE)
BiocManager::install("edgeR", force = TRUE, update = TRUE, ask = FALSE)
library(limma)
library(edgeR)

#library(dplyr)
#library(ggplot2)
#library(pheatmap)
#library(RColorBrewer)

rm(list = ls()) # deletes all user-defined variables
gc() # Clear Unused Memory (Garbage Collection)


# 1. Load raw count matrix -----------#
counts <- read.table("GSE113218_raw_counts_GRCh38.p13_NCBI.tsv", 
                     header = TRUE, row.names = 1, sep = "\t")

# 2. Define condition labels (Control vs BFT) -----------#
# Samples: 3 Control, 3 BFT, 3 Control, 3 BFT
condition <- factor(c(rep("Control", 3), rep("BFT", 3), 
                      rep("Control", 3), rep("BFT", 3)))

# OPTIONAL - Sample Reordering
# Explicitly set factor levels so that Controls first, then BFT
condition <- factor(condition, levels = c("Control", "BFT"))
# Reorder columns in 'counts' based on the condition
counts <- counts[, order(condition)]
# Reorder the condition vector to match new column order
condition <- condition[order(condition)]
# Check result
head(counts)  # preview reordered columns
condition            # confirm new condition order


# 3. Create DGEList object -----------#
#    DGEList() id from edger
dge <- DGEList(counts = counts, group = condition)

# Plot
#  Library Sizes (Before Filtering)
barplot(dge$samples$lib.size/1e6,
        names = condition,
        las = 2,
        col = condition,
        main = "Library Sizes (Million Reads)", ylab = "Library Size (millions)")
#add legend
legend("topleft", legend = levels(factor(condition)), fill = unique(condition), cex = 0.8)


# 4. Filter lowly expressed genes -----------#
keep <- filterByExpr(dge)
dge.f <- dge[keep, , keep.lib.sizes = FALSE]

#  compare the number of genes before and after filtering
nrow(dge)      # Total genes before filtering
nrow(dge.f)    # Genes retained after filtering

# boxplot for raw expression profiles
# Calculate logCPM BEFORE normalization
logCPM_raw <- cpm(dge.f, log=TRUE)

# Plot boxplot of logCPM (pre-normalization)
boxplot(logCPM_raw, las=2, col=condition, 
        main="Boxplot of logCPM Before Normalization")


# 5. Normalize # TMM Normalization -----------#
dge.n <- calcNormFactors(dge.f)
norm_counts <- cpm(dge.n, log=TRUE)
boxplot(norm_counts, las=2, col=condition, 
        main="Boxplot of logCPM After Normalization")

# Compare Before vs After Normalization
par(mfrow = c(1, 2))  # 1 row, 2 columns
boxplot(logCPM_raw, las=2, names = condition, col=condition, main="Before Normalization", medcol = "white", lwd = 0.5, ylab=c(-5, 10)) # Dotted line for the median
boxplot(norm_counts, las=2,  names = condition, col=condition, main="After Normalization", medcol = "white", lwd = 0.5) # Dotted line for the median
par(mfrow = c(1, 1))  # reset layout



#---------------------------------#
# ANNOTATION
# Annotating Normalized DGE Data
#--------------------------------#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(AnnotationDbi) 

gene_ids <- rownames(dge.n)
gene_anno <- select(org.Hs.eg.db,
                    keys = gene_ids,
                    keytype = "ENTREZID",
                    columns = c("SYMBOL", "GENENAME"))
head(gene_anno)
# Convert to data frame and keep ENTREZID as a column
norm_counts_df <- as.data.frame(norm_counts)
norm_counts_df$ENTREZID <- rownames(norm_counts_df)

# Merge with gene annotation (from org.Hs.eg.db)
library(dplyr)
annotated_norm_counts <- norm_counts_df %>%
  left_join(gene_anno, by = "ENTREZID") %>%
  dplyr::select(SYMBOL, everything()) %>%
  dplyr::relocate(ENTREZID, .after = last_col())  # SYMBOL first, ENTREZID last

write.csv(annotated_norm_counts, "annotated_norm_counts.csv", row.names = F)
#---------------------------------# ANNOTATION DONE




#-----------------------#
# PCA Plot
#-----------------------#
# Transpose logCPM matrix
norm_counts_t <- t(norm_counts)
# Compute Principal Components
pca_result <- prcomp(norm_counts_t, scale. = TRUE)


library(ggplot2)
# Calculate percentage of variance explained by each PC
percent_variance <- 100 * (pca_result$sdev^2) / sum(pca_result$sdev^2)
# Create a data frame for plotting
variance_df <- data.frame(
  PC = factor(1:length(percent_variance)),
  Percentage = percent_variance)
# Plot the percentage of variance explained by each PC
ggplot(variance_df, aes(x = PC, y = Percentage)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), vjust = -0.3, size = 3) +
  labs(
    x = "Principal Component",
    y = "Percentage of Variance Explained",
    title = "Percentage of Variance Explained by Each Principal Component"
  ) +
  theme_minimal()


# PCA - PC1 & PC2
pca_df <- data.frame(PC1 = pca_result$x[,1],PC2 = pca_result$x[,2],
                     Condition = condition, Sample = colnames(norm_counts))
# Load ggplot2
library(ggplot2)
# Plot with ellipses for condition groups
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, fill = Condition)) +
  #stat_ellipse(geom = "polygon", level = 0.5, alpha = 0.2, show.legend = FALSE) +  # 50% confidence
  geom_point(size = 3) +
  geom_text(aes(label = Sample), vjust = -1, size = 3, show.legend = FALSE) +
  theme_minimal() +
  labs(title = "PCA with Separated Condition")



# 6. Matric design & voom transformation -----------#
design <- model.matrix(~condition) 
v <- voom(dge.n, design, plot = TRUE)

# 7. Fit linear model & Empirical Bayes moderation -----------#
fit <- lmFit(v, design)
fit <- eBayes(fit)

# 8. DE genes (BFT vs Control) -----------#
de_results <- topTable(fit, coef = 2, number = Inf, sort.by = "logFC")
head(de_results)

#-------------------------#
# ANNOTATE DEGs
#-------------------------#
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)

# Add gene ID as a column for merging
de_results$ENTREZID <- rownames(de_results)

gene_anno <- select(org.Hs.eg.db,
                    keys = de_results$ENTREZID,  # Use ENTREZIDs from de_results
                    keytype = "ENTREZID",
                    columns = c("SYMBOL", "GENENAME"))

# Merge annotation with DE results
de_annotated <- merge(de_results, gene_anno, by = "ENTREZID", all.x = TRUE)

# Keep relevant columns and order them
de_annotated <- de_annotated[, c("SYMBOL", "logFC", "P.Value", "AveExpr", "t", "adj.P.Val", "B", "GENENAME", "ENTREZID")]
head(de_annotated)


#-------------------------#
# VOLCANO plot
#-------------------------#
# Create a status column for coloring
de_annotated$status <- "Not Significant"
de_annotated$status[de_annotated$logFC > 0.2 & de_annotated$P.Value < 0.05] <- "Upregulated"
de_annotated$status[de_annotated$logFC < -0.2 & de_annotated$P.Value < 0.05] <- "Downregulated"
de_annotated$status <- factor(de_annotated$status, levels = c("Upregulated", "Downregulated", "Not Significant"))
# Volcano plot
ggplot(de_annotated, aes(x=logFC, y=-log10(P.Value), color=status)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "darkgrey")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(P-value)", color = "Gene Status")



#-----------------------------#
#  HEATMAP for filtered DEGs
#----------------------------#
# Filter DEGs: adj.P.Val < 0.05 and |logFC| > 1
library(dplyr)
deg_filtered <- filter(de_annotated, P.Value < 0.05 & abs(logFC) > 0.6 & SYMBOL != "NA")
cat("Number of DEGs:", nrow(deg_filtered), "\n")
write.csv(deg_filtered, "deg_filtered.csv", row.names = FALSE)

deg_filtered

# HEATMAP matrix preparation #-------------------------#
# Subset normalized counts (assumes annotated_norm_counts includes SYMBOL column)
heatmap_data <- annotated_norm_counts[annotated_norm_counts$SYMBOL %in% deg_filtered$SYMBOL, ]
# Remove rows with NA in SYMBOL
heatmap_data <- heatmap_data[!is.na(heatmap_data$SYMBOL), ]
# Set SYMBOL as rownames and remove the SYMBOL column
rownames(heatmap_data) <- heatmap_data$SYMBOL
heatmap_data <- heatmap_data[, -which(colnames(heatmap_data) == "SYMBOL")]
# Convert to matrix (adjust columns as per sample count)
heatmap_matrix <- as.matrix(heatmap_data[, 1:12])
# SAMPLE ANNOTATION #-------------------------#
annotation_col <- data.frame(Condition = factor(condition, levels = unique(condition)))
rownames(annotation_col) <- colnames(heatmap_matrix)

ann_colors <- list(Condition = c("Control" = "darkred", "BFT" = "#1f78b4"))

# PLOT HEATMAP for filtered DEGs #-------------------------#
library(pheatmap)
library(RColorBrewer)
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
         main = "Heatmap of Filtered DEGs (|logFC| > 0.5, p < 0.05)")




#-------------------------#
# Correlation Plot
#-------------------------#
# Gene-Gene Correlation Plot for Top 30 DE Genes
library(corrplot)
# Subset the top 30 genes (rows = genes, columns = samples)
DEGs_expression <- heatmap_matrix[1:30, ]
# Ensure rownames have gene symbols
gene_names <- rownames(DEGs_expression)
# Transpose to get genes as columns for cor()
gene_correlation_matrix <- cor(t(DEGs_expression))
# Create the correlation plot
corrplot(gene_correlation_matrix,
         method = "circle",
         type = "full",
         order = "hclust", # Cluster by correlation
         col = colorRampPalette(c("blue", "white", "red"))(200),
         tl.col = "black",
         tl.cex = 0.8,
         cl.cex = 1,
         title = "Gene-Gene Correlation Plot for Top 30 DE Genes",
         mar = c(0,0,2,0),
         addrect = 3)  # Draw rectangles around 3 clusters (adjust as needed))  # Add margin for title
