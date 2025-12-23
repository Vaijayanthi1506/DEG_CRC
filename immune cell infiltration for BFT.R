# Load without row.names
expr <- read.csv("annotated_bft_norm.csv", check.names = FALSE)

# First column = gene names
gene_names <- expr[,1]

# Keep only rows where gene name is not NA or blank
valid_rows <- which(!is.na(gene_names) & gene_names != "")
expr <- expr[valid_rows, ]
gene_names <- gene_names[valid_rows]

# Remove first column (now only numeric expression values left)
expr <- expr[,-1]

# Assign unique row names
rownames(expr) <- make.unique(gene_names)

# Convert to numeric matrix
expr_mat <- as.matrix(expr)
# 2. Load immune gene sets
library(GSEABase)
geneSets <- getGmt("immune_signatures.gmt")

# 3. Run ssGSEA (new API style)
library(GSVA)
params <- ssgseaParam(expr_mat, geneSets)
gsva_res <- gsva(params)

# 4. Heatmap (using ComplexHeatmap for bold labels)
library(ComplexHeatmap)
Heatmap(gsva_res,
        name = "ssGSEA Score",
        row_names_gp = gpar(fontface = "bold", fontsize = 10),
        cluster_rows = TRUE,
        cluster_columns = TRUE)
dev.off()



# Create heatmap object
ht <- Heatmap(gsva_res,
              name = "ssGSEA Score",
              row_names_gp = gpar(fontface = "bold", fontsize = 10),
              cluster_rows = TRUE,
              cluster_columns = TRUE)

# 1. Show interactively
draw(ht)

# 2. Save as PDF
pdf("ssGSEA_heatmap.pdf", width = 8, height = 6)
draw(ht)
dev.off()

# 3. Save as PNG
png("ssGSEA_heatmap.png", width = 2000, height = 1500, res = 300)
draw(ht)
dev.off()







library(pheatmap)

# Define red-white-blue color palette
my_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Heatmap with scaling (row-wise)
pheatmap(gsva_res,
         scale = "row",                       # normalize across samples
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = my_colors,                   # custom colors
         fontsize_row = 10,
         fontsize_col = 10,
         main = "Immune Cell Enrichment BFT (ssGSEA)")
