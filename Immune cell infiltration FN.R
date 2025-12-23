# Load libraries
library(GSVA)
library(pheatmap)

# 1. Read your normalized expression matrix
# Make sure the first column is gene names, rest are samples
expr <- read.csv("annotated_fn_norm.csv", check.names = FALSE)

# If the first column is gene names
rownames(expr) <- expr[,1]
expr <- expr[,-1]

# If duplicates exist, GSVA will handle them better if you collapse them
expr <- expr[!duplicated(rownames(expr)), ]   # removes duplicate genes
expr_mat <- as.matrix(expr)

dim(expr_mat)        # should show genes x samples
head(rownames(expr_mat))  # confirm gene names look fine
head(colnames(expr_mat))  # confirm sample names look fine
# Read with first column as rownames
expr <- read.csv("annotated_fn_norm.csv", check.names = FALSE)

# First column = gene names
gene_names <- expr[,1]

# Remove the first column (now only samples remain)
expr <- expr[,-1]

# Assign rownames
rownames(expr) <- make.unique(gene_names)   # <-- this ensures duplicates get ".1", ".2", etc.

# Convert to matrix
expr_mat <- as.matrix(expr)

# Check
dim(expr_mat)
head(rownames(expr_mat))
head(colnames(expr_mat))










expr <- read.csv("annotated_fn_norm.csv", row.names = 1, check.names = FALSE)

# Optional: remove duplicates and junk rows
expr <- expr[!duplicated(rownames(expr)), ]
expr <- expr[grepl("^[A-Za-z0-9]+$", rownames(expr)), ]

# Convert to matrix
expr_mat <- as.matrix(expr)

expr_mat <- expr_mat[!rownames(expr_mat) %in% c("---"), ]
# 2. Load immune gene sets (GMT file)
geneSets <- getGmt("immune_signatures.gmt")

# 3. Run ssGSEA
params <- ssgseaParam(expr_mat, geneSets)
gsva_res <- gsva(params)

# 4. Heatmap of enrichment scores
pheatmap(gsva_res,
         scale = "row",                       # scale across samples
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize_row = 10,
         face = "bold"
         main = "Immune Cell Enrichment (ssGSEA)")




# Load required packages
library(pheatmap)
library(grid)  # needed for gpar()

# Generate heatmap with bold row labels
pheatmap(gsva_res,
         scale = "row",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize_row = 10,
         main = "Immune Cell Enrichment (ssGSEA)",
         labels_row = rownames(gsva_res),
         annotation_names_gp = gpar(fontface = "bold")  # bold labels
)



# Read CSV without assigning row names
expr <- read.csv("annotated_fn_norm.csv", check.names = FALSE)

# Make the first column (gene names) unique and set as rownames
gene_names <- make.unique(expr[,1])  # this adds .1, .2, etc. for duplicates
rownames(expr) <- gene_names

# Remove the first column (gene names) now that they are rownames
expr_mat <- as.matrix(expr[,-1])

# Ensure numeric
mode(expr_mat) <- "numeric"

# Check first few rows
head(expr_mat)

# Load GSVA package
BiocManager::install("GSVA", version = "3.17")
library(GSVA)

# Define single-gene gene sets
gene_list <- list(
  FOS = "FOS",
  ENPP2 = "ENPP2",
  FERMT2 = "FERMT2",
  EGR1 = "EGR1"
)

# Run GSVA using ssGSEA
gsva_res <- gsva(expr_mat, gene_list, method = "ssgsea", verbose = TRUE)

# Check the results
head(gsva_res)





# Install/load packages
if(!require(GSVA)) BiocManager::install("GSVA")
if(!require(GSEABase)) BiocManager::install("GSEABase")
if(!require(ComplexHeatmap)) BiocManager::install("ComplexHeatmap")
library(GSVA)
library(GSEABase)
library(ComplexHeatmap)


# Load without row.names first
expr <- read.csv("annotated_bft_norm.csv", check.names = FALSE)

# First column = gene names
gene_names <- expr[,1]

# Remove first column from expression matrix
expr <- expr[,-1]

# Handle duplicate gene names
rownames(expr) <- make.unique(gene_names)   # <-- ensures duplicates get ".1", ".2", etc.

# Convert to numeric matrix
expr_mat <- as.matrix(expr)

# 1. Load RNA-seq expression matrix
expr <- read.csv("annotated_bft_norm.csv", row.names = 1, check.names = FALSE)
expr <- expr[!duplicated(rownames(expr)), ]
expr <- expr[grepl("^[A-Za-z0-9]+$", rownames(expr)), ]
expr_mat <- as.matrix(expr)

# Remove junk rows
expr_mat <- expr_mat[!rownames(expr_mat) %in% c("---"), ]

# 2. Load immune gene sets (GMT)
geneSets <- getGmt("immune_signatures.gmt")


# 3. Run ssGSEA
gsva_res <- gsva(expr_mat, geneSets, method="ssgsea", kcdf="Gaussian", abs.ranking=TRUE)

gsva_res <- gsva(expr = expr_mat,
                 gset.idx.list = geneSets,
                 method = "ssgsea",
                 kcdf = "Gaussian",
                 abs.ranking = TRUE)

# 4. Heatmap with bold row labels
Heatmap(gsva_res,
        name = "ssGSEA Score",
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontface = "bold", fontsize = 10),
        cluster_rows = TRUE,
        cluster_columns = TRUE)





