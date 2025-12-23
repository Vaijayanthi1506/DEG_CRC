#----------------------#
# FUNCTIONAL ANNOTATION
#----------------------#
# Install packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "enrichplot", "ggplot2", "DOSE"))
BiocManager::install("ReactomePA") # optional

# Load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)
library(DOSE)


# Prepare Gene List -------#
# A list of DE genes in "filtered_deg.csv" which contains ENTREZID column (or SYMBOL, which we can convert).
filtered_deg <- read.csv("DEGs_filtered.csv", h=T)

# Convert gene SYMBOLS to Entrez IDs # if required
# gene_symbols <- filtered_deg$SYMBOL
# gene_entrez <- bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# GO Enrichment Analysis ----------#
go_result <- enrichGO(gene         = filtered_deg$SYMBOL,
                      OrgDb        = org.Hs.eg.db,
                      keyType      = "SYMBOL",
                      ont          = "ALL",  # other Options: BP, CC, MF
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      readable      = TRUE)
write.csv(as.data.frame(go_result), "GO_enrichment_results.csv")

# Barplot
barplot(go_result, split = "ONTOLOGY", width = 0.9) +     # adjust if needed # bars thicker/wider
  facet_wrap(~ONTOLOGY, scales = "free") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),        # smaller y-axis text
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),      # y-axis label size
        #strip.text = element_text(size = 9), # facet title text
        legend.text = element_text(size = 8), legend.title = element_text(size = 9),       # legend title
        plot.title = element_text(size = 10, face = "bold"))


#-----------------------------#
# KEGG Pathways Analysis
#----------------------------#
kegg_result <- enrichKEGG(gene         = filtered_deg$ENTREZID,
                          organism     = 'hsa',     # human
                          pvalueCutoff = 0.05)
write.csv(as.data.frame(kegg_result), "KEGG_enrichment_results.csv")

# Convert Entrez to symbols for readability
kegg_result <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
barplot(kegg_result, showCategory = 20) + ggtitle("KEGG Pathway Enrichment") + scale_fill_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 8, face = "bold"),        # smaller y-axis text
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"),      # y-axis label size
        strip.text = element_text(size = 9, face = "bold"), # facet title text
        legend.text = element_text(size = 8, face = "bold"), legend.title = element_text(size = 9, face = "bold"),       # legend title
        plot.title = element_text(size = 10, face = "bold"))


#-----------------------------#
# REACTOME Pathways Analysis
#----------------------------#
reactome_result <- enrichPathway(gene         = filtered_deg$ENTREZID,
                                 organism     = "human",
                                 pvalueCutoff = 0.05,
                                 readable     = TRUE)
write.csv(as.data.frame(reactome_result), "Reactome_pathway_enrichment.csv")

barplot(reactome_result, showCategory = 15) + ggtitle("Reactome Pathway Barplot") + scale_fill_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),        # smaller y-axis text
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),      # y-axis label size
        strip.text = element_text(size = 9), # facet title text
        legend.text = element_text(size = 8), legend.title = element_text(size = 9),       # legend title
        plot.title = element_text(size = 10, face = "bold"))

