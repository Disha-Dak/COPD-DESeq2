# Install required packages if you haven't done so
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "tximport", "airway", "tidyverse", "pheatmap", "EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db"),force = TRUE)
BiocManager::install(c("GenomicRanges", force = TRUE))
BiocManager::install("GenomicInfoDb", version = "3.20")
BiocManager::install("enrichplot")

# Load the necessary libraries
library(DESeq2)
library(tximport)
library(airway)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)
library(ggrepel)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GenomicRanges)
library(enrichplot)

# Load sample data
data("airway")

# This dataset contains a SummarizedExperiment object with RNA-seq counts for airway epithelial cells
se <- airway

# Convert the SummarizedExperiment to a DESeqDataSet for downstream analysis
dds <- DESeqDataSet(se, design = ~ cell + dex)

# Prefiltering: remove rows with very few reads
dds <- dds[rowSums(counts(dds)) > 1, ]

# Run the DESeq function for normalization and differential expression
dds <- DESeq(dds)

# Results: specify comparison and get differentially expressed genes
res <- results(dds, contrast = c("dex", "trt", "untrt"))

# Summary of differential expression results
summary(res)

# View the top differentially expressed genes
resOrdered <- res[order(res$padj),]
head(resOrdered)

# Visualize results: MA-plot
plotMA(res, ylim = c(-5, 5))

# Transform counts for clustering/heatmap
dds_rlog <- rlog(dds, blind = TRUE)

# Heatmap of top 20 differentially expressed genes
topVarGenes <- head(order(rowVars(assay(dds_rlog)), decreasing = TRUE), 20)
mat <- assay(dds_rlog)[topVarGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = as.data.frame(colData(dds_rlog)[, c("cell", "dex")]))

# Volcano plot for differential expression using EnhancedVolcano
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano Plot',
                pCutoff = 0.05,
                FCcutoff = 1.5)

# Gene Ontology (GO) enrichment analysis
# Extract the list of significantly upregulated genes
sigGenes <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange > 1)]

# Use keytype = "ENSEMBL" to map Ensembl IDs to Entrez IDs
sigGenesEntrez <- mapIds(org.Hs.eg.db,
                         keys = sigGenes,
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "list") # Alternatively, use "first" if you want first mapping

# Filter to keep only unique (single) mappings
sigGenesEntrez_filtered <- sapply(sigGenesEntrez, function(x) if(length(x) == 1) x else NA)
sigGenesEntrez_filtered <- na.omit(sigGenesEntrez_filtered)

# Remove NA values (cases where no mapping was found)
sigGenesEntrez <- na.omit(sigGenesEntrez)

# Convert to character and remove NA or non-numeric values
sigGenesEntrez <- as.character(sigGenesEntrez)
sigGenesEntrez <- sigGenesEntrez[!is.na(sigGenesEntrez) & !is.na(as.numeric(sigGenesEntrez))]

# Print a few values to verify
print(head(sigGenesEntrez))

# Perform GO enrichment analysis
ego <- enrichGO(gene = sigGenesEntrez,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.1,  # Less stringent cutoff to get some results
                readable = TRUE)

# Check if there are any results and visualize
if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  # Plot GO enrichment if results are available
  barplot(ego, showCategory = 10, title = "GO Enrichment Analysis - Biological Processes")
} else {
  cat("No significant GO terms were found for the given gene list.\n")
}

# Principal Component Analysis (PCA) plot
plotPCA(dds_rlog, intgroup = c("cell", "dex"))
