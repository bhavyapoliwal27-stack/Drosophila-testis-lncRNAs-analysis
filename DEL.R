#Load Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(DESeq2)
library(circlize)
library(ComplexHeatmap)

#Input file
lnc <- read.table("Drosophila_melanogaster_Lnc.TSV", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Prepare data
selected_samples <- c("SRR4245530","SRR4245531","SRR4245532",
                      "SRR4245533","SRR4245534","SRR4245535",
                      "SRR4245536","SRR4245537","SRR4245538")

lnc_counts <- lnc[, c("Geneid", selected_samples)]

rownames(lnc_counts) <- lnc_counts$Geneid
lnc_counts <- lnc_counts[, -1]
lnc_counts <- as.matrix(lnc_counts)

conditions <- c("WT", "WT", "WT", "TUT", "TUT", "BAM", "BAM", "BGCN", "BGCN")

coldata <- data.frame(row.names = selected_samples,
                      condition = conditions)

lnc_counts <- apply(lnc_counts, 2, as.numeric)
rownames(lnc_counts) <- lnc$Geneid

#Iniate DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = lnc_counts,
  colData = coldata,
  design = ~ condition
)

dds <- DESeq(dds)

#threshold = 0.05
res <- results(dds, alpha = 0.05)
res_ordered <- res[order(res$padj), ]
head(res_ordered)

vsd <- vst(dds, blind = FALSE)

# adjust p value
res_ordered <- res[!is.na(res$padj), ]

# Order by smallest adjusted p-value
res_ordered <- res_ordered[order(res_ordered$padj), ]

# Extract the top 50 gene IDs (not indices)
top_genes <- rownames(res_ordered)[1:50]


mat <- assay(vsd)[top_genes, ]
mat <- t(scale(t(mat)))

col_fun <- colorRamp2(c(min(mat), 0, max(mat)), c("pink", "white", "purple"))

#Assign colors
condition_colors <- c("WT" = "red", "TUT" = "green",
                      "BAM" = "brown", "BGCN" = "blue")


ha <- HeatmapAnnotation(
  Condition = coldata$condition,
  col = list(Condition = condition_colors)
)

#Heatmap Visualization
png('DEL.png', width = 10, height = 7)
Heatmap(
  mat,
  name = "Expression",
  col = col_fun,
  top_annotation = ha, 
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE
)
dev.off()

#Similar script of DEG
#Input file
genes <- read.table("Drosophila_melanogaster_PCG.TSV", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

selected_samples <- c("SRR4245530","SRR4245531","SRR4245532",
                      "SRR4245533","SRR4245534","SRR4245535",
                      "SRR4245536","SRR4245537","SRR4245538")

gene_counts <- genes[, c("Geneid", selected_samples)]

rownames(gene_counts) <- gene_counts$Geneid
gene_counts <- gene_counts[, -1]
gene_counts <- as.matrix(gene_counts)

conditions1 <- c("WT", "WT", "WT", "TUT", "TUT", "BAM", "BAM", "BGCN", "BGCN")

coldata1 <- data.frame(row.names = selected_samples,
                      condition = conditions)

gene_counts <- apply(gene_counts, 2, as.numeric)
rownames(gene_counts) <- genes$Geneid


dds <- DESeqDataSetFromMatrix(
  countData = gene_counts,
  colData = coldata1,
  design = ~ condition
)

dds <- DESeq(dds)

res <- results(dds, alpha = 0.05)
res_ordered <- res[order(res$padj), ]
head(res_ordered)

vsd <- vst(dds, blind = FALSE)
# Remove NA values in padj
res_ordered <- res_ordered[!is.na(res_ordered$padj), ]

# Select top 50 most significant genes
top_genes <- rownames(res_ordered)[1:50]


mat1 <- assay(vsd)[top_genes, ]
mat1 <- t(scale(t(mat1)))

col_fun <- colorRamp2(c(min(mat1), 0, max(mat1)), c("lightblue", "white", "blue"))

condition_colors <- c("WT" = "green", "TUT" = "pink",
                      "BAM" = "brown", "BGCN" = "purple")


ha <- HeatmapAnnotation(
  Condition = coldata$condition,
  col = list(Condition = condition_colors)
)

png('DEG.png', width = 10, height = 7)
Heatmap(
  mat1,
  name = "Expression",
  col = col_fun,
  top_annotation = ha, 
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE
)
dev.off()