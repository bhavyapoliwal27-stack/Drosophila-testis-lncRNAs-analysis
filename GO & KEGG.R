# Load libraries
library(dplyr)
library(readr)
library(clusterProfiler)
library(org.Dm.eg.db)
library(enrichplot)

#Input files
lnc_counts <- read.table("Drosophila_melanogaster_Lnc.TSV", row.names = 1, header = TRUE)
pcg_counts <- read.table("Drosophila_melanogaster_PCG.TSV", row.names = 1, header = TRUE)

#Check files
colnames(lnc_counts)
colnames(pcg_counts)
ncol(lnc_counts)
ncol(pcg_counts)
colnames(pcg_counts) <- colnames(lnc_counts)

#merged files
combined_counts <- rbind(lnc_counts, pcg_counts)
colnames(combined_counts)
rownames(combined_counts)
combined_counts <- combined_counts[, 6:ncol(combined_counts)]

#Select top 100 genes
top100_genes <- combined_counts %>%
  mutate(mean_expr = rowMeans(.)) %>%
  arrange(desc(mean_expr)) %>%
  head(100) %>%
  rownames()

#Prepare GO/KEGG data
gene_df <- bitr(top100_genes,
                fromType = "FLYBASE",   # or "SYMBOL", depending on your IDs
                toType = "ENTREZID",
                OrgDb = org.Dm.eg.db)

#Gene Ontology
ego <- enrichGO(gene          = gene_df$ENTREZID,
                OrgDb         = org.Dm.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

#Visualize GO
barplot(ego, showCategory = 20, title = "GO Enrichment")
cnetplot(ego, categorySize="pvalue", foldChange=NULL)
ggsave("GO_enrichment_bar_CC.pdf", plot = last_plot(), width = 10, height = 8)

#Rename genes according to KEGG database
gene_kegg <- bitr_kegg(gene_df$ENTREZID, fromType = "ncbi-geneid",
                       toType = "kegg", organism = "dme")

#KEGG
ekegg <- enrichKEGG(
  gene         = gene_kegg$kegg,
  organism     = "dme",
  pvalueCutoff = 0.05
)

#Visualize KEGG
barplot(ekegg, showCategory = 20, title = "KEGG Enrichment")
cnetplot(ekegg, categorySize="pvalue", foldChange=NULL)
ggsave("KEGG_pathway_bar.pdf", plot = last_plot(), width = 10, height = 8)

