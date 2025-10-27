#Load libraries
library(reshape2)
library(pheatmap)
library(tidyverse)
library(ComplexHeatmap)
library(colorRamp2)
library(circlize)
library(stringr)
library(gtools)
library(WGCNA)
library(dplyr)
library(RCy3)
library(tibble)
library(magrittr)
library(RColorBrewer)

allowWGCNAThreads(nThreads = 12)

#Input and filter data
data <- read.table("Drosophila_melanogaster_counts.TSV", header = T, 
                   row.names = 1, sep = '\t')
data <- round(data, 2)
phenoData <- read.table("Drosophila_phenodata.TSV", header = T, sep = '\t')
colnames(phenoData)[2] <- 'Condition'
phenoData <- phenoData[,c(-3:-5)]
phenoData <- phenoData %>% column_to_rownames("Run")
trait <- phenoData %>% mutate(Condition=case_when(Condition=='wild type' ~ 0, Condition=='tut mutant' ~ 1,
                                                  Condition=='bam mutant' ~ 2, TRUE ~ 3))

data <- as.data.frame(data[, rownames(phenoData)])
data <- t(data)
goodgenes <- goodGenes(data)
data <- data[,goodgenes]
powers <- c(1:10, seq(12, 40, by = 2))
sft <- pickSoftThreshold(data, powerVector = powers, verbose = 5)

#pick soft threshold value
pdf('get_powers.pdf', width = 10, height = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit", type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, col = "red")
abline(h = 0.85, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

picked_power <- 14

#Create network
temp_cor <- cor
cor <- WGCNA::cor

netwk <- blockwiseModules(data,
                          power = picked_power,
                          networkType = "signed",
                          deepSplit = 2,
                          pamRespectsDendro = FALSE,
                          minModuleSize = 30,
                          maxBlockSize = 14000,
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "ER",
                          numericLabels = TRUE,
                          verbose = 3)

cor <- temp_cor

# Visualize dendrogram
mergedColors <- labels2colors(netwk$colors)
unmergedColors <- labels2colors(netwk$unmergedColors)
pdf('cluster_dendrogram.pdf')
plotDendroAndColors(netwk$dendrograms[[1]], 
                    cbind(mergedColors[netwk$blockGenes[[1]]], unmergedColors[netwk$blockGenes[[1]]]),
                    c("Merged","Unmerged"), dendroLabels = FALSE, addGuide = TRUE, hang = 0.03)
dev.off()

#Generate modules
module_df <- data.frame(gene_id = names(netwk$colors),
                        colors = mergedColors)

MEs0 <- moduleEigengenes(data, mergedColors)$eigengenes
MEs0 <- orderMEs(MEs0)
MEs0$treatment <- rownames(MEs0)

MEs0 <- MEs0 %>%
  left_join(phenoData %>% rownames_to_column("treatment"), by = "treatment")

module_order <- names(MEs0) %>% grep("ME", ., value = TRUE) %>% gsub("ME", "", .)
mME <- MEs0 %>%
  pivot_longer(starts_with("ME")) %>%
  mutate(name = gsub("ME", "", name),
         name = factor(name, levels = module_order))

write.table(mME, 'Module_trait.TSV', row.names = FALSE, quote = FALSE, sep = '\t')

MEs <- moduleEigengenes(data, mergedColors)$eigengenes
moduleTraitCor <- cor(MEs, trait[, c('Condition')], use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(data))

pdf('Module-trait_relationships.pdf')
mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
dev.off()

#Select modules of interest
modules_of_interest <- c("black","brown","yellow","magenta","greenyellow","red")
submod <- module_df %>% filter(colors %in% modules_of_interest)
row.names(module_df) <- module_df$gene_id
subexpr <- t(data)[submod$gene_id, ]

submod_df <- as.data.frame(subexpr) %>%
  mutate(gene_id = row.names(.)) %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "expression") %>%
  mutate(module = module_df[gene_id, ]$colors)

#Create tsv files
write.table(submod %>% filter(colors == "magenta"), 'magenta_module.TSV', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(submod %>% filter(colors == "yellow"), 'yellow_module.TSV', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(submod %>% filter(colors == "greenyellow"), 'greenyellow_module.TSV', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(submod %>% filter(colors == "black"), 'black_module.TSV', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(submod %>% filter(colors == "brown"), 'brown_module.TSV', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(submod %>% filter(colors == "red"), 'red_module.TSV', row.names = FALSE, sep = '\t', quote = FALSE)

#Visualize normalized expression 
pdf('normalized_expression.pdf', width = 10, height = 6)
submod_df %>%
  ggplot(aes(x = sample, y = expression, group = gene_id, color = module)) +
  geom_line(alpha = 0.3) +
  facet_grid(rows = vars(module), scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text.y = element_text(angle = 0)
  ) +
  labs(x = "Sample", y = "Normalized Expression")
dev.off()

#Prepare module-trait data
data <- read.table("Module_trait.TSV", header = TRUE, sep = '\t')
data$value <- as.numeric(data$value)
data <- dcast(data, treatment ~ name)
data <- data[mixedorder(data$treatment),]
rownames(data) <- NULL
data <- as.data.frame(t(data))
colnames(data) <- data[1, ]
data <- data[-1, ]
ann <- read.table("Drosophila_phenodata.TSV", header = T, sep = '\t')
ann <- ann[,c(1,2)]
colnames(ann)[2] <- 'Condition'
ann <- ann %>% column_to_rownames('Run')
data <- data[,rownames(ann)]
data <- data %>% arrange(SRR4245530)
data <- as.matrix(data)
rownamex <- rownames(data)
data1 <- apply(data, 2, as.numeric)
rownames(data1) <- rownamex

#colMain <- colorRampPalette(colors = c("royalblue", "white", "indianred3"))(15)
colMain <- colorRamp2(
  breaks = c(-1, 0,1),
  colors = c("royalblue", "white", "indianred3")
)


ann$Condition <- factor(ann$Condition, 
                        levels = c("wild type", "tut mutant", "bam mutant", "bcgn mutant"))
colours <- list('Condition' = c("wild type"= "#6aa84f", "tut mutant"="#3d85c6",
                                "bam mutant"= "#674ea7", "bcgn mutant" = "#f1c232" ))

colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

#Visualize module-trait relationship
pdf('Module_trait_relationship.pdf', width = 8, height = 6.5)
Heatmap(data1,
        col = colMain,
        cluster_columns = FALSE,
        cluster_rows = T,
        show_column_names = TRUE,
        column_title_side = 'bottom',
        column_title = 'Condition',show_row_dend = F,
        row_title = 'Modules', 
        show_row_names = TRUE,
        column_names_centered = TRUE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 10), 
        column_names_gp = gpar(fontsize = 8),
        #column_names_rot = 0,
        border = FALSE, 
        bottom_annotation = colAnn,
        heatmap_legend_param = list(
          title = "Correlation", 
          title_gp = gpar(fontsize = 10),
          labels_gp = gpar(fontsize = 8))
)
dev.off()

#Perform TOM 
TOM = TOMsimilarityFromExpr(t(subexpr),
                           power = picked_power)

row.names(TOM) = row.names(subexpr)
colnames(TOM) = row.names(subexpr)

#Create edgelist for WGCNA 
edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

head(edge_list)
write_delim(edge_list,
            file = "edgelist.tsv",
            delim = "\t")


#Load edgelist
edges <- read.delim("edgelist.tsv", header = TRUE, sep = "\t")

#Filter data
top_edges <- edges %>% arrange(desc(correlation)) %>% head(100)

nodes <- data.frame(
  id = unique(c(top_edges$gene1, top_edges$gene2)),
  stringsAsFactors = FALSE
)

# Collect module info for both gene1 and gene2 columns
modules1 <- top_edges %>% select(id = gene1, module = module1)
modules2 <- top_edges %>% select(id = gene2, module = module2)

# Combine and keep the first module per gene
modules <- bind_rows(modules1, modules2) %>%
  group_by(id) %>%
  summarise(module = first(module), .groups = "drop")

# Join back
nodes <- nodes %>% left_join(modules, by = "id")

#Connect to cytoscape
cytoscapePing()

#Prepare data for Cytoscape 
edges_cyto <- top_edges %>%
  rename(
    source = gene1,
    target = gene2
  ) %>%
  mutate(
    interaction = "correlated" 
  )

createNetworkFromDataFrames(
  nodes = nodes,
  edges = edges_cyto,
  title = "Gene Correlation Network",
  collection = "MyNetworks"
)


#Apply visual mappings
modules <- unique(nodes$module)

analyzeNetwork()

nodeTable <- getTableColumns("node", c("name", "Degree"))
head(nodeTable[order(-nodeTable$Degree), ])

# Use a better layout 
layoutNetwork("force-directed")


# Adjust node size or edge thickness by correlation
setEdgeLineWidthMapping(
  table.column = "correlation",
  mapping.type = "continuous",
  widths = c(1, 10),
  style.name = "default"
)

range(top_edges$correlation)

setEdgeColorMapping(
  table.column = "correlation",
  table.column.values = c(min(top_edges$correlation),
                          mean(top_edges$correlation),
                          max(top_edges$correlation)),
  colors = c("#d3d3d3", "#377eb8", "#e41a1c"),
  style.name = "default"
)

# Identify top 5 by degree
topHubs <- head(nodeTable[order(-nodeTable$Degree), "name"], 5)

# Make them larger and red
setNodeSizeBypass(topHubs, 80)
setNodeColorBypass(topHubs, "red")

fitContent()

exportImage("WGCNA Network", type = "PDF", height = 1200, width = 1200, resolution = 300)
write.csv(nodeTable[, c("name", "Degree")], "hub_list.csv", row.names = FALSE)

#Perform the same with SFP data
#Load and filter SFP data
sfp_list <- read.delim("Drosophila_melanogaster_SFPs.TSV", header = TRUE, sep = "\t")

sfp_genes <- sfp_list$Gene_ID
sfp_edges <- edges %>%
  filter(gene1 %in% sfp_genes | gene2 %in% sfp_genes)
  
degree_counts <- sfp_edges %>%
  select(gene1, gene2) %>%
  tidyr:: pivot_longer(cols = c(gene1, gene2), values_to = "gene") %>%
  count(gene, name = "degree") %>%
  arrange(desc(degree))


sfp_top_edges <- sfp_edges %>% arrange(desc(correlation)) %>% head(50)

sfp_nodes <- data.frame(
  id = unique(c(sfp_top_edges$gene1, sfp_top_edges$gene2)),
  stringsAsFactors = FALSE
)

# Collect module info for both gene1 and gene2 columns
sfp_modules1 <- sfp_top_edges %>% select(id = gene1, module = module1)
sfp_modules2 <- sfp_top_edges %>% select(id = gene2, module = module2)

# Combine and keep the first module per gene
sfp_modules <- bind_rows(sfp_modules1, sfp_modules2) %>%
  group_by(id) %>%
  summarise(module = first(module), .groups = "drop")

# Join back
sfp_nodes <- sfp_nodes %>% left_join(sfp_modules, by = "id")

#Prepare data for Cytoscape
sfp_edges_cyto <- sfp_top_edges %>%
  rename(
    source = gene1,
    target = gene2
  ) %>%
  mutate(
    interaction = "correlated"  
  )

createNetworkFromDataFrames(
  nodes = sfp_nodes,
  edges = sfp_edges_cyto,
  title = "WGCNA_SFP",
  collection = "MyNetworks"
)


#Apply visual mappings
sfp_modules <- unique(sfp_nodes$module)

analyzeNetwork()

sfp_nodeTable <- getTableColumns("node", c("name", "Degree"))
head(sfp_nodeTable[order(-sfp_nodeTable$Degree), ])

# Use a better layout (organic or force-directed)
layoutNetwork("force-directed")


# Adjust node size or edge thickness by correlation
setEdgeLineWidthMapping(
  table.column = "correlation",
  mapping.type = "continuous",
  widths = c(1, 10),
  style.name = "default"
)

range(sfp_top_edges$correlation)

setEdgeColorMapping(
  table.column = "correlation",
  table.column.values = c(min(sfp_top_edges$correlation),
                          mean(sfp_top_edges$correlation),
                          max(sfp_top_edges$correlation)),
  colors = c("#d3d3d3", "#377eb8", "#e41a1c"),
  style.name = "default"
)

# Identify top 5 by degree
sfp_topHubs <- head(sfp_nodeTable[order(-sfp_nodeTable$Degree), "name"], 5)

# Make them larger and red
setNodeSizeBypass(sfp_topHubs, 80)
setNodeColorBypass(sfp_topHubs, "red")

fitContent()

exportImage("SFP_WGCNA Network", type = "PDF", height = 1200, width = 1200, resolution = 300)
write.csv(nodeTable[, c("name", "Degree")], "SFP_hub_list.csv", row.names = FALSE)

  
