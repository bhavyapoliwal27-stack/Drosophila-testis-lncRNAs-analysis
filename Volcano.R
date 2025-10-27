#Install packages
BiocManager::install("EnhancedVolcano")
BiocManager::install("apeglm")

#Load libraries
library(DESeq2)
library(EnhancedVolcano)
library(apeglm)

#For lnc
#Input lnc files
counts <- read.table("C:\\Users\\91982\\Downloads\\Drosophila_melanogaster_Lnc.TSV", header = T, 
                     stringsAsFactors=FALSE, sep = '\t')

#Sample IDs
samples <- c("SRR4245530","SRR4245531","SRR4245532",
             "SRR4245533","SRR4245534",
             "SRR4245535","SRR4245536",
             "SRR4245537","SRR4245538")

#Prepare the data
countData <- counts[, c("Geneid", samples)]
rownames(countData) <- countData$Geneid
countData <- countData[, -1]

#Set Conditions
colData <- data.frame(row.names = samples,
                      condition = c("WT","WT","WT",
                                    "tut","tut",
                                    "bam","bam",
                                    "bgcn","bgcn")
)


#Prepare individual data
samples_bam_wt <- c("SRR4245530","SRR4245531","SRR4245532","SRR4245535","SRR4245536")
sub_counts <- countData[, samples_bam_wt]
colData_sub <- data.frame(row.names = samples_bam_wt,
                          condition = c("WT","WT","WT","bam","bam"))

#Establish DESeq2
dds_sub <- DESeqDataSetFromMatrix(countData = sub_counts,
                                  colData = colData_sub,
                                  design = ~ condition)

#Filter low counts
keep <- rowSums(counts(dds_sub)) >= 10
dds_sub <- dds_sub[keep, ]

#Perform DESEq2
dds_sub <- DESeq(dds_sub)
resultsNames(dds_sub)

#Reduce extreme log-fold change values
resLFC <- lfcShrink(dds_sub, coef=2, type="apeglm")  
res_df <- as.data.frame(resLFC)
res_df$gene <- rownames(res_df)

#Replace padj zero values
min_nonzero <- min(res_df$padj[res_df$padj > 0], na.rm = TRUE)
res_df$padj[res_df$padj == 0] <- min_nonzero

#Volcano plot
pdf('bam-WT lnc Volcano.pdf', width = 10, height = 7)
EnhancedVolcano(res_df,
                lab = res_df$gene,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = 'bam vs WT')
dev.off()

#For genes
#Input lnc files
counts1 <- read.table("C:\\Users\\91982\\Downloads\\Drosophila_melanogaster_PCG.TSV", header = T, 
                     stringsAsFactors=FALSE, sep = '\t')

#Sample IDs
samples1 <- c("SRR4245530","SRR4245531","SRR4245532",
             "SRR4245533","SRR4245534",
             "SRR4245535","SRR4245536",
             "SRR4245537","SRR4245538")

#Prepare the data
countData1 <- counts1[, c("Geneid", samples1)]
rownames(countData1) <- countData1$Geneid
countData1 <- countData1[, -1]

#Set Conditions
colData1 <- data.frame(row.names = samples1,
                      condition = c("WT","WT","WT",
                                    "tut","tut",
                                    "bam","bam",
                                    "bgcn","bgcn")
)

#Prepare individual data
samples_bgcn_wt1 <- c("SRR4245530","SRR4245531","SRR4245532","SRR4245537","SRR4245538")
sub_counts1 <- countData1[, samples_bgcn_wt1]
colData_sub1 <- data.frame(row.names = samples_bgcn_wt1,
                          condition = c("WT","WT","WT","bgcn","bgcn"))

#Establish DESeq2
dds_sub1 <- DESeqDataSetFromMatrix(countData = sub_counts1,
                                  colData = colData_sub1,
                                  design = ~ condition)

#Filter low counts
keep1 <- rowSums(counts(dds_sub1)) >= 10
dds_sub1 <- dds_sub1[keep1, ]

#Perform DESEq2
dds_sub1 <- DESeq(dds_sub1)
resultsNames(dds_sub1)

#Reduce extreme log-fold change values
resLFC1 <- lfcShrink(dds_sub1, coef=2, type="apeglm")  
res_df1 <- as.data.frame(resLFC1)
res_df1$gene <- rownames(res_df1)
 
#Replace padj zero values
min_nonzero1 <- min(res_df1$padj[res_df1$padj > 0], na.rm = TRUE)
res_df1$padj[res_df1$padj == 0] <- min_nonzero1

#Volcano plot
pdf('bgcn-WT Gene Volcano.pdf', width = 10, height = 7)
EnhancedVolcano(res_df1,
                lab = res_df1$gene,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = 'bgcn vs WT')
dev.off()


#Similarly for all mutants with lnc/Gene
colData2 <- data.frame(
  row.names = samples,
  condition = c("WT","WT","WT",   
                "mutant","mutant", 
                "mutant","mutant", 
                "mutant","mutant") 
)

keep2 <- rowSums(countData) >= 10
countData <- countData[keep2, ]

colData2$condition <- factor(colData2$condition, levels = c("WT","mutant"))


dds2 <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = colData2,
  design    = ~ condition
)

dds2 <- DESeq(dds2)

resultsNames(dds2)
resLFC2 <- lfcShrink(dds2, coef="condition_mutant_vs_WT", type="apeglm")


min_nonzero <- min(resLFC2$padj[resLFC2$padj > 0], na.rm=TRUE)
resLFC2$padj[resLFC2$padj == 0] <- min_nonzero

pdf("WT-Mutant lnc Volcano.pdf", width = 10, height = 7)
EnhancedVolcano(
  resLFC2, 
  lab = rownames(resLFC2),
  x = "log2FoldChange",  
  y = "padj",
  xlab = bquote(~Log[2]~ 'fold change'),
  ylab = bquote(~-Log[10]~ 'adj. p-value'),
  title = "Volcano Plot: Mutants vs WT",
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 2.0,
  labSize = 3.0
)

dev.off()