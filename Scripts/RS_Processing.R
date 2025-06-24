library(org.Hs.eg.db)
library(dplyr)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(scatterplot3d)
library(plotly)
library(biomaRt)
library(corrplot)
library(rgl)

exp.data = apply(exp.data, 2, as.numeric)
mode(exp.data)

exp.data.agg = aggregate(exp.data, by = list(data$GeneSymbol), FUN = mean)
dim(exp.data.agg)
names(exp.data.agg)
rownames(exp.data.agg) = exp.data.agg$Group.1
exp.data.agg = exp.data.agg %>% dplyr::select(-Group.1)

is.numeric(exp.data.agg)
exp.data.agg1 <- as.matrix(exp.data.agg)
exp.data.agg1 <- apply(exp.data.agg1, 2, as.numeric)
rownames(exp.data.agg1) <- rownames(exp.data.agg)
dim(exp.data.agg1)

varRow <- apply(exp.data.agg1, 1, var, na.rm = TRUE)
constRow <- (varRow == 0 | is.na(varRow))
exp.data.agg1 <- exp.data.agg1[!constRow, ]
varRow <- varRow[!constRow]
dim(exp.data.agg1)
exp = exp.data.agg1

library(DESeq2)
exp <- round(exp)

metadata$Condition <- factor(metadata$Condition, levels = c("Normal", "Glioblastoma", "Recurrent_Glioblastoma"))

dds <- DESeqDataSetFromMatrix(
  countData = exp,
  colData = metadata,
  design = ~ Condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ]

dds.run <- DESeq(dds)

res_glioblastoma_vs_normal <- results(dds.run, contrast = c("Condition", "Glioblastoma", "Normal"), alpha = 0.05)
res_glioblastoma_vs_normal <- res_glioblastoma_vs_normal[complete.cases(res_glioblastoma_vs_normal), ]

res_recurrent_vs_normal <- results(dds.run, contrast = c("Condition", "Recurrent_Glioblastoma", "Normal"), alpha = 0.05)
res_recurrent_vs_normal <- res_recurrent_vs_normal[complete.cases(res_recurrent_vs_normal), ]

summary(res_glioblastoma_vs_normal)
summary(res_recurrent_vs_normal)

sig_glioblastoma <- as.data.frame(res_glioblastoma_vs_normal)
sig_glioblastoma <- sig_glioblastoma[sig_glioblastoma$padj < 0.05, ]

sig_recurrent <- as.data.frame(res_recurrent_vs_normal)
sig_recurrent <- sig_recurrent[sig_recurrent$padj < 0.05, ]


res_glioblastoma_df <- as.data.frame(res_glioblastoma_vs_normal)
res_recurrent_df <- as.data.frame(res_recurrent_vs_normal)

plotMA(res_glioblastoma_vs_normal, main = "Glioblastoma vs Normal (MA Plot)", alpha = 0.05)
plotMA(res_recurrent_vs_normal, main = "Recurrent_Glioblastoma vs Normal (MA Plot)", alpha = 0.05)
