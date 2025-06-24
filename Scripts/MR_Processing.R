library(affy)
library(affyPLM) 
library(GEOquery)
library(limma)
library(readxl)
library(dplyr)
library(tidyverse)
library(pheatmap)
library(affy)
library(affyPLM)
library(hgu133plus2.db)
library(ggplot2)


celFilesID = "GSE50161"
cels = list.files(pattern = "CEL")
cels

affyData = ReadAffy(celfile.path = celFilesDirectory,  cdfname = "hgu133plus2cdf")
affyData

eset = threestep(affyData,
                 background.method = "IdealMM",  
                 normalize.method = "quantile",  
                 summary.method = "average.log") 

range(exprs(eset))    
exp.data = exprs(eset)

exp.data = apply(exp.data, 2, as.numeric)
mode(exp.data)

exp.data.agg = aggregate(exp.data, by = list(data2$symbol), FUN = mean)
dim(exp.data.agg)

names(exp.data.agg)
rownames(exp.data.agg) = exp.data.agg$Group.1

exp.data.agg = exp.data.agg %>% 
  select(-Group.1)

names(exp.data.agg) <- sub("_.*", "", names(exp.data.agg))

dim(exp.data.agg)

varRow <- apply(exp.data.agg, 1, var, na.rm = TRUE)
constRow <- (varRow == 0 | is.na(varRow))
sum(constRow)
exp.data.agg <- exp.data.agg[!constRow, ]
dim(exp.data.agg)

exp <- exp.data.agg
all(colnames(exp) == metadata$SampleID)

metadata = metadata[match(colnames(exp_scaled), metadata$SampleID), ]
unique(metadata$Condition)

design = model.matrix(~ 0 + Condition, data = metadata)
colnames(design) <- gsub("Condition", "", colnames(design))
design

contrast.matrix <- makeContrasts(
  Glioblastoma_vs_Normal = Glioblastoma - Normal,
  Ependymoma_vs_Normal = Ependymoma - Normal,
  Medulloblastoma_vs_Normal = Medulloblastoma - Normal,
  Astrocytoma_vs_Normal = Astrocytoma - Normal,  
  levels = design
)

fit <- lmFit(exp, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results_Glioblastoma_vs_Normal <- topTable(fit2, coef = "Glioblastoma_vs_Normal", number = Inf,adjust = "BH")
results_Ependymoma_vs_Normal <- topTable(fit2, coef = "Ependymoma_vs_Normal", number = Inf,adjust = "BH")
results_Medulloblastoma_vs_Normal <- topTable(fit2, coef = "Medulloblastoma_vs_Normal", number = Inf,adjust = "BH")
results_Astrocytoma_vs_Normal <- topTable(fit2, coef = "Astrocytoma_vs_Normal", number = Inf,adjust = "BH")
