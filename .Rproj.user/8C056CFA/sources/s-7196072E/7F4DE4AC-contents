---
date: "2020-10-12"
diagram: true
math: true
title: 7 - SCT Preprocessing
---

```
## Loading R packages
library(Seurat)
library(SeuratWrappers)
library(future)

## Setting up the threads for 12 workers
options(future.globals.maxSize = +Inf)
plan("multiprocess", workers = 12)

## Loading data Cellranger output format
data <- Read10X("/Users/biagi/PhD/AdipoSNAP/data/10x/")

## Reading metadata annotation
anno <- read.table("/Users/biagi/PhD/AdipoSNAP/data/GSE133486_10XAdiposeNuclei.metaData.tsv")

## Creating Seurat object
data <- CreateSeuratObject(counts = data,
                           project = "10x", 
                           meta.data = anno, 
                           min.cells = 3,
                           min.features = 200)

## Subset remaining only N condition
data <- subset(data, cells = rownames(data@meta.data)[which(data@meta.data$condition == "N")])

## Cell cycle genes
m.s.genes <- c("Gmnn", "Rad51", "Prim1", "Dscc1", "Cdca7", "Slbp", "Mcm7", "Cenpu", "Pold3", 
               "Ccne2", "Mcm4", "Polr1b", "Fen1", "Rad51ap1", "Tyms", "Usp1", "Rrm2", "Wdr76", 
               "Dtl", "Rrm1", "Gins2", "Tipin", "Hells", "Ubr7", "Chaf1b", "Clspn", "E2f8", "Mcm5", 
               "Nasp", "Pcna", "Mrpl36", "Rfc2", "Cdc45", "Casp8ap2", "Mcm6", "Exo1", "Pola1", "Cdc6", 
               "Ung", "Uhrf1", "Blm", "Msh2")
m.g2m.genes <- c("Cdk1", "Tmpo", "Smc4", "Tacc3", "Mki67", "Ckap2l", "Cks2", "Cdc25c", "Nusap1", "Kif11", 
                 "Top2a", "Cdca3", "Cks1b", "Ect2", "Ckap5", "Ckap2", "Cenpa", "Cdca2", "Ncapd2", "Aurkb", 
                 "Cenpf", "Gtse1", "Birc5", "Bub1", "Cdca8", "Anp32e", "Rangap1", "Tpx2", "Hjurp", "Lbr", 
                 "Dlgap5", "Psrc1", "Ndc80", "Nek2", "Cbx5", "Ube2c", "Gas2l3", "G2e3", "Cdc20", "Hmgb2", 
                 "Cenpe", "Nuf2", "Anln", "Ttk", "Kif2c", "Kif20b", "Aurka", "Hmmr", "Pimreg", "Cks1brt", 
                 "Tubb4b", "Kif23", "Ccnb2", "Ctcf")

## nFeature_RNA based cell thresholding
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

## Scoring cell cycle phases
data <- CellCycleScoring(data, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)

## Using regularized negative binomial regression to normalize UMI count data regressing out nCount_RNA, S.Score and G2M.Score variables
data <- SCTransform(data, vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score"), verbose = TRUE)

## Running PCA
data <- RunPCA(data)

## Defining the best number of PC that has a explained variance greater than 80%
dp <- data@reductions$pca@stdev
dp <- dp^2
for (z in 1:length(dp)) {
  soma <- sum(dp[1:z])/sum(dp)
  if (soma >= 0.8) {
    best_pc <- z
    break()
  }
}

## Clustering
data <- FindNeighbors(data, dims = 1:best_pc)
data <- FindClusters(data)

## RunUMAP
data <- RunUMAP(data, dims = 1:best_pc)

## RunTSNE
data <- RunTSNE(data, dims = 1:best_pc, max_iter = 2000, perplexity = 30, verbose = T)

## Saving RDS file
saveRDS(data, "/Users/biagi/PhD/AdipoSNAP/output/10x/10x_SCT_Processed.rds")

## Running Adaptively-thresholded Low Rank Approximation (ALRA)
data <- SeuratWrappers::RunALRA(data)
saveRDS(data, "/Users/biagi/PhD/AdipoSNAP/output/10x/10x_SCT_Processed_ALRA.rds")
```