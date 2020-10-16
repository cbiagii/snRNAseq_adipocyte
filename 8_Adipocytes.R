# DESCRIPTION ####
# The following code ASSUMES all dependencies in R have been installed See file(s): 
# 1_environment_setup.R
# The purpose of this code is to generate de Seurat object reproducing the analysis performed by Rajbhandari et al., 2020

# LOAD LIBRARIES ####
# Run the following code once you have Seurat installed
library(Seurat)

## Loading metacell results
marks_colors <- NULL
marks_colors <- rbind(marks_colors, c("Adipocyte_1", "Acsl1", "#0000b3", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Adipocyte_2", "Plin4", "#0000cc", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Adipocyte_3", "Mlxipl", "#0000e6", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Adipocyte_4", "Pck1", "#0000ff", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Adipocyte_5", "Adrb3", "#1a1aff", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Endothelial_1", "Btnl9", "#00cd00", 1, 1.8))
marks_colors <- rbind(marks_colors, c("Endothelial_2", "Ushbp1", "#00b300", 1, 1.8))
marks_colors <- rbind(marks_colors, c("Endothelial_3", "Egfl7", "#009a00", 1, 1.8))
marks_colors <- rbind(marks_colors, c("Endothelial_4", "Mcf2l", "#008000", 1, 1.8))
marks_colors <- rbind(marks_colors, c("Endothelial_5", "Ptprb", "#006700", 1, 1.8))
marks_colors <- rbind(marks_colors, c("Immune_1", "Zeb2", "#ff7f7f", 1, 0.9))
marks_colors <- rbind(marks_colors, c("Immune_2", "Trps1", "#ff6666", 1, 0.9))
marks_colors <- rbind(marks_colors, c("Immune_3", "Runx1", "#ff4c4c", 1, 0.9))
marks_colors <- rbind(marks_colors, c("Immune_4", "Ptprc", "#ff3232", 1, 0.9))
marks_colors <- rbind(marks_colors, c("Immune_5", "Adap2", "#ff1919", 1, 0.9))
marks_colors <- rbind(marks_colors, c("Progenitor_1", "Dcn", "#ffff4d", 1, 2.4))
marks_colors <- rbind(marks_colors, c("Progenitor_2", "Celf2", "#ffff33", 1, 2.4))
marks_colors <- rbind(marks_colors, c("Progenitor_3", "Meg3", "#ffff1a", 1, 2.4))
marks_colors <- rbind(marks_colors, c("Progenitor_4", "Col1a2", "#ffff00", 1, 2.4))
marks_colors <- rbind(marks_colors, c("Progenitor_5", "Col3a1", "#e6e600", 1, 2.4))
marks_colors <- as.data.frame(marks_colors)
colnames(marks_colors) <- c("group", "gene", "color", "priority", "T_fold")
marks_colors$priority <- as.integer(marks_colors$priority)
marks_colors$T_fold <- as.numeric(marks_colors$T_fold)

load("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell_SCT/db/mc2d.test_2dproj.Rda")
dims <- data.frame(x = object@sc_x, 
                   y = object@sc_y)

load("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell_SCT/db/mc.test_mc_f.Rda")
tmp1 <- data.frame(cells = names(object@mc), cols = object@mc)
tmp2 <- data.frame(cols = object@colors)
teste <- merge(tmp1, tmp2, by.x = "cols", by.y = "row.names")
teste$cols <- NULL

teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Adipocyte", marks_colors$group)], "Adipocytes", "Unknown")
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Progenitor", marks_colors$group)], "Progenitors", teste$cellType)
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Immune", marks_colors$group)], "Immunes", teste$cellType)
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Endothelial", marks_colors$group)], "Endothelials", teste$cellType)
tab <- merge(dims, teste, by.x = "row.names", by.y = "cells")

## Assign sccaf clusters to metacell results
data <- readRDS("/Users/biagi/PhD/AdipoSNAP/output/10x/10x_SCT_Processed.rds")
infos <- read.table("/Users/biagi/PhD/AdipoSNAP/SCCAF/Adipocytes/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- as.factor(new_cluster)
data$clusters_sccaf <- new_cluster

infos <- data@meta.data
infos <- infos[tab$Row.names, ]

final <- merge(infos, tab, by.x = "row.names", by.y = "Row.names")

Idents(data) <- data$clusters_sccaf
new.cluster.ids <- c("E1", "P1", "P2", "A1", "A2", "P3", "P4", "I1", "A3", "A4", "E2", "I2", "I3", "P5")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
Idents(data) <- factor(Idents(data), levels = c('A1', 'A2', 'A3', 'A4', 'E1', 'E2', 'I1', 'I2', 'I3', 'P1', 'P2', 'P3', 'P4', 'P5'))

## Selecting only the 4 adipocytes clusters
cells <- names(Idents(data))[which(Idents(data) == "A1" | Idents(data) == "A2" | Idents(data) == "A3" | Idents(data) == "A4")]

## Subsetting adipocyte clusters from all
data <- readRDS("/Users/biagi/PhD/AdipoSNAP/output/10x/10x_Processed.rds")
data <- subset(data, cells = cells)

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
data <- subset(data, subset = nFeature_RNA < 2000)

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

## RunUMAP
data <- RunUMAP(data, dims = 1:best_pc)

## RunTSNE
data <- RunTSNE(data, dims = 1:best_pc, max_iter = 2000, perplexity = 30, verbose = T)

## Saving RDS file
saveRDS(data, "/Users/biagi/PhD/AdipoSNAP/output/10x/Adipocytes.rds")