convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  return(humanx)
}


library(Seurat)
library(SeuratWrappers)


data <- Read10X("/projects/cangen/coliveir/Miguel/data/10x/")
anno <- read.table("/projects/cangen/coliveir/Miguel/data/GSE133486_10XAdiposeNuclei.metaData.tsv")


data <- CreateSeuratObject(counts = data,
                           project = "10x", 
                           meta.data = anno, 
                           min.cells = 3,
                           min.features = 200)

data <- subset(data, cells = rownames(data@meta.data)[which(data@meta.data$condition == "N")])

m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)

#FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

data <- NormalizeData(data)

data <- FindVariableFeatures(data)

data <- CellCycleScoring(data, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)

data <- ScaleData(data, vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score"))

data <- RunPCA(data)

dp <- data@reductions$pca@stdev
dp <- dp^2
for (z in 1:length(dp)) {
  soma <- sum(dp[1:z])/sum(dp)
  if (soma >= 0.8) {
    best_pc <- z
    break()
  }
}

data <- FindNeighbors(data, dims = 1:best_pc)
data <- FindClusters(data)

data <- RunUMAP(data, dims = 1:best_pc)
data <- RunTSNE(data, dims = 1:best_pc, max_iter = 5000, perplexity = 50, verbose = T)

saveRDS(data, "/projects/cangen/coliveir/Miguel/output/10x/10x_Processed.rds")


##########################################################
##########################################################


library(Seurat)

data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/10x_Processed.rds")


### Figure 1B
# Adipocito: Adrb3
# Endoteliais: Pecam1
# Imune: Ptprc
# Progenitoras: Cd34
png("/projects/cangen/coliveir/Miguel/output/10x/tSNE_CellType_Markers.png", width=7, height=6, res = 500, units = "in")
FeaturePlot(data, c("Adrb3", "Pecam1", "Ptprc", "Cd34"), cols = c("grey", 'red'), reduction = 'tsne')
dev.off()

png("/projects/cangen/coliveir/Miguel/output/10x/tSNE_CellType_Paper.png", width=7, height=6, res = 500, units = "in")
TSNEPlot(data, group.by = 'CellTypeRefined', label = T)
dev.off()


## Figure S1A
library(clustree)
data <- FindNeighbors(data, features = VariableFeatures(data))
data <- FindClusters(data, resolution=seq(0,1,0.1))

png("/projects/cangen/coliveir/Miguel/output/10x/clustree.png", width=10, height=12, res = 500, units = "in")
clustree(data)
dev.off()

TSNEPlot(data, group.by = 'RNA_snn_res.1')


## Figure 2A
library(ggplot2)
library(ggpubr)

pt <- FeaturePlot(data, c('Lipe', 'Adrb3', 'Pnpla2', 'Plin1', 'Fasn', 'Dgat1', 'Pck1', 'Adipoq', 'Retn', 'Cidec'), cols = c("#D0CACA", '#C72901'), reduction = 'tsne', combine = F)
pt <- lapply(pt, function(x){x + xlab('') + ylab('') + NoLegend()})
figure <- ggarrange(plotlist = pt, ncol = 5, nrow = 2, legend = 'none')
png("/projects/cangen/coliveir/Miguel/output/10x/Fig2A.png", width=20, height=9, res = 500, units = "in")
annotate_figure(figure,
                left = text_grob("Expression level", color = "black", rot = 90, size = 16))
dev.off()


### Figure 1C
library(Seurat)

marks_colors <- NULL
marks_colors <- rbind(marks_colors, c("Adipocyte", "Adrb3", "blue", 1, 2))
marks_colors <- rbind(marks_colors, c("Endothelial", "Pecam1", "green", 1, 1))
marks_colors <- rbind(marks_colors, c("Immune_1", "Ptprc", "#ff748c", 1, 0.5))
marks_colors <- rbind(marks_colors, c("Immune_2", "Cd19", "#ff8fa3", 1, 0.5))
marks_colors <- rbind(marks_colors, c("Progenitor_1", "Cd34", "#ffa500", 1, 2))
marks_colors <- rbind(marks_colors, c("Progenitor_2", "Pdgfra", "#ffb732", 1, 2))
marks_colors <- as.data.frame(marks_colors)
colnames(marks_colors) <- c("group", "gene", "color", "priority", "T_fold")
marks_colors$priority <- as.integer(marks_colors$priority)
marks_colors$T_fold <- as.numeric(marks_colors$T_fold)

library(ggplot2)
load("/projects/cangen/coliveir/Miguel/output/10x/metacell/db/mc2d.test_2dproj.Rda")
dims <- data.frame(x = object@sc_x, 
                   y = object@sc_y)

load("/projects/cangen/coliveir/Miguel/output/10x/metacell/db/mc.test_mc_f.Rda")
tmp1 <- data.frame(cells = names(object@mc), cols = object@mc)
tmp2 <- data.frame(cols = object@colors)
teste <- merge(tmp1, tmp2, by.x = "cols", by.y = "row.names")
teste$cols <- NULL

teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Adipocyte", marks_colors$group)], "Adipocytes", "Unknown")
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Progenitor", marks_colors$group)], "Progenitors", teste$cellType)
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Immune", marks_colors$group)], "Immunes", teste$cellType)
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Endothelial", marks_colors$group)], "Endothelials", teste$cellType)
tab <- merge(dims, teste, by.x = "row.names", by.y = "cells")

data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/10x_Processed.rds")
Idents(data) <- factor(data$CellTypeRefined)
infos <- data@meta.data
infos <- infos[tab$Row.names, ]

final <- merge(infos, tab, by.x = "row.names", by.y = "Row.names")

new.cluster.ids <- c("P1", "P2", "P3", "E1", "I1", "I2", "I3", "I4", "A", "I5", "E2", "I6", "E3", "I7", "E4", "E5", "E6")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
Idents(data) <- factor(Idents(data), levels = c("A", "E1", "E2", "E3", "E4", "E5", "E6", "I1", "I2", "I3", "I4", "I5", "I6", "I7", "P1", "P2", "P3"))

cls <- c("#FFA500", "#329932", "#4ca64c", "#66b266", "#7fbf7f", "#99cc99", "#b2d8b2", "#ff9999", "#ffa3a3", "#ffadad", "#ffb7b7", "#ffc1c1", "#ffcccc", "#ffd6d6", "#6666ff", "#7f7fff", "#9999ff")

png("/projects/cangen/coliveir/Miguel/output/10x/Fig1C.png", width=12, height=15, res = 500, units = "in")
VlnPlot(data, c("Adrb3", "Lipe", "Pecam1", "Ptprc", "Cd34"), ncol = 1, cols = cls)
dev.off()

png("/projects/cangen/coliveir/Miguel/output/10x/Fig1C_2.png", width=7, height=6, res = 500, units = "in")
DotPlot(data, features = c("Adrb3", "Lipe", "Pecam1", "Ptprc", "Cd34"))
dev.off()



### Boxplot
library(reshape2)
library(cowplot)
tmp <- data@assays$RNA@data["Lipe", ]
tmp <- tmp[rownames(data@meta.data)]
names(tmp) <- data$timpoint
tmp <- tmp[tmp > 0]
tmp <- data.frame(values = tmp, cond = names(tmp))
pt <- ggplot(tmp, aes(x=cond, y=values, fill=cond)) +
  geom_boxplot() + theme_cowplot() + xlab("") + ylab("Score") + 
  guides(fill=guide_legend(title="Samples"))




#### SCCAF results
data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/10x_Processed.rds")

infos <- read.table("/projects/cangen/coliveir/scRNA_output/SCCAF/Adipocytes/results/obs.csv")

#res2
new_cluster <- infos$res2
names(new_cluster) <- rownames(infos)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster
png("/projects/cangen/coliveir/Miguel/output/10x/Fig3A_1.png", width=7, height=6, res = 500, units = "in")
DimPlot(data, reduction='umap')
dev.off()

#result
new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster

png("/projects/cangen/coliveir/Miguel/output/10x/Fig3A_2.png", width=7, height=6, res = 500, units = "in")
DimPlot(data, reduction='tsne')
dev.off()



library(scCATCH)
library(Seurat)
library(future)

clu_markers <- findmarkergenes(object = data,
                               species = 'Mouse',
                               cluster = 'All',
                               match_CellMatch = FALSE,
                               cancer = NULL,
                               tissue = NULL,
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05)

clu_ann <- scCATCH(object = clu_markers$clu_markers,
                   species = 'Mouse',
                   tissue = "Adipose tissue")

convertSeurat <- function(seurat_object, scCATCH_anno) {
  tmp1 <- data.frame(cluster = levels(Idents(seurat_object)))
  tmp <- merge(tmp1, scCATCH_anno, by = 'cluster', all = T)
  tmp$cell_type[which(is.na(tmp$cell_type))] <- "Unclassified"
  
  new.cluster.ids <- tmp$cell_type
  names(new.cluster.ids) <- levels(seurat_object)
  seurat_object <- RenameIdents(seurat_object, new.cluster.ids)
  
  return(seurat_object)
}

data <- convertSeurat(data, clu_ann)

png("/projects/cangen/coliveir/Miguel/output/10x/Fig2A_2.png", width=12, height=6, res = 500, units = "in")
DimPlot(data, reduction='tsne', pt.size = 1, cols = c("#f08e14", "#41ae74", "#b00cb0", "#e5e50b"))
dev.off()
















### FindAllMArkers
library(Seurat)
library(future)

#options(future.globals.maxSize = 20 * 1024 ^ 3)
plan("multiprocess", workers = 8)

data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/10x_Processed.rds")
Idents(data) <- as.factor(data$CellTypeRefined)
markers <- FindAllMarkers(data, logfc.threshold = 0)
write.csv(markers, "/projects/cangen/coliveir/Miguel/output/10x/FindAllMarkers_CellTypeRefined_2.csv", row.names = F, quote = F)

# Correlation
Lipe_exp <- FetchData(object = data, vars = "Lipe", slot = "data")
matrix_mod <- as.matrix(data@assays$RNA@data)[, rownames(Lipe_exp)[which(Lipe_exp$Lipe > 3.5)]]
gene <- as.numeric(matrix_mod["Lipe",])
correlations <- t(apply(matrix_mod, 1, function(x){c(cor = cor.test(gene, x)[["estimate"]], pvalue = cor.test(gene, x)[["p.value"]])}))
correlations <- as.data.frame(correlations)
correlations <- correlations[order(correlations$cor.cor, decreasing = T), ]