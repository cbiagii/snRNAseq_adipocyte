library(Seurat)
library(SeuratWrappers)
library(future)

options(future.globals.maxSize = 20 * 1024 ^ 3)
plan("multiprocess", workers = 6)

data <- Read10X("/projects/cangen/coliveir/Miguel/data/10x/")
anno <- read.table("/projects/cangen/coliveir/Miguel/data/GSE133486_10XAdiposeNuclei.metaData.tsv")


data <- CreateSeuratObject(counts = data,
                           project = "10x", 
                           meta.data = anno, 
                           min.cells = 3,
                           min.features = 200)

data <- subset(data, cells = rownames(data@meta.data)[which(data@meta.data$condition == "N")])

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

#FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

data <- CellCycleScoring(data, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)

data <- SCTransform(data, vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score"), verbose = TRUE)

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
data <- RunTSNE(data, dims = 1:best_pc, max_iter = 2000, perplexity = 30, verbose = T)

saveRDS(data, "/projects/cangen/coliveir/Miguel/output/10x/10x_SCT_Processed.rds")


##########################################################
##########################################################


library(Seurat)

data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/10x_SCT_Processed.rds")

png("/projects/cangen/coliveir/Miguel/output/10x/tSNE_CellType_Markers.png", width=7, height=6, res = 500, units = "in")
DimPlot(data, group.by = 'CellTypeRefined', reduction='tsne', pt.size = 0.1)
dev.off()



#### SCCAF results
infos <- read.table("/projects/cangen/coliveir/scRNA_output/SCCAF/Adipocytes/results/obs.csv")

#result
new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster

library(scCATCH)
library(Seurat)

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