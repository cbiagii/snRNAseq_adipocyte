convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  return(humanx)
}


library(Seurat)

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

load("/Users/biagi/cangen/coliveir/Miguel/output/10x/metacell_SCT/db/mc2d.test_2dproj.Rda")
dims <- data.frame(x = object@sc_x, 
                   y = object@sc_y)

load("/Users/biagi/cangen/coliveir/Miguel/output/10x/metacell_SCT/db/mc.test_mc_f.Rda")
tmp1 <- data.frame(cells = names(object@mc), cols = object@mc)
tmp2 <- data.frame(cols = object@colors)
teste <- merge(tmp1, tmp2, by.x = "cols", by.y = "row.names")
teste$cols <- NULL

teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Adipocyte", marks_colors$group)], "Adipocytes", "Unknown")
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Progenitor", marks_colors$group)], "Progenitors", teste$cellType)
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Immune", marks_colors$group)], "Immunes", teste$cellType)
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Endothelial", marks_colors$group)], "Endothelials", teste$cellType)
tab <- merge(dims, teste, by.x = "row.names", by.y = "cells")

data <- readRDS("/Users/biagi/cangen/coliveir/Miguel/output/10x/10x_SCT_Processed.rds")
infos <- read.table("/Users/biagi/cangen/coliveir/scRNA_output/SCCAF/Adipocytes/results/obs.csv")

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

cells <- names(Idents(data))[which(Idents(data) == "A1" | Idents(data) == "A2" | Idents(data) == "A3" | Idents(data) == "A4")]

data <- readRDS("/Users/biagi/cangen/coliveir/Miguel/output/10x/10x_Processed.rds")
data <- subset(data, cells = cells)

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
data <- subset(data, subset = nFeature_RNA < 2000)

#data <- NormalizeData(data)

#data <- FindVariableFeatures(data)

data <- CellCycleScoring(data, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)

data <- SCTransform(data, vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score"), verbose = TRUE)
#data <- ScaleData(data, vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score"), features = rownames(data))

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

#data <- FindNeighbors(data, dims = 1:best_pc)
#data <- FindClusters(data)

data <- RunUMAP(data, dims = 1:best_pc)
data <- RunTSNE(data, dims = 1:best_pc, max_iter = 2000, perplexity = 30, verbose = T)

saveRDS(data, "/Users/biagi/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")


####################################################
####################################################
## SCCAF results
library(Seurat)

data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")

infos <- read.table("/projects/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster

png("/projects/cangen/coliveir/Miguel/output/10x/Fig2A_1.png", width=7, height=6, res = 500, units = "in")
DimPlot(data, reduction='tsne', pt.size = 1, cols = c("#ff5c5c", "purple", "#85ff00", "#43b4f9", "orange", "#ffa3b3"))
dev.off()

png("/projects/cangen/coliveir/Miguel/output/10x/Fig2B_1.png", width=15, height=9, res = 500, units = "in")
VlnPlot(data, c('Lipe', 'Adrb3', 'Pnpla2', 'Plin1', 'Fasn', 'Dgat1', 'Pck1', 'Adipoq', 'Retn', 'Cidec'))
dev.off()

## scCATCH
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


png("/projects/cangen/coliveir/Miguel/output/10x/Fig2B_2.png", width=15, height=9, res = 500, units = "in")
FeaturePlot(data, c('Lipe', 'Adrb3', 'Pnpla2', 'Plin1', 'Fasn', 'Dgat1', 'Pck1', 'Adipoq', 'Retn', 'Cidec'), cols = c("#D0CACA", '#C72901'), reduction = 'tsne')
dev.off()

new_idents <- Idents(data)
new_idents <- gsub("\\<Schwalie et al.Nature.G2, Schwalie et al.Nature.P2\\>", "SchwalieG2_P2", new_idents)
new_idents <- gsub("\\<Fat Cell\\>", "FatCell", new_idents)
new_idents <- gsub("\\<Schwalie et al.Nature.P2\\>", "SchwalieP2", new_idents)
new_idents <- gsub("\\<Schwalie et al.Nature.G4\\>", "SchwalieG4", new_idents)
Idents(data) <- factor(new_idents, levels = c("FatCell", "SchwalieP2", "SchwalieG4", "SchwalieG2_P2"))

png("/projects/cangen/coliveir/Miguel/output/10x/Fig2B_3.png", width=15, height=9, res = 500, units = "in")
VlnPlot(data, c('Lipe', 'Adrb3', 'Pnpla2', 'Plin1', 'Fasn', 'Dgat1', 'Pck1', 'Adipoq', 'Retn', 'Cidec'), )
dev.off()



library(ggpubr)
library(dplyr)
markers <- FindAllMarkers(data)
mks <- markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_logFC)

genes <- c("Nnat", "Gas6", "Hmgcs1", "Col15a1", "Flt1", "Btnl9", "Vldlr", "Grk3")

a1 <- FeaturePlot(data, features = genes, cols = c("#D0CACA", '#C72901'), reduction = 'tsne', ncol = 2, combine = F)
a2 <- VlnPlot(data, features = genes, cols = c("#f08e14", "#41ae74", "#b00cb0", "#e5e50b"), ncol = 2, combine = F)

a1 <- lapply(a1, function(x){x + xlab('') + ylab('') + NoLegend() + theme_void()})
a2 <- lapply(a2, function(x){x + xlab('') + ylab('') + NoLegend() + theme_void()})

b <- c(a1, a2)
b <- b[c(1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16)]

png("/projects/cangen/coliveir/Miguel/output/10x/Fig2C_1.png", width=10, height=9, res = 500, units = "in")
ggarrange(plotlist = b, ncol = 4, nrow = 4, legend = "none")
dev.off()





#Figure 3
library(Seurat)

data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")

infos <- read.table("/projects/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster

png("/projects/cangen/coliveir/Miguel/output/10x/Fig3A_1.png", width=15, height=9, res = 500, units = "in")
FeaturePlot(data, c('Ppara', 'Ucp1', 'Dio2', 'Prdm16', 'Elovl3', 'Cidea', 'Cpt1b', 'Cox7a1', 'Cox8b', 'Fgf21', 'Cidec', 'Tbx1'), cols = c("#D0CACA", '#C72901'), reduction = 'tsne', ncol = 4)
dev.off()

png("/projects/cangen/coliveir/Miguel/output/10x/Fig3A_2.png", width=15, height=9, res = 500, units = "in")
VlnPlot(data, c('Ppara', 'Ucp1', 'Dio2', 'Prdm16', 'Elovl3', 'Cidea', 'Cpt1b', 'Cox7a1', 'Cox8b', 'Fgf21', 'Cidec', 'Tbx1'), ncol = 4)
dev.off()


library(clusterProfiler)
library(org.Mm.eg.db)

types <- names(table(mks$cluster))
kegg_results <- list()
for (i in seq_along(types)) {
  tmp1 <- subset(mks, cluster == types[i])
  
  tmp2 <- bitr(geneID = tmp1$gene, 
               fromType = "SYMBOL", toType = "ENTREZID", 
               OrgDb = org.Mm.eg.db, drop = TRUE)
  
  kegg <- enrichKEGG(gene = tmp2$ENTREZID, 
                     organism = "mmu")
  
  kegg_results[[i]] <- as.data.frame(kegg)
}
names(kegg_results) <- types



clu_ann[c(2,6),]
clu_ann[3,]
clu_ann[c(4,5),]



tmp2 <- bitr(geneID = gsub(" ", "", strsplit(clu_ann[1,2], ",")[[1]]), 
             fromType = "SYMBOL", toType = "ENTREZID", 
             OrgDb = org.Mm.eg.db, drop = TRUE)

kegg <- enrichKEGG(gene = tmp2$ENTREZID, 
                   organism = "mmu")

go <- enrichGO(gene = tmp2$SYMBOL, 
               OrgDb = org.Mm.eg.db, 
               keyType = "SYMBOL")


#Fat Cell
"Slc7a10"

#Schwalie et al.Nature.P2
"Col15a1, Noct, Sh3pxd2a"
"Adamts12, Auts2, Btg2, Csad, Ddx3x, Fosb, H3f3b, Irs2, Klf6, Map4k4, Neat1, Slc5a3, Zfand5"

#Schwalie et al.Nature.G4
"Kdr, Timp3"

#Schwalie et al.Nature.G2, Schwalie et al.Nature.P2
"Sgpl1, Col27a1, Psap"
"Lpl"

VlnPlot(data, features = "Lpl")