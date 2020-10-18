## Loading R packages
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(fgsea)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)


#################################
########### Figure 4A ###########
#################################
data <- readRDS("/Users/biagi/PhD/AdipoSNAP/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/AdipoSNAP/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
data$clusters <- new_cluster

Idents(data) <- as.factor(data$timpoint)
data <- subset(data, cells = names(which(data$timpoint == '4day' | data$timpoint == 'RT' | data$timpoint == 'CL')))

Idents(data) <- data$clusters

data$UCP1 <- factor(ifelse(data@assays$SCT@data["Ucp1", ] > 0, "High", "Low"))

cols = c("#11c78b", "#800080", "#e57400", "#0000FF", "#dfdf0d")

cid <- as.numeric(1:length(unique(data$clusters)))[data$clusters]
dic_cid <- data.frame(grp = levels(data$clusters))

dims <- data.frame(sample_id = rownames(Embeddings(data, 'tsne')), 
                   cluster_id = data$clusters, 
                   x = Embeddings(data, 'tsne')[,1], 
                   y = Embeddings(data, 'tsne')[,2])

ggplot() +
  geom_point(data = dims,
             aes(x = x,
                 y = y,
                 color = factor(cluster_id)), 
             size = 1) + 
  scale_color_manual(values = rep("#F5F5F5", length(cols))) + 
  guides(color = guide_legend(override.aes = list(size=7))) + 
  labs(color = '', x = 't-SNE1', y = 't-SNE2') +
  theme_classic() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size=8)) + 
  annotate("point",
           dims[["x"]][which(dims[["cluster_id"]] == 'Ad1')],
           dims[["y"]][which(dims[["cluster_id"]] == 'Ad1')], 
           colour = cols[1], size = 1) + 
  annotate("point",
           dims[["x"]][which(dims[["cluster_id"]] == 'Ad2')],
           dims[["y"]][which(dims[["cluster_id"]] == 'Ad2')], 
           colour = cols[2], size = 1) +
  annotate("point",
           dims[["x"]][which(dims[["cluster_id"]] == 'Ad5')],
           dims[["y"]][which(dims[["cluster_id"]] == 'Ad5')], 
           colour = cols[5], size = 1) + 
  annotate("point",
           dims[names(data$UCP1[which(data$UCP1 == 'High')]), 3],
           dims[names(data$UCP1[which(data$UCP1 == 'High')]), 4], 
           colour = 'red', size = 1.5)

#pieplot UCP1 High
data <- readRDS("/Users/biagi/PhD/AdipoSNAP/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/AdipoSNAP/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster
data$cluster <- Idents(data)

Idents(data) <- as.factor(data$timpoint)
data <- subset(data, cells = names(which(data$timpoint == '4day' | data$timpoint == 'RT' | data$timpoint == 'CL')))

Idents(data) <- data$cluster
data <- subset(data, idents = list('Ad1', 'Ad2', 'Ad5'))

data$gene <- factor(ifelse(data@assays$SCT@data['Ucp1', ] > 0, "High", "Low"))

df <- data.frame(
  class = c("Ad1", "Ad2", 'Ad5'),
  n = c(sum(Idents(data) == 'Ad1'), sum(Idents(data) == 'Ad2'), sum(Idents(data) == 'Ad5')),
  value = c(round((length(names(data$cluster[data$cluster == 'Ad1'])[names(data$cluster[data$cluster == 'Ad1']) %in% names(data$gene[which(data$gene == 'High')])])/sum(Idents(data) == 'Ad1'))*100, 2), round((length(names(data$cluster[data$cluster == 'Ad2'])[names(data$cluster[data$cluster == 'Ad2']) %in% names(data$gene[which(data$gene == 'High')])])/sum(Idents(data) == 'Ad2'))*100, 2), round((length(names(data$cluster[data$cluster == 'Ad5'])[names(data$cluster[data$cluster == 'Ad5']) %in% names(data$gene[which(data$gene == 'High')])])/sum(Idents(data) == 'Ad5'))*100, 2)))
df$class <- factor(df$class, levels = c("Ad1", "Ad2", "Ad5"))

df <- df %>%
  arrange(desc(class)) %>%
  mutate(text_y = cumsum(value) - value/2)
df$pos = (cumsum(c(0, df$value)) + c(df$value / 2, .01))[1:nrow(df)]

ggplot(df, aes(x = "", y = value, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = c("#11c78b", "#800080", "#dfdf0d")) +
  theme_void() + labs(fill = "Clusters") +
  geom_text_repel(aes(x = 1.5, y = pos, label = paste0(value, "%")), nudge_x = .1, segment.size = .7, show.legend = FALSE) +
  ggtitle('Ucp1')



#################################
########### Figure 4B ###########
#################################
data <- readRDS("/Users/biagi/PhD/AdipoSNAP/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/AdipoSNAP/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster
data$cluster <- Idents(data)

data <- subset(data, idents = list('Ad1'))
data <- ScaleData(data, features = rownames(data))

data$UCP1 <- factor(ifelse(data@assays$SCT@data["Ucp1", ] > 0, "High", "Low"))
Idents(data) <- data$UCP1

markers <- FindAllMarkers(data, logfc.threshold = 0, only.pos = F)

markers <- subset(markers, p_val_adj < 0.05 & cluster == 'High')
markers <- markers[order(markers$avg_logFC, decreasing = T), ]

genes1 <- c(head(markers$gene, 10), tail(markers$gene, 10))

mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)

DoHeatmap(data, features = genes1, angle = 0, size = 5, label = F) +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression')



#################################
########### Figure 4C ###########
#################################
data <- readRDS("/Users/biagi/PhD/AdipoSNAP/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/AdipoSNAP/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster
data$cluster <- Idents(data)

data <- subset(data, idents = list('Ad1'))
data <- ScaleData(data, features = rownames(data))

data$UCP1 <- factor(ifelse(data@assays$SCT@data["Ucp1", ] > 0, "High", "Low"))
Idents(data) <- data$UCP1

markers <- FindAllMarkers(data, logfc.threshold = 0, only.pos = F)

markers <- subset(markers, p_val_adj < 0.05 & cluster == 'High')
markers <- markers[order(markers$avg_logFC, decreasing = T), ]

ort <- read.table("/Users/biagi/PhD/AdipoSNAP/Orthologs_human_mouse.txt", sep = ",", header = T)
pathways <- gmtPathways('/Users/biagi/PhD/AdipoSNAP/msigdb.v7.1.symbols.gmt')

comp1 <- rbind(head(markers, 50), tail(markers, 35))

comp1 <- merge(ort, comp1, by.x = "Mouse.gene.name", by.y = "gene")
comp1$Mouse.gene.name <- comp1$Gene.stable.ID <- comp1$Mouse.gene.stable.ID <- NULL
colnames(comp1)[1] <- "Gene"
ranks <- structure(comp1$avg_logFC, names = comp1$Gene)
ranks <- sort(ranks)
fgseaRes1 <- fgsea(pathways, ranks, minSize=5)
fgseaRes1 <- subset(fgseaRes1, pval < 0.05)
fgseaRes1 <- fgseaRes1[order(fgseaRes1$NES, decreasing = T), ]
genes1 <- fgseaRes1$leadingEdge
for (i in 1:length(genes1)) {
  fgseaRes1$genes[i] <- paste(genes1[[i]], collapse = ';')
}
fgseaRes1$leadingEdge <- NULL

pathways <- c('GO_POSITIVE_REGULATION_OF_COLD_INDUCED_THERMOGENESIS', 'GO_CELLULAR_CARBOHYDRATE_METABOLIC_PROCESS', 'GO_RESPONSE_TO_EXTRACELLULAR_STIMULUS', 'GO_MITOCHONDRION', 'GO_LIPID_OXIDATION', 'GO_CELLULAR_KETONE_METABOLIC_PROCESS', 'GO_INTRINSIC_COMPONENT_OF_PLASMA_MEMBRANE', 'GO_PROTEIN_HOMODIMERIZATION_ACTIVITY', 'GO_RESPONSE_TO_CYTOKINE', 'GO_REGULATION_OF_KINASE_ACTIVITY', 'RUAN_RESPONSE_TO_TROGLITAZONE_UP', 'GERHOLD_ADIPOGENESIS_UP')
fgseaRes1 <- fgseaRes1[which(fgseaRes1$pathway %in% pathways), ]
fgseaRes1$class <- ifelse(fgseaRes1$NES > 0, 'Up-regulated', 'Down-regulated')
fgseaRes1 <- fgseaRes1[order(fgseaRes1$class, fgseaRes1$NES), ]
fgseaRes1$order <- seq_len(nrow(fgseaRes1))

ggplot(fgseaRes1, aes(order, NES, fill = class)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  xlab(NULL) + ylab("Normalized Enrichment Score (NES)") +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "right", legend.title = element_text(size = 10), plot.title = element_text(hjust = 0.5)) +
  scale_fill_discrete(name = "Enrichment") +
  scale_x_continuous(breaks = fgseaRes1$order, labels = fgseaRes1$pathway, expand = c(0, 0)) +
  ggtitle('GSEA for UCP1 High/Low for AD1 cluster')



#################################
########### Figure 4D ###########
#################################
data <- readRDS("/Users/biagi/PhD/AdipoSNAP/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/AdipoSNAP/SCCAF/AdipocytesOnly/results/obs.csv")
new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster
data$cluster <- Idents(data)
data <- subset(data, idents = list('Ad1'))
tmp1 <- factor(ifelse(data@assays$SCT@data["Ucp1", ] > 0, "Ad1_High", "Ad1_Low"))


data <- readRDS("/Users/biagi/PhD/AdipoSNAP/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/AdipoSNAP/SCCAF/AdipocytesOnly/results/obs.csv")
new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster
data$cluster <- Idents(data)
data <- subset(data, idents = list('Ad1'))

tmp2 <- as.character(Idents(data))
tmp3 <- names(Idents(data))
tmp2[which(tmp3 %in% names(tmp1[tmp1 == 'Ad1_High']))] <- 'Ad1_High'
tmp2[which(tmp3 %in% names(tmp1[tmp1 == 'Ad1_Low']))] <- 'Ad1_Low'

tmp <- as.factor(structure(tmp2, names = tmp3))
Idents(data) <- tmp
data$cluster <- Idents(data)

data <- SeuratWrappers::RunALRA(data)
data <- ScaleData(data, features = rownames(data))

genes <- c('Acaa2', 'Ehhadh', 'Slc27a2', 'Acadm', 'Scp2', 'Acadvl', 'Hadha', 'Hadhb', 'Acaa1b', 'Abcd3', 
           'Sdhaf2', 'Pdha1', 'Pdhb', 'Idh3g', 'Mdh1', 'Sdha', 'Ndufs4', 'Sdhc', 'Suclg1', 'Mdh2', 
           'Akt1', 'Kiss1r', 'Mif', 'Cpt1b', 'Fabp3', 'Acsl5', 'Mfsd2a', 'Abcc4', 'Acsl3', 'Pla2g12a', 
           'Aldoa', 'Gpi1', 'Pgk1', 'Gapdh', 'Pfkp', 'Pkm', 'Ogdh', 'Ppara', 'Prkaa1', 'Gale', 'Entpd5', 'Stat3', 'Foxk2', 'Igf1', 
           "Daglb", "Slc2a4", "Gpd1", "Pcx", "Abhd5", "Adrb3", "Lipe", "Pck1", "Pnpla2")
genes <- genes[genes %in% rownames(data)]

mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)

DoHeatmap(data, features = genes, angle = 0, size = 5, label = F) +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression')



#################################
########### Figure 4F ###########
#################################
expressionInput <- c(Set1 = 962, Set2 = 244, Set3 = 73, Set4 = 35, `Set1&Set2` = 315, `Set1&Set3` = 16, `Set2&Set3` = 1, `Set1&Set2&Set3` = 46)

upset(fromExpression(expressionInput), set_size.show = T, shade.alpha	= 1)



#################################
########### Figure 4G ###########
#################################
data <- readRDS("/Users/biagi/PhD/AdipoSNAP/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/AdipoSNAP/SCCAF/AdipocytesOnly/results/obs.csv")

ort <- read.table("/Users/biagi/PhD/AdipoSNAP/Orthologs_human_mouse.txt", sep = ",", header = T)

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster
data$cluster <- Idents(data)

data <- subset(data, idents = list('Ad1'))
data <- ScaleData(data, features = rownames(data))

data$UCP1 <- factor(ifelse(data@assays$SCT@data["Ucp1", ] > 0, "High", "Low"))
Idents(data) <- data$UCP1

genes_High <- c('PPARD', 'HNF4A', 'ESR1', 'PPARG', 'SALL4', 'CEBPD', 'EGR1', 'NANOG', 'STAT3', 'BHLHE40')
genes_High <- ort$Mouse.gene.name[which(ort$Gene.name %in% genes_High)]
genes_High <- genes_High[genes_High %in% rownames(data)]
tabHigh <- as.matrix(data@assays$SCT@data[genes_High, names(data$UCP1)[which(data$UCP1 == 'High')]])
tabLow <- as.matrix(data@assays$SCT@data[genes_High, names(data$UCP1)[which(data$UCP1 == 'Low')]])

df <- data.frame(High = rowSums(tabHigh), Low = rowSums(tabLow))
df <- t(df)

x <- as.matrix(df)
m = apply(x, 1, mean, na.rm = T)
s = apply(x, 1, sd, na.rm = T)
res <- (x - m)/s
cn = colnames(res)

ba <- HeatmapAnnotation(
  text = anno_text(cn, rot = 0, location = unit(0.9, "npc"), just = "centre", gp = gpar(fontsize = 8)),
  annotation_height = max_text_width(cn)
)

breaks <- seq(-2,2, by= 0.1)
Heatmap(as.matrix(res), 
        bottom_annotation = ba, 
        name = "zscore", column_title = "Ad1 UCP1 High", width = 1, 
        show_row_names = T, show_column_names = F,
        cluster_rows = F, cluster_columns = T,
        col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)), 
        heatmap_height = unit(6, "cm"), row_names_gp = gpar(fontsize = 8))


genes_Low <- c('PPARG', 'SALL4', 'ESR1', 'TP63', 'AR', 'GATA2')
genes_Low <- ort$Mouse.gene.name[which(ort$Gene.name %in% genes_Low)]
genes_Low <- genes_Low[genes_Low %in% rownames(data)]
tabHigh <- as.matrix(data@assays$SCT@data[genes_Low, names(data$UCP1)[which(data$UCP1 == 'High')]])
tabLow <- as.matrix(data@assays$SCT@data[genes_Low, names(data$UCP1)[which(data$UCP1 == 'Low')]])

df <- data.frame(High = rowSums(tabHigh), Low = rowSums(tabLow))
df <- t(df)

x <- as.matrix(df)
m = apply(x, 1, mean, na.rm = T)
s = apply(x, 1, sd, na.rm = T)
res <- (x - m)/s
cn = colnames(res)

ba <- HeatmapAnnotation(
  text = anno_text(cn, rot = 0, location = unit(0.9, "npc"), just = "centre", gp = gpar(fontsize = 8)),
  annotation_height = max_text_width(cn)
)

breaks <- seq(-2,2, by= 0.1)
Heatmap(as.matrix(res), 
        bottom_annotation = ba, 
        name = "zscore", column_title = "Ad1 UCP1 Low", width = 1, 
        show_row_names = T, show_column_names = F,
        cluster_rows = F, cluster_columns = T,
        col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)), 
        heatmap_height = unit(6, "cm"), row_names_gp = gpar(fontsize = 8))