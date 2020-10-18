## Loading R packages
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)


##################################
########### Figure S4A ###########
##################################
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

plotList <- list()
genes <- c('Ppara', 'Dio2', 'Prdm16', 'Elovl3', 'Cox8b')
for (i in 1:length(genes)) {
  data$gene <- factor(ifelse(data@assays$SCT@data[genes[i], ] > 0, "High", "Low"))
  
  df <- data.frame(
    class = c("Ad1", "Ad2", 'Ad5'),
    n = c(sum(Idents(data) == 'Ad1'), sum(Idents(data) == 'Ad2'), sum(Idents(data) == 'Ad5')),
    value = c(round((length(names(data$cluster[data$cluster == 'Ad1'])[names(data$cluster[data$cluster == 'Ad1']) %in% names(data$gene[which(data$gene == 'High')])])/sum(Idents(data) == 'Ad1'))*100, 2), round((length(names(data$cluster[data$cluster == 'Ad2'])[names(data$cluster[data$cluster == 'Ad2']) %in% names(data$gene[which(data$gene == 'High')])])/sum(Idents(data) == 'Ad2'))*100, 2), round((length(names(data$cluster[data$cluster == 'Ad5'])[names(data$cluster[data$cluster == 'Ad5']) %in% names(data$gene[which(data$gene == 'High')])])/sum(Idents(data) == 'Ad5'))*100, 2)))
  df$class <- factor(df$class, levels = c("Ad1", "Ad2", "Ad5"))
  
  df <- df %>%
    arrange(desc(class)) %>%
    mutate(text_y = cumsum(value) - value/2)
  df$pos = (cumsum(c(0, df$value)) + c(df$value / 2, .01))[1:nrow(df)]
  
  plotList[[i]] <- ggplot(df, aes(x = "", y = value, fill = class)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0)+
    scale_fill_manual(values = c("#11c78b", "#800080", "#dfdf0d")) +
    theme_void() + labs(fill = "Clusters") +
    geom_text_repel(aes(x = 1.5, y = pos, label = paste0(value, "%")), nudge_x = .1, segment.size = .7, show.legend = FALSE) +
    ggtitle(genes[i])
}

ggarrange(plotlist = plotList)



##################################
########### Figure S4B ###########
##################################
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

mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)

DoHeatmap(data, features = markers$gene, angle = 0, size = 5, label = F) +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=3)) +
  labs(color='UCP1 Expression')



##################################
########### Figure S4C ###########
##################################
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

# SERCA2
genes <- c('Atp2a1', 'Atp2a2', 'Tmtc4', 'Adra1a', 'Arpc2')
genes <- genes[genes %in% rownames(data)]
mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)
ht1 <- DoHeatmap(data, features = genes, angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression')
y <- ht1$data %>% drop_na()
x <- y %>% group_by(Identity) %>% select(Feature, Cell, Identity, Expression) %>%
  spread(key = Feature, value = Expression)
w <- y %>% select(Feature, Cell, Expression) %>%
  spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()
pt <- Heatmap(w, cluster_columns = F)
DoHeatmap(data, features = rownames(w)[row_order(pt)], angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression') + ggtitle('SERCA2') + 
  theme(plot.title = element_text(hjust = 0.5))

# Arginine/creatine and proline metabolism v
genes <- c('Slc6a8', 'Gatm', 'Gamt', 'Ckmt1', 'Ckmt2')
genes <- genes[genes %in% rownames(data)]
mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)
ht1 <- DoHeatmap(data, features = genes, angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression')
y <- ht1$data %>% drop_na()
x <- y %>% group_by(Identity) %>% select(Feature, Cell, Identity, Expression) %>%
  spread(key = Feature, value = Expression)
w <- y %>% select(Feature, Cell, Expression) %>%
  spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()
pt <- Heatmap(w, cluster_columns = F)
DoHeatmap(data, features = rownames(w)[row_order(pt)], angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression') + ggtitle('Arginine/creatine and proline metabolism v') + 
  theme(plot.title = element_text(hjust = 0.5))

# Glycolytic process
genes <- c('Actn3', 'Adpgk', 'Aldoa', 'Aldoart1', 'Aldoart2', 'Aldob', 'Aldoc', 'App', 'Bpgm', 'Cbfa2t3', 'Ddit4', 'Dhtkd1', 'Eif6', 'Eno1', 'Eno1b', 'Eno2', 'Eno3', 'Eno4', 'Entpd5', 'Ep300', 'Esrrb', 'Fbp1', 'Foxk1', 'Foxk2', 'Gale', 'Galk1', 'Galt', 'Gapdh', 'Gapdhs', 'Gck', 'Gm3839', 'Gm10358', 'Gm11214', 'Gm12117', 'Gm15294', 'Gpd1', 'Gpi1', 'Hdac4', 'Hif1a', 'Hk1', 'Hk2', 'Hk3', 'Hkdc1', 'Htr2a', 'Ier3', 'Ifng', 'Igf1', 'Il3', 'Ins2', 'Insr', 'Jmjd8', 'Khk', 'Mif', 'Mlxipl', 'Mpi', 'Mtch2', 'Myc', 'Myog', 'Ncor1', 'Nupr1', 'Ogdh', 'Ogt', 'P2rx7', 'Pfkfb2', 'Pfkl', 'Pfkm', 'Pfkp', 'Pgam1', 'Pgam2', 'Pgk1', 'Pgk2', 'Pklr', 'Pkm', 'Ppara', 'Ppargc1a', 'Prkaa1', 'Prkaa2', 'Prkag2', 'Prkag3', 'Prxl2c', 'Psen1', 'Sirt6', 'Slc2a6', 'Slc4a1', 'Slc4a4', 'Stat3', 'Tigar', 'Tkfc', 'Tpi1', 'Trex1', 'Zbtb7a', 'Zbtb20')
genes <- genes[genes %in% rownames(data)]
mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)
ht1 <- DoHeatmap(data, features = genes, angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression')
y <- ht1$data %>% drop_na()
x <- y %>% group_by(Identity) %>% select(Feature, Cell, Identity, Expression) %>%
  spread(key = Feature, value = Expression)
w <- y %>% select(Feature, Cell, Expression) %>%
  spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()
pt <- Heatmap(w, cluster_columns = F)
DoHeatmap(data, features = rownames(w)[row_order(pt)], angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression') + ggtitle('Glycolytic Process') + 
  theme(plot.title = element_text(hjust = 0.5))



##################################
########### Figure S4H ###########
##################################
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

genes <- unique(c('SLC4A4', 'TXNIP', 'HADHB', 'HADHA', 'ITPR1', 'ATP6V1H', 'ESRRA', 'ETFDH', 'EGLN1', 'SLC25A42', 'NEAT1', 'ATXN2', 'SLC25A20', 'SUCLA2', 'HADHB', 'HADHA', 'SF3B2', 'UHRF1BP1L', 'COL27A1', 'PANK1', 'TECPR1', 'TOMM40', 'SLC20A2', 'AKAP1', 'ABCD3', 'EHHADH', 'TOB2', 'DEPTOR', 'PDK4', 'ITPR1', 'CNTNAP1', 'SIRT3', 'WDR91', 'CS', 'LRRC39', 'UHRF1BP1L', 'ESRRA', 'RBPMS', 'ACO2', 'KCNK3', 'GPD2', 'ATP1A2', 'ANK2', 'SDC4', 'FGD4', 'PKN1', 'RBPMS', 'HK2', 'PKM'))
genes <- ort$Mouse.gene.name[which(ort$Gene.name %in% genes)]
genes <- genes[genes %in% rownames(data)]

tabHigh <- as.matrix(data@assays$SCT@data[genes, names(data$UCP1)[which(data$UCP1 == 'High')]])
tabLow <- as.matrix(data@assays$SCT@data[genes, names(data$UCP1)[which(data$UCP1 == 'Low')]])

df <- data.frame(High = rowSums(tabHigh), Low = rowSums(tabLow))
df <- t(df)
x <- as.matrix(df)
m = apply(x, 1, mean, na.rm = T)
s = apply(x, 1, sd, na.rm = T)
res <- (x - m)/s
cn = colnames(res)

ba <- HeatmapAnnotation(
  text = anno_text(cn, rot = 0, location = unit(0.9, "npc"), just = "centre", gp = gpar(fontsize = 3)),
  annotation_height = max_text_width(cn)
)

breaks <- seq(-2,2, by= 0.1)
Heatmap(as.matrix(res), 
        bottom_annotation = ba, 
        name = "zscore", column_title = "Ad1 UCP1 Targets TF High", width = 1, 
        show_row_names = T, show_column_names = F,
        cluster_rows = F, cluster_columns = T,
        col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)), 
        heatmap_height = unit(6, "cm"), row_names_gp = gpar(fontsize = 8))


genes <- unique(c('ACSL1', 'CIDEC', 'SLC1A5', 'RETN', 'FASN', 'ADRB3', 'ABCC5', 'SH3PXD2A', 'NRIP1', 'FASN', 'SLC1A5', 'XIST', 'ACSL1', 'SH3PXD2A', 'ABCC5', 'ACACA', 'GHR', 'SH3PXD2A', 'SORBS1', 'MAPK6'))
genes <- ort$Mouse.gene.name[which(ort$Gene.name %in% genes)]
genes <- genes[genes %in% rownames(data)]

tabHigh <- as.matrix(data@assays$SCT@data[genes, names(data$UCP1)[which(data$UCP1 == 'High')]])
tabLow <- as.matrix(data@assays$SCT@data[genes, names(data$UCP1)[which(data$UCP1 == 'Low')]])

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
        name = "zscore", column_title = "Ad1 UCP1 Targets TF Low", width = 1, 
        show_row_names = T, show_column_names = F,
        cluster_rows = F, cluster_columns = T,
        col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)), 
        heatmap_height = unit(6, "cm"), row_names_gp = gpar(fontsize = 8))