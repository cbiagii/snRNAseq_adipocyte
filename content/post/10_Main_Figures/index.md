---
date: "2020-10-09"
diagram: true
math: true
title: 10 - Main Figures
---

## Figure 1
```
## Loading R packages
library(dplyr)
library(enrichR)
library(ggplot2)
library(ggpubr)
library(ggtext)
library(glue)
library(metacell)
library(Nebulosa)
library(RColorBrewer)
library(Seurat)


#################################
########### Figure 1B ###########
#################################
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

load("/Users/biagi/PhD/Adipocyte/output/10x/metacell_SCT/db/mc2d.test_2dproj.Rda")
dims <- data.frame(x = object@sc_x,
                   y = object@sc_y)

load("/Users/biagi/PhD/Adipocyte/output/10x/metacell_SCT/db/mc.test_mc_f.Rda")
tmp1 <- data.frame(cells = names(object@mc), cols = object@mc)
tmp2 <- data.frame(cols = object@colors)
teste <- merge(tmp1, tmp2, by.x = "cols", by.y = "row.names")
teste$cols <- NULL

teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Adipocyte", marks_colors$group)], "Adipocytes", "Unknown")
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Progenitor", marks_colors$group)], "Progenitors", teste$cellType)
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Immune", marks_colors$group)], "Immunes", teste$cellType)
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Endothelial", marks_colors$group)], "Endothelials", teste$cellType)
tab <- merge(dims, teste, by.x = "row.names", by.y = "cells")

data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/10x_SCT_Processed.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/Adipocytes/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- as.factor(new_cluster)
data$clusters_sccaf <- new_cluster

infos <- data@meta.data
infos <- infos[tab$Row.names, ]

final <- merge(infos, tab, by.x = "row.names", by.y = "Row.names")

A <- subset(final, cellType == "Adipocytes")
E <- subset(final, cellType == "Endothelials")
I <- subset(final, cellType == "Immunes")
P <- subset(final, cellType == "Progenitors")

Idents(data) <- data$clusters_sccaf
new.cluster.ids <- c("EN1", "PG1", "PG2", "AD1", "AD2", "PG3", "PG4", "IM1", "AD3", "AD4", "EN2", "IM2", "IM3", "PG5")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
Idents(data) <- factor(Idents(data), levels = c('AD1', 'AD2', 'AD3', 'AD4', 'EN1', 'EN2', 'IM1', 'IM2', 'IM3', 'PG1', 'PG2', 'PG3', 'PG4', 'PG5'))

cls <- c("#FFA500", "#ffae19", "#ffb732", "#ffc04c", "#66b266", "#7fbf7f", "#ffa3a3", "#ffadad", "#ffb7b7", "#4043ff", "#6666ff", "#7f7fff", "#9999ff", "#adadff")

TSNEPlot(data, label = FALSE, pt.size = 0.3, cols = cls) +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  theme_classic() + labs(color = "Cell Type") +
  theme(legend.position = "bottom")


a <- round((sum(Idents(data) == "AD1" | Idents(data) == "AD2" | Idents(data) == "AD3" | Idents(data) == "AD4")/ncol(data))*100, 2); a1 <- sum(Idents(data) == "AD1" | Idents(data) == "AD2" | Idents(data) == "AD3" | Idents(data) == "AD4")
b <- round((sum(Idents(data) == "EN1" | Idents(data) == "EN2")/ncol(data))*100, 2); b1 <- sum(Idents(data) == "EN1" | Idents(data) == "EN2")
c <- round((sum(Idents(data) == "IM1" | Idents(data) == "IM2" | Idents(data) == "IM3" | Idents(data) == "IM4")/ncol(data))*100, 2); c1 <- sum(Idents(data) == "IM1" | Idents(data) == "IM2" | Idents(data) == "IM3" | Idents(data) == "IM4")
d <- round((sum(Idents(data) == "PG1" | Idents(data) == "PG2" | Idents(data) == "PG3" | Idents(data) == "PG4")/ncol(data))*100, 2); d1 <- sum(Idents(data) == "PG1" | Idents(data) == "PG2" | Idents(data) == "PG3" | Idents(data) == "PG4")

df <- data.frame(
  class = c("Adipocyte", "Endothelial", "Immune", "Progenitor"),
  n = c(a1, b1, c1, d1),
  value = c(a, b, c, d)
)

df$class <- factor(df$class, levels = c("Adipocyte", "Endothelial", "Immune", "Progenitor"))

df <- df %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(value) - 0.5*value)

ggplot(df, aes(x = "", y = value, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = paste0(value, "%")), color = "white", fontface = "bold")+
  scale_fill_manual(values = c("#FFA500", "#329932", "#ff9999", "#6666ff")) +
  theme_void() + labs(fill = "Cell Type")



#################################
########### Figure 1C ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/10x_SCT_Processed_ALRA.rds")

pt <- plot_density(data, c("Adrb3", "Pecam1", "Ptprc", "Cd34", "Pdgfra", "Itgb1"), reduction = "tsne", pal = "inferno", combine = FALSE)
pt <- lapply(pt, function(x) {
  x + theme_bw() + xlab("t-SNE 1") + ylab("t-SNE 2") +
    theme(plot.title = element_text(hjust = 0.5, face = "italic"), 
          panel.background = element_rect(colour = "black", size = 0.1), 
          axis.ticks.length = unit(.2, "cm"), 
          axis.text = element_text(size = 11), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.position = "none")
})

ggarrange(plotlist = pt)



#################################
########### Figure 1D ###########
#################################
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

load("/Users/biagi/PhD/Adipocyte/output/10x/metacell_SCT/db/mc2d.test_2dproj.Rda")
dims <- data.frame(x = object@sc_x,
                   y = object@sc_y)

load("/Users/biagi/PhD/Adipocyte/output/10x/metacell_SCT/db/mc.test_mc_f.Rda")
tmp1 <- data.frame(cells = names(object@mc), cols = object@mc)
tmp2 <- data.frame(cols = object@colors)
teste <- merge(tmp1, tmp2, by.x = "cols", by.y = "row.names")
teste$cols <- NULL

teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Adipocyte", marks_colors$group)], "Adipocytes", "Unknown")
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Progenitor", marks_colors$group)], "Progenitors", teste$cellType)
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Immune", marks_colors$group)], "Immunes", teste$cellType)
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Endothelial", marks_colors$group)], "Endothelials", teste$cellType)
tab <- merge(dims, teste, by.x = "row.names", by.y = "cells")

data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/10x_SCT_Processed.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/Adipocytes/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- as.factor(new_cluster)
data$clusters_sccaf <- new_cluster

infos <- data@meta.data
infos <- infos[tab$Row.names, ]

final <- merge(infos, tab, by.x = "row.names", by.y = "Row.names")

Idents(data) <- data$clusters_sccaf
new.cluster.ids <- c("EN", "PG", "PG", "AD", "AD", "PG", "PG", "IM", "AD", "AD", "EN", "IM", "IM", "PG")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
Idents(data) <- factor(Idents(data), levels = c('AD', 'EN', 'IM', 'PG'))

markers <- FindAllMarkers(data, only.pos = TRUE)

top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

data2 <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/10x_SCT_Processed_ALRA.rds")

Idents(data2) <- Idents(data)
rm(data)
data2 <- ScaleData(data2)

cls <- c("#FFA500", "#329932", "#ff9999", "#6666ff")
mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)

DoHeatmap(data2, features = top20$gene, group.colors = cls, angle = 0, size = 5, label = FALSE) +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text = element_text(size = 5))



#################################
########### Figure 1E ###########
#################################
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

load("/Users/biagi/PhD/Adipocyte/output/10x/metacell/db/mc2d.test_2dproj.Rda")
dims <- data.frame(x = object@sc_x,
                   y = object@sc_y)

load("/Users/biagi/PhD/Adipocyte/output/10x/metacell/db/mc.test_mc_f.Rda")
tmp1 <- data.frame(cells = names(object@mc), cols = object@mc)
tmp2 <- data.frame(cols = object@colors)
teste <- merge(tmp1, tmp2, by.x = "cols", by.y = "row.names")
teste$cols <- NULL

teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Adipocyte", marks_colors$group)], "Adipocytes", "Unknown")
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Progenitor", marks_colors$group)], "Progenitors", teste$cellType)
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Immune", marks_colors$group)], "Immunes", teste$cellType)
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Endothelial", marks_colors$group)], "Endothelials", teste$cellType)
tab <- merge(dims, teste, by.x = "row.names", by.y = "cells")

data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/10x_SCT_Processed_ALRA.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/Adipocytes/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- as.factor(new_cluster)
data$clusters_sccaf <- new_cluster

infos <- data@meta.data
infos <- infos[tab$Row.names, ]

final <- merge(infos, tab, by.x = "row.names", by.y = "Row.names")

A <- subset(final, cellType == "Adipocytes")
E <- subset(final, cellType == "Endothelials")
I <- subset(final, cellType == "Immunes")
P <- subset(final, cellType == "Progenitors")

Idents(data) <- data$clusters_sccaf
new.cluster.ids <- c("EN1", "PG1", "PG2", "AD1", "AD2", "PG3", "PG4", "IM1", "AD3", "AD4", "EN2", "IM2", "IM3", "PG5")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
Idents(data) <- factor(Idents(data), levels = rev(c('AD1', 'AD2', 'AD3', 'AD4', 'EN1', 'EN2', 'IM1', 'IM2', 'IM3', 'PG1', 'PG2', 'PG3', 'PG4', 'PG5')))

DotPlot(data, features = c("Acsl1", "Plin4", "Mlxipl", "Pck1", "Adrb3",
                           "Btnl9", "Ushbp1", "Egfl7", "Mcf2l", "Ptprb",
                           "Trps1", "Runx1", "Ptprc", "Adap2",
                           "Dcn", "Celf2", "Meg3", "Col1a2", "Col3a1"), cols = c("grey", 'red')) +
  xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Nebulosa expression plot
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/10x_SCT_Processed_ALRA.rds")

pt <- plot_density(data, c("Mlxipl", "Egfl7", "Runx1", "Celf2"), reduction = "tsne", pal = "inferno", combine = FALSE)
pt <- lapply(pt, function(x) {
  x + theme_void() + xlab("") + ylab("") +
    theme(plot.title = element_text(hjust = 0.5, face = "italic"),
          legend.position = "none")
})

ggarrange(plotlist = pt)



#################################
########### Figure 1F ###########
#################################
dbs <- c("KEGG_2019_Mouse", "WikiPathways_2019_Mouse", "Jensen_TISSUES", "GO_Biological_Process_2018")

genes_A <- readLines("/Users/biagi/PhD/Adipocyte/Figures/update/Fig1D_2_Markers_A.txt")
genes_E <- readLines("/Users/biagi/PhD/Adipocyte/Figures/update/Fig1D_2_Markers_E.txt")
genes_I <- readLines("/Users/biagi/PhD/Adipocyte/Figures/update/Fig1D_2_Markers_I.txt")
genes_P <- readLines("/Users/biagi/PhD/Adipocyte/Figures/update/Fig1D_2_Markers_P.txt")

genes <- list(Adipocyte = genes_A,
              Endothelial = genes_E,
              Immune = genes_I,
              Progenitor = genes_P)

results <- lapply(genes, enrichr, dbs)

plotlist <- list()
pathlist <- NULL
tmplist <- NULL
paths <- NULL
for (i in 1:length(dbs)) {
  tmp <- lapply(results, `[[`, i)
  tmp <- mapply(cbind, tmp, "type" = names(tmp), SIMPLIFY = FALSE)
  tmp <- do.call("rbind", tmp)
  
  tmp <- subset(tmp, Adjusted.P.value < 0.05)
  
  tmp$Adjusted.P.value <- -log10(tmp$Adjusted.P.value)
  tmp$Combined.Score <- log2(tmp$Combined.Score)
  
  tmplist <- rbind(tmplist, data.frame(tmp, class = dbs[i]))
  
  if (dbs[i] == "GO_Biological_Process_2018") {
    tmp$Term <- gsub("\\s*\\([^\\)]+\\)","",as.character(tmp$Term))
  }
  tmp <- tmp %>% group_by(type) %>% top_n(n = 5, wt = Combined.Score)
  tmp$Term <- factor(tmp$Term, levels = rev(unique(tmp$Term)))
  pathlist <- rbind(pathlist, data.frame(tmp %>% group_by(type) %>% top_n(n = 2, wt = Combined.Score), class = dbs[i]))
  paths <- rbind(paths, data.frame(tmp, class = dbs[i]))
  
  plotlist[[i]] <- ggplot(tmp, aes(x = type, y = Term)) +
    geom_point(aes(size = Adjusted.P.value, color = Combined.Score)) +
    theme_bw(base_size = 9) +
    scale_colour_gradient(limits = c(0, max(tmp$Combined.Score) + 0.5), high = "#2b9348", low = "#eeef20") +
    xlab(NULL) + ylab(NULL) +
    ggtitle(dbs[i]) + labs(color = "Combined Score", size = "-log10(padj)")
}
tmplist <- tmplist %>% group_by(type) %>% arrange(desc(Combined.Score), .by_group = TRUE)

paths <- paths[which(paths$Term %in% c("PPAR signaling pathway", "AMPK signaling pathway", "Fatty acid biosynthesis", "Focal adhesion", "ECM-receptor interaction", "Cell adhesion molecules (CAMs)", "Adherens junction", "3T3-L1 cell", "Abdominal adipose tissue", "Subcutaneous adipose tissue", "Vasculature", "Umbilical vein", "Vascular cell", "Common cardinal vein", "Mesenchyme", "triglyceride biosynthetic process", "fatty acid transmembrane transport")), ]

paths <- paths %>%
  mutate(
    class_with_color = ifelse(class == "KEGG_2019_Mouse", glue("<strong><span style='color:#0000FF'>{Term}</span></strong>"),
                              ifelse(class == "GO_Biological_Process_2018", glue("<strong><span style='color:#9900CC'>{Term}</span></strong>"),
                                     ifelse(class == "Jensen_TISSUES", glue("<strong><span style='color:#CC0000'>{Term}</span></strong>"),
                                            glue("<strong><span style='color:'>{Term}</span></strong>"))))
  )
paths <- paths %>% group_by(type) %>% arrange(desc(Combined.Score), .by_group = TRUE)
paths$class_with_color <- factor(paths$class_with_color, levels = rev(unique(paths$class_with_color)))

ggplot(paths, aes(x = type, y = class_with_color)) +
  geom_point(aes(size = Adjusted.P.value, color = Combined.Score)) +
  theme_bw(base_size = 8) +
  scale_colour_gradient(limits = c(0, max(paths$Combined.Score) + 0.5), high = "#2b9348", low = "#eeef20") +
  xlab(NULL) + ylab(NULL) +
  ggtitle(NULL) + labs(color = "Combined Score", size = "-log10(padj)") +
  theme(axis.text.y = element_markdown(), axis.text = element_text(size = 12))
```


## Figure 2
```
## Loading R packages
library(dplyr)
library(enrichR)
library(ggplot2)
library(ggpubr)
library(ggtext)
library(glue)
library(RColorBrewer)
library(Seurat)
library(SeuratWrappers)


#################################
########### Figure 2A ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster

TSNEPlot(data, pt.size = 1, cols = c("#11c78b", "#800080", "#e57400", "#0000FF", "#dfdf0d")) +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  theme_classic() + labs(color = "Cluster") +
  theme(legend.position = "bottom")



#################################
########### Figure 2B ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster

markers <- FindAllMarkers(data, logfc.threshold = 0, only.pos = FALSE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

data <- SeuratWrappers::RunALRA(data)
data <- ScaleData(data)

cls <- c("#11c78b", "#800080", "#e57400", "#0000FF", "#dfdf0d")
mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)

DoHeatmap(data, features = top10$gene, group.colors = cls, angle = 0, size = 5, label = FALSE) +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text = element_text(size = 6))



#################################
########### Figure 2C ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster

DotPlot(data, features = c('Adipoq', 'Adrb3', 'Cidec', 'Dgat1', 'Fasn', 'Lipe', 'Pck1', 'Plin1', 'Pnpla2', 'Retn'), cols = c("grey", 'red'), assay = "RNA") +
  xlab(NULL) + ylab(NULL) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#################################
########### Figure 2D ###########
#################################
dbs <- c("KEGG_2019_Mouse", "WikiPathways_2019_Mouse", "Jensen_TISSUES", "GO_Biological_Process_2018", "Mouse_Gene_Atlas")

genes_1 <- readLines("/Users/biagi/PhD/Adipocyte/Figures/update/Fig2C_Markers_1.txt")
genes_2 <- readLines("/Users/biagi/PhD/Adipocyte/Figures/update/Fig2C_Markers_2.txt")
genes_3 <- readLines("/Users/biagi/PhD/Adipocyte/Figures/update/Fig2C_Markers_3.txt")
genes_4 <- readLines("/Users/biagi/PhD/Adipocyte/Figures/update/Fig2C_Markers_4.txt")
genes_5 <- readLines("/Users/biagi/PhD/Adipocyte/Figures/update/Fig2C_Markers_5.txt")

genes <- list(Ad1 = genes_1,
              Ad2 = genes_2,
              Ad3 = genes_3,
              Ad4 = genes_4,
              Ad5 = genes_5)

results <- lapply(genes, enrichr, dbs)

plotlist <- list()
pathlist <- NULL
tmplist <- NULL
paths <- NULL
for (i in 1:length(dbs)) {
  tmp <- lapply(results, `[[`, i)
  tmp <- mapply(cbind, tmp, "type" = names(tmp), SIMPLIFY = FALSE)
  tmp <- do.call("rbind", tmp)
  
  tmp <- subset(tmp, Adjusted.P.value < 0.05)
  
  tmp$Adjusted.P.value <- -log10(tmp$Adjusted.P.value)
  tmp$Combined.Score <- log2(tmp$Combined.Score)
  
  tmplist <- rbind(tmplist, data.frame(tmp, class = dbs[i]))
  
  if (dbs[i] == "WikiPathways_2019_Mouse" | dbs[i] == "GO_Biological_Process_2018") {
    tmp$Term <- gsub("\\s*\\([^\\)]+\\)","",as.character(tmp$Term))
  }
  tmp <- tmp %>% group_by(type) %>% top_n(n = 5, wt = Combined.Score)
  tmp$Term <- factor(tmp$Term, levels = rev(unique(tmp$Term)))
  pathlist <- rbind(pathlist, data.frame(tmp %>% group_by(type) %>% top_n(n = 2, wt = Combined.Score), class = dbs[i]))
  paths <- rbind(paths, data.frame(tmp, class = dbs[i]))
  
  plotlist[[i]] <- ggplot(tmp, aes(x = type, y = Term)) +
    geom_point(aes(size = Adjusted.P.value, color = Combined.Score)) +
    theme_bw(base_size = 9) +
    scale_colour_gradient(limits = c(0, max(tmp$Combined.Score) + 0.5), high = "#2b9348", low = "#eeef20") +
    xlab(NULL) + ylab(NULL) +
    ggtitle(dbs[i]) + labs(color = "Combined Score", size = "-log10(padj)")
}
tmplist <- tmplist %>% group_by(type) %>% arrange(desc(Combined.Score), .by_group = TRUE)

paths <- paths[which(paths$Term %in% c("PPAR signaling pathway", "Fatty acid degradation", "Glyoxylate and dicarboxylate metabolism", "Insulin resistance", "Biosynthesis of unsaturated fatty acids", "Fatty acid biosynthesis", "Adipocytokine signaling pathway", "Propanoate metabolism", "IL-17 signaling pathway", "Fat digestion and absorption", "Mitochondrial Gene Expression WP1263", "Fatty acid oxidation WP2318", "Fatty Acid Biosynthesis WP336", "TCA Cycle WP434", "Cholesterol metabolism WP4346", "Subcutaneous adipose tissue", "3T3-L1 cell", "triglyceride biosynthetic process", "triglyceride metabolic process", "long-chain fatty acid transport", "neutral amino acid transport", "regulation of sequestering of triglyceride", "cardiolipin acyl-chain remodeling", "acylglycerol acyl-chain remodeling", "adipose brown")), ]

paths <- paths %>%
  mutate(
    class_with_color = ifelse(class == "KEGG_2019_Mouse", glue("<strong><span style='color:#0000FF'>{Term}</span></strong>"),
                              ifelse(class == "WikiPathways_2019_Mouse", glue("<strong><span style='color:#f28804'>{Term}</span></strong>"),
                                     ifelse(class == "Jensen_TISSUES", glue("<strong><span style='color:#CC0000'>{Term}</span></strong>"),
                                            ifelse(class == "GO_Biological_Process_2018", glue("<strong><span style='color:#9900CC'>{Term}</span></strong>"),
                                                   glue("<strong><span style='color:#18a997'>{Term}</span></strong>")))))
  )

paths <- paths %>% group_by(type) %>% arrange(desc(Combined.Score), .by_group = TRUE)
paths$class_with_color <- factor(paths$class_with_color, levels = rev(unique(paths$class_with_color)))

ggplot(paths, aes(x = type, y = class_with_color)) +
  geom_point(aes(size = Adjusted.P.value, color = Combined.Score)) +
  theme_bw(base_size = 8) +
  scale_colour_gradient(limits = c(0, max(paths$Combined.Score) + 0.5), high = "#2b9348", low = "#eeef20") +
  xlab(NULL) + ylab(NULL) +
  ggtitle(NULL) + labs(color = "Combined Score", size = "-log10(padj)") +
  theme(axis.text.y = element_markdown(), axis.text = element_text(size = 12))



#################################
########### Figure 2E ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster
data$clusters <- Idents(data)

data <- SeuratWrappers::RunALRA(data)
data <- ScaleData(data)

genes <- c('Adipoq', 'Pnpla2', 'Fasn', 'Lep', 'Cebpa', 'Pparg', 'Cidec', 'Car3', 'Cd36', 'Gadd45g', 'Serpine1', 'Vegfa')
genes <- genes[genes %in% rownames(data)]

df <- NULL
for (i in 1:length(genes)) {
  c3 <- data.frame(value = data@assays$alra@data[genes[i], names(Idents(data))[which(Idents(data) == 'Ad3')]], cluster = 'Ad3', gene = genes[i])
  c4 <- data.frame(value = data@assays$alra@data[genes[i], names(Idents(data))[which(Idents(data) == 'Ad4')]], cluster = 'Ad4', gene = genes[i])
  df <- rbind(df, rbind(c3, c4))
}
df$cluster <- as.factor(df$cluster)
df$gene <- factor(df$gene, levels = genes)

intr <- aggregate(df[, 1], list(df$gene), mean)
colnames(intr) <- c('gene', 'Z')

ggboxplot(df, x = "cluster", y = "value",
                 color = "cluster", add = "jitter") +
  facet_wrap(~ gene) + geom_hline(data = intr, aes(yintercept = Z), linetype = 'dashed') +
  stat_compare_means(label = 'p.signif', method = 't.test', label.x = 1.5, label.y = 3.4) +
  xlab(NULL) + ylab('Normalized Expression') + theme_bw() +
  scale_color_manual(values = c('Ad3' = "#e57400", 'Ad4' = "#0000FF"))
```


## Figure 3
```
## Loading R packages
library(dplyr)
library(fgsea)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(monocle)
library(Nebulosa)
library(RColorBrewer)
library(scales)
library(Seurat)
library(SeuratWrappers)
library(viridis)

source('/Users/biagi/PhD/Adipocyte/2_Functions.R')


#################################
########### Figure 3A ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
Idents(data) <- as.factor(data$timpoint)
data <- subset(data, cells = names(which(data$timpoint == '4day' | data$timpoint == 'RT' | data$timpoint == 'CL')))
data$timpoint <- gsub('4day', 'Cold', data$timpoint)

TSNEPlot(data, pt.size = 1, group.by = 'timpoint') +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  theme_classic() + labs(color = "Cluster") +
  theme(legend.position = "bottom")



#################################
########### Figure 3B ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
Idents(data) <- as.factor(data$timpoint)
data <- subset(data, cells = names(which(data$timpoint == '4day' | data$timpoint == 'RT' | data$timpoint == 'CL')))

data <- SeuratWrappers::RunALRA(data)
data <- ScaleData(data)

genes <- c('Lipe', 'Adrb3', 'Pnpla2', 'Plin1', 'Fasn', 'Acadm', 'Gk', 'Adipoq', 'Retn', 'Cidec', 'Ppara', 'Ucp1', 'Dio2', 'Prdm16', 'Elovl3')
genes <- genes[genes %in% rownames(data)]

## Nebulosa expression plot
pt <- plot_density(data, genes, reduction = "tsne", pal = "inferno", size = 0.5, combine = FALSE)

pt <- lapply(pt, function(x){ x + theme_bw() + 
    theme(panel.background = element_rect(colour = "black", size = 0.1), 
          axis.ticks.length = unit(.2, "cm"), axis.text = element_text(size = 11), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5, face = "italic", size = 10), 
          legend.position = "none")})

ggarrange(plotlist = pt, ncol = 5, nrow = 3, common.legend = TRUE)



#################################
########### Figure 3C ###########
#################################
ort <- read.table("/projects/cangen/coliveir/cellphonedb/Orthologs_human_mouse.txt", sep = ",", header = TRUE)
pathways <- gmtPathways('/Users/biagi/PhD/Adipocyte/msigdb.v7.1.symbols.gmt')

## ColdxRT
comp1 <- read.table("/Users/biagi/PhD/Adipocyte/Figures/update/Diff_Adipocytes_ColdxRT_ALL.txt")
comp1 <- subset(comp1, p_val_adj < 0.05)
comp1 <- comp1[order(comp1$avg_logFC, decreasing = TRUE), ]
comp1 <- rbind(head(comp1, 50), tail(comp1, 50))
comp1 <- merge(ort, comp1, by.x = "Mouse.gene.name", by.y = "row.names")
comp1$Mouse.gene.name <- comp1$Gene.stable.ID <- comp1$Mouse.gene.stable.ID <- NULL
colnames(comp1)[1] <- "Gene"
ranks <- structure(comp1$avg_logFC, names = comp1$Gene)
ranks <- sort(ranks)
fgseaRes1 <- fgsea(pathways, ranks, minSize = 5)
fgseaRes1 <- subset(fgseaRes1, padj < 0.05)
fgseaRes1 <- fgseaRes1[order(fgseaRes1$NES, decreasing = TRUE), ]
genes1 <- fgseaRes1$leadingEdge
for (i in 1:length(genes1)) {
  fgseaRes1$genes[i] <- paste(genes1[[i]], collapse = ';')
}
fgseaRes1$leadingEdge <- NULL
enrich_ColdxRT <- read.table('/Users/biagi/PhD/Adipocyte/Figures/update/fgsea_ColdxRT.txt', header = TRUE)
pathways_ColdxRT <- c('GO_ORGANOPHOSPHATE_BIOSYNTHETIC_PROCESS', 'GO_FATTY_ACID_TRANSPORT', 'GO_LONG_CHAIN_FATTY_ACID_TRANSPORT', 'GO_MONOCARBOXYLIC_ACID_TRANSPORT', 'GO_RESPONSE_TO_FATTY_ACID', 'REACTOME_NGF_STIMULATED_TRANSCRIPTION', 'REACTOME_NUCLEAR_EVENTS_KINASE_AND_TRANSCRIPTION_FACTOR_ACTIVATION', 'OSWALD_HEMATOPOIETIC_STEM_CELL_IN_COLLAGEN_GEL_UP', 'NAGASHIMA_EGF_SIGNALING_UP', 'NAGASHIMA_NRG1_SIGNALING_UP')
enrich_ColdxRT <- enrich_ColdxRT[which(enrich_ColdxRT$pathway %in% pathways_ColdxRT), ]
enrich_ColdxRT$class <- ifelse(enrich_ColdxRT$NES > 0, 'Up-regulated', 'Down-regulated')
enrich_ColdxRT <- enrich_ColdxRT[order(enrich_ColdxRT$class, enrich_ColdxRT$NES), ]
enrich_ColdxRT$order <- seq_len(nrow(enrich_ColdxRT))

ggplot(enrich_ColdxRT, aes(order, NES, fill = class)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  xlab(NULL) + ylab("Normalized Enrichment Score (NES)") +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "right", legend.title = element_text(size = 10), plot.title = element_text(hjust = 0.5)) +
  scale_fill_discrete(name = "Enrichment") +
  scale_x_continuous(breaks = enrich_ColdxRT$order, labels = enrich_ColdxRT$pathway, expand = c(0, 0)) +
  ggtitle('GSEA for Cold x RT')

## CLxRT
comp2 <- read.table("/Users/biagi/PhD/Adipocyte/Figures/update/Diff_Adipocytes_CLxRT_ALL.txt")
comp2 <- subset(comp2, p_val_adj < 0.05)
comp2 <- comp2[order(comp2$avg_logFC, decreasing = TRUE), ]
comp2 <- rbind(head(comp2, 50), tail(comp2, 50))
comp2 <- merge(ort, comp2, by.x = "Mouse.gene.name", by.y = "row.names")
comp2$Mouse.gene.name <- comp2$Gene.stable.ID <- comp2$Mouse.gene.stable.ID <- NULL
colnames(comp2)[1] <- "Gene"
ranks <- structure(comp2$avg_logFC, names = comp2$Gene)
ranks <- sort(ranks)
fgseaRes2 <- fgsea(pathways, ranks, minSize = 5)
fgseaRes2 <- subset(fgseaRes2, padj < 0.05)
fgseaRes2 <- fgseaRes2[order(fgseaRes2$NES, decreasing = TRUE), ]
genes2 <- fgseaRes2$leadingEdge
for (i in 1:length(genes2)) {
  fgseaRes2$genes[i] <- paste(genes2[[i]], collapse = ';')
}
fgseaRes2$leadingEdge <- NULL
enrich_CLxRT <- read.table('/Users/biagi/PhD/Adipocyte/Figures/update/fgsea_CLxRT.txt', header = TRUE)
pathways_CLxRT <- c('GO_ORGANOPHOSPHATE_METABOLIC_PROCESS', 'GO_MITOCHONDRIAL_ENVELOPE', 'GO_ATP_METABOLIC_PROCESS', 'GO_OXIDOREDUCTASE_ACTIVITY', 'GO_OXIDATIVE_PHOSPHORYLATION', 'GO_NEGATIVE_REGULATION_OF_PROTEIN_METABOLIC_PROCESS', 'GO_RESPONSE_TO_MECHANICAL_STIMULUS', 'GO_NEGATIVE_REGULATION_OF_BIOSYNTHETIC_PROCESS', 'NAGASHIMA_EGF_SIGNALING_UP', 'NAGASHIMA_NRG1_SIGNALING_UP')
enrich_CLxRT <- enrich_CLxRT[which(enrich_CLxRT$pathway %in% pathways_CLxRT), ]
enrich_CLxRT$class <- ifelse(enrich_CLxRT$NES > 0, 'Up-regulated', 'Down-regulated')
enrich_CLxRT <- enrich_CLxRT[order(enrich_CLxRT$class, enrich_CLxRT$NES), ]
enrich_CLxRT$order <- seq_len(nrow(enrich_CLxRT))

ggplot(enrich_CLxRT, aes(order, NES, fill = class)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  xlab(NULL) + ylab("Normalized Enrichment Score (NES)") +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "right", legend.title = element_text(size = 10), plot.title = element_text(hjust = 0.5)) +
  scale_fill_discrete(name = "Enrichment") +
  scale_x_continuous(breaks = enrich_CLxRT$order, labels = enrich_CLxRT$pathway, expand = c(0, 0)) +
  ggtitle('GSEA for CL x RT')



#################################
########### Figure 3D ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
data$timpoint <- gsub('4day', 'Cold', data$timpoint)

infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster
data$clusters <- Idents(data)

data <- SeuratWrappers::RunALRA(data)
data <- ScaleData(data)

anno_CL <- subset(data@meta.data, timpoint == 'CL')
anno_RT <- subset(data@meta.data, timpoint == 'RT')
anno_24hr <- subset(data@meta.data, timpoint == '24hr')
anno_48hr <- subset(data@meta.data, timpoint == '48hr')
anno_Cold <- subset(data@meta.data, timpoint == 'Cold')

tsne <- as.data.frame(Embeddings(data, 'tsne'))
meta <- data@meta.data
df <- merge(tsne, meta, by = 'row.names')
rownames(df) <- df$Row.names
df$Row.names <- NULL

## CL
ggplot() +
  geom_point(data = df, aes(x = tSNE_1, y = tSNE_2, color = timpoint), size = 0.5) +
  scale_color_manual(values = rep("gray75", 5)) +
  guides(color = guide_legend(override.aes = list(size = 7))) +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  labs(color = '') + ggtitle(paste0('CL (', nrow(anno_CL), ' cells)')) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", plot.title = element_text(hjust = 0.5, size = 8)) +
  annotate("point", df[rownames(anno_CL)[anno_CL$clusters == 'Ad1'], 'tSNE_1'],
           df[rownames(anno_CL)[anno_CL$clusters == 'Ad1'], 'tSNE_2'], colour = "#11c78b", size = 1) +
  annotate("point", df[rownames(anno_CL)[anno_CL$clusters == 'Ad2'], 'tSNE_1'],
           df[rownames(anno_CL)[anno_CL$clusters == 'Ad2'], 'tSNE_2'], colour = "#800080", size = 1) +
  annotate("point", df[rownames(anno_CL)[anno_CL$clusters == 'Ad3'], 'tSNE_1'],
           df[rownames(anno_CL)[anno_CL$clusters == 'Ad3'], 'tSNE_2'], colour = "#e57400", size = 1) +
  annotate("point", df[rownames(anno_CL)[anno_CL$clusters == 'Ad4'], 'tSNE_1'],
           df[rownames(anno_CL)[anno_CL$clusters == 'Ad4'], 'tSNE_2'], colour = "#0000FF", size = 1) +
  annotate("point", df[rownames(anno_CL)[anno_CL$clusters == 'Ad5'], 'tSNE_1'],
           df[rownames(anno_CL)[anno_CL$clusters == 'Ad5'], 'tSNE_2'], colour = "#dfdf0d", size = 1)

df2 <- data.frame(
  class = c("Ad1", "Ad2", "Ad3", "Ad4", 'Ad5'),
  n = c(sum(anno_CL == "Ad1"), sum(anno_CL == "Ad2"), sum(anno_CL == "Ad3"), sum(anno_CL == "Ad4"), sum(anno_CL == "Ad5")),
  value = c(round((sum(anno_CL == "Ad1")/nrow(anno_CL))*100, 2), round((sum(anno_CL == "Ad2")/nrow(anno_CL))*100, 2), round((sum(anno_CL == "Ad3")/nrow(anno_CL))*100, 2), round((sum(anno_CL == "Ad4")/nrow(anno_CL))*100, 2), round((sum(anno_CL == "Ad5")/nrow(anno_CL))*100, 2))
)
df2$class <- factor(df2$class, levels = c("Ad1", "Ad2", "Ad3", "Ad4", "Ad5"))

df2 <- df2 %>%
  arrange(desc(class)) %>%
  mutate(text_y = cumsum(value) - value/2)
df2$pos = (cumsum(c(0, df2$value)) + c(df2$value / 2, .01))[1:nrow(df2)]

ggplot(df2, aes(x = "", y = value, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = c("#11c78b", "#800080", "#e57400", "#0000FF", "#dfdf0d")) +
  theme_void() + labs(fill = "Clusters") +
  geom_text_repel(aes(x = 1.5, y = pos, label = paste0(value, "%")), nudge_x = .1,
                  segment.size = .7, show.legend = FALSE)

## RT
ggplot() +
  geom_point(data = df, aes(x = tSNE_1, y = tSNE_2, color = timpoint), size = 0.5) +
  scale_color_manual(values = rep("gray75", 5)) +
  guides(color = guide_legend(override.aes = list(size = 7))) +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  labs(color = '') + ggtitle(paste0('RT (', nrow(anno_RT), ' cells)')) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", plot.title = element_text(hjust = 0.5, size = 8)) +
  annotate("point", df[rownames(anno_RT)[anno_RT$clusters == 'Ad1'], 'tSNE_1'],
           df[rownames(anno_RT)[anno_RT$clusters == 'Ad1'], 'tSNE_2'], colour = "#11c78b", size = 1) +
  annotate("point", df[rownames(anno_RT)[anno_RT$clusters == 'Ad2'], 'tSNE_1'],
           df[rownames(anno_RT)[anno_RT$clusters == 'Ad2'], 'tSNE_2'], colour = "#800080", size = 1) +
  annotate("point", df[rownames(anno_RT)[anno_RT$clusters == 'Ad3'], 'tSNE_1'],
           df[rownames(anno_RT)[anno_RT$clusters == 'Ad3'], 'tSNE_2'], colour = "#e57400", size = 1) +
  annotate("point", df[rownames(anno_RT)[anno_RT$clusters == 'Ad4'], 'tSNE_1'],
           df[rownames(anno_RT)[anno_RT$clusters == 'Ad4'], 'tSNE_2'], colour = "#0000FF", size = 1) +
  annotate("point", df[rownames(anno_RT)[anno_RT$clusters == 'Ad5'], 'tSNE_1'],
           df[rownames(anno_RT)[anno_RT$clusters == 'Ad5'], 'tSNE_2'], colour = "#dfdf0d", size = 1)

df2 <- data.frame(
  class = c("Ad1", "Ad2", "Ad3", "Ad4", 'Ad5'),
  n = c(sum(anno_RT == "Ad1"), sum(anno_RT == "Ad2"), sum(anno_RT == "Ad3"), sum(anno_RT == "Ad4"), sum(anno_RT == "Ad5")),
  value = c(round((sum(anno_RT == "Ad1")/nrow(anno_RT))*100, 2), round((sum(anno_RT == "Ad2")/nrow(anno_RT))*100, 2), round((sum(anno_RT == "Ad3")/nrow(anno_RT))*100, 2), round((sum(anno_RT == "Ad4")/nrow(anno_RT))*100, 2), round((sum(anno_RT == "Ad5")/nrow(anno_RT))*100, 2))
)
df2$class <- factor(df2$class, levels = c("Ad1", "Ad2", "Ad3", "Ad4", "Ad5"))

df2 <- df2 %>%
  arrange(desc(class)) %>%
  mutate(text_y = cumsum(value) - value/2)
df2$pos = (cumsum(c(0, df2$value)) + c(df2$value / 2, .01))[1:nrow(df2)]

ggplot(df2, aes(x = "", y = value, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = c("#11c78b", "#800080", "#e57400", "#0000FF", "#dfdf0d")) +
  theme_void() + labs(fill = "Clusters") +
  geom_text_repel(aes(x = 1.5, y = pos, label = paste0(value, "%")), nudge_x = .1,
                  segment.size = .7, show.legend = FALSE)

## Cold
ggplot() +
  geom_point(data = df, aes(x = tSNE_1, y = tSNE_2, color = timpoint), size = 0.5) +
  scale_color_manual(values = rep("gray75", 5)) +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  guides(color = guide_legend(override.aes = list(size = 7))) +
  labs(color = '') + ggtitle(paste0('Cold (', nrow(anno_Cold), ' cells)')) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", plot.title = element_text(hjust = 0.5, size = 8)) +
  annotate("point", df[rownames(anno_Cold)[anno_Cold$clusters == 'Ad1'], 'tSNE_1'],
           df[rownames(anno_Cold)[anno_Cold$clusters == 'Ad1'], 'tSNE_2'], colour = "#11c78b", size = 1) +
  annotate("point", df[rownames(anno_Cold)[anno_Cold$clusters == 'Ad2'], 'tSNE_1'],
           df[rownames(anno_Cold)[anno_Cold$clusters == 'Ad2'], 'tSNE_2'], colour = "#800080", size = 1) +
  annotate("point", df[rownames(anno_Cold)[anno_Cold$clusters == 'Ad3'], 'tSNE_1'],
           df[rownames(anno_Cold)[anno_Cold$clusters == 'Ad3'], 'tSNE_2'], colour = "#e57400", size = 1) +
  annotate("point", df[rownames(anno_Cold)[anno_Cold$clusters == 'Ad4'], 'tSNE_1'],
           df[rownames(anno_Cold)[anno_Cold$clusters == 'Ad4'], 'tSNE_2'], colour = "#0000FF", size = 1) +
  annotate("point", df[rownames(anno_Cold)[anno_Cold$clusters == 'Ad5'], 'tSNE_1'],
           df[rownames(anno_Cold)[anno_Cold$clusters == 'Ad5'], 'tSNE_2'], colour = "#dfdf0d", size = 1)

df2 <- data.frame(
  class = c("Ad1", "Ad2", "Ad3", "Ad4", 'Ad5'),
  n = c(sum(anno_Cold == "Ad1"), sum(anno_Cold == "Ad2"), sum(anno_Cold == "Ad3"), sum(anno_Cold == "Ad4"), sum(anno_Cold == "Ad5")),
  value = c(round((sum(anno_Cold == "Ad1")/nrow(anno_Cold))*100, 2), round((sum(anno_Cold == "Ad2")/nrow(anno_Cold))*100, 2), round((sum(anno_Cold == "Ad3")/nrow(anno_Cold))*100, 2), round((sum(anno_Cold == "Ad4")/nrow(anno_Cold))*100, 2), round((sum(anno_Cold == "Ad5")/nrow(anno_Cold))*100, 2))
)
df2$class <- factor(df2$class, levels = c("Ad1", "Ad2", "Ad3", "Ad4", "Ad5"))

df2 <- df2 %>%
  arrange(desc(class)) %>%
  mutate(text_y = cumsum(value) - value/2)
df2$pos = (cumsum(c(0, df2$value)) + c(df2$value / 2, .01))[1:nrow(df2)]

ggplot(df2, aes(x = "", y = value, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = c("#11c78b", "#800080", "#e57400", "#0000FF", "#dfdf0d")) +
  theme_void() + labs(fill = "Clusters") +
  geom_text_repel(aes(x = 1.5, y = pos, label = paste0(value, "%")), nudge_x = .1,
                  segment.size = .7, show.legend = FALSE)



#################################
########### Figure 3E ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- factor(data$timpoint)
data$cluster <- factor(new_cluster)

data <-SeuratWrappers:: RunALRA(data)
data <- ScaleData(data, features = rownames(data))

genes <- c('Acaa2', 'Ehhadh', 'Slc27a2', 'Acadm', 'Scp2', 'Acadvl', 'Hadha', 'Hadhb', 'Acaa1b', 'Abcd3', 
           'Sdhaf2', 'Pdha1', 'Pdhb', 'Idh3g', 'Mdh1', 'Sdha', 'Ndufs4', 'Sdhc', 'Suclg1', 'Mdh2', 
           'Akt1', 'Kiss1r', 'Mif', 'Cpt1b', 'Fabp3', 'Acsl5', 'Mfsd2a', 'Abcc4', 'Acsl3', 'Pla2g12a', 
           'Aldoa', 'Gpi1', 'Pgk1', 'Gapdh', 'Pfkp', 'Pkm', 'Ogdh', 'Ppara', 'Prkaa1', 'Gale', 'Entpd5', 'Stat3', 'Foxk2', 'Igf1', 
           "Daglb", "Slc2a4", "Gpd1", "Pcx", "Abhd5", "Adrb3", "Lipe", "Pck1", "Pnpla2", 
           'Acaca', 'Acly', 'Fasn', 'Srebf1', 'Dgat2', 'Insig1', 'Hadh', 'Elovl6', 'Lpl')
genes <- genes[genes %in% rownames(data)]

mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)
DoHeatmap(data, features = genes, angle = 0, size = 5, label = FALSE, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text = element_text(size = 6)) +
  labs(color = 'Ucp1 Expression')



#################################
########### Figure 3G ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster
data$clusters <- Idents(data)

data <- SeuratWrappers::RunALRA(data)
data <- ScaleData(data)

genes <- c('Acadm', 'Adipoq', 'Adrb3', 'Lep', 'Fasn', 
           'Ucp1', 'Cidea', 'Dio2', 'Elovl3', 'Cox8b', 
           'Cpt1b', 'Gk', 'Plin5', 'Pnpla2', 'Lpl', 
           'Acaca', 'Gadd45g', 'Fabp3', 'Ppara', 'Slc2a4')
genes <- genes[genes %in% rownames(data)]

df <- NULL
for (i in 1:length(genes)) {
  c1 <- data.frame(value = data@assays$alra@data[genes[i], names(Idents(data))[which(Idents(data) == 'Ad1')]], cluster = 'Ad1', gene = genes[i])
  c2 <- data.frame(value = data@assays$alra@data[genes[i], names(Idents(data))[which(Idents(data) == 'Ad2')]], cluster = 'Ad2', gene = genes[i])
  c3 <- data.frame(value = data@assays$alra@data[genes[i], names(Idents(data))[which(Idents(data) == 'Ad3')]], cluster = 'Ad3', gene = genes[i])
  c4 <- data.frame(value = data@assays$alra@data[genes[i], names(Idents(data))[which(Idents(data) == 'Ad4')]], cluster = 'Ad4', gene = genes[i])
  c5 <- data.frame(value = data@assays$alra@data[genes[i], names(Idents(data))[which(Idents(data) == 'Ad5')]], cluster = 'Ad5', gene = genes[i])
  
  df <- rbind(df, rbind(c1, c2, c3, c4, c5))
}
df$cluster <- as.factor(df$cluster)
df$gene <- factor(df$gene, levels = genes)

intr <- aggregate(df[, 1], list(df$gene), mean)
colnames(intr) <- c('gene', 'Z')

ggboxplot(df, x = "cluster", y = "value",
                 color = "cluster", add = "jitter") +
  facet_wrap(~ gene) + geom_hline(data = intr, aes(yintercept = Z), linetype = 'dashed') +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", label.x = 1.5, label.y = 3.5) +
  xlab(NULL) + ylab('Normalized Expression') + theme_bw() +
  scale_color_manual(values = c('Ad1' = "#11c78b", 'Ad2' = "#800080", 'Ad3' = "#e57400", 'Ad4' = "#0000FF", 'Ad5' = "#dfdf0d"))



#################################
########### Figure 3H ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
data$timpoint <- gsub('4day', 'Cold', data$timpoint)

infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")
new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster
data$cluster <- Idents(data)

Idents(data) <- as.factor(data$timpoint)
data <- subset(data, cells = names(which(data$timpoint == 'CL' | data$timpoint == 'RT')))

DefaultAssay(data) <- "RNA"
mat <- data@assays$RNA@data

fdata <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
pdata <- data@meta.data

pd <- new("AnnotatedDataFrame", data = pdata)
fd <- new("AnnotatedDataFrame", data = fdata)
cds <- newCellDataSet(as(mat, "dgCMatrix"),
                      phenoData = pd,
                      featureData = fd, 
                      expressionFamily = negbinomial.size())

rm(data, mat, fdata, pdata, fd, pd)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

diff_test_res <- differentialGeneTest(cds[expressed_genes, ], fullModelFormulaStr = "~cluster", cores = 6, verbose = TRUE)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree', verbose = TRUE)
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 3)

plot_cell_trajectory(cds, color_by = "Pseudotime")

plot_cell_trajectory(cds, color_by = "cluster") + 
  scale_color_manual(values = c("Ad1" = "#11c78b", "Ad2" = "#800080", "Ad3" = "#e57400", "Ad4" = "#0000FF", "Ad5" = "#dfdf0d"))

plotMonocle(cds, c('Ppara', 'Prdm16', 'Ppargc1a', 'Egr2', 'Dbp', 'Xbp1', 'Clock', 'Zbtb43', 'Zbtb7a'))



#################################
########### Figure 3I ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
data$timpoint <- gsub('4day', 'Cold', data$timpoint)

infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")
new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster
data$cluster <- Idents(data)

Idents(data) <- as.factor(data$timpoint)
data <- subset(data, cells = names(which(data$timpoint == 'Cold' | data$timpoint == 'RT')))

DefaultAssay(data) <- "RNA"
mat <- data@assays$RNA@data

fdata <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
pdata <- data@meta.data

pd <- new("AnnotatedDataFrame", data = pdata)
fd <- new("AnnotatedDataFrame", data = fdata)
cds <- newCellDataSet(as(mat, "dgCMatrix"),
                      phenoData = pd,
                      featureData = fd, 
                      expressionFamily = negbinomial.size())

rm(data, mat, fdata, pdata, fd, pd)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

diff_test_res <- differentialGeneTest(cds[expressed_genes, ], fullModelFormulaStr = "~cluster", cores = 6, verbose = TRUE)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree', verbose = TRUE)
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 4)

plot_cell_trajectory(cds, color_by = "Pseudotime")

plot_cell_trajectory(cds, color_by = "cluster") + 
  scale_color_manual(values = c("Ad1" = "#11c78b", "Ad2" = "#800080", "Ad3" = "#e57400", "Ad4" = "#0000FF", "Ad5" = "#dfdf0d"))

plotMonocle(cds, c('Ppara', 'Prdm16', 'Ppargc1a', 'Srebf1', 'Cebpa', 'Egr1', 'Foxn3', 'Esrra', 'Gtf2ird1'))
```


## Figure 4
```
## Loading R packages
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(fgsea)
library(ggplot2)
library(ggrepel)
library(Nebulosa)
library(RColorBrewer)
library(Seurat)
library(SeuratWrappers)


#################################
########### Figure 4A ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")

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
  guides(color = guide_legend(override.aes = list(size = 7))) + 
  labs(color = '', x = 't-SNE1', y = 't-SNE2') +
  theme_classic() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 8)) + 
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

## pieplot UCP1 High
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")

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
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")

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

markers <- FindAllMarkers(data, logfc.threshold = 0, only.pos = FALSE)

markers <- subset(markers, p_val_adj < 0.05 & cluster == 'High')
markers <- markers[order(markers$avg_logFC, decreasing = TRUE), ]

genes1 <- c(head(markers$gene, 10), tail(markers$gene, 10))

mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)

DoHeatmap(data, features = genes1, angle = 0, size = 5, label = FALSE) +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text = element_text(size = 6)) +
  labs(color = 'UCP1 Expression')



#################################
########### Figure 4C ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")

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

markers <- FindAllMarkers(data, logfc.threshold = 0, only.pos = FALSE)

markers <- subset(markers, p_val_adj < 0.05 & cluster == 'High')
markers <- markers[order(markers$avg_logFC, decreasing = TRUE), ]

ort <- read.table("/Users/biagi/PhD/Adipocyte/Orthologs_human_mouse.txt", sep = ",", header = TRUE)
pathways <- gmtPathways('/Users/biagi/PhD/Adipocyte/msigdb.v7.1.symbols.gmt')

comp1 <- rbind(head(markers, 50), tail(markers, 35))

comp1 <- merge(ort, comp1, by.x = "Mouse.gene.name", by.y = "gene")
comp1$Mouse.gene.name <- comp1$Gene.stable.ID <- comp1$Mouse.gene.stable.ID <- NULL
colnames(comp1)[1] <- "Gene"
ranks <- structure(comp1$avg_logFC, names = comp1$Gene)
ranks <- sort(ranks)
fgseaRes1 <- fgsea(pathways, ranks, minSize = 5)
fgseaRes1 <- subset(fgseaRes1, pval < 0.05)
fgseaRes1 <- fgseaRes1[order(fgseaRes1$NES, decreasing = TRUE), ]
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
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")
new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster
data$cluster <- Idents(data)
data <- subset(data, idents = list('Ad1'))
tmp1 <- factor(ifelse(data@assays$SCT@data["Ucp1", ] > 0, "Ad1_High", "Ad1_Low"))


data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")
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
           "Daglb", "Slc2a4", "Gpd1", "Pcx", "Abhd5", "Adrb3", "Lipe", "Pck1", "Pnpla2", 
           'Acaca', 'Acly', 'Fasn', 'Srebf1', 'Dgat2', 'Insig1', 'Hadh', 'Elovl6', 'Lpl')
genes <- genes[genes %in% rownames(data)]

mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)

DoHeatmap(data, features = genes, angle = 0, size = 5, label = FALSE) +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text = element_text(size = 6)) +
  labs(color = 'Ucp1 Expression')



#################################
########### Figure 4F ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")

ort <- read.table("/Users/biagi/PhD/Adipocyte/Orthologs_human_mouse.txt", sep = ",", header = TRUE)

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
m = apply(x, 1, mean, na.rm = TRUE)
s = apply(x, 1, sd, na.rm = TRUE)
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
        show_row_names = TRUE, show_column_names = FALSE,
        cluster_rows = FALSE, cluster_columns = TRUE,
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
m = apply(x, 1, mean, na.rm = TRUE)
s = apply(x, 1, sd, na.rm = TRUE)
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
        show_row_names = TRUE, show_column_names = FALSE,
        cluster_rows = FALSE, cluster_columns = TRUE,
        col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)), 
        heatmap_height = unit(6, "cm"), row_names_gp = gpar(fontsize = 8))

        
        
#################################
########### Figure 4G ###########
#################################
data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
data$clusters <- new_cluster

Idents(data) <- as.factor(data$timpoint)
data <- subset(data, cells = names(which(data$timpoint == '4day' | data$timpoint == 'RT' | data$timpoint == 'CL')))

pt <- plot_density(data, c('Slc36a2', 'Acadm', 'Slc27a1', 'Fasn'), reduction = "tsne", pal = "inferno", combine = FALSE)

pt <- lapply(pt, function(x){ x + theme_bw() + 
               theme(panel.background = element_rect(colour = "black", size = 0.1), 
                     axis.ticks.length = unit(.2, "cm"), 
                     axis.text = element_text(size = 11), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())})

ggarrange(plotlist = pt)
```