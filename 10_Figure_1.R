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