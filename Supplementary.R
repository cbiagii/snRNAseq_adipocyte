### Loading librarys
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(metacell)
library(pheatmap)
library(enrichR)
library(forcats)
library(ggtext)
library(RColorBrewer)
library(monocle)
library(slingshot)
library(SingleCellExperiment)
library(ggbeeswarm)
library(tradeSeq)
library(tidyr)
library(tibble)
library(pbapply)
library(gam)
library(tidyverse)


##################################
########### Figure S1A ###########
##################################
load("/projects/cangen/coliveir/Miguel/output/10x/metacell/db/mc.test_mc_f.Rda")
lfp = log2(object@mc_fp)

marks_colors <- NULL
marks_colors <- rbind(marks_colors, c("Adipocyte_1", "Acsl1", "#0000b3", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Adipocyte_2", "Plin4", "#0000cc", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Adipocyte_3", "Mlxipl", "#0000e6", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Adipocyte_4", "Pck1", "#0000ff", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Adipocyte_5", "Adrb3", "#1a1aff", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Endothelial_1", "Btnl9", "#00cd00", 1, 1))
marks_colors <- rbind(marks_colors, c("Endothelial_2", "Flt1", "#00b300", 1, 1))
marks_colors <- rbind(marks_colors, c("Endothelial_3", "Kdr", "#009a00", 1, 1))
marks_colors <- rbind(marks_colors, c("Endothelial_4", "Cdh13", "#008000", 1, 1))
marks_colors <- rbind(marks_colors, c("Endothelial_5", "Cyyr1", "#006700", 1, 1))
marks_colors <- rbind(marks_colors, c("Immune_1", "Zeb2", "#ff7f7f", 1, 1))
marks_colors <- rbind(marks_colors, c("Immune_2", "Trps1", "#ff6666", 1, 1))
marks_colors <- rbind(marks_colors, c("Immune_3", "Runx1", "#ff4c4c", 1, 1))
marks_colors <- rbind(marks_colors, c("Immune_4", "Ptprc", "#ff3232", 1, 1))
marks_colors <- rbind(marks_colors, c("Immune_5", "Adap2", "#ff1919", 1, 1))
marks_colors <- rbind(marks_colors, c("Progenitor_1", "Dcn", "#ffff4d", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Progenitor_2", "Celf2", "#ffff33", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Progenitor_3", "Meg3", "#ffff1a", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Progenitor_4", "Col1a2", "#ffff00", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Progenitor_5", "Col3a1", "#e6e600", 1, 2.5))
marks_colors <- as.data.frame(marks_colors)
colnames(marks_colors) <- c("group", "gene", "color", "priority", "T_fold")
marks_colors$priority <- as.integer(marks_colors$priority)
marks_colors$T_fold <- as.numeric(marks_colors$T_fold)


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

A <- subset(final, cellType == "Adipocytes")
E <- subset(final, cellType == "Endothelials")
I <- subset(final, cellType == "Immunes")
P <- subset(final, cellType == "Progenitors")

new.cluster.ids <- c("PG", "PG", "PG", "EN", "IM", "PG", "PG", "IM", "AD", "IM", "EN", "IM", "IM", "IM", "EN", "EN", "PG")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
Idents(data) <- factor(Idents(data), levels = c("AD", "EN", "IM", "PG"))

cls <- c("#FFA500", "#329932", "#ff9999", "#6666ff")

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S1A_1.eps", width = 7, height = 6)
TSNEPlot(data, label = FALSE, pt.size = 0.3, cols = cls) +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  theme_classic() + labs(color = "Cluster") +
  theme(legend.position="bottom")
dev.off()

a <- round((sum(Idents(data) == "AD")/ncol(data))*100, 2); a1 <- sum(Idents(data) == "AD")
b <- round((sum(Idents(data) == "EN")/ncol(data))*100, 2); b1 <- sum(Idents(data) == "EN")
c <- round((sum(Idents(data) == "IM")/ncol(data))*100, 2); c1 <- sum(Idents(data) == "IM")
d <- round((sum(Idents(data) == "PG")/ncol(data))*100, 2); d1 <- sum(Idents(data) == "PG")

df <- data.frame(
  class = c("AD", "EN", "IM", "PG"),
  n = c(a1, b1, c1, d1),
  value = c(a, b, c, d)
)

df$class <- factor(df$class, levels = c("AD", "EN", "IM", "PG"))

df <- df %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(value) - 0.5*value)

pt <- ggplot(df, aes(x = "", y = value, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = paste0(value, "%")), color = "white", fontface = "bold")+
  scale_fill_manual(values = c("#FFA500", "#329932", "#ff9999", "#6666ff")) +
  theme_void()

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S1A_2.eps", width = 7, height = 7)
print(pt)
dev.off()



##################################
########### Figure S1B ###########
##################################
data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/10x_Processed.rds")
Idents(data) <- factor(data$CellTypeRefined)

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S1B.eps", width = 7, height = 6)
TSNEPlot(data, label = FALSE, pt.size = 0.3) +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  theme_classic() + labs(color = "Cluster") +
  theme(legend.position="bottom")
dev.off()



##################################
########### Figure S1C ###########
##################################
filelist = list.files(path = "/projects/cangen/coliveir/scRNA_output/SCCAF/Adipocytes",
                      pattern = "sccaf_assess", recursive = T, full.names = T)
fnames <- gsub("sccaf_assess_", "", basename(filelist))
fnames <- gsub(".txt", "", fnames)

datalist <- lapply(filelist, function(x)read.csv(x))
names(datalist) <- fnames
datafr <- do.call("rbind", datalist)

pt <- ggplot(datafr, aes(x = Round, y = Accuracy, fill = Type)) +
  geom_violin(
    trim = FALSE,
    position = position_dodge(0.9)
  ) +
  scale_fill_manual(values = c("#d20c0c", "#1111C9")) +
  geom_boxplot(
    width = 0.05,
    position = position_dodge(0.9), colour = "white", outlier.colour = NA, show.legend = FALSE
  ) +
  xlab("") + ggtitle("Test and cross validation accurancies per SCCAF round.") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.9, 0.15), legend.title = element_blank())

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S1C.eps", width = 6, height = 6)
print(pt)
dev.off()



##################################
########### Figure S1D ###########
##################################
import warnings
warnings.filterwarnings("ignore")
from SCCAF import *
  
ad = sc.read("/projects/cangen/coliveir/scRNA_output/SCCAF/Adipocytes/results.h5ad")

y_prob, y_pred, y_test, clf, cvsm, acc = SCCAF_assessment(ad.X, ad.obs['L1_result'],n_jobs=8)
aucs = plot_roc(y_prob, y_test, clf, cvsm=cvsm, acc=acc)
plt.savefig('/projects/cangen/coliveir/Miguel/paper/Figure_SD.eps')



##################################
########### Figure S1E ###########
##################################
data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/10x_SCT_Processed_ALRA.rds")

pt <- FeaturePlot(data, c("Adrb3", "Pecam1", "Ptprc", "Cd34", "Pdgfra", "Itgb1"),
                  cols = c("grey", 'red'), reduction = 'tsne', pt.size = 0.1, combine = F)
pt <- lapply(pt, function(x) {
  x + theme_classic() + xlab("t-SNE 1") + ylab("t-SNE 2") +
    theme(plot.title = element_text(hjust = 0.5, face = "italic"))
})

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S1D.eps", width = 9, height = 6)
ggarrange(plotlist = pt, common.legend = T)
dev.off()



##################################
########### Figure S1F ###########
##################################
scdb_init("/Users/biagi/cangen/coliveir/Miguel/output/10x/metacell_SCT2/db", force_reinit=T)
scfigs_init("/Users/biagi/cangen/coliveir/Miguel/output/10x/metacell_SCT2/figs")

load("/Users/biagi/cangen/coliveir/Miguel/output/10x/metacell_SCT2/db/mc.test_mc_f.Rda")
lfp = log2(object@mc_fp)

mc = scdb_mc("test_mc_f")
gset = scdb_gset("test_markers")

gene_folds = mc@mc_fp

good_marks = intersect(names(gset@gene_set), rownames(mc@mc_fp))
mc_ord = 1:ncol(mc@mc_fp)

mat = log2(gene_folds[good_marks, mc_ord])
mat = pmax(pmin(mat, 3), -3)


mat_A <- mat[, which(mc@colors == "blue")]
mat_A <- mat_A[rowSums(mat_A) > quantile(rowSums(mat_A), 0.9), ]

mat_E <- mat[, which(mc@colors == "green")]
mat_E <- mat_E[rowSums(mat_E) > quantile(rowSums(mat_E), 0.9), ]

mat_I <- mat[, which(mc@colors %in% c("#ff748c", "#ff8fa3"))]
mat_I <- mat_I[rowSums(mat_I) > quantile(rowSums(mat_I), 0.9), ]

mat_P <- mat[, which(mc@colors %in% c("#ffa500", "#ffb732"))]
mat_P <- mat_P[rowSums(mat_P) > quantile(rowSums(mat_P), 0.9), ]

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S1E_1.eps", width = 8, height = 6)
pheatmap(mat_A, fontsize = 8, main = 'Adipocytes', legend = TRUE, treeheight_row = 0, treeheight_col = 0)
dev.off()

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S1E_2.eps", width = 8, height = 6)
pheatmap(mat_E, fontsize = 8, main = 'Endothelials', legend = TRUE, treeheight_row = 0, treeheight_col = 0)
dev.off()

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S1E_3.eps", width = 8, height = 6)
pheatmap(mat_I, fontsize = 8, main = 'Immunes', legend = TRUE, treeheight_row = 0, treeheight_col = 0)
dev.off()

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S1E_4.eps", width = 8, height = 6)
pheatmap(mat_P, fontsize = 8, main = 'Progenitors', legend = TRUE, treeheight_row = 0, treeheight_col = 0)
dev.off()



##################################
########### Figure S1G ###########
##################################
dbs <- listEnrichrDbs()
dbs <- c("KEGG_2019_Mouse", "WikiPathways_2019_Mouse", "Jensen_TISSUES", "GO_Biological_Process_2018")

genes_A <- readLines("/Users/biagi/cangen/coliveir/Miguel/Figures/update/Fig1D_2_Markers_A.txt")
genes_E <- readLines("/Users/biagi/cangen/coliveir/Miguel/Figures/update/Fig1D_2_Markers_E.txt")
genes_I <- readLines("/Users/biagi/cangen/coliveir/Miguel/Figures/update/Fig1D_2_Markers_I.txt")
genes_P <- readLines("/Users/biagi/cangen/coliveir/Miguel/Figures/update/Fig1D_2_Markers_P.txt")

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
  tmp <- mapply(cbind, tmp, "type" = names(tmp), SIMPLIFY=F)
  tmp <- do.call("rbind", tmp)

  tmp <- subset(tmp, Adjusted.P.value < 0.05)

  tmp$Adjusted.P.value <- -log10(tmp$Adjusted.P.value)
  tmp$Combined.Score <- log2(tmp$Combined.Score)

  tmplist <- rbind(tmplist, data.frame(tmp, class = dbs[i]))

  tmp <- tmp %>% group_by(type) %>% top_n(n = 5, wt = Combined.Score)
  tmp$Term <- factor(tmp$Term, levels = rev(unique(tmp$Term)))
  pathlist <- rbind(pathlist, data.frame(tmp %>% group_by(type) %>% top_n(n = 2, wt = Combined.Score), class = dbs[i]))
  paths <- rbind(paths, data.frame(tmp, class = dbs[i]))

  plotlist[[i]] <- ggplot(tmp, aes(x = type, y = Term)) +
    geom_point(aes(size = Adjusted.P.value, color = Combined.Score)) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) +
    scale_colour_gradient(limits=c(0, max(tmp$Combined.Score)+0.5), high = "#2b9348", low = "#eeef20") +
    xlab(NULL) + ylab(NULL) +
    ggtitle(dbs[i]) + labs(color = "Combined Score", size = "-log10(padj)")
}
tmplist <- tmplist %>% group_by(type) %>% arrange(desc(Combined.Score), .by_group = TRUE)
#write.table(tmplist, file.path("/projects/cangen/coliveir/Miguel/Figures/update/Enrichment/Figure_1D_Table.txt"), col.names = T, row.names = F, sep = "\t", quote = F)

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S1F_1.eps", width = 8, height = 6)
plotlist[[1]] + theme(legend.position = "right",
                      axis.title = element_blank(),
                      axis.text = element_blank(),
                      panel.background = element_blank())
dev.off()

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S1F_2.eps", width = 8, height = 6)
plotlist[[2]] + theme(legend.position = "right",
                      axis.title = element_blank(),
                      axis.text = element_blank(),
                      panel.background = element_blank())
dev.off()

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S1F_3.eps", width = 8, height = 6)
plotlist[[3]] + theme(legend.position = "right",
                      axis.title = element_blank(),
                      axis.text = element_blank(),
                      panel.background = element_blank())
dev.off()

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S1F_4.eps", width = 8, height = 6)
plotlist[[4]] + theme(legend.position = "right",
                      axis.title = element_blank(),
                      axis.text = element_blank(),
                      panel.background = element_blank())
dev.off()



##################################
########### Figure S2A ###########
##################################
filelist = list.files(path = "/projects/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly",
                      pattern = "sccaf_assess", recursive = T, full.names = T)
fnames <- gsub("sccaf_assess_", "", basename(filelist))
fnames <- gsub(".txt", "", fnames)

datalist <- lapply(filelist, function(x)read.csv(x))
names(datalist) <- fnames
datafr <- do.call("rbind", datalist)

pt <- ggplot(datafr, aes(x = Round, y = Accuracy, fill = Type)) +
  geom_violin(
    trim = FALSE,
    position = position_dodge(0.9)
  ) +
  scale_fill_manual(values = c("#d20c0c", "#1111C9")) +
  geom_boxplot(
    width = 0.05,
    position = position_dodge(0.9), colour = "white", outlier.colour = NA, show.legend = FALSE
  ) +
  xlab("") + ggtitle("Test and cross validation accurancies per SCCAF round.") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.9, 0.15), legend.title = element_blank())

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_2A_2.eps", width = 6, height = 6)
print(pt)
dev.off()



##################################
########### Figure S2B ###########
##################################
import warnings
warnings.filterwarnings("ignore")
from SCCAF import *

  ad = sc.read("/projects/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results.h5ad")

y_prob, y_pred, y_test, clf, cvsm, acc = SCCAF_assessment(ad.X, ad.obs['L1_result'],n_jobs=8)
aucs = plot_roc(y_prob, y_test, clf, cvsm=cvsm, acc=acc)
plt.savefig('/projects/cangen/coliveir/Miguel/paper/Figure_2A_3.eps')



##################################
########### Figure S2C ###########
##################################
data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")

infos <- read.table("/projects/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster

pt <- VlnPlot(data, c('Adipoq', 'Adrb3', 'Cidec', 'Dgat1', 'Fasn', 'Lipe', 'Pck1', 'Plin1', 'Pnpla2', 'Retn'),
              cols = c("#11c78b", "#800080", "#e57400", "#0000FF", "#dfdf0d"), assay = "RNA", combine = F)
pt <- lapply(pt, function(x) {
  x + theme_bw(base_size = 9) + xlab(NULL) + ylab(NULL) +
    theme(plot.title = element_text(hjust = 0.5))
})

figure <- ggarrange(plotlist = pt, nrow = 2, ncol = 5, common.legend = T, legend = "bottom")

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S2B.eps", width = 12, height = 6)
annotate_figure(figure, left = text_grob("Expression Level", rot = 90))
dev.off()



##################################
########### Figure S2D ###########
##################################
volcano.plot = function(res, upGenes = NULL, downGenes = NULL){
  mut <- as.data.frame(res)
  mut <- na.omit(mut)
  mutateddf <- dplyr::mutate(mut, sig=ifelse(mut$gene %in% upGenes,"Up_regulated", ifelse(mut$gene %in% downGenes , "Down_regulated", "Not_different")))
  rownames(mutateddf) <- rownames(mut)
  input <- cbind(gene=rownames(mutateddf), mutateddf)
  colnames(input)[which(colnames(input)=="sig")] <- "Significance"
  input[,1] <- NULL
  input[which(input[["p_val_adj"]] == 0), "p_val_adj"] <- min(input[which(input[["p_val_adj"]] != 0), "p_val_adj"], na.rm = TRUE) * 10^-1

  p <- ggplot(input, aes(avg_logFC, -log10(p_val_adj))) +
    geom_point(colour="white") +
    ggtitle("") +
    theme_bw() +
    scale_y_continuous(limits = c(0, -log10(input$p_val_adj)))
  p <- p + geom_point(data=subset(input, input$Significance == 'Not_different'), aes(avg_logFC, -log10(p_val_adj)), colour="gray70") +
    geom_point(data=subset(input, input$Significance == 'Up_regulated'), aes(avg_logFC, -log10(p_val_adj)), colour="firebrick4") +
    geom_point(data=subset(input, input$Significance == 'Down_regulated'), aes(avg_logFC, -log10(p_val_adj)), colour="dodgerblue") +
    xlab("logFC") + ylab("-log10(padj)")
  p <- p + geom_text_repel(data=input[c(head(upGenes, 5), head(downGenes, 5)), ], aes(label=gene))

  return(p)
}

data <- readRDS("/Users/biagi/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")

infos <- read.table("/Users/biagi/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster

markers <- FindAllMarkers(data, logfc.threshold = 0, only.pos = F)

markers_1 <- subset(markers, cluster == "Ad1")
mks_1 <- subset(markers_1, p_val_adj < 0.01)
mks_1 <- mks_1[order(mks_1$p_val_adj), ]
upGenes <- head(subset(mks_1, avg_logFC > 0)$gene, 50)
downGenes <- head(subset(mks_1, avg_logFC < 0)$gene, 50)
pt_1 <- volcano.plot(res = markers_1, upGenes = upGenes, downGenes = downGenes) + ggtitle("Ad1") + theme(plot.title = element_text(hjust = 0.5))
#write.table(c(upGenes, downGenes), "/projects/cangen/coliveir/Miguel/Figures/update/Fig2C_Markers_1.txt", row.names = F, col.names = F, quote = F, sep = "\t")


markers_2 <- subset(markers, cluster == "Ad2")
mks_2 <- subset(markers_2, p_val_adj < 0.01)
mks_2 <- mks_2[order(mks_2$p_val_adj), ]
upGenes <- head(subset(mks_2, avg_logFC > 0)$gene, 50)
downGenes <- head(subset(mks_2, avg_logFC < 0)$gene, 50)
pt_2 <- volcano.plot(res = markers_2, upGenes = upGenes, downGenes = downGenes) + ggtitle("Ad2") + theme(plot.title = element_text(hjust = 0.5))
#write.table(c(upGenes, downGenes), "/projects/cangen/coliveir/Miguel/Figures/update/Fig2C_Markers_2.txt", row.names = F, col.names = F, quote = F, sep = "\t")


markers_3 <- subset(markers, cluster == "Ad3")
mks_3 <- subset(markers_3, p_val_adj < 0.01)
mks_3 <- mks_3[order(mks_3$p_val_adj), ]
upGenes <- head(subset(mks_3, avg_logFC > 0)$gene, 50)
downGenes <- head(subset(mks_3, avg_logFC < 0)$gene, 50)
pt_3 <- volcano.plot(res = markers_3, upGenes = upGenes, downGenes = downGenes) + ggtitle("Ad3") + theme(plot.title = element_text(hjust = 0.5))
#write.table(c(upGenes, downGenes), "/projects/cangen/coliveir/Miguel/Figures/update/Fig2C_Markers_3.txt", row.names = F, col.names = F, quote = F, sep = "\t")


markers_4 <- subset(markers, cluster == "Ad4")
mks_4 <- subset(markers_4, p_val_adj < 0.01)
mks_4 <- mks_4[order(mks_4$p_val_adj), ]
upGenes <- head(subset(mks_4, avg_logFC > 0)$gene, 50)
downGenes <- head(subset(mks_4, avg_logFC < 0)$gene, 50)
pt_4 <- volcano.plot(res = markers_4, upGenes = upGenes, downGenes = downGenes) + ggtitle("Ad4") + theme(plot.title = element_text(hjust = 0.5))
#write.table(c(upGenes, downGenes), "/projects/cangen/coliveir/Miguel/Figures/update/Fig2C_Markers_4.txt", row.names = F, col.names = F, quote = F, sep = "\t")


markers_5 <- subset(markers, cluster == "Ad5")
mks_5 <- subset(markers_5, p_val_adj < 0.01)
mks_5 <- mks_5[order(mks_5$p_val_adj), ]
upGenes <- head(subset(mks_5, avg_logFC > 0)$gene, 50)
downGenes <- head(subset(mks_5, avg_logFC < 0)$gene, 50)
pt_5 <- volcano.plot(res = markers_5, upGenes = upGenes, downGenes = downGenes) + ggtitle("Ad5") + theme(plot.title = element_text(hjust = 0.5))
#write.table(c(upGenes, downGenes), "/projects/cangen/coliveir/Miguel/Figures/update/Fig2C_Markers_5.txt", row.names = F, col.names = F, quote = F, sep = "\t")

postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S2C.eps", width = 10, height = 8)
ggarrange(pt_1, pt_2, pt_3, pt_4, pt_5, nrow = 2, ncol = 3, common.legend = T)
dev.off()



##################################
########### Figure S2E ###########
##################################
dbs <- listEnrichrDbs()
dbs <- c("Jensen_TISSUES", "Mouse_Gene_Atlas")

genes_1 <- readLines("/Users/biagi/cangen/coliveir/Miguel/Figures/update/Fig2C_Markers_1.txt")
genes_2 <- readLines("/Users/biagi/cangen/coliveir/Miguel/Figures/update/Fig2C_Markers_2.txt")
genes_3 <- readLines("/Users/biagi/cangen/coliveir/Miguel/Figures/update/Fig2C_Markers_3.txt")
genes_4 <- readLines("/Users/biagi/cangen/coliveir/Miguel/Figures/update/Fig2C_Markers_4.txt")
genes_5 <- readLines("/Users/biagi/cangen/coliveir/Miguel/Figures/update/Fig2C_Markers_5.txt")

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
  tmp <- mapply(cbind, tmp, "type" = names(tmp), SIMPLIFY=F)
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
    theme_bw(base_size = 12) +
    scale_colour_gradient(limits=c(0, max(tmp$Combined.Score)+0.5), high = "#2b9348", low = "#eeef20") +
    xlab(NULL) + ylab(NULL) +
    ggtitle(dbs[i]) + labs(color = "Combined Score", size = "-log10(padj)")
}
tmplist <- tmplist %>% group_by(type) %>% arrange(desc(Combined.Score), .by_group = TRUE)
#write.table(tmplist, "/projects/cangen/coliveir/Miguel/Figures/update/Enrichment/Figure_2C_2_Table.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## Figuew S2C
postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S2D.eps", width = 8, height = 6)
print(plotlist[[1]])
dev.off()

## Figuew S2D
postscript("/projects/cangen/coliveir/Miguel/paper/Figure_S2E.eps", width = 8, height = 6)
print(plotlist[[2]])
dev.off()



##################################
########### Figure S3A ###########
##################################
data <- readRDS("/Users/biagi/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
data$timpoint <- gsub('4day', 'Cold', data$timpoint)

Idents(data) <- as.factor(data$timpoint)

markers_ColdxRT <- FindMarkers(data, ident.1 = "Cold", ident.2 = "RT", logfc.threshold = 0)
markers_ColdxRT$gene <- rownames(markers_ColdxRT)

markers_CLxRT <- FindMarkers(data, ident.1 = "CL", ident.2 = "RT", logfc.threshold = 0)
markers_CLxRT$gene <- rownames(markers_CLxRT)

data <- RunALRA(data)
data <- ScaleData(data)

data2 <- subset(data, idents = c("Cold", "RT", "CL"))

top1 <- markers_ColdxRT %>% top_n(n = 50, wt = avg_logFC)
top2 <- markers_CLxRT %>% top_n(n = 50, wt = avg_logFC)
top10 <- rbind(top1, top2)

cls <- c("blue", "red", "green")
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

ht1 <- DoHeatmap(data2, features = top10$gene, group.colors = cls, angle = 0, size = 5, label = F) +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='Timepoint')

ht2 <- ht1 + theme(legend.position = "none",
                   panel.grid = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   panel.background = element_blank())

postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3A_1.eps", width = 8, height = 6)
print(ht1)
dev.off()

tiff("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3A_2.tiff", width = 8, height = 6, units = 'in', res = 600)
print(ht2)
dev.off()



##################################
########### Figure S3B ###########
##################################
dbs <- listEnrichrDbs()
dbs <- c("KEGG_2019_Mouse", "WikiPathways_2019_Mouse", "Jensen_TISSUES", "GO_Biological_Process_2018")

genes_ColdxRT <- readLines("/Users/biagi/cangen/coliveir/Miguel/Figures/update/Diff_Adipocytes_ColdxRT.txt")
genes_CLxRT <- readLines("/Users/biagi/cangen/coliveir/Miguel/Figures/update/Diff_Adipocytes_CLxRT.txt")

genes <- list(ColdxRT = genes_ColdxRT,
              CLxRT = genes_CLxRT)

results <- lapply(genes, enrichr, dbs)

plotlist <- list()
pathlist <- NULL
tmplist <- NULL
paths <- NULL
for (i in 1:length(dbs)) {
  tmp <- lapply(results, `[[`, i)
  tmp <- mapply(cbind, tmp, "type" = names(tmp), SIMPLIFY=F)
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
    scale_colour_gradient(limits=c(0, max(tmp$Combined.Score)+0.5), high = "#2b9348", low = "#eeef20") +
    xlab(NULL) + ylab(NULL) +
    ggtitle(dbs[i]) + labs(color = "Combined Score", size = "-log10(padj)")
}
tmplist <- tmplist %>% group_by(type) %>% arrange(desc(Combined.Score), .by_group = TRUE)
#write.table(tmplist, "/Users/biagi/cangen/coliveir/Miguel/Figures/update/Enrichment/Diff_Table.txt", col.names = T, row.names = F, sep = "\t", quote = F)

paths <- paths[which(paths$Term %in% c("PPAR signaling pathway", "Propanoate metabolism", "Regulation of lipolysis in adipocytes",
                                       "Fatty acid biosynthesis", "Fatty acid degradation", "Valine, leucine and isoleucine degradation",
                                       "PPAR signaling pathway WP2316", "Triacylglyceride Synthesis WP386", "Subcutaneous adipose tissue",
                                       "3T3-L1 cell", "Adipocyte", "long-chain fatty acid transport", "fatty acid transmembrane transport",
                                       "regulation of sequestering of triglyceride", "intracellular lipid transport")), ]

paths <- paths %>%
  mutate(
    class_with_color = ifelse(class == "KEGG_2019_Mouse", glue::glue("<strong><span style='color:#0000FF'>{Term}</span></strong>"),
                              ifelse(class == "WikiPathways_2019_Mouse", glue::glue("<strong><span style='color:#f28804'>{Term}</span></strong>"),
                                     ifelse(class == "Jensen_TISSUES", glue::glue("<strong><span style='color:#CC0000'>{Term}</span></strong>"),
                                            ifelse(class == "GO_Biological_Process_2018", glue::glue("<strong><span style='color:#9900CC'>{Term}</span></strong>"),
                                                   glue::glue("<strong><span style='color:#18a997'>{Term}</span></strong>")))))
  )


paths <- paths %>% group_by(type) %>% arrange(desc(Combined.Score), .by_group = TRUE)
paths$class_with_color <- factor(paths$class_with_color, levels = rev(unique(paths$class_with_color)))

pt <- ggplot(paths, aes(x = type, y = class_with_color)) +
  geom_point(aes(size = Adjusted.P.value, color = Combined.Score)) +
  theme_bw(base_size = 8) +
  scale_colour_gradient(limits=c(0, max(paths$Combined.Score)+0.5), high = "#2b9348", low = "#eeef20") +
  xlab(NULL) + ylab(NULL) +
  ggtitle(NULL) + labs(color = "Combined Score", size = "-log10(padj)") +
  theme(axis.text.y = element_markdown(), axis.text=element_text(size=12))

postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3B.eps", width = 8, height = 6)
print(pt)
dev.off()



##################################
########### Figure S3C ###########
##################################
data <- readRDS("/Users/biagi/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")

infos <- read.table("/Users/biagi/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- factor(data$timpoint)
data$cluster <- factor(new_cluster)

data <- RunALRA(data)
data <- ScaleData(data, features = rownames(data))


##Fatty Acid Oxidation
genes1 <- c('Abcb11', 'Abcd1', 'Abcd2', 'Abcd3', 'Abcd4', 'Acaa1a', 'Acaa1b', 'Acaa2', 'Acacb', 'Acad11', 'Acadl', 'Acadm', 'Acads', 'Acadvl', 'Acat1', 'Acat2', 'Acat3', 'Acox1', 'Acox2', 'Acox3', 'Acoxl', 'Acsbg2', 'Acsl5', 'Adh4', 'Adh5', 'Adh7', 'Adipoq', 'Adipor1', 'Adipor2', 'Akt1', 'Akt2', 'Alox12', 'Appl2', 'Auh', 'Bdh2', 'C1qtnf2', 'C1qtnf9', 'Cd36', 'Cnr1', 'Cpt1a', 'Cpt2', 'Crat', 'Crot', 'Cygb', 'Cyp4v3', 'Cyp24a1', 'Dbi', 'Decr1', 'Dgat1', 'Dgat2', 'Echdc1', 'Echdc2', 'Echs1', 'Eci1', 'Eci2', 'Eci3', 'Ehhadh', 'Etfa', 'Etfb', 'Etfbkmt', 'Etfdh', 'Fabp1', 'Fabp3', 'Gcdh', 'Gm45753', 'Gm49387', 'Hacl1', 'Hadh', 'Hadha', 'Hadhb', 'Hao1', 'Hao2', 'Hsd17b4', 'Hsd17b10', 'Ilvbl', 'Irs1', 'Irs2', 'Ivd', 'Lep', 'Lonp2', 'Mapk14', 'Mfsd2a', 'Mir199a-2', 'Mir214', 'Mir696', 'Mlycd', 'Mtor', 'Nr4a3', 'Nucb2', 'Pdk4', 'Pex2', 'Pex5', 'Pex7', 'Pex13', 'Phyh', 'Plin5', 'Por', 'Ppara', 'Ppard', 'Pparg', 'Ppargc1a', 'Prkaa1', 'Scp2', 'Sesn2', 'Sirt4', 'Slc25a17', 'Slc27a2', 'Sox9', 'Twist1', 'Tysnd1')
genes1 <- genes1[genes1 %in% rownames(data)]
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
ht1 <- DoHeatmap(data, features = genes1, angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression')
y <- ht1$data %>% drop_na()
x <- y %>% group_by(Identity) %>% dplyr::select(Feature, Cell, Identity, Expression) %>%
  tidyr::spread(key = Feature, value = Expression)
w <- y %>% dplyr::select(Feature, Cell, Expression) %>%
  tidyr::spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()
pt <- Heatmap(w, cluster_columns = F)
ht1 <- DoHeatmap(data, features = rownames(w)[row_order(pt)], angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression') + ggtitle('Fatty Acid Oxidation') + 
  theme(plot.title = element_text(hjust = 0.5))
ht2 <- ht1 + theme(legend.position = "none",
                   panel.grid = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   panel.background = element_blank())
postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3D_1.eps", width = 8, height = 6)
print(ht1)
dev.off()

tiff("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3D_1.tiff", width = 8, height = 6, units = 'in', res = 600)
print(ht2)
dev.off()


##Tricarboxylic Acid Cycle
genes2 <- c('4933405O20Rik', 'Aco1', 'Aco2', 'Cs', 'Csl', 'Dhtkd1', 'Dlat', 'Dlst', 'Fh1', 'Idh1', 'Idh2', 'Idh3a', 'Idh3b', 'Idh3g', 'Ireb2', 'Mdh1', 'Mdh1b', 'Mdh2', 'Ndufs4', 'Ogdh', 'Ogdhl', 'Pdha1', 'Pdha2', 'Pdhb', 'Sdha', 'Sdhaf2', 'Sdhb', 'Sdhc', 'Sdhd', 'Sucla2', 'Suclg1', 'Suclg2')
genes2 <- genes2[genes2 %in% rownames(data)]
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
ht1 <- DoHeatmap(data, features = genes2, angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression')
y <- ht1$data %>% drop_na()
x <- y %>% group_by(Identity) %>% dplyr::select(Feature, Cell, Identity, Expression) %>%
  tidyr::spread(key = Feature, value = Expression)
w <- y %>% dplyr::select(Feature, Cell, Expression) %>%
  tidyr::spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()
pt <- Heatmap(w, cluster_columns = F)
ht1 <- DoHeatmap(data, features = rownames(w)[row_order(pt)], angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression') + ggtitle('Tricarboxylic Acid Cycle') + 
  theme(plot.title = element_text(hjust = 0.5))
ht2 <- ht1 + theme(legend.position = "none",
                   panel.grid = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   panel.background = element_blank())
postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3D_2.eps", width = 8, height = 6)
print(ht1)
dev.off()

tiff("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3D_2.tiff", width = 8, height = 6, units = 'in', res = 600)
print(ht2)
dev.off()


##Fat Acid Transport
genes3 <- c('Abcc1', 'Abcc2', 'Abcc4', 'Abcd1', 'Abcd2', 'Abcd3', 'Abcd4', 'Ace', 'Acsl1', 'Acsl3', 'Acsl4', 'Acsl5', 'Acsl6', 'Agtr2', 'Akt1', 'Akt2', 'Anxa1', 'Apoe', 'Atp5j', 'Avpr1b', 'Bdkrb2', 'Cd36', 'Cpt1b', 'Crot', 'Cyp4f18', 'Drd2', 'Drd3', 'Drd4', 'Edn1', 'Eprs', 'Erfe', 'Fabp1', 'Fabp2', 'Fabp3', 'Fabp4', 'Fabp5', 'Hnf1a', 'Hrh2', 'Il1a', 'Il1b', 'Irs2', 'Kiss1r', 'Lep', 'Lhcgr', 'Map2k6', 'Mapk9', 'Mfsd2a', 'Mif', 'Nmb', 'Nmur2', 'Nos2', 'Ntsr1', 'Oc90', 'Oxt', 'P2rx7', 'P2ry2', 'Pla2g1b', 'Pla2g2c', 'Pla2g2d', 'Pla2g2e', 'Pla2g2f', 'Pla2g3', 'Pla2g4a', 'Pla2g4f', 'Pla2g5', 'Pla2g6', 'Pla2g10', 'Pla2g12a', 'Pla2g12b', 'Pla2r1', 'Plin2', 'Pnpla8', 'Pparg', 'Ptges', 'Repin1', 'Rps6kb1', 'Slc2a1', 'Slc5a8', 'Slc22a22', 'Slc25a17', 'Slc27a1', 'Slc27a2', 'Slc27a3', 'Slc27a4', 'Slc27a5', 'Slc27a6', 'Slco2a1', 'Slco3a1', 'Spx', 'Sstr4', 'Syk', 'Thbs1', 'Tnfrsf11a', 'Tnfsf11')
genes3 <- genes3[genes3 %in% rownames(data)]
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
ht1 <- DoHeatmap(data, features = genes3, angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression')
y <- ht1$data %>% drop_na()
x <- y %>% group_by(Identity) %>% dplyr::select(Feature, Cell, Identity, Expression) %>%
  tidyr::spread(key = Feature, value = Expression)
w <- y %>% dplyr::select(Feature, Cell, Expression) %>%
  tidyr::spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()
pt <- Heatmap(w, cluster_columns = F)
ht1 <- DoHeatmap(data, features = rownames(w)[row_order(pt)], angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression') + ggtitle('Fat Acid Transport') + 
  theme(plot.title = element_text(hjust = 0.5))
ht2 <- ht1 + theme(legend.position = "none",
                   panel.grid = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   panel.background = element_blank())
postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3D_3.eps", width = 8, height = 6)
print(ht1)
dev.off()

tiff("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3D_3.tiff", width = 8, height = 6, units = 'in', res = 600)
print(ht2)
dev.off()


##Triglyceride/Fatty Acid Cycle
genes4 <- c('Aadac', 'Abhd5', 'Daglb', 'Pnpla2', 'Pck1', 'Pcx' , 'Gpd1', 'Slc2a4', 'Lipe', 'Adrb3')
genes4 <- genes4[genes4 %in% rownames(data)]
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
ht1 <- DoHeatmap(data, features = genes4, angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression')
y <- ht1$data %>% drop_na()
x <- y %>% group_by(Identity) %>% dplyr::select(Feature, Cell, Identity, Expression) %>%
  tidyr::spread(key = Feature, value = Expression)
w <- y %>% dplyr::select(Feature, Cell, Expression) %>%
  tidyr::spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()
pt <- Heatmap(w, cluster_columns = F)
ht1 <- DoHeatmap(data, features = rownames(w)[row_order(pt)], angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression') + ggtitle('Triglyceride/Fatty Acid Cycle') + 
  theme(plot.title = element_text(hjust = 0.5))
ht2 <- ht1 + theme(legend.position = "none",
                   panel.grid = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   panel.background = element_blank())
postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3D_4.eps", width = 8, height = 6)
print(ht1)
dev.off()

tiff("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3D_4.tiff", width = 8, height = 6, units = 'in', res = 600)
print(ht2)
dev.off()



##################################
########### Figure S3D ###########
##################################
## CL
data <- readRDS("/Users/biagi/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
data$timpoint <- gsub('4day', 'Cold', data$timpoint)

infos <- read.table("/Users/biagi/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")
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
                      expressionFamily=negbinomial.size())

rm(data, mat, fdata, pdata, fd, pd)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

diff_test_res <- differentialGeneTest(cds[expressed_genes, ], fullModelFormulaStr = "~cluster", cores = 6, verbose = T)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree', verbose = TRUE)
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 3)

postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3D_1.eps", width = 8, height = 6)
plot_cell_trajectory(cds, color_by = "timpoint")
dev.off()

postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3D_2.eps", width = 8, height = 6)
plot_cell_trajectory(cds, color_by="State", state_number_size = 1) + facet_wrap(~State)
dev.off()

source('/Users/biagi/cangen/coliveir/plotMonocle.R')
postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3D_3.eps", width = 8, height = 15)
print(plotMonocle(cds, c('Ucp1', 'Ppara', 'Dio2', 'Chst1', 'Plppr3', 'Nnat', 'Pim1', 'Adcy3', 'Cs')))
dev.off()


## Cold
data <- readRDS("/Users/biagi/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
data$timpoint <- gsub('4day', 'Cold', data$timpoint)

infos <- read.table("/Users/biagi/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")
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
                      expressionFamily=negbinomial.size())

rm(data, mat, fdata, pdata, fd, pd)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

diff_test_res <- differentialGeneTest(cds[expressed_genes, ], fullModelFormulaStr = "~cluster", cores = 6, verbose = T)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree', verbose = TRUE)
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 4)

postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3D_4.eps", width = 8, height = 6)
plot_cell_trajectory(cds, color_by = "timpoint")
dev.off()

postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3D_5.eps", width = 8, height = 6)
plot_cell_trajectory(cds, color_by="State", state_number_size = 1) + facet_wrap(~State)
dev.off()

source('/Users/biagi/cangen/coliveir/plotMonocle.R')
postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S3D_6.eps", width = 8, height = 15)
print(plotMonocle(cds, c('Ucp1', 'Ppara', 'Dio2', 'Ccdc80', 'Slc7a10', 'Tmem43', 'Adcy3', 'Perm1', 'Fabp3')))
dev.off()



##################################
########### Figure S4A ###########
##################################
data <- readRDS("/Users/biagi/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")

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

  piePlot <- ggplot(df, aes(x = "", y = value, fill = class)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0)+
    scale_fill_manual(values = c("#11c78b", "#800080", "#dfdf0d")) +
    theme_void() + labs(fill = "Clusters") +
    geom_text_repel(aes(x = 1.5, y = pos, label = paste0(value, "%")), nudge_x = .1, segment.size = .7, show.legend = FALSE) +
    ggtitle(genes[i])

  postscript(paste0("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S4A_", i, ".eps"), width = 7, height = 7)
  print(piePlot)
  dev.off()
}



##################################
########### Figure S4B ###########
##################################
data <- readRDS("/Users/biagi/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")

infos <- read.table("/Users/biagi/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")

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

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

ht1 <- DoHeatmap(data, features = markers$gene, angle = 0, size = 5, label = F) +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=3)) +
  labs(color='UCP1 Expression')

ht2 <- ht1 + theme(legend.position = "none",
                   panel.grid = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   panel.background = element_blank())

postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S4B_1.eps", width = 8, height = 6)
print(ht1)
dev.off()

tiff("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S4B_2.tiff", width = 8, height = 6, units = 'in', res = 600)
print(ht2)
dev.off()



##################################
########### Figure S4C ###########
##################################
data <- readRDS("/Users/biagi/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")
new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster
data$cluster <- Idents(data)
data <- subset(data, idents = list('Ad1'))
tmp1 <- factor(ifelse(data@assays$SCT@data["Ucp1", ] > 0, "Ad1_High", "Ad1_Low"))


data <- readRDS("/Users/biagi/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
infos <- read.table("/Users/biagi/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")
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

data <- RunALRA(data)
data <- ScaleData(data, features = rownames(data))


genes8 <- c('Atp2a1', 'Atp2a2', 'Tmtc4', 'Adra1a', 'Arpc2')
genes8 <- genes8[genes8 %in% rownames(data)]
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
ht1 <- DoHeatmap(data, features = genes8, angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression')
y <- ht1$data %>% drop_na()
x <- y %>% group_by(Identity) %>% dplyr::select(Feature, Cell, Identity, Expression) %>%
  tidyr::spread(key = Feature, value = Expression)
w <- y %>% dplyr::select(Feature, Cell, Expression) %>%
  tidyr::spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()
pt <- Heatmap(w, cluster_columns = F)
ht1 <- DoHeatmap(data, features = rownames(w)[row_order(pt)], angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression') + ggtitle('SERCA2') + 
  theme(plot.title = element_text(hjust = 0.5))
ht2 <- ht1 + theme(legend.position = "none",
                   panel.grid = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   panel.background = element_blank())
postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S4C_1.eps", width = 8, height = 6)
print(ht1)
dev.off()

tiff("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S4C_1.tiff", width = 8, height = 6, units = 'in', res = 600)
print(ht2)
dev.off()


genes6 <- c('Slc6a8', 'Gatm', 'Gamt', 'Ckmt1', 'Ckmt2')
genes6 <- genes6[genes6 %in% rownames(data)]
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
ht1 <- DoHeatmap(data, features = genes6, angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression')
y <- ht1$data %>% drop_na()
x <- y %>% group_by(Identity) %>% dplyr::select(Feature, Cell, Identity, Expression) %>%
  tidyr::spread(key = Feature, value = Expression)
w <- y %>% dplyr::select(Feature, Cell, Expression) %>%
  tidyr::spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()
pt <- Heatmap(w, cluster_columns = F)
ht1 <- DoHeatmap(data, features = rownames(w)[row_order(pt)], angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression') + ggtitle('Arginine/creatine and proline metabolism v') + 
  theme(plot.title = element_text(hjust = 0.5))
ht2 <- ht1 + theme(legend.position = "none",
                   panel.grid = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   panel.background = element_blank())
postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S4C_2.eps", width = 8, height = 6)
print(ht1)
dev.off()

tiff("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S4C_2.tiff", width = 8, height = 6, units = 'in', res = 600)
print(ht2)
dev.off()


genes <- c('Actn3', 'Adpgk', 'Aldoa', 'Aldoart1', 'Aldoart2', 'Aldob', 'Aldoc', 'App', 'Bpgm', 'Cbfa2t3', 'Ddit4', 'Dhtkd1', 'Eif6', 'Eno1', 'Eno1b', 'Eno2', 'Eno3', 'Eno4', 'Entpd5', 'Ep300', 'Esrrb', 'Fbp1', 'Foxk1', 'Foxk2', 'Gale', 'Galk1', 'Galt', 'Gapdh', 'Gapdhs', 'Gck', 'Gm3839', 'Gm10358', 'Gm11214', 'Gm12117', 'Gm15294', 'Gpd1', 'Gpi1', 'Hdac4', 'Hif1a', 'Hk1', 'Hk2', 'Hk3', 'Hkdc1', 'Htr2a', 'Ier3', 'Ifng', 'Igf1', 'Il3', 'Ins2', 'Insr', 'Jmjd8', 'Khk', 'Mif', 'Mlxipl', 'Mpi', 'Mtch2', 'Myc', 'Myog', 'Ncor1', 'Nupr1', 'Ogdh', 'Ogt', 'P2rx7', 'Pfkfb2', 'Pfkl', 'Pfkm', 'Pfkp', 'Pgam1', 'Pgam2', 'Pgk1', 'Pgk2', 'Pklr', 'Pkm', 'Ppara', 'Ppargc1a', 'Prkaa1', 'Prkaa2', 'Prkag2', 'Prkag3', 'Prxl2c', 'Psen1', 'Sirt6', 'Slc2a6', 'Slc4a1', 'Slc4a4', 'Stat3', 'Tigar', 'Tkfc', 'Tpi1', 'Trex1', 'Zbtb7a', 'Zbtb20')
genes <- genes[genes %in% rownames(data)]
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
ht1 <- DoHeatmap(data, features = genes, angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression')
y <- ht1$data %>% drop_na()
x <- y %>% group_by(Identity) %>% dplyr::select(Feature, Cell, Identity, Expression) %>%
  tidyr::spread(key = Feature, value = Expression)
w <- y %>% dplyr::select(Feature, Cell, Expression) %>%
  tidyr::spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()
pt <- Heatmap(w, cluster_columns = F)
ht1 <- DoHeatmap(data, features = rownames(w)[row_order(pt)], angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression') + ggtitle('Glycolytic Process') + 
  theme(plot.title = element_text(hjust = 0.5))
ht2 <- ht1 + theme(legend.position = "none",
                   panel.grid = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   panel.background = element_blank())
postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S4C_3.eps", width = 8, height = 6)
print(ht1)
dev.off()

tiff("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S4C_3.tiff", width = 8, height = 6, units = 'in', res = 600)
print(ht2)
dev.off()



##################################
########### Figure S4H ###########
##################################
data <- readRDS("/Users/biagi/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")

infos <- read.table("/Users/biagi/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")

ort <- read.table("/Users/biagi/cangen/coliveir/cellphonedb/Orthologs_human_mouse.txt", sep = ",", header = T)

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
ht <- Heatmap(as.matrix(res), 
              bottom_annotation = ba, 
              name = "zscore", column_title = "Ad1 UCP1 Targets TF High", width = 1, 
              show_row_names = T, show_column_names = F,
              cluster_rows = F, cluster_columns = T,
              col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)), 
              heatmap_height = unit(6, "cm"), row_names_gp = gpar(fontsize = 8))

postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S4D_1.eps", width = 8, height = 6)
print(ht)
dev.off()


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
ht <- Heatmap(as.matrix(res), 
              bottom_annotation = ba, 
              name = "zscore", column_title = "Ad1 UCP1 Targets TF Low", width = 1, 
              show_row_names = T, show_column_names = F,
              cluster_rows = F, cluster_columns = T,
              col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)), 
              heatmap_height = unit(6, "cm"), row_names_gp = gpar(fontsize = 8))

postscript("/Users/biagi/cangen/coliveir/Miguel/paper/Figure_S4D_2.eps", width = 8, height = 6)
print(ht)
dev.off()