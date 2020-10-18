## Loading R packages
library(Seurat)
library(ggplot2)
library(metacell)
library(dplyr)
library(ggpubr)
library(pheatmap)
library(enrichR)


##################################
########### Figure S1A ###########
##################################
data <- readRDS("/Users/biagi/PhD/AdipoSNAP/output/10x/10x_Processed.rds")
Idents(data) <- factor(data$CellTypeRefined)

TSNEPlot(data, label = FALSE, pt.size = 0.3) +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  theme_classic() + labs(color = "Cluster") +
  theme(legend.position="bottom")



##################################
########### Figure S1B ###########
##################################
load("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell/db/mc.test_mc_f.Rda")
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

load("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell/db/mc2d.test_2dproj.Rda")
dims <- data.frame(x = object@sc_x,
                   y = object@sc_y)

load("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell/db/mc.test_mc_f.Rda")
tmp1 <- data.frame(cells = names(object@mc), cols = object@mc)
tmp2 <- data.frame(cols = object@colors)
teste <- merge(tmp1, tmp2, by.x = "cols", by.y = "row.names")
teste$cols <- NULL

teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Adipocyte", marks_colors$group)], "Adipocytes", "Unknown")
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Progenitor", marks_colors$group)], "Progenitors", teste$cellType)
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Immune", marks_colors$group)], "Immunes", teste$cellType)
teste$cellType <- ifelse(teste$cols.y %in% marks_colors$color[grep("Endothelial", marks_colors$group)], "Endothelials", teste$cellType)
tab <- merge(dims, teste, by.x = "row.names", by.y = "cells")

data <- readRDS("/Users/biagi/PhD/AdipoSNAP/output/10x/10x_Processed.rds")
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

TSNEPlot(data, label = FALSE, pt.size = 0.3, cols = cls) +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  theme_classic() + labs(color = "Cluster") +
  theme(legend.position="bottom")


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

ggplot(df, aes(x = "", y = value, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = paste0(value, "%")), color = "white", fontface = "bold")+
  scale_fill_manual(values = c("#FFA500", "#329932", "#ff9999", "#6666ff")) +
  theme_void()



##################################
########### Figure S1C ###########
##################################
filelist = list.files(path = "/Users/biagi/PhD/AdipoSNAP/SCCAF/Adipocytes",
                      pattern = "sccaf_assess", recursive = T, full.names = T)
fnames <- gsub("sccaf_assess_", "", basename(filelist))
fnames <- gsub(".txt", "", fnames)

datalist <- lapply(filelist, function(x)read.csv(x))
names(datalist) <- fnames
datafr <- do.call("rbind", datalist)

ggplot(datafr, aes(x = Round, y = Accuracy, fill = Type)) +
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



##################################
########### Figure S1D ###########
##################################
import warnings
warnings.filterwarnings("ignore")
from SCCAF import *
  
ad = sc.read("/Users/biagi/PhD/AdipoSNAP/SCCAF/Adipocytes/results.h5ad")

y_prob, y_pred, y_test, clf, cvsm, acc = SCCAF_assessment(ad.X, ad.obs['L1_result'],n_jobs=8)
aucs = plot_roc(y_prob, y_test, clf, cvsm=cvsm, acc=acc)
plt.savefig('/Users/biagi/PhD/AdipoSNAP/paper/Figure_SD.eps')



##################################
########### Figure S1E ###########
##################################
data <- readRDS("/Users/biagi/PhD/AdipoSNAP/output/10x/10x_SCT_Processed_ALRA.rds")

pt <- FeaturePlot(data, c("Adrb3", "Pecam1", "Ptprc", "Cd34", "Pdgfra", "Itgb1"),
                  cols = c("grey", 'red'), reduction = 'tsne', pt.size = 0.1, combine = F)
pt <- lapply(pt, function(x) {
  x + theme_classic() + xlab("t-SNE 1") + ylab("t-SNE 2") +
    theme(plot.title = element_text(hjust = 0.5, face = "italic"))
})

ggarrange(plotlist = pt, common.legend = T)



##################################
########### Figure S1F ###########
##################################
scdb_init("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell_SCT2/db", force_reinit=T)
scfigs_init("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell_SCT2/figs")

load("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell_SCT2/db/mc.test_mc_f.Rda")
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

pheatmap(mat_A, fontsize = 8, main = 'Adipocytes', legend = TRUE, treeheight_row = 0, treeheight_col = 0)

pheatmap(mat_E, fontsize = 8, main = 'Endothelials', legend = TRUE, treeheight_row = 0, treeheight_col = 0)

pheatmap(mat_I, fontsize = 8, main = 'Immunes', legend = TRUE, treeheight_row = 0, treeheight_col = 0)

pheatmap(mat_P, fontsize = 8, main = 'Progenitors', legend = TRUE, treeheight_row = 0, treeheight_col = 0)



##################################
########### Figure S1G ###########
##################################
dbs <- c("KEGG_2019_Mouse", "WikiPathways_2019_Mouse", "Jensen_TISSUES", "GO_Biological_Process_2018")

genes_A <- readLines("/Users/biagi/PhD/AdipoSNAP/Figures/update/Fig1D_2_Markers_A.txt")
genes_E <- readLines("/Users/biagi/PhD/AdipoSNAP/Figures/update/Fig1D_2_Markers_E.txt")
genes_I <- readLines("/Users/biagi/PhD/AdipoSNAP/Figures/update/Fig1D_2_Markers_I.txt")
genes_P <- readLines("/Users/biagi/PhD/AdipoSNAP/Figures/update/Fig1D_2_Markers_P.txt")

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

ggarrange(plotlist = plotlist)