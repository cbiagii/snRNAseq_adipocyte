## Loading R packages
library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(enrichR)
library(dplyr)

source('0_Function_volcanoPlot.R') #####ARRUMAR CAMINHO**********


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
annotate_figure(figure, left = text_grob("Expression Level", rot = 90))



##################################
########### Figure S2D ###########
##################################
data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
infos <- read.table("/projects/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")

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

markers_2 <- subset(markers, cluster == "Ad2")
mks_2 <- subset(markers_2, p_val_adj < 0.01)
mks_2 <- mks_2[order(mks_2$p_val_adj), ]
upGenes <- head(subset(mks_2, avg_logFC > 0)$gene, 50)
downGenes <- head(subset(mks_2, avg_logFC < 0)$gene, 50)
pt_2 <- volcano.plot(res = markers_2, upGenes = upGenes, downGenes = downGenes) + ggtitle("Ad2") + theme(plot.title = element_text(hjust = 0.5))

markers_3 <- subset(markers, cluster == "Ad3")
mks_3 <- subset(markers_3, p_val_adj < 0.01)
mks_3 <- mks_3[order(mks_3$p_val_adj), ]
upGenes <- head(subset(mks_3, avg_logFC > 0)$gene, 50)
downGenes <- head(subset(mks_3, avg_logFC < 0)$gene, 50)
pt_3 <- volcano.plot(res = markers_3, upGenes = upGenes, downGenes = downGenes) + ggtitle("Ad3") + theme(plot.title = element_text(hjust = 0.5))

markers_4 <- subset(markers, cluster == "Ad4")
mks_4 <- subset(markers_4, p_val_adj < 0.01)
mks_4 <- mks_4[order(mks_4$p_val_adj), ]
upGenes <- head(subset(mks_4, avg_logFC > 0)$gene, 50)
downGenes <- head(subset(mks_4, avg_logFC < 0)$gene, 50)
pt_4 <- volcano.plot(res = markers_4, upGenes = upGenes, downGenes = downGenes) + ggtitle("Ad4") + theme(plot.title = element_text(hjust = 0.5))

markers_5 <- subset(markers, cluster == "Ad5")
mks_5 <- subset(markers_5, p_val_adj < 0.01)
mks_5 <- mks_5[order(mks_5$p_val_adj), ]
upGenes <- head(subset(mks_5, avg_logFC > 0)$gene, 50)
downGenes <- head(subset(mks_5, avg_logFC < 0)$gene, 50)
pt_5 <- volcano.plot(res = markers_5, upGenes = upGenes, downGenes = downGenes) + ggtitle("Ad5") + theme(plot.title = element_text(hjust = 0.5))

ggarrange(pt_1, pt_2, pt_3, pt_4, pt_5, nrow = 2, ncol = 3, common.legend = T)



##################################
########### Figure S2E ###########
##################################
dbs <- c("Jensen_TISSUES", "Mouse_Gene_Atlas")

genes_1 <- readLines("/projects/cangen/coliveir/Miguel/Figures/update/Fig2C_Markers_1.txt")
genes_2 <- readLines("/projects/cangen/coliveir/Miguel/Figures/update/Fig2C_Markers_2.txt")
genes_3 <- readLines("/projects/cangen/coliveir/Miguel/Figures/update/Fig2C_Markers_3.txt")
genes_4 <- readLines("/projects/cangen/coliveir/Miguel/Figures/update/Fig2C_Markers_4.txt")
genes_5 <- readLines("/projects/cangen/coliveir/Miguel/Figures/update/Fig2C_Markers_5.txt")

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

ggarrange(plotlist = plotlist)



##################################
########### Figure S2H ###########
##################################
data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
infos <- read.table("/projects/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- new_cluster

count_raw <- data@assays$SCT@counts[, rownames(data@meta.data)]
count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)

ort <- read.table("/projects/cangen/coliveir/cellphonedb/Orthologs_human_mouse.txt", sep = ",", header = T)
mat <- merge(ort, count_norm, by.x = "Mouse.gene.name", by.y = "row.names")
mat$Mouse.gene.name <- mat$Gene.stable.ID <- mat$Mouse.gene.stable.ID <- NULL
colnames(mat)[1] <- "Gene"

anno <- data.frame(Cell = names(Idents(data)), 
                   cluster = Idents(data), row.names = NULL)

write.table(mat, "/projects/cangen/coliveir/cellphonedb/counts.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(anno, "/projects/cangen/coliveir/cellphonedb/meta.txt", row.names = F, col.names = T, sep = "\t", quote = F)


#Python
source /home/coliveir/cpdb-venv/bin/activate
cd /projects/cangen/coliveir/cellphonedb

cellphonedb method statistical_analysis meta.txt counts.txt --threads=16 --counts-data=gene_name --pvalue=0.05 --iterations=1000