## Loading R packages
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(ggpubr)
library(fgsea)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(monocle)
library(scales)
library(viridis)

source('/Users/biagi/cangen/coliveir/plotMonocle.R')


#################################
########### Figure 3A ###########
#################################
data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
Idents(data) <- as.factor(data$timpoint)
data <- subset(data, cells = names(which(data$timpoint == '4day' | data$timpoint == 'RT' | data$timpoint == 'CL')))
data$timpoint <- gsub('4day', 'Cold', data$timpoint)

TSNEPlot(data, pt.size = 1, group.by = 'timpoint') +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  theme_classic() + labs(color = "Cluster") +
  theme(legend.position="bottom")



#################################
########### Figure 3B ###########
#################################
data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
Idents(data) <- as.factor(data$timpoint)
data <- subset(data, cells = names(which(data$timpoint == '4day' | data$timpoint == 'RT' | data$timpoint == 'CL')))

data <- SeuratWrappers::RunALRA(data)
data <- ScaleData(data)

genes <- c('Lipe', 'Adrb3', 'Pnpla2', 'Plin1', 'Fasn', 'Acadm', 'Gk', 'Adipoq', 'Retn', 'Cidec', 'Ppara', 'Ucp1', 'Dio2', 'Prdm16', 'Elovl3')
genes <- genes[genes %in% rownames(data)]

pt <- FeaturePlot(data, genes,
                  cols = c("grey", 'red'), reduction = 'tsne', pt.size = 0.1, combine = F)
pt <- lapply(pt, function(x) {
  x + theme_classic() + xlab("t-SNE 1") + ylab("t-SNE 2") +
    theme(plot.title = element_text(hjust = 0.5, face = "italic", size = 10), legend.position = "none")
})

ggarrange(plotlist = pt, ncol = 5, nrow = 3, common.legend = T)



#################################
########### Figure 3C ###########
#################################
ort <- read.table("/projects/cangen/coliveir/cellphonedb/Orthologs_human_mouse.txt", sep = ",", header = T)
pathways <- gmtPathways('/projects/cangen/coliveir/Miguel/msigdb.v7.1.symbols.gmt')

#ColdxRT
comp1 <- read.table("/projects/cangen/coliveir/Miguel/Figures/update/Diff_Adipocytes_ColdxRT_ALL.txt")
comp1 <- subset(comp1, p_val_adj < 0.05)
comp1 <- comp1[order(comp1$avg_logFC, decreasing = T), ]
comp1 <- rbind(head(comp1, 50), tail(comp1, 50))
comp1 <- merge(ort, comp1, by.x = "Mouse.gene.name", by.y = "row.names")
comp1$Mouse.gene.name <- comp1$Gene.stable.ID <- comp1$Mouse.gene.stable.ID <- NULL
colnames(comp1)[1] <- "Gene"
ranks <- structure(comp1$avg_logFC, names = comp1$Gene)
ranks <- sort(ranks)
fgseaRes1 <- fgsea(pathways, ranks, minSize=5)
fgseaRes1 <- subset(fgseaRes1, padj < 0.05)
fgseaRes1 <- fgseaRes1[order(fgseaRes1$NES, decreasing = T), ]
genes1 <- fgseaRes1$leadingEdge
for (i in 1:length(genes1)) {
  fgseaRes1$genes[i] <- paste(genes1[[i]], collapse = ';')
}
fgseaRes1$leadingEdge <- NULL
enrich_ColdxRT <- read.table('/projects/cangen/coliveir/Miguel/Figures/update/fgsea_ColdxRT.txt', header = T)
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

#CLxRT
comp2 <- read.table("/projects/cangen/coliveir/Miguel/Figures/update/Diff_Adipocytes_CLxRT_ALL.txt")
comp2 <- subset(comp2, p_val_adj < 0.05)
comp2 <- comp2[order(comp2$avg_logFC, decreasing = T), ]
comp2 <- rbind(head(comp2, 50), tail(comp2, 50))
comp2 <- merge(ort, comp2, by.x = "Mouse.gene.name", by.y = "row.names")
comp2$Mouse.gene.name <- comp2$Gene.stable.ID <- comp2$Mouse.gene.stable.ID <- NULL
colnames(comp2)[1] <- "Gene"
ranks <- structure(comp2$avg_logFC, names = comp2$Gene)
ranks <- sort(ranks)
fgseaRes2 <- fgsea(pathways, ranks, minSize=5)
fgseaRes2 <- subset(fgseaRes2, padj < 0.05)
fgseaRes2 <- fgseaRes2[order(fgseaRes2$NES, decreasing = T), ]
genes2 <- fgseaRes2$leadingEdge
for (i in 1:length(genes2)) {
  fgseaRes2$genes[i] <- paste(genes2[[i]], collapse = ';')
}
fgseaRes2$leadingEdge <- NULL
enrich_CLxRT <- read.table('/projects/cangen/coliveir/Miguel/Figures/update/fgsea_CLxRT.txt', header = T)
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
data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
data$timpoint <- gsub('4day', 'Cold', data$timpoint)

infos <- read.table("/projects/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")

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

### CL
ggplot() +
  geom_point(data = df, aes(x = tSNE_1, y = tSNE_2, color = timpoint), size = 0.5) +
  scale_color_manual(values = rep("gray75", 5)) +
  guides(color = guide_legend(override.aes = list(size=7))) +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  labs(color = '') + ggtitle(paste0('CL (', nrow(anno_CL), ' cells)')) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none", plot.title = element_text(hjust = 0.5, size=8)) +
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

### RT
ggplot() +
  geom_point(data = df, aes(x = tSNE_1, y = tSNE_2, color = timpoint), size = 0.5) +
  scale_color_manual(values = rep("gray75", 5)) +
  guides(color = guide_legend(override.aes = list(size=7))) +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  labs(color = '') + ggtitle(paste0('RT (', nrow(anno_RT), ' cells)')) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none", plot.title = element_text(hjust = 0.5, size=8)) +
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

### Cold
ggplot() +
  geom_point(data = df, aes(x = tSNE_1, y = tSNE_2, color = timpoint), size = 0.5) +
  scale_color_manual(values = rep("gray75", 5)) +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  guides(color = guide_legend(override.aes = list(size=7))) +
  labs(color = '') + ggtitle(paste0('Cold (', nrow(anno_Cold), ' cells)')) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", plot.title = element_text(hjust = 0.5, size=8)) +
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
data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
infos <- read.table("/projects/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")

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
           "Daglb", "Slc2a4", "Gpd1", "Pcx", "Abhd5", "Adrb3", "Lipe", "Pck1", "Pnpla2")
genes <- genes[genes %in% rownames(data)]

mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)
DoHeatmap(data, features = genes, angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='UCP1 Expression')



#################################
########### Figure 3G ###########
#################################
data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
infos <- read.table("/projects/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")

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
           'Acaa2', 'Sdhc', 'Fabp3', 'Ppara', 'Slc2a4')
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



##################################
######## Figure 3H and 3I ########
##################################
## CL
data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
data$timpoint <- gsub('4day', 'Cold', data$timpoint)

infos <- read.table("/projects/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")
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

plot_cell_trajectory(cds, color_by = "Pseudotime")

plot_cell_trajectory(cds, color_by = "cluster") + 
  scale_color_manual(values = c("Ad1" = "#11c78b", "Ad2" = "#800080", "Ad3" = "#e57400", "Ad4" = "#0000FF", "Ad5" = "#dfdf0d"))

plotMonocle(cds, c('Pparg', 'Creb1', 'Atf2', 'Egr2', 'Dbp', 'Xbp1', 'Clock', 'Zbtb43', 'Zbtb7a'))

## Cold
data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
data$timpoint <- gsub('4day', 'Cold', data$timpoint)

infos <- read.table("/projects/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")
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

plot_cell_trajectory(cds, color_by = "Pseudotime")

plot_cell_trajectory(cds, color_by = "cluster") + 
  scale_color_manual(values = c("Ad1" = "#11c78b", "Ad2" = "#800080", "Ad3" = "#e57400", "Ad4" = "#0000FF", "Ad5" = "#dfdf0d"))

plotMonocle(cds, c('Pparg', 'Creb1', 'Atf2', 'Srebf1', 'Cebpa', 'Egr1', 'Foxn3', 'Esrra', 'Gtf2ird1'))