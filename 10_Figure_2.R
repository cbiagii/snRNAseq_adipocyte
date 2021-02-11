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