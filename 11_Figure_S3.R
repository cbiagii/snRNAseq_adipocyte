## Loading R packages
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(enrichR)
library(glue)
library(ggtext)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(monocle)
library(scales)
library(ggpubr)
library(viridis)

source('/Users/biagi/cangen/coliveir/plotMonocle.R')


##################################
########### Figure S3A ###########
##################################
data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
data$timpoint <- gsub('4day', 'Cold', data$timpoint)

Idents(data) <- as.factor(data$timpoint)

markers_ColdxRT <- FindMarkers(data, ident.1 = "Cold", ident.2 = "RT", logfc.threshold = 0)
markers_ColdxRT$gene <- rownames(markers_ColdxRT)

markers_CLxRT <- FindMarkers(data, ident.1 = "CL", ident.2 = "RT", logfc.threshold = 0)
markers_CLxRT$gene <- rownames(markers_CLxRT)

data <- SeuratWrappers::RunALRA(data)
data <- ScaleData(data)

data2 <- subset(data, idents = c("Cold", "RT", "CL"))

top1 <- markers_ColdxRT %>% top_n(n = 50, wt = avg_logFC)
top2 <- markers_CLxRT %>% top_n(n = 50, wt = avg_logFC)
topGenes <- rbind(top1, top2)

cls <- c("blue", "red", "green")
mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)

DoHeatmap(data2, features = topGenes$gene, group.colors = cls, angle = 0, size = 5, label = F) +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='Timepoint')



##################################
########### Figure S3B ###########
##################################
dbs <- c("KEGG_2019_Mouse", "WikiPathways_2019_Mouse", "Jensen_TISSUES", "GO_Biological_Process_2018")

genes_ColdxRT <- readLines("/projects/cangen/coliveir/Miguel/Figures/update/Diff_Adipocytes_ColdxRT.txt")
genes_CLxRT <- readLines("/projects/cangen/coliveir/Miguel/Figures/update/Diff_Adipocytes_CLxRT.txt")

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

paths <- paths[which(paths$Term %in% c("PPAR signaling pathway", "Propanoate metabolism", "Regulation of lipolysis in adipocytes", "Fatty acid biosynthesis", "Fatty acid degradation", "Valine, leucine and isoleucine degradation", "PPAR signaling pathway WP2316", "Triacylglyceride Synthesis WP386", "Subcutaneous adipose tissue", "3T3-L1 cell", "Adipocyte", "long-chain fatty acid transport", "fatty acid transmembrane transport", "regulation of sequestering of triglyceride", "intracellular lipid transport")), ]

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
  scale_colour_gradient(limits=c(0, max(paths$Combined.Score)+0.5), high = "#2b9348", low = "#eeef20") +
  xlab(NULL) + ylab(NULL) +
  ggtitle(NULL) + labs(color = "Combined Score", size = "-log10(padj)") +
  theme(axis.text.y = element_markdown(), axis.text=element_text(size=12))



##################################
########### Figure S3C ###########
##################################
data <- readRDS("/projects/cangen/coliveir/Miguel/output/10x/Adipocytes.rds")
infos <- read.table("/projects/cangen/coliveir/scRNA_output/SCCAF/AdipocytesOnly/results/obs.csv")

new_cluster <- infos$L1_result
names(new_cluster) <- rownames(infos)
new_cluster <- new_cluster + 1
new_cluster <- paste0("Ad", new_cluster)
new_cluster <- as.factor(new_cluster)
Idents(data) <- factor(data$timpoint)
data$cluster <- factor(new_cluster)

data <- SeuratWrappers::RunALRA(data)
data <- ScaleData(data, features = rownames(data))

##Fatty Acid Oxidation
genes1 <- c('Abcb11', 'Abcd1', 'Abcd2', 'Abcd3', 'Abcd4', 'Acaa1a', 'Acaa1b', 'Acaa2', 'Acacb', 'Acad11', 'Acadl', 'Acadm', 'Acads', 'Acadvl', 'Acat1', 'Acat2', 'Acat3', 'Acox1', 'Acox2', 'Acox3', 'Acoxl', 'Acsbg2', 'Acsl5', 'Adh4', 'Adh5', 'Adh7', 'Adipoq', 'Adipor1', 'Adipor2', 'Akt1', 'Akt2', 'Alox12', 'Appl2', 'Auh', 'Bdh2', 'C1qtnf2', 'C1qtnf9', 'Cd36', 'Cnr1', 'Cpt1a', 'Cpt2', 'Crat', 'Crot', 'Cygb', 'Cyp4v3', 'Cyp24a1', 'Dbi', 'Decr1', 'Dgat1', 'Dgat2', 'Echdc1', 'Echdc2', 'Echs1', 'Eci1', 'Eci2', 'Eci3', 'Ehhadh', 'Etfa', 'Etfb', 'Etfbkmt', 'Etfdh', 'Fabp1', 'Fabp3', 'Gcdh', 'Gm45753', 'Gm49387', 'Hacl1', 'Hadh', 'Hadha', 'Hadhb', 'Hao1', 'Hao2', 'Hsd17b4', 'Hsd17b10', 'Ilvbl', 'Irs1', 'Irs2', 'Ivd', 'Lep', 'Lonp2', 'Mapk14', 'Mfsd2a', 'Mir199a-2', 'Mir214', 'Mir696', 'Mlycd', 'Mtor', 'Nr4a3', 'Nucb2', 'Pdk4', 'Pex2', 'Pex5', 'Pex7', 'Pex13', 'Phyh', 'Plin5', 'Por', 'Ppara', 'Ppard', 'Pparg', 'Ppargc1a', 'Prkaa1', 'Scp2', 'Sesn2', 'Sirt4', 'Slc25a17', 'Slc27a2', 'Sox9', 'Twist1', 'Tysnd1')
genes1 <- genes1[genes1 %in% rownames(data)]
mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)
ht1 <- DoHeatmap(data, features = genes1, angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='Expression')
y <- ht1$data %>% drop_na()
x <- y %>% group_by(Identity) %>% select(Feature, Cell, Identity, Expression) %>%
  spread(key = Feature, value = Expression)
w <- y %>% select(Feature, Cell, Expression) %>%
  spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()
pt <- Heatmap(w, cluster_columns = F)
DoHeatmap(data, features = rownames(w)[row_order(pt)], angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='Expression') + ggtitle('Fatty Acid Oxidation') + 
  theme(plot.title = element_text(hjust = 0.5))

##Tricarboxylic Acid Cycle
genes2 <- c('4933405O20Rik', 'Aco1', 'Aco2', 'Cs', 'Csl', 'Dhtkd1', 'Dlat', 'Dlst', 'Fh1', 'Idh1', 'Idh2', 'Idh3a', 'Idh3b', 'Idh3g', 'Ireb2', 'Mdh1', 'Mdh1b', 'Mdh2', 'Ndufs4', 'Ogdh', 'Ogdhl', 'Pdha1', 'Pdha2', 'Pdhb', 'Sdha', 'Sdhaf2', 'Sdhb', 'Sdhc', 'Sdhd', 'Sucla2', 'Suclg1', 'Suclg2')
genes2 <- genes2[genes2 %in% rownames(data)]
mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)
ht1 <- DoHeatmap(data, features = genes2, angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='Expression')
y <- ht1$data %>% drop_na()
x <- y %>% group_by(Identity) %>% select(Feature, Cell, Identity, Expression) %>%
  spread(key = Feature, value = Expression)
w <- y %>% select(Feature, Cell, Expression) %>%
  spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()
pt <- Heatmap(w, cluster_columns = F)
DoHeatmap(data, features = rownames(w)[row_order(pt)], angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='Expression') + ggtitle('Tricarboxylic Acid Cycle') + 
  theme(plot.title = element_text(hjust = 0.5))

##Fat Acid Transport
genes3 <- c('Abcc1', 'Abcc2', 'Abcc4', 'Abcd1', 'Abcd2', 'Abcd3', 'Abcd4', 'Ace', 'Acsl1', 'Acsl3', 'Acsl4', 'Acsl5', 'Acsl6', 'Agtr2', 'Akt1', 'Akt2', 'Anxa1', 'Apoe', 'Atp5j', 'Avpr1b', 'Bdkrb2', 'Cd36', 'Cpt1b', 'Crot', 'Cyp4f18', 'Drd2', 'Drd3', 'Drd4', 'Edn1', 'Eprs', 'Erfe', 'Fabp1', 'Fabp2', 'Fabp3', 'Fabp4', 'Fabp5', 'Hnf1a', 'Hrh2', 'Il1a', 'Il1b', 'Irs2', 'Kiss1r', 'Lep', 'Lhcgr', 'Map2k6', 'Mapk9', 'Mfsd2a', 'Mif', 'Nmb', 'Nmur2', 'Nos2', 'Ntsr1', 'Oc90', 'Oxt', 'P2rx7', 'P2ry2', 'Pla2g1b', 'Pla2g2c', 'Pla2g2d', 'Pla2g2e', 'Pla2g2f', 'Pla2g3', 'Pla2g4a', 'Pla2g4f', 'Pla2g5', 'Pla2g6', 'Pla2g10', 'Pla2g12a', 'Pla2g12b', 'Pla2r1', 'Plin2', 'Pnpla8', 'Pparg', 'Ptges', 'Repin1', 'Rps6kb1', 'Slc2a1', 'Slc5a8', 'Slc22a22', 'Slc25a17', 'Slc27a1', 'Slc27a2', 'Slc27a3', 'Slc27a4', 'Slc27a5', 'Slc27a6', 'Slco2a1', 'Slco3a1', 'Spx', 'Sstr4', 'Syk', 'Thbs1', 'Tnfrsf11a', 'Tnfsf11')
genes3 <- genes3[genes3 %in% rownames(data)]
mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)
ht1 <- DoHeatmap(data, features = genes3, angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='Expression')
y <- ht1$data %>% drop_na()
x <- y %>% group_by(Identity) %>% select(Feature, Cell, Identity, Expression) %>%
  spread(key = Feature, value = Expression)
w <- y %>% select(Feature, Cell, Expression) %>%
  spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()
pt <- Heatmap(w, cluster_columns = F)
DoHeatmap(data, features = rownames(w)[row_order(pt)], angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='Expression') + ggtitle('Fat Acid Transport') + 
  theme(plot.title = element_text(hjust = 0.5))

##Triglyceride/Fatty Acid Cycle
genes4 <- c('Aadac', 'Abhd5', 'Daglb', 'Pnpla2', 'Pck1', 'Pcx' , 'Gpd1', 'Slc2a4', 'Lipe', 'Adrb3')
genes4 <- genes4[genes4 %in% rownames(data)]
mapal <- colorRampPalette(brewer.pal(11,"RdBu"))(256)
ht1 <- DoHeatmap(data, features = genes4, angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='Expression')
y <- ht1$data %>% drop_na()
x <- y %>% group_by(Identity) %>% select(Feature, Cell, Identity, Expression) %>%
  spread(key = Feature, value = Expression)
w <- y %>% select(Feature, Cell, Expression) %>%
  spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()
pt <- Heatmap(w, cluster_columns = F)
DoHeatmap(data, features = rownames(w)[row_order(pt)], angle = 0, size = 5, label = F, group.by = 'cluster') +
  scale_fill_gradientn(colours = rev(mapal)) +
  theme(axis.text=element_text(size=6)) +
  labs(color='Expression') + ggtitle('Triglyceride/Fatty Acid Cycle') + 
  theme(plot.title = element_text(hjust = 0.5))



##################################
########### Figure S3H ###########
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



##################################
########### Figure S3I ###########
##################################
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

plot_cell_trajectory(cds, color_by = "timpoint")

plot_cell_trajectory(cds, color_by="State", state_number_size = 1) + facet_wrap(~State)

plotMonocle(cds, c('Ucp1', 'Ppara', 'Dio2', 'Chst1', 'Plppr3', 'Nnat', 'Pim1', 'Adcy3', 'Cs'))



##################################
########### Figure S3J ###########
##################################
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

plot_cell_trajectory(cds, color_by = "timpoint")

plot_cell_trajectory(cds, color_by="State", state_number_size = 1) + facet_wrap(~State)

plotMonocle(cds, c('Ucp1', 'Ppara', 'Dio2', 'Ccdc80', 'Slc7a10', 'Tmem43', 'Adcy3', 'Perm1', 'Fabp3'))