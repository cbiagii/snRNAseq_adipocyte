# DESCRIPTION ####
# The following code ASSUMES all dependencies in R have been installed See file(s): 
# 1_environment_setup.R
# The purpose of this code is to generate de Seurat object reproducing the analysis performed by Rajbhandari et al., 2020

# LOAD LIBRARIES ####
# Run the following code once you have Seurat installed
library(Seurat)
library(future)
library(metacell)
library(ggpubr)

## Convert an Seurat object to 10x format
library(Seurat)
library(DropletUtils)
data <- readRDS("/Users/biagi/PhD/AdipoSNAP/output/10x/10x_Processed.rds")
write10xCounts(x = data@assays$RNA@counts, path = "/Users/biagi/PhD/AdipoSNAP/output/10x/metacell/data")


if(!dir.exists("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell/db")) dir.create("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell/db")
scdb_init("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell/db", force_reinit=T)

mcell_import_scmat_10x("test1", base_dir="/Users/biagi/PhD/AdipoSNAP/output/10x/metacell/data")

mat = scdb_mat("test1")

if(!dir.exists("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell/figs")) dir.create("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell/figs")
scfigs_init("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell/figs")

mcell_plot_umis_per_cell("test1")
mat = scdb_mat("test1")
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
ig_genes = c(grep("^Igj", nms, v=T), 
             grep("^Igh",nms,v=T),
             grep("^Igk", nms, v=T), 
             grep("^Igl", nms, v=T))
bad_genes = unique(c(grep("^mt-", nms, v=T), grep("^mtmr", nms, v=T), grep("^Mtnd", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))
#bad_genes

mcell_mat_ignore_genes(new_mat_id="test1", mat_id="test1", bad_genes, reverse=F)


mcell_mat_ignore_small_cells("test1", "test1", 700)

mcell_add_gene_stat(gstat_id="test1", mat_id="test1", force=T)

mcell_gset_filter_varmean(gset_id="test_feats", gstat_id="test1", T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id = "test_feats", gstat_id="test1", T_tot=100, T_top3=2)

mcell_plot_gstats(gstat_id="test1", gset_id="test_feats")

mcell_add_cgraph_from_mat_bknn(mat_id="test1", 
                               gset_id = "test_feats", 
                               graph_id="test_graph",
                               K=150,
                               dsamp=F)

mcell_coclust_from_graph_resamp(
  coc_id="test_coc500", 
  graph_id="test_graph",
  min_mc_size=20, 
  p_resamp=0.75, n_resamp=1000)

mcell_mc_from_coclust_balanced(
  coc_id="test_coc500", 
  mat_id= "test1",
  mc_id= "test_mc", 
  K=20, min_mc_size=20, alpha=2)

mcell_plot_outlier_heatmap(mc_id="test_mc", mat_id = "test1", T_lfc=3)
mcell_mc_split_filt(new_mc_id="test_mc_f", 
                    mc_id="test_mc", 
                    mat_id="test1",
                    T_lfc=3, plot_mats=F)

mcell_gset_from_mc_markers(gset_id="test_markers", mc_id="test_mc_f")

load("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell/db/mc.test_mc_f.Rda")
lfp = log2(object@mc_fp)


########### 1st Round ###########
marks_colors <- NULL
marks_colors <- rbind(marks_colors, c("Adipocyte", "Adrb3", "blue", 1, 2))
marks_colors <- rbind(marks_colors, c("Endothelial", "Pecam1", "green", 1, 1))
marks_colors <- rbind(marks_colors, c("Immune_1", "Ptprc", "#ff748c", 1, 0.5))
marks_colors <- rbind(marks_colors, c("Immune_2", "Cd19", "#ff8fa3", 1, 0.5))
marks_colors <- rbind(marks_colors, c("Progenitor_1", "Cd34", "#ffa500", 1, 2))
marks_colors <- rbind(marks_colors, c("Progenitor_2", "Pdgfra", "#ffb732", 1, 2))
marks_colors <- as.data.frame(marks_colors)
colnames(marks_colors) <- c("group", "gene", "color", "priority", "T_fold")
marks_colors$priority <- as.integer(marks_colors$priority)
marks_colors$T_fold <- as.numeric(marks_colors$T_fold)

mc = scdb_mc("test_mc_f")
gene_folds = mc@mc_fp

load("/Users/biagi/PhD/AdipoSNAP/output/10x/metacell/db/gset.test_markers.Rda")
gset <- object
good_marks = intersect(names(gset@gene_set), rownames(mc@mc_fp))
mc_ord = 1:ncol(mc@mc_fp)

mat = log2(gene_folds[good_marks, mc_ord])
mat = pmax(pmin(mat, 3), -3)

mat_A <- mat[, which(mc@colors == "blue")]
mat_A <- mat_A[rowSums(mat_A) > quantile(rowSums(mat_A), 0.9), ]
rowMeans(mat_A[names(head(sort(rowSums(mat_A), decreasing = T), 5)), ])

mat_E <- mat[, which(mc@colors == "green")]
mat_E <- mat_E[rowSums(mat_E) > quantile(rowSums(mat_E), 0.9), ]
rowMeans(mat_E[names(head(sort(rowSums(mat_E), decreasing = T), 5)), ])

mat_I <- mat[, which(mc@colors %in% c("#ff748c", "#ff8fa3"))]
mat_I <- mat_I[rowSums(mat_I) > quantile(rowSums(mat_I), 0.9), ]
rowMeans(mat_I[names(head(sort(rowSums(mat_I), decreasing = T), 5)), ])

mat_P <- mat[, which(mc@colors %in% c("#ffa500", "#ffb732"))]
mat_P <- mat_P[rowSums(mat_P) > quantile(rowSums(mat_P), 0.9), ]
rowMeans(mat_P[names(head(sort(rowSums(mat_P), decreasing = T), 5)), ])


items <- list(Adipocytes = mat_A, Endothelials = mat_E, Immunes = mat_I, Progenitors = mat_P)
plot_list=list()
for (i in 1:length(items)){
  x <- pheatmap(items[[i]], fontsize = 6, main = names(items)[i], legend = FALSE, treeheight_row = 0, treeheight_col = 0)
  plot_list[[i]] = x[[4]]     ##to save each plot into a list. note the [[4]]
}

gggpubr(plotlist = plot_list, ncol = 2)

########### 2st Round ###########
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
marks_colors <- rbind(marks_colors, c("Immune_1", "Zeb2", "#ff7f7f", 1, 0.85))
marks_colors <- rbind(marks_colors, c("Immune_2", "Trps1", "#ff6666", 1, 0.85))
marks_colors <- rbind(marks_colors, c("Immune_3", "Runx1", "#ff4c4c", 1, 0.85))
marks_colors <- rbind(marks_colors, c("Immune_4", "Ptprc", "#ff3232", 1, 0.85))
marks_colors <- rbind(marks_colors, c("Immune_5", "Adap2", "#ff1919", 1, 0.85))
marks_colors <- rbind(marks_colors, c("Progenitor_1", "Dcn", "#ffff4d", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Progenitor_2", "Celf2", "#ffff33", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Progenitor_3", "Meg3", "#ffff1a", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Progenitor_4", "Col1a2", "#ffff00", 1, 2.5))
marks_colors <- rbind(marks_colors, c("Progenitor_5", "Col3a1", "#e6e600", 1, 2.5))
marks_colors <- as.data.frame(marks_colors)
colnames(marks_colors) <- c("group", "gene", "color", "priority", "T_fold")
marks_colors$priority <- as.integer(marks_colors$priority)
marks_colors$T_fold <- as.numeric(marks_colors$T_fold)



mc_colorize("test_mc_f", marker_colors=marks_colors)

mc = scdb_mc("test_mc_f")

mcell_mc_plot_marks(mc_id="test_mc_f", gset_id="test_markers", mat_id="test1", plot_cells = F)

lfp = log2(mc@mc_fp)

mcell_mc2d_force_knn(mc2d_id="test_2dproj",mc_id="test_mc_f", graph_id="test_graph")
tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id="test_2dproj")

