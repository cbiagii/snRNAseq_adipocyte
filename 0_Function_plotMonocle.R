plotMonocle <- function(cds, gene) {
  if (sum(gene %in% rownames(cds)) == 0) {
    stop('None gene found in dataset')
  }
  
  if (sum(gene %in% rownames(cds)) != length(gene)) {
    gene <- gene[gene %in% rownames(cds)]
  }
  
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(monocle))
  suppressPackageStartupMessages(library(scales))
  suppressPackageStartupMessages(library(ggpubr))
  suppressPackageStartupMessages(library(viridis))
  
  return_rotation_mat <- function(theta) {
    theta <- theta/180 * pi
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 
           nrow = 2)
  }
  monocle_theme_opts <- function() {
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
      theme(panel.border = element_blank()) +
      theme(axis.line.x = element_line(size=0.25, color="black")) +
      theme(axis.line.y = element_line(size=0.25, color="black")) +
      theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
      theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
      theme(panel.background = element_rect(fill='white')) +
      theme(legend.key=element_blank())
  }
  
  
  tmp <- cds@assayData$exprs[gene, ]
  if (length(gene) == 1) {
    tmp <- rescale(tmp, to = c(-2,2))
    cds[[gene]] <- tmp
  } else {
    tmp <- apply(tmp, 1, function(x) rescale(x, to = c(-2,2)))
    for (i in 1:ncol(tmp)) {
      cds[[colnames(tmp)[i]]] <- tmp[,i]
    }
  }
  
  pt <- plot_cell_trajectory(cds, color_by = gene[1])
  
  reduced_dim_coords <- reducedDimK(cds)
  ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>% 
    select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>% 
    mutate(sample_name = rownames(.), sample_state = rownames(.))
  dp_mst <- minSpanningTree(cds)
  
  edge_df <- dp_mst %>% igraph::as_data_frame() %>% select_(source = "from", target = "to") %>% left_join(ica_space_df %>% select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1", source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>% left_join(ica_space_df %>% select_(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1", target_prin_graph_dim_2 = "prin_graph_dim_2"), by = "target")
  
  rot_mat <- return_rotation_mat(0)
  cn2 <- c("source_prin_graph_dim_1", "source_prin_graph_dim_2")
  cn3 <- c("target_prin_graph_dim_1", "target_prin_graph_dim_2")
  edge_df[, cn2] <- as.matrix(edge_df[, cn2]) %*% t(rot_mat)
  edge_df[, cn3] <- as.matrix(edge_df[, cn3]) %*% t(rot_mat)
  
  
  data_df <- pt$data
  
  g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
  
  g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                   y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                   yend = "target_prin_graph_dim_2"), size = 0.75, 
                        linetype = "solid", na.rm = TRUE, data = edge_df)
  
  mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
  branch_point_df <- ica_space_df %>% slice(match(mst_branch_nodes, sample_name)) %>% mutate(branch_point_idx = seq_len(n()))
  g <- g + geom_point(aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, branch_point_df) + 
    geom_text(aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2", label = "branch_point_idx"), size = 4, color = "white", na.rm = TRUE, branch_point_df)
  
  g <- g + monocle_theme_opts() + xlab("Component 1") + ylab("Component 2") + 
    theme(legend.position = "top", legend.key.height = grid::unit(0.35, "in")) + theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white"))
  
  plotlist <- list()
  for (i in 1:length(gene)) {
    plotlist[[i]] <- g + geom_point(data = data_df[which(data_df[[gene[i]]] < 0), ], aes_string(color = paste0('`', gene[i], '`')), size = I(1), na.rm = TRUE) + 
      geom_point(data = data_df[which(data_df[[gene[i]]] > 0), ], aes_string(color = paste0('`', gene[i], '`')), size = I(1.5), na.rm = TRUE) + 
      scale_color_viridis(option = 'C', discrete = F, end = 0.9) + ggtitle(gene[i]) + 
      theme(plot.title = element_text(hjust = 0.5)) + labs(color = "")
  }
  
  pt2 <- ggarrange(plotlist = plotlist, common.legend = T)
  
  return(pt2)
}