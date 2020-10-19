---
title: "1 - Setting up R environment"
diagram: yes
date: '2020-10-18'
math: yes
---

In this section we will show you how to configure your R environment by installing all the necessary dependencies to run the following codes. The dependencies. The dependencies include packages from the CRAN and Bioconductor repositories. To install, just use the following commands:

```r
# Install additional CRAN packages
CRAN_packages <- c('BiocManager', 'circlize', 'dplyr', 'enrichR', 'future', 'ggplot2', 
                   'ggpubr', 'ggrepel', 'ggtext', 'glue', 'pheatmap', 'RColorBrewer', 
                   'scales', 'tibble', 'tidyr', 'UpSetR', 'viridis')
new.packages <- CRAN_packages[!(CRAN_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages)
}
  
# Install additional Bioconductor packages
Bioconductor_packages <- c('ComplexHeatmap', 'DropletUtils', 'fgsea', 'monocle')
new.packages <- Bioconductor_packages[!(Bioconductor_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  BiocManager::install(new.packages)
}
  
# Install reticulate package
install.packages('reticulate')

# Install Seurat and SeuratWrappers packages
install.packages('Seurat')
remotes::install_github('satijalab/seurat-wrappers')

# Install metacell package
if(length('metacell'[!('metacell' %in% installed.packages()[,"Package"])])) {
  BiocManager::install("tanaylab/metacell")
}
```