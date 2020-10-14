# DESCRIPTION ####
# The following code will install all dependencies in R allowing for analysis to be performed on AdipoSNAP paper

# INSTALLATION ####
# Install additional CRAN packages
CRAN_packages <- c('BiocManager', 'future', 'ggpubr', 'dplyr', 'ggplot2', 
                   'enrichR', 'forcats', 'ggtext', 'ggrepel', 'RColorBrewer', 'ggbeeswarm', 'circlize', 
                   'pheatmap', 'pbapply', 'gam', 'UpSetR', 'tidyr', 'tibble', 'tidyverse')
new.packages <- CRAN_packages[!(CRAN_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Install additional Bioconductor packages
Bioconductor_packages <- c('slingshot', 'SingleCellExperiment', 'fgsea', 'ComplexHeatmap', 'monocle', 'tradeSeq', 'SummarizedExperiment')
new.packages <- Bioconductor_packages[!(Bioconductor_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

# Install reticulate package
install.packages('reticulate')

# Install Seurat and SeuratWrappers packages
install.packages('Seurat')
install.packages('SeuratWrappers')

# Install metacell package
if(length('metacell'[!('metacell' %in% installed.packages()[,"Package"])])) BiocManager::install("tanaylab/metacell")


# LOAD LIBRARIES ####
# Restart Rstudio or R

library(BiocManager)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(enrichR)
library(fgsea)
library(forcats)
library(future)
library(gam)
library(ggbeeswarm)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggtext)
library(metacell)
library(monocle)
library(pbapply)
library(pheatmap)
library(RColorBrewer)
library(Seurat)
library(SeuratWrappers)
library(SingleCellExperiment)
library(slingshot)
library(SummarizedExperiment)
library(tibble)
library(tidyr)
library(tidyverse)
library(tradeSeq)
library(UpSetR)

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle")
packageVersion("metacell")