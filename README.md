# AdipoSNAP: Adipose Single-Nuclei Analysis Pipeline <img src="Logo.png" align="right" width="300" />


# Table of Contents
- [What is this?](#what-is-this)
- [Workflow](#workflow)
- [How can I use this data, and where can I find it?](#how-can-i-use-this-data-and-where-can-i-find-it)
	- [Downloading Data files](#downloading-data-files)
- [Analysis and visualization programs](#analysis-and-visualization-programs)
	- [R and R's integrated developmental environment RStudio:](#r-and-rs-integrated-developmental-environment-rstudio)
	- [scRNAseq analysis pipeline SEURAT developed by the Satija lab:](#scrnaseq-analysis-pipeline-seurat-developed-by-the-satija-lab)
	- [Pseudotemporal gene expression analysis using Monocle developed by the Trapnell Lab:](#pseudotemporal-gene-expression-analysis-using-monocle-developed-by-the-trapnell-lab)
	- [Cell type classification using Metacell:](#cell-type-classification-using-metacell)
	- [Cellular component prediction using cellphonedb:](#cellular-component-prediction-using-cellphonedb)
	- [Finding the optimal number of clusters using SCCAF:](#finding-the-optimal-number-of-clusters-using-sccaf)
- [Setting up the right environment](#setting-up-the-right-environment)
- [Computing environment](#computing-environment)
	- [Hardware and OS environment for running Seurat and Monocle](#hardware-and-os-environment-for-running-seurat-and-monocle)
- [Citation](#citation)
- [Contributors](#contributors)
- [Acknowledgements](#acknowledgements)



## What is this?
This repository contains coding scripts utilized for the analysis performed in the "Single-nuclei reconstruction of the adipose tissue using AdipoSNAP (Adipose Single-Nuclei Analysis Pipeline) reveals the mature adipocyte landscape underlying thermogenic response" publication [(XXX)](XXX). The purpose of providing the code here is to allow for transparency and robust data-analysis reproducibility. The methodology has already been described extensively in the manuscript. However, this analysis relies heavily on powerful scRNAseq analysis algorithms like [Seurat](https://satijalab.org/seurat/) [(Butler et al., 2018: Nature Biotechnology;](https://www.nature.com/articles/nbt.4096) [Stuart et al., 2018: Cell)](https://www.sciencedirect.com/science/article/pii/S0092867419305598?via%3Dihub), [SCCAF](https://github.com/SCCAF/sccaf), [Metacell](https://tanaylab.github.io/metacell/) [(Baran et al., 2019: Genome Biology)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1812-2) and [cellphonedb](https://www.cellphonedb.org) [(Efremova et al., 2020: Nature;](https://www.nature.com/articles/s41596-020-0292-x) [Vento-Tormo et al., 2018: Nature)](https://www.nature.com/articles/s41586-018-0698-6) (for a complete list of dependencies and code utilized see analysis & visualization programs).


## Workflow
AdipoSNAP contains three main steps:
1. **Dataset:**
2. **Clustering:**
    1. Overclustering
    2. Optimal number of clusters
    3. Cell type identification
    4. Markers expression
3. **Main analysis:**
    1. Differential Expression
    2. Functional Enrichment
    3. Transdifferentiation
    4. Cellular component prediction
<img src="Workflow.png" align="center">


## How can I use this data, and where can I find it?
### Downloading Data files
Public data files utilized in this analysis have been downloaded from [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/), gene expression data repository at the NIH. Data are part of the GSE133486 high-thoroughput sequencing repository and can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133486). The Cellranger output files were renamed to 'matrix.mtx.gz', 'barcodes.tsv.gz' and 'features.tsv.gz' to allow Seurat to read these files.


## Analysis and visualization programs
#### R and R's integrated developmental environment RStudio:
1. [R v4.0.2 (x64 bit)](https://cran.r-project.org/bin/macosx/base/)
2. [RStudio v1.3.1073 (x64 bit)](https://www.rstudio.com/products/rstudio/download/)
4. [Tutorial for R](https://cran.r-project.org/doc/manuals/r-release/R-intro.html)
5. [Tutorial for RStudio](https://resources.rstudio.com/)
#### scRNAseq analysis pipeline SEURAT developed by the Satija lab:
1. [Source code for Seurat v3.2.2](https://cran.r-project.org/web/packages/Seurat/index.html)
2. [Tutorials for Seurat](https://satijalab.org/seurat/)
#### Pseudotemporal gene expression analysis using Monocle developed by the Trapnell Lab:
1. [Source code for Monocle v2.16.0](https://bioconductor.org/packages/release/bioc/html/monocle.html)
2. [Tutorial for Monocle](http://cole-trapnell-lab.github.io/monocle-release/docs/#constructing-single-cell-trajectories)
#### Cell type classification using Metacell:
1. [Source code for metacell v0.3.41](https://github.com/tanaylab/metacell/releases/tag/v0.3.41)
2. [Tutorial for metacell](https://tanaylab.github.io/metacell/)
#### Cellular component prediction using cellphonedb:
1. [Source code for cellphonedb v2.1.4](https://github.com/Teichlab/cellphonedb/releases/tag/v2.1.4)
2. [Tutorial for cellphonedb](https://github.com/Teichlab/cellphonedb)
#### Finding the optimal number of clusters using SCCAF:
1. [Source code for SCCAF v0.09](https://github.com/SCCAF/sccaf/releases/tag/0.09)
2. [Tutorial for SCCAF](https://github.com/SCCAF/sccaf)


## Setting up the right environment
1. Install R and Rstudio
2. Once you have installed R and RStudio, run the [**1_environmetn_setup.R**](https://github.com/cbiagii/AdipoSNAP/blob/master/1_environment_setup.R) script.
3. Together with *Seurat*, a conda environment called ***r-reticulate*** will be installed. We will install the *SCCAF* and *cellphonedb* modules within this environment so that we can run the Python code inside R using the ***reticulate*** package previously installed. So, to check the installed environment just type the following commands in a new terminal:
```
conda env list
```
4. After checking the full name of the environment mentioned above, we will load it:
```
conda activate /Users/biagi/Library/r-miniconda/envs/r-reticulate
```
5. Finally, we will install the modules with the following commands:
```
pip install sccaf
pip install cellphonedb
```
6. If you run into problems, please open a new issue, you can do this by going to 'issues' and clicking on the 'new issue' icon. We will help you replicate our analysis! Do not fear single cell analysis!


## Computing environment
### Hardware and OS environment for running Seurat and Monocle
#### Environment 1
1. Processor: Intel Sandy Bridge E5-2670 (16cores x 16 threads)
2. RAM: 25GB
3. OS: CentOS 6.5
#### Environment 2
1. Hardware integrated into the [Pegasus Supercomputing array](http://ccs.miami.edu/ac/service/pegasus/) at the [University of Miami](https://welcome.miami.edu/)


# Citation
Qadir, M.M.F., Alvarez-Cubela, S., Klein, D., Van Dijk, J., Anquela, R.M., Lanzoni, G., Sadiq, S., Moreno-Hernandez, Y.B., Navarro-Rubio, B., Garcia, M.T., Diaz, A., Johnson, K., Sant, D., Ricordi, C., Griswold, T., Pastori, R.L., Dominguez-bendala, J. (2020) Proceedings of the National Academy of Sciences. Single cell resolution analysis of the human pancreatic ductal progenitor cell niche. Apr 2020, 201918314; DOI: 10.1073/pnas.1918314117


# Contributors
1. Carlos Alberto Oliveira de Biagi Junior - [Github](https://github.com/cbiagii) - [University of Sao Paulo](https://www.fmrp.usp.br/en/) - to contact please [Email](mailto:biagi@usp.br)
2. Sarah Santiloni Cury - [Sao Paulo State University](https://www.international.unesp.br/) - to contact please [Email](mailto:santiloni.cury@unesp.br)
3. Miguel Luiz Batista Junior [Boston University School of Medicine](https://www.bumc.bu.edu/busm/) - to contact please [Email](mailto:migueljr4@me.com)


# Acknowledgements
1. Diabetes Research Institute Foundation (DRIF)
2. The Inserra family
3. The Fred and Mabel R. Parks Foundation
4. The Tonkinson Foundation
5. ADA Grant #1-19-ICTS-078
6. NIH Grant #1R43DK105655-01
7. NIH Grant #2R44DK105655-02
8. NIH/NIDDK HIRN Grant #U01DK120393 (These studies are part of this grant)
9. IIE Fulbright

## If this was useful to you please dont forget to cite, star and fork this repository.