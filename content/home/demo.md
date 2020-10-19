+++
# A Demo section created with the Blank widget.
# Any elements can be added in the body: https://sourcethemes.com/academic/docs/writing-markdown-latex/
# Add more sections by duplicating this file and customizing to your requirements.

widget = "blank"  # See https://sourcethemes.com/academic/docs/page-builder/
headless = true  # This file represents a page section.
active = true  # Activate this widget? true/false
weight = 15  # Order that this section will appear.

title = ""
subtitle = ""

[design]
  # Choose how many columns the section has. Valid values: 1 or 2.
  columns = "1"

[design.background]
  # Apply a background color, gradient, or image.
  #   Uncomment (by removing `#`) an option to apply it.
  #   Choose a light or dark text color by setting `text_color_light`.
  #   Any HTML color name or Hex value is valid.

  # Background color.
   #color = "#f3fdf3"
  
  # Background gradient.
   #gradient_start = "DeepSkyBlue"
   #gradient_end = "SkyBlue"
  
  # Background image.
  #image = "bubbles-wide.jpg"  # Name of image in `static/media/`.
  #image_darken = 0.3  # Darken the image? Range 0-1 where 0 is transparent and 1 is opaque.
  #image_size = "cover"  #  Options are `cover` (default), `contain`, or `actual` size.
  #image_position = "center"  # Options include `left`, `center` (default), or `right`.
  #image_parallax = true  # Use a fun parallax-like fixed background effect? true/false

  # Text color (true=light or false=dark).
  text_color_light = false

[design.spacing]
  # Customize the section spacing. Order is top, right, bottom, left.
  padding = ["20px", "0", "20px", "0"]

[advanced]
 # Custom CSS. 
 css_style = ""
 
 # CSS class.
 css_class = ""
+++


# AdipoSNAP: Adipose Single-Nuclei Analysis Pipeline <img src="Logo.png" align="right" width="300" />
## Table of Contents
- [What is this?](#what-is-this)
- [Workflow](#workflow)
- [How can I use this data, and where can I find it?](#how-can-i-use-this-data-and-where-can-i-find-it)
	- [Downloading Data files](#downloading-data-files)
- [Analysis and visualization programs](#analysis-and-visualization-programs)
	- [R and R's integrated developmental environment RStudio](#r-and-rs-integrated-developmental-environment-rstudio)
	- [scRNAseq analysis pipeline Seurat](#scrnaseq-analysis-pipeline-seurat)
	- [Pseudotemporal gene expression analysis using Monocle](#pseudotemporal-gene-expression-analysis-using-monocle)
	- [Cell type classification using Metacell](#cell-type-classification-using-metacell)
	- [Cellular component prediction using cellphonedb](#cellular-component-prediction-using-cellphonedb)
	- [Finding the optimal number of clusters using SCCAF](#finding-the-optimal-number-of-clusters-using-sccaf)
- [Setting up the right environment](#setting-up-the-right-environment)
- [Citation](#citation)
- [Contributors](#contributors)
- [Acknowledgements](#acknowledgements)


### What is this?
This repository contains coding scripts utilized for the analysis performed in the "Single-nuclei reconstruction of the adipose tissue using AdipoSNAP (Adipose Single-Nuclei Analysis Pipeline) reveals the mature adipocyte landscape underlying thermogenic response" publication [(XXX)](XXX). The purpose of providing the code here is to allow for transparency and robust data-analysis reproducibility. The methodology has already been described extensively in the manuscript. However, this analysis relies heavily on powerful scRNAseq analysis algorithms like [`Seurat`](https://satijalab.org/seurat/) [(Butler et al., 2018: Nature Biotechnology;](https://www.nature.com/articles/nbt.4096) [Stuart et al., 2018: Cell)](https://www.sciencedirect.com/science/article/pii/S0092867419305598?via%3Dihub), [`SCCAF`](https://github.com/SCCAF/sccaf), [`Metacell`](https://tanaylab.github.io/metacell/) [(Baran et al., 2019: Genome Biology)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1812-2) and [`cellphonedb`](https://www.cellphonedb.org) [(Efremova et al., 2020: Nature;](https://www.nature.com/articles/s41596-020-0292-x) [Vento-Tormo et al., 2018: Nature)](https://www.nature.com/articles/s41586-018-0698-6) (for a complete list of dependencies and code utilized see analysis & visualization programs).


### Workflow
AdipoSNAP contains three main steps:
1. **`Dataset`**
2. **`Clustering`**
    1. Overclustering;
    2. Optimal number of clusters;
    3. Cell type identification;
    4. Markers expression.
3. **`Main analysis`**
    1. Differential Expression;
    2. Functional Enrichment;
    3. Transdifferentiation;
    4. Cellular component prediction.

<img src="Workflow.png" align="center">


### How can I use this data, and where can I find it?
#### Downloading Data files
Public data files utilized in this analysis have been downloaded from [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/), gene expression data repository at the NIH. Data are part of the GSE133486 high-thoroughput sequencing repository and can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133486). The Cellranger output files were renamed to 'matrix.mtx.gz', 'barcodes.tsv.gz' and 'features.tsv.gz' to allow Seurat to read these files.


### Analysis and visualization programs
#### R and R's integrated developmental environment RStudio:
1. [R v4.0.2 (x64 bit)](https://cran.r-project.org/bin/macosx/base/)
2. [RStudio v1.3.1073 (x64 bit)](https://www.rstudio.com/products/rstudio/download/)
4. [Tutorial for R](https://cran.r-project.org/doc/manuals/r-release/R-intro.html)
5. [Tutorial for RStudio](https://resources.rstudio.com/)
#### scRNAseq analysis pipeline Seurat:
1. [Source code for Seurat v3.2.2](https://cran.r-project.org/web/packages/Seurat/index.html)
2. [Tutorials for Seurat](https://satijalab.org/seurat/)
#### Pseudotemporal gene expression analysis using Monocle:
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


### Setting up the right environment
1. Install R and Rstudio.
2. Once you have installed R and RStudio, run the [**setup script**](https://cbiagii.github.io/AdipoSNAP/post/1_setup/) to install R packages dependencies.
3. Together with *Seurat*, a conda environment called ***r-reticulate*** will be installed. We will install the `anndata`, `scanpy`, `igraph` and `louvain` modules within this environment so that we can run the Python code inside R using the ***reticulate*** package previously installed. So, to check the installed environment full name just type the following commands in a new terminal:
```
conda env list
```
4. After checking the full name of the environment mentioned above, we will load it (replace the path below with the similar one shown on your terminal):
```
conda activate /Users/biagi/Library/r-miniconda/envs/r-reticulate
```
5. Then, we will install the modules with the following commands:
```
pip install anndata
pip install scanpy
pip install python-igraph
pip install louvain
```
6. Next step is to deactivate the ***r-reticulate*** environment typing:
```
conda deactivate /Users/biagi/Library/r-miniconda/envs/r-reticulate
```
7. Finally, we will install `SCCAF` and `cellphonedb` modules:
```
pip install sccaf
pip install cellphonedb
```
8. Any problem installing the softwares, please click in the following links: [`anndata`](https://anndata.readthedocs.io/en/latest/), [`scanpy`](https://scanpy.readthedocs.io/en/stable/installation.html), [`igraph`](https://igraph.org/python/), [`louvain`](https://louvain-igraph.readthedocs.io/en/latest/install.html), [`SCCAF`](https://github.com/SCCAF/sccaf) and [`cellphonedb`](https://github.com/Teichlab/cellphonedb).
9. If you run into problems, please open a new issue, you can do this by going to *'issues'* and clicking on the *'new issue'* icon. We will help you replicate our analysis! Do not fear single cell analysis!


### Citation
`xxxxx`
<!---
Qadir, M.M.F., Alvarez-Cubela, S., Klein, D., Van Dijk, J., Anquela, R.M., Lanzoni, G., Sadiq, S., Moreno-Hernandez, Y.B., Navarro-Rubio, B., Garcia, M.T., Diaz, A., Johnson, K., Sant, D., Ricordi, C., Griswold, T., Pastori, R.L., Dominguez-bendala, J. (2020) Proceedings of the National Academy of Sciences. Single cell resolution analysis of the human pancreatic ductal progenitor cell niche. Apr 2020, 201918314; DOI: 10.1073/pnas.1918314117
--->


### Contributors
1. Carlos Alberto Oliveira de Biagi Junior - [Github](https://github.com/cbiagii) - [University of Sao Paulo](https://www.fmrp.usp.br/en/) - to contact please [Email](mailto:biagi@usp.br)
2. Sarah Santiloni Cury - [Sao Paulo State University](https://www.international.unesp.br/) - to contact please [Email](mailto:santiloni.cury@unesp.br)
3. Miguel Luiz Batista Junior - [Boston University School of Medicine](https://www.bumc.bu.edu/busm/) - to contact please [Email](mailto:migueljr4@me.com)


### Acknowledgements
`xxxxx`
<!---
1. Diabetes Research Institute Foundation (DRIF)
2. The Inserra family
3. The Fred and Mabel R. Parks Foundation
4. The Tonkinson Foundation
5. ADA Grant #1-19-ICTS-078
6. NIH Grant #1R43DK105655-01
7. NIH Grant #2R44DK105655-02
8. NIH/NIDDK HIRN Grant #U01DK120393 (These studies are part of this grant)
9. IIE Fulbright
-->


{{% alert note %}}
***If this was useful to you please don't forget to cite, star and fork this repository!***
{{% /alert %}}

