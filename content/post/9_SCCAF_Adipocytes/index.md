---
date: "2020-10-10"
diagram: true
math: true
title: 9 - Finding optimal number of clusters in adipocytes cells using SCCAF
---

Creating `h5ad` file to use as input in SCCAF:
```r
library(Seurat)

source('/Users/biagi/PhD/Adipocyte/2_Functions.R')

data <- readRDS("/Users/biagi/PhD/Adipocyte/output/10x/Adipocytes.rds")
SeuratToH5ad(data, "/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/data.h5ad", "RNA", 1)
```


Defining the optimal accuracy number:
```python
# python to define the optimal accuracy number
import warnings
warnings.filterwarnings("ignore")
from SCCAF import *

ad = sc.read("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/data.h5ad")

y_prob, y_pred, y_test, clf, cvsm, acc = SCCAF_assessment(ad.X, ad.obs['L1_Round0'],n = 100)
aucs = plot_roc(y_prob, y_test, clf, cvsm = cvsm, acc = acc)
plt.show()
```


Optimisation and general purpose usage:
```python
##### Optimisation and general purpose usage
cd /Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/
sccaf -i /Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/data.h5ad --optimise --skip-assessment -s L1_Round0 -a 0.669 -c 10 --produce-rounds-summary -o /Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results.h5ad --optimisation-plots-output /Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results.pdf
```


Parallel run of assessments:
```python
sccaf-assess -i /Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results.h5ad -o /Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/sccaf_assess_L1_Round0.txt --slot-for-existing-clustering L1_Round0 --iterations 20 --cores 16
sccaf-assess -i /Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results.h5ad -o /Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/sccaf_assess_L1_Round1.txt --slot-for-existing-clustering L1_Round1 --iterations 20 --cores 16
sccaf-assess -i /Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results.h5ad -o /Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/sccaf_assess_L1_Round2.txt --slot-for-existing-clustering L1_Round2 --iterations 20 --cores 16
```


Merging parallel runs to produce plot:
```python
sccaf-assess-merger -i /Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly -r /Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/rounds.txt -o /Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/rounds-acc-comparison-plot.png
```


Exporting SCCAF results
```python
import scanpy as sc

adata = sc.read("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results.h5ad")
adata.write_csvs("/Users/biagi/PhD/Adipocyte/SCCAF/AdipocytesOnly/results", sep = '\t', skip_data = True)
```