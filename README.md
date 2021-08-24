# Prep4DeepDEP: an R package to prepare data for DeepDEP

## Introduction

*Prep4DeepDEP* is an R package to prepare the input and output data files required by the Python package of [*DeepDEP*](https://codeocean.com/capsule/7914207/tree/). *DeepDEP* is a deep learning model that uses the baseline genomic profiles of cancer cell lines (CCLs) or tumors to predict their gene dependencies (i.e., the degree to which knocking out a gene inhibits the survival of cancer cells).

*Prep4DeepDEP* generates the genomic and gene fingerprint data tables from user's datasets. *Prep4DeepDEP* has two main modes:
- The *‘Prediction’* mode generates data for *DeepDEP* to predict gene dependency scores of unscreened CCLs or tumors. It extracts and orders the required genomic features from user’s genome-wide datasets, and generates the functional fingerprints of gene dependencies of interest (DepOIs) from a user-provided list or the default 1,298 genes we studied in the paper. For the copy number alteration (CNA) data, an embedded R function (*PrepCNA*) converts copy-number segments to bins (every 10k bases in the genome) and calculate per-bin CNA scores.
- The *‘Training’* mode generates data to train a new *DeepDEP* model using user's genome-wide genomic data and gene dependency scores from an in-house CRISPR screening experiment. It creates data tables of genomics and gene dependencies for all CCL-DepOI pairs (number of samples = number of CCLs x number of DepOIs). Functional fingerprints are generated based on the list of genes available in the gene dependency dataset.

Please refer to the paper and [*DeepDEP* package](https://codeocean.com/capsule/7914207/tree/) about how to use the generated data tables for DeepDEP model training and prediction.

## Installation from GitHub ##
```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("ChenLabGCCRI/Prep4DeepDEP")
```

## Run Prep4DeepDEP ##
```R
Prep4DeepDEP <- function(exp.data = NULL, mut.data = NULL,
                         meth.data = NULL,cna.data  = NULL,
                         dep.data = NULL, mode = c("Training","Prediction"),
                         filename.out = "data_out")
```

## Examples ##
Please see the examples in the example.r file:
```R
path <- system.file("examples/",package = "Prep4DeepDEP")
file.edit(paste0(path,"example.r"))
```

## Inputs/outputs
#### Inputs:
- Single or multi-omic data (gene mutation, gene expression, DNA methylation, and CNA) of CCLs or tumors. Dimension: #genomic features by #CCLs/tumors.
- List of DepOIs of interest with or without corresponding gene dependency data (*training*: required; *prediction*: optional). Dimension: #DepOIs by #CCLs if screening data are available for the *training* mode, or #DepOIs by 1 (gene symbols of DepOIs) for the *prediction* mode.
#### Outputs:
*Prediction* mode
- *.txt* file for each genomic data. Dimension: #genomic features by #CCLs/tumors.
- *.txt* file for gene fingerprints. Dimension: #fingerprint features (3,115 chemical and genetic perturbation [CGP] gene sets) by #DepOIs.

*Training* mode
- *.txt* file for each genomic data. Dimension: #genomic features by #CCL-DepOI pairs.
- *.txt* file for gene fingerprints. Dimension: #fingerprint features (3115) by #CCL-DepOIs.

Please refer to the manual of detailed descriptions of each input and output parameter:
```R
?Prep4DeepDEP
```

## Flowchart
<img align="center" src="./sketch/Prep4DeepDEP.png?raw=true">

Note: examples of the output data (paths labeled in blue) are available at the [*DeepDEP* repository](https://codeocean.com/capsule/7914207/tree/).

## R Version
The functionality of this package was developed and tested using R version 3.6.1 (2019-07-05) on an x86_64-pc-linux-gnu platform.

## Reference
Chiu YC, Zheng S, Wang LJ, Iskra BS, Rao MK, Houghton PJ, Huang Y, Chen Y. 
**"Predicting and characterizing a cancer dependency map of tumors with deep learning."** *Science Advances*. 2021;7(34). Epub 2021/08/22. doi: 10.1126/sciadv.abh1275.
