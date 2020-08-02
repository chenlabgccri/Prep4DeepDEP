# Prep4DeepDEP: an R package to prepare data for DeepDEP

## Introduction

*Prep4DeepDEP* is an R package to prepare the input and output data files required by the Python package of [*DeepDEP*](https://codeocean.com/capsule/3348251/tree/v1) (will be pubished with our manuscript). DeepDEP is a deep learning model to predict the gene dependency profile of an unscreened cancer cell line (CCL) or impracticable-to-screen tumors based on the baseline genomic profiles.

*Prep4DeepDEP* generates the genomic and gene fingerprint data tables from datasets provided by the user. Prep4DeepDEP has two main modes:
- The ‘Prediction’ mode generates data for *DeepDEP* to predict gene dependency scores of unscreened CCLs or tumors. It extracts and orders the required genomic features from user’s genome-wide genomic datasets, and generates functional fingerprints from a list of user-provided genes or the default 1,298 genes we studied in the manuscript. For the copy number alteration (CNA) data, an embedded R function (*PrepCNA*) is embedded to covert copy-number segments to bins and calculate per-bin CNA scores.
- The ‘Training’ mode generates data to train a new *DeepDEP* model using user's genomic and gene dependency data. It creates txt tables of genomic and gene dependency data for all CCL-DepOI pairs (number of samples = number of CCLs x number of DepOIs pair) based on user’s dataset. Functional fingerprints are generated based on the list of genes in the gene dependency data.

Please refer to the paper and [*DeepDEP* package](https://codeocean.com/capsule/3348251/tree/v1) about how to use the genertated data tables for DeepDEP model training and prediction.

## Installation from GitHub ##
```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("chenlabgccri/Prep4DeepDEP")
```
## Run Prep4DeepDEP ##
```R
Prep4DeepDEP()
```
Please read the manual of each input and output parameter:
```R
?Prep4DeepDEP
```

## Flowchart
<img align="center" src="./sketch/Prep4DeepDEP.png?raw=true">

## Reference
Chiu YC, Zheng S, Wang LJ, Iskra BS, Rao MK, Houghton PJ, Huang Y, Chen Y.
**"DeepDEP: deep learning of a cancer dependency map using cancer genomics."**
*Nature Communications.* In revision.
