# Bulk short-read RNA-seq data

## Overview

- [Creating count, TPM and metadata](make_matrix.ipynb)
    - Creating TPM and count matrix by reading rsem files
    - Creating metadata based on the file name
    - Filtering TPM and count matrix by only using PolyA genes
    - Check sex and genotype specific genes such as Thy1, Xsit to avoid any mislabeling samples

- differential expression gene (DEG) analysis
    - using [PyDESeq2](PyDESEq2.ipynb)
    - using [EdgR](DEG.ipynb)

- Weighted correlation network analysis (WGCNA)
    - Done by [PyWGCNA](PyWGCNA.ipynb)

## Directory structure

- main matrices
    - located in the data folder
    - countMatrix contains the count value
    - expressionList contains TPM
    - ployA means only has polyA genes instead of all genes
    - sorted means that the columns are sorted based on the time/sex/genotype
    - experimentList contains metadata
    
- DEG
    - located in the DEG folder
    - The name of each file contains the criteria of age/sex unless we used all ages or both sexes
    - The first genotype is the first group and the second genotype is the second group which means that the first genotype vs. second genotype was done
    
- PyWGCNA
    - figures located in the figures folder
