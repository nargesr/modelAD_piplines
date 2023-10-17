# modelAD_piplines for bulk short read RNA-seq data

## Overview

- [Creating count, TPM and metadata](make_matrix.ipynb)
    - Creating TPM and count matrix by reading rsem files
    - Creating metadata based on the file name
    - Filtering TPM and count matrix by only using PolyA genes
    - Check sex and genotype specific genes such as Thy1, Xsit to avoid any mislabling samples

- differential expression gene (DEG) analysis
    - using [PyDESeq2](PyDESEq2.ipynb)
    - using [EdgR](DEG.ipynb)

- Weighted correlation network analysis (WGCNA)
    - Done by PyWGCNA(PyWGCNA.ipynb)

## Directory structure

- main matrices
    - located in data folder
    - countMatrix contains count value
    - expressionList contains TPM
    - ployA means only has polyA genes instead of all genes
    - sorted means that the columns sorted based on the time/sex/genotype
    - experimentList contains metadata
    
- DEG
    - located in DEG folder
    - name of each files contains the criteria of age/sex unless we used all ages or both sexes
    - first genotype is first group and second genotype is second group which means that first genotype vs. second genotype were done
    
- PyWGCNA
    - figures located in figures folder
