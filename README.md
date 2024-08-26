# modelAD_piplines

Pipelines for quantification and analyzing model AD data in the Mortazavi lab

## Table of Content
- Fetching information of mice from [CLIMB](https://climb.bio)

- [Bulk Short Read RNA_seq](short_read/): Analysing bulk short read RNA_seq using RSEM output
    - aggregating TPM and count of all mice in the same study
    - Filter metadata and genes (optional)
    - Differential Expression gene (DEG) analysis ([EdgR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) and [PyDESeq2](https://pydeseq2.readthedocs.io/en/latest/))
    - [PyWGCNA](https://github.com/mortazavilab/PyWGCNA)

- [Bulk Long Read RNA_seq](long_read/): Under construction!
 



