# modelAD_piplines

Pipelines for quantification and analyzing model AD data in the Mortazavi lab

## Table of Content
- Fetching information of mice from [CLIMB](https://climb.bio)

- [Bulk Short Read RNA_seq Analysis](short_read/): Analysing bulk short read RNA_seq using RSEM output
    - Aggregating TPM and count of all mice in the same study
    - Filter metadata and genes (optional)
    - Differential Expression gene (DEG) analysis ([EdgR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) and [PyDESeq2](https://pydeseq2.readthedocs.io/en/latest/))
    - Weighted correlation network analysis ([PyWGCNA](https://github.com/mortazavilab/PyWGCNA))

- [Bulk Long Read RNA_seq Analysis](long_read/): Analysing bulk long read RNA_seq using Kallisto
  - Quantify reads using Kallisto long read
  - Aggregating TPM and count of all mice in the same study
  - Differential Gene Expression (DGE) analysis ([PyDESeq2](https://pydeseq2.readthedocs.io/en/latest/))
  - Differential Expression Isoform (DIE) analysis ([PyDESeq2](https://pydeseq2.readthedocs.io/en/latest/))
  - SWAN report (_under construction_)
  
- Validating new mouse model using bulk long read RNA_seq (_under construction_)
 



