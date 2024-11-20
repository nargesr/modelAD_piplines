# modelAD_piplines


Pipelines for quantifying, analyzing, and submitting model-AD data in the Mortazavi lab

## Table of Content
- [Prepare metadata to upload the data into Synapse](CLIMB/):
  - [Fetching information on mice from CLIMB](CLIMB/example/)
  - prepare metadata table per each assay

- [Bulk Short Read RNA_seq Analysis](short_read/): Analysing bulk short read RNA_seq using RSEM output
  - Aggregating TPM and count of all mice in the same study
  - Filter metadata and genes (optional)
  - Differential Expression gene (DEG) analysis ([EdgR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) and [PyDESeq2](https://pydeseq2.readthedocs.io/en/latest/))
  - Weighted correlation network analysis ([PyWGCNA](https://github.com/mortazavilab/PyWGCNA))

- [Bulk Long Read RNA_seq Analysis](long_read): Analysing bulk long read RNA_seq using Kallisto
  - Quantify reads using Kallisto long read
  - Aggregating TPM and count of all mice in the same study
  - Differential Gene Expression (DGE) analysis ([PyDESeq2](https://pydeseq2.readthedocs.io/en/latest/))
  - Differential Expression Isoform (DIE) analysis ([PyDESeq2](https://pydeseq2.readthedocs.io/en/latest/))
  - SWAN report ([SWAN_vis](https://github.com/mortazavilab/swan_vis))
  
- [Validating new mouse model using bulk long read RNA_seq](variant_validation): Validate new variant splicing through bulk long read RNA_seq using minimap and IGV
  - Prepare custom GTF and Fasta files if we have a chimeric gene or humanized gene
  - Map reads using [minimap2](https://github.com/lh3/minimap2)
  - Visualize the coverage tracks using [IGV](https://igv.org/)
 



