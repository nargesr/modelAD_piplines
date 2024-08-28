# Pipeline for analyzing bulk long read RNA_seq data using long read Kallisto

## Setting up and installing necessary packages

### Install kallisto (>= v0.51.0)

kallisto is located in [this repo](https://github.com/pachterlab/kallisto).

Follow commands for installation

```bash
git clone https://github.com/pachterlab/kallisto
cd kallisto
mkdir build
cd build
cmake -DMAX_KMER_SIZE=64 ..
make 
```
### Install other necessary packages

1. [bustools](https://github.com/BUStools/bustools?tab=readme-ov-file#installation)
2. [kb-python](https://github.com/pachterlab/kb_python?tab=readme-ov-file#installation)
3. python packages
    - pandas
    - numpy
    - anndata
    - PyDESeq2

## overview of pipeline

- Run Kallisto
    - make indices
    - Quantify reads with kallisto
    - More info: [Run Kallisto](kallisto-lr.md)

- Concat count and TPM matrices
    - After preparing transcript-level count and TPM matrices for each sample, it's time to concat them.
    - more information on [this notebook](analysis_pipeline.ipynb)

- differential isoform expression (DIE) analysis
    - using [PyDESeq2](https://pydeseq2.readthedocs.io/en/latest/)
    - more information on [this notebook](analysis_pipeline.ipynb)

- differential gene expression (DGE) analysis
    - using [PyDESeq2](https://pydeseq2.readthedocs.io/en/latest/)
    - more information on [this notebook](analysis_pipeline.ipynb)

- SWAN reports
    - under construction

## Directory structure

- data
    - sample_name.txt: text file contains sample name and corresponding fastq files ([more info](kallisto-lr.md))
    - Kallisto outputs per sample ([more info](kallisto-lr.md))
    - sample_metadata.csv: CSV table contains information about each sample. sample name should match the sample name in sample_name.txt
    - transcript_exp_tpm.h5ad: h5ad file contains TPM values along with sample information as obs and transcript information as var.
    - transcript_exp_count.h5ad: h5ad file contains Count values along with sample information as obs and transcript information as var.
    
- Differential Isoform Expression (DIE)
    - located in the DIE folder
    - The name of each file contains the criteria of age/sex unless we used all ages or both sexes
    - The first genotype is the first group and the second genotype is the second group which means that the first genotype vs. second genotype was done
 
- Differential Gene Expression (DGE)
    - located in the DEG folder
    - The name of each file contains the criteria of age/sex unless we used all ages or both sexes
    - The first genotype is the first group and the second genotype is the second group which means that the first genotype vs. second genotype was done
    
- SWAN
    - under construction
