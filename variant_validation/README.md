# Pipeline for validating new variant splicing

## Setting up and installing necessary packages

1. [minimap2](https://github.com/lh3/minimap2)
2. [samtools](https://www.htslib.org/)
3. python packages
    - [pandas](https://pandas.pydata.org/docs/index.html)
    - [AGEpy](https://github.com/mpg-age-bioinformatics/AGEpy)

If you are running it on the HPC, you can use `module load` to load necessary packages.

For python packages, you can use `pip install [pythonPackage]` or check their installation section.

## Overview of pipeline

- Create custom GTF and FASTA file if necessary
    - If we have chimeric gene (hClu) or humanized gene (hMapt)
    - Subtract GENCODE human GTF and human FASTA files for given genes and treat it as a new chromosome
    - More info: [Making custom GTF and FASTA](make_custom_GTF_FASTA.ipynb)

- Map and Visualize reads
    - Use minimap to map reads using mouse/customized Fasta and GTF files
    - Use IGV to look at the reads to validate the new mouse model

## Directory structure

- ref
    - located in the ref folder
    - Contain GTF and Fasta references
    
- minimap
    - located in the minimap folder
    - Contain one folder for each sample

- sample_name.txt: text file contains sample name and corresponding fastq files ([more info](minimap_IGV.md))

