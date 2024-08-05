# Pipeline for analysing lond read data using long read Kallisto

## Setting up and install neccessary packages

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
### Install other neccessary pakcages

1. [bustools](https://github.com/BUStools/bustools?tab=readme-ov-file#installation)
2. [kb-python](https://github.com/pachterlab/kb_python?tab=readme-ov-file#installation)
3. python packages
    - pandas
    - numpy
    - anndata
    - scanpy


## Run Kallisto

- make indices
- Quantify reads with kallisto

More info: [Run Kallisto](kallisto-lr.md)

## Concate count and tpm matrices

Once you have transcript-level count and tpm matrices for each sample, it's time to concat them. Please look at the [this notebook](analysis_pipeline.ipynb) for more information.
