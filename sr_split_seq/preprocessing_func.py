## import library
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import anndata as ad

import warnings
warnings.filterwarnings("ignore")

def add_well(adata, bc1):
    adata.obs['bc1_well'] = [bcs[16:] for bcs in adata.obs.index.tolist()]
    adata.obs['bc1_well'].replace(bc1.sequence.values, bc1.well.values, inplace=True)
    
    return adata

def run_scrublet(adata):
    wells = adata.obs.bc1_well.unique().tolist()
    for well in wells:
        print(f"Calculating doublet for well {well}")
        tmp = adata[adata.obs.bc1_well == well]
        scrub = scr.Scrublet(tmp.X)
        adata.obs.loc[tmp.obs.index, 'doublet_scores'],_ = scrub.scrub_doublets(min_counts=1, 
                                                                                min_cells=1, 
                                                                                min_gene_variability_pctl=85, 
                                                                                n_prin_comps=30)
        
    return adata


def add_cell_information(adata, sample_info, on='bc1_well'):
    adata.obs['index'] = adata.obs.index.tolist()
    adata.obs = adata.obs.merge(sample_info, 
                                on=on)
    adata.obs.index = adata.obs['index'].values + "_" + adata.obs.Mouse_Tissue_ID.values
    adata.obs.drop(['index'], axis=1, inplace=True)
    return adata


def preprocessing_per_sublibraries(experiment_id,
                                   sublibraries,
                                   kallisto_output_path,
                                   bc1):
    ## Read data
    annData_name = f"{kallisto_output_path}/{experiment_id}_{sublibraries}/counts_unfiltered_modified/adata.h5ad"
    adata = sc.read_h5ad(annData_name)

    adata = add_well(adata, bc1)
    adata = run_scrublet(adata)
    
    return adata

def aggregate_subliberaries(experiment_id, sublibraries, anndata_path):
    adata = None
    for sublibrary in sublibraries:
        adata_sublibrary = sc.read_h5ad(f"{anndata_path}/{experiment_id}_{sublibrary}/adata.h5ad")
        if adata is None:
            adata = adata_sublibrary
        else:
            adata = ad.concat([adata, adata_sublibrary])
    return adata