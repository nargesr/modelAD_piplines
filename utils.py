import warnings
import os
import sys
import itertools
import copy


import pandas as pd
import numpy as np
import scanpy as sc

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

import gseapy as gp
from gseapy import dotplot

import seaborn as sns
from adjustText import adjust_text
import matplotlib.pylab as plt
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['pdf.fonttype'] = 42
sns.set_context('paper')

def plot_volcano(df, 
                 ofile, 
                 obs_filtering, 
                 obs_condition,
                 l2fc_thresh=0, 
                 adj_p_thresh=0.05,
                 label='transcript_name'):
    

    df = pd.read_csv(df, sep='\t', index_col=0)

    df['label'] = df[label].tolist()

    df.loc[df.DE == "No", 'label'] = ""
    
    # Calculate counts
    num_up = df[df.DE == "Up"].shape[0]
    num_down = df[df.DE == "Down"].shape[0]
    num_not_significant = df[df.DE == "No"].shape[0]
    
    # Plotting
    plt.scatter(x=df['log2FoldChange'], y=df['padj'].apply(lambda x: -np.log10(x)), s=1,
                label=f"Not significant (n={num_not_significant})")
    down = df[df.DE == "Down"]
    down.sort_values(["padj"], inplace=True)
    plt.scatter(x=down['log2FoldChange'], y=down['padj'].apply(lambda x: -np.log10(x)), s=3,
                label=f"Up-regulated in {obs_filtering[obs_condition][0]} (n={num_down})", color="blue")
    up = df[df.DE == "Up"]
    up.sort_values(["padj"], inplace=True)
    plt.scatter(x=up['log2FoldChange'], y=up['padj'].apply(lambda x: -np.log10(x)), s=3,
                label=f"Up-regulated in {obs_filtering[obs_condition][1]} (n={num_up})", color="red")
    texts = []
    for i in range(min(10, up.shape[0])):
        texts.append(plt.text(x=up.iloc[i]['log2FoldChange'],
                              y=-np.log10(up.iloc[i]['padj']),
                              s=up.iloc[i][label]))
    for i in range(min(10, down.shape[0])):
        texts.append(plt.text(x=down.iloc[i]['log2FoldChange'],
                              y=-np.log10(down.iloc[i]['padj']),
                              s=down.iloc[i][label]))
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    plt.xlabel("logFC")
    plt.ylabel("-log10(adj p-value)")
    plt.axvline(l2fc_thresh, color="grey", linestyle="--")
    plt.axhline(-np.log10(adj_p_thresh), color="grey", linestyle="--")
    # Adjust the legend with a numerical font size
    plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.35))  # Change the font size here
    
    title = ""
    for key in obs_filtering.keys():
        if key != obs_condition:
            title = f"{title} {'_'.join(obs_filtering[key])}"
    title = f"{title} samples\n{obs_filtering[obs_condition][0]} Vs. {obs_filtering[obs_condition][1]}"
    plt.title(title)
    plt.tight_layout()

    plt.savefig(ofile, bbox_inches='tight', dpi=500)
    plt.show()

def run_deseq2(adata,
               obs_filtering,
               obs_condition,
               l2fc_thresh=0, 
               adj_p_thresh=0.05,
               label='transcript_name',
               ofile="test"):
    """
    Parameters:
        fname (str): Name of input AnnData
        obs_filtering (dict): group of filtering we used
        obs_conditions (str): condition for comparing. must be a key in obs_filtering
        how (str): {'gene', 'iso'}
        ofile (str): Name of output file
        threads (int): Number of threads to run on
    """

    # wrong number of conditions
    if len(obs_filtering[obs_condition]) != 2:
        raise ValueError(f'{obs_filtering[obs_condition]} not valid. Please provide exactly two conditions to compare')

    for col in obs_filtering.keys():
        if col not in adata.obs.columns:
            raise ValueError(f'{col} not found in adata')
        adata = adata[adata.obs[col].isin(obs_filtering[col])]


    # remove unexpressed stuff
    sc.pp.filter_genes(adata, min_counts=1)

    # convert to int
    adata.X = np.round(adata.X).astype(int)

    # run test
    dds = DeseqDataSet(adata=adata,
                   design_factors=obs_condition,
                   refit_cooks=True)
    dds.deseq2()
    
    stat_res = DeseqStats(dds,
                          contrast=[obs_condition] + obs_filtering[obs_condition])
    stat_res.summary()

    # save output
    df = stat_res.results_df

    df = pd.concat([df, adata.var], axis=1)
    
    # call things as upregulated or downregulated
    df['DE'] = 'No'
    df.loc[(df.log2FoldChange >= l2fc_thresh)&\
           (df.padj <= adj_p_thresh), 'DE'] = 'Up'
    df.loc[(df.log2FoldChange <= -1*l2fc_thresh)&\
           (df.padj <= adj_p_thresh), 'DE'] = 'Down'

    df.to_csv(f"{ofile}.tsv", sep='\t')

    plot_volcano(df=f"{ofile}.tsv", 
                 ofile=f"{ofile}.pdf", 
                 obs_filtering=obs_filtering,
                 obs_condition=obs_condition, 
                 l2fc_thresh=l2fc_thresh, 
                 adj_p_thresh=adj_p_thresh,
                 label=label)


def Gene_ontology_analysis(df, 
                                   l2fc_thresh=0, 
                                   adj_p_thresh=0.05, 
                                   sets=['GO_Biological_Process_2023'], 
                                   organism='Mouse',
                                   background=None,
                                   p_value=1, 
                                   file_name=None, 
                                   **kwargs):
    """
        Doing GO analysis

        :param df: file path
        :type df: str
        :param sets: str, list, tuple of Enrichr Library name(s). or custom defined gene_sets (dict, or gmt file) (you can add any Enrichr Libraries from here: https://maayanlab.cloud/Enrichr/#stats) only need to fill if the type is GO or KEGG
        :type sets: str, list, tuple
        :param p_value: Defines the pValue threshold. (default: 0.05)
        :type p_value: float
        :param file_name: name of the file you want to use to save plot (default is up/down)
        :type file_name: str
        :param kwargs: Other keyword arguments are passed through to the underlying gseapy.enrichr() finction
        :type kwargs: key, value pairings
    """
    df = pd.read_csv(df, sep='\t')

    # call things as upregulated or downregulated
    df['tmp'] = 'No'
    df.loc[(df.log2FoldChange >= l2fc_thresh)&\
           (df.padj <= adj_p_thresh), 'tmp'] = 'Up'
    df.loc[(df.log2FoldChange <= -1*l2fc_thresh)&\
           (df.padj <= adj_p_thresh), 'tmp'] = 'Down'

    # Calculate counts
    up = df[df.DE == "Up"].gname.fillna("").values.tolist()

    try:
        enr = gp.enrichr(gene_list=up,
                                 gene_sets=sets,
                                 organism=organism,
                                 outdir=f"{file_name}_up",
                                 cutoff=p_value,
                                 **kwargs)
        dotplot(enr.res2d,
                        title=f"Gene ontology in up-regulated genes",
                        cmap='viridis_r',
                        cutoff=p_value,
                        ofname=f"{file_name}_up.pdf")
    except:
        print(f"No enrich terms when cutoff = {p_value} for up-regulated genes")

    
    down = df[df.DE == "Down"].gname.fillna("").values.tolist()

    try:
        enr = gp.enrichr(gene_list=down,
                                 gene_sets=sets,
                                 organism=organism,
                                 outdir=f"{file_name}_down",
                                 cutoff=p_value,
                                 **kwargs)
        dotplot(enr.res2d,
                        title=f"Gene ontology in down-regulated genes",
                        cmap='viridis_r',
                        cutoff=p_value,
                        ofname=f"{file_name}_down.pdf")
    except:
        print(f"No enrich terms when cutoff = {p_value} for up-regulated genes")

