import warnings
import os
import sys
import itertools
import copy


import pandas as pd
import numpy as np
import scanpy as sc

import swan_vis as swan
import cerberus

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

def create_swan_obj(study, 
                    gtf_cerberus,
                    ab_cerberus,
                    annot,
                    meta,
                    swan_output):
    "Create a swan object from cerberus output"

    sg = swan.SwanGraph()

    sg.add_annotation(annot)

    for gtf in gtf_cerberus:
        sg.add_transcriptome(f"cerberus/{study}/{gtf}", include_isms=True)

    for ab in ab_cerberus:
        sg.add_abundance(f"cerberus/{study}/{ab}")

    for ab in ab_cerberus:
        sg.add_abundance(f"cerberus/{study}/{ab}", how='gene')

    sg.add_metadata(meta)
    sg.save_graph(swan_output)

def make_reports(gname):
  sg.gen_report(gname,
                'figures/'+gname,
                metadata_cols=['cell_type'],
                cmap='viridis',
                transcript_col='tname',
                novelty=True,
                indicate_novel=True,
                layer='tpm')

  sg.gen_report(gname,
                'figures/'+gname,
                metadata_cols=['cell_type'],
                cmap='magma',
                transcript_col='tname',
                novelty=True,
                layer='pi',
                browser=True)

def plot_volcano(df, 
                 ofile, 
                 obs_filtering, 
                 obs_condition,
                 kind='gene',
                 l2fc_thresh=0, 
                 adj_p_thresh=0.05):
    

    df = pd.read_csv(df, sep='\t')

    if kind == 'gene':
        df['label'] = df.gname
    elif kind == 'transcript':
        df['label'] = df.tname

    if kind == 'gene':
        label = 'gname'
    elif kind == 'transcript':
        label = 'tname'

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

def run_deseq2(sg,
               how,
               obs_filtering,
               obs_condition,
               l2fc_thresh=0, 
               adj_p_thresh=0.05,
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

    if how == 'gene':
        adata = sg.gene_adata
    elif how == 'iso':
        adata = sg.adata

    # wrong number of conditions
    if len(obs_filtering[obs_condition]) != 2:
        raise ValueError(f'{obs_filtering[obs_condition]} not valid. Please provide exactly two conditions to compare')

    for col in obs_filtering.keys():
        if col not in adata.obs.columns:
            raise ValueError(f'{obs_col} not found in adata')
        adata = adata[adata.obs[col].isin(obs_filtering[col])]

    
    # remove novel genes
    if how == 'gene':
        adata = adata[:, adata.var.loc[adata.var.index.str.contains('ENSMUS')].index].copy()

    # remove unexpressed stuff
    sc.pp.filter_genes(adata, min_counts=1)

    # densify matrix
    adata.X = np.array(adata.X.todense())

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
    df.reset_index(inplace=True)

    #return df

    if how == 'iso':
        g_df = sg.t_df[['tid', 'tname']].drop_duplicates().reset_index()
        merge_thing = 'tid'
    elif how == 'gene':
        g_df = sg.t_df[['gid', 'gname']].drop_duplicates().reset_index()
        df['gid_stable'] = cerberus.get_stable_gid(df, 'gid')
        df.drop('gid', axis=1, inplace=True)
        df.rename({'gid_stable':'gid'}, axis=1, inplace=True)
        merge_thing = 'gid'

    df = df.merge(g_df, how='left', on=merge_thing)

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
                 kind=how,
                 l2fc_thresh=l2fc_thresh, 
                 adj_p_thresh=adj_p_thresh)


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



def parse_config_file_analysis(fname,
                      meta_fname,
                      p_meta_fname,
                      geno_fname,
                      an_meta_fname,
                      auto_dedupe=True):

    """
    Parameters:
        fname (str): Path to config file fname. One line per input fastq.
        meta_fname (str): Path to file with metadata information.
        p_meta_fname (str): Path to pseudochromosome metadata information
        geno_fname (str): Path to genotype metadata information
        an_meta_fname (str): Path to analysis metadata information
        datasets_per_run (int): Number of datasets to process in each TALON run
        auto_dedupe (bool): Automatically deduplicate duplicate fastqs that result from
            successive Porechop rounds
        include_pseudochrom (bool): Include models with pseudochrom loci, which need
            some preprocessing / different treatment

    Returns:
        df (pandas DataFrame): DF w/ pipeline information; one line per fastq
        dataset_df (pandas DataFrame): DF w/ dataset information; one line per mouse
    """

    df, p_df = parse_config_common(fname,
                          meta_fname,
                          p_meta_fname,
                          geno_fname,
                          auto_dedupe=True)


    # limit to just the studies and genotypes requested
    an_df = pd.read_csv(an_meta_fname, sep='\t')
    i = len(an_df[['genotype', 'study']].drop_duplicates().index)
    an_df['genotype_study'] = an_df['genotype']+' '+an_df['study']
    df['genotype_study'] = df['genotype']+' '+df['study']
    p_df['genotype_study'] = p_df['genotype']+' '+p_df['study']
    # import pdb; pdb.set_trace()
    # p_df = p_df.loc[(p_df.genotype.isin(genotypes))&\
    #                 (p_df.study.isin(studies))]
    p_df = p_df.loc[p_df.genotype_study.isin(an_df.genotype_study.tolist())]
    p_df.drop('genotype_study', axis=1, inplace=True)
    i2 = len(p_df[['genotype', 'study']].drop_duplicates().index)
    df = df.loc[df.genotype_study.isin(an_df.genotype_study.tolist())]
    df.drop('genotype_study', axis=1, inplace=True)


    # df = df.loc[(df.genotype.isin(genotypes))&\
    #             (df.study.isin(studies))]
    i3 = len(p_df[['genotype', 'study']].drop_duplicates().index)

    genotypes = an_df.genotype.unique().tolist()
    studies = an_df.study.unique().tolist()

    if not (i==i2==i3):
        genotypes = list(set(genotypes)-\
                         set(df.genotype.unique().tolist()))
        studies = list(set(studies)-\
                         set(df.study.unique().tolist()))
        warnings.warn(f'Genotypes {genotypes} and studies {studies} not found. Is this expected?')

    # assign a cerberus run to each "sample" (study+genotype+sex+age+tissue)
    # but first sort on study and sample such that they will always be ordered in the same way
    # this should freeze our results
    gb_cols = ['study', 'genotype', 'sex', 'age', 'tissue']
    df = df.sort_values(by=gb_cols, ascending=True)
    temp = df.copy(deep=True)
    temp = temp[gb_cols].groupby(gb_cols).count().reset_index()
    temp['cerberus_run'] = [i+1 for i in range(len(temp.index))]
    df = df.merge(temp, how='left', on=gb_cols)
    df['cerberus_run'] = df.cerberus_run.astype(str)

    # add in analysis stuff
    an_df = pd.read_csv(an_meta_fname, sep='\t')
    p_df = p_df.merge(an_df, how='left',
                  on=['genotype', 'study'])

    # add cerberus run info
    p_df = p_df.merge(df[gb_cols+['cerberus_run']].drop_duplicates(), how='left',
                 on=gb_cols)

    # sanitize genotype alias internally (int) w/ characters better for file names
    exp = '[^0-9a-zA-Z-_]+'
    p_df.loc[p_df.genotype_alias.isnull(), 'genotype_alias'] = p_df.loc[p_df.genotype_alias.isnull(), 'genotype']
    p_df['genotype_alias_int'] = p_df.genotype_alias.str.replace(exp, '_', regex=True)

    # add in columns for comparisons
    p_df['genotype_sex'] = p_df['genotype_alias_int']+'_'+p_df['sex']

    # sanitize analysis
    exp = '[^0-9a-zA-Z-_]+'
    p_df['analysis'] = p_df.analysis.str.replace(exp, '_', regex=True)

    return df, p_df

