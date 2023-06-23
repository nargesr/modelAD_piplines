import swan_vis as swan
import matplotlib.pyplot as plt

annot = '/home/ubuntu/ECCB2022/data/reference/gencode.v40.annotation.encode_pilot_regions.gtf'
gtf = 'h1_talon.gtf'
ab = 'h1_talon_abundance.tsv'

sg = swan.SwanGraph()
sg.add_annotation(annot)
sg.add_transcriptome(gtf)
sg.add_abundance(ab)

sg.save_graph('swan')

# add metadata
m = {'PacBio_cDNA_H1_DE_1_alignments_ENCFF362CPC_subset': 'h1_de_1',
     'PacBio_cDNA_H1_DE_2_alignments_ENCFF168AZV_subset': 'h1_de_2',
     'PacBio_cDNA_H1_DE_3_alignments_ENCFF940WVU_subset': 'h1_de_3',
     'PacBio_cDNA_H1_ESC_1_alignments_ENCFF213XDA_subset': 'h1_1',
     'PacBio_cDNA_H1_ESC_2_alignments_ENCFF281VKZ_subset': 'h1_2',
     'PacBio_cDNA_H1_ESC_3_alignments_ENCFF250BDM_subset': 'h1_3'}
sg.adata.obs['sample'] = sg.adata.obs.dataset.map(m)
# sg.adata.obs.head()
meta = sg.adata.obs.copy(deep=True)
meta['cell_type'] = meta['sample'].str.rsplit('_', n=1, expand=True)[0]

meta = meta[['dataset', 'sample', 'cell_type']]
meta.to_csv('swan_metadata.tsv', sep='\t', index=False)
meta = 'swan_metadata.tsv'
sg.add_metadata(meta)

# add colors for h1 and h1_de
c_dict = {'h1': 'darkorchid', 'h1_de': 'darkgoldenrod'}
sg.set_metadata_colors('cell_type', c_dict)

sg.save_graph('swan')

# intron retention and exon skipping
es_df = sg.find_es_genes(verbose=False)
ir_df = sg.find_ir_genes(verbose=False)

ir_df.head()

ir_tid = ir_df.tid.tolist()[0]
ir_gid = ir_df.gid.tolist()[0]

# inspect a gene with IR and ES
sg.plot_graph(ir_gid, indicate_novel=True)

sg.plot_transcript_path(ir_tid, indicate_novel=True)

es_df.head()

es_gid = es_df.gid.tolist()[0]
es_tid = es_df.tid.tolist()[0]

sg.plot_graph(es_gid, indicate_novel=True)

sg.plot_transcript_path(es_tid, indicate_novel=True)

sg.plot_transcript_path(es_tid, browser=True)

# differential expression tests
sg = swan.read('swan.p')

obs_col = 'cell_type'
# degs = sg.de_gene_test(obs_col)

# dets = sg.de_transcript_test(obs_col)

die, results = sg.die_gene_test(obs_col=obs_col, verbose=True)

die_df = sg.get_die_genes(obs_col=obs_col,
                          obs_conditions=['h1', 'h1_de'],
                          p=0.05, dpi=10)

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

plt.clf()
make_reports('DES')
make_reports('POGZ')
make_reports('PI4KB')

tid = sg.t_df.loc[(sg.t_df.gname == 'PI4KB')&(sg.t_df.novelty!='Known')].tid.tolist()[0]
sg.plot_transcript_path(tid, indicate_novel=True)

sg.plot_each_transcript_in_gene('PI4KB', 'figures/pi4kb', indicate_novel=True)
