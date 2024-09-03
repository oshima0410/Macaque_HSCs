#!/usr/bin/env python
# coding: utf-8

# Cell type annotation and marker genes

import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
from igraph import *
from fa2 import ForceAtlas2
import networkx as nx
import seaborn as sns
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=120, dpi_save=300)
sc.settings.n_jobs = 30
pd.options.display.max_columns = None
warnings.simplefilter(action = "ignore", category = FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="scanpy")
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")

## data import
adata = sc.read_h5ad('./anndataobjs/Mf.HSPC.scVI.processed.h5ad')

## create a dictionary to map cluster to annotation label
## cell types are manually annotated based on differentially expressed genes (wilcoxon or scVI-DE).
cluster2annotation = {
     '0': 'HSC/MPP',
     '1': 'EP.1',
     '2': 'GMP',
     '3': 'EP.2',
     '4': 'GrP',
     '5': 'MoP.2',
     '6': 'Pre Pro-B',
     '7': 'pDC',
     '8': 'Pro-B',
     '9': 'MEMP',
     '10': 'MkP',
     '11': 'Mast',
     '12': 'Prolif-MPP.1',
     '13': 'Prolif-unspecified',
     '14': 'Prolif-GrP',
     '15': 'MoP.1',
     '16': 'HSC/CLP',
     '17': 'Prolif-MPP.2',
     '18': 'EoP',
     '19': 'Macrophage',
     '20': 'Endothelial',
     '21': 'Fibroblast',
     '22': 'B',
     '23': 'B', #(both annotated 'B')
     '24': 'Neut',
     '25': 'NK',
     '26': 'Erythroid',
}
adata.obs['cell_type'] = adata.obs['leiden_sub_R2'].map(cluster2annotation).astype('category')

## visualization
combined_palette = sns.color_palette("Set1") + sns.color_palette("Set2") + sns.color_palette("Set3")
sns.set_palette(combined_palette)
color_dict={
     'HSC/MPP':combined_palette[0],
     'EP.1':combined_palette[1],
     'GMP':combined_palette[2],
     'EP.2':combined_palette[3],
     'GrP':combined_palette[4],
     'MoP.2':combined_palette[5],
     'Pre Pro-B':combined_palette[6],
     'pDC':combined_palette[7],
     'Pro-B':combined_palette[8],
     'MEMP':combined_palette[9],
     'MkP':combined_palette[10],
     'Mast':combined_palette[11],
     'Prolif-MPP.1':combined_palette[12],
     'Prolif-unspecified':combined_palette[13],
     'Prolif-GrP':combined_palette[14],
     'MoP.1':combined_palette[15],
     'HSC/CLP':combined_palette[16],
     'Prolif-MPP.2':combined_palette[17],
     'EoP':combined_palette[18],
     'Macrophage':combined_palette[19],
     'Endothelial':combined_palette[20],
     'Fibroblast':combined_palette[21],
     'B':combined_palette[22],
     'Neut':combined_palette[24],
     'NK':combined_palette[25],
     'Erythroid':combined_palette[26]}
sc.pl.draw_graph(adata, color='cell_type', legend_loc='on data',size=10,
           frameon=False, add_outline=True, palette = color_dict,
           legend_fontsize=20, legend_fontweight='bold', legend_fontoutline= True, title='')

sc.pl.umap(adata, color='cell_type',size=10,
           frameon=True, add_outline=True, legend_loc= 'on data',  palette = color_dict,
           legend_fontsize=25, legend_fontweight='bold', legend_fontoutline= True, title='')


## Find marker genes
sc.tl.rank_genes_groups(adata, groupby = 'cell_type', method = 'wilcoxon', tie_correct = True, pts = True, key_added = 'wilcoxon_celltype')
celltype_marker = pd.DataFrame(adata.uns['wilcoxon_celltype']['names'])

sc.tl.dendrogram(adata, groupby="cell_type", use_rep="X_scVI")

markers = sc.get.rank_genes_groups_df(adata, None, key='wilcoxon_celltype')
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > 1.0) & (markers.pct_nz_group > 0.05)]
markers = markers[~markers['names'].str.startswith('LOC')]

### extract cell type specific genes
markers[markers['group']=='EP.1'][:5].names.tolist()

## Dotplot
genes_to_show = {
    'HSC/MPP':['HOPX', 'IFITM1', 'HLF'],
    'Prolif-MPP.1':['MXD3', 'TOP2A', 'NDC80'],
    'Prolif-MPP.2':['PIMREG', 'CDC25B', 'CENPA'],
    'GMP':['RHEX', 'CSF3R', 'CEBPA'],
    'MoP.1':['ARHGAP24', 'LYST','VCAN'],
    'MoP.2':['SERPINB10','FCGR1A','RETN'],
    'GrP':['CTSG', 'AZU1', 'PRTN3'],
    'EoP':['EPX', 'PRG2', 'IL5RA'],
    'Prolif-GrP':['ELANE', 'CTSG', 'AZU1'],
    'NK':['GZMK', 'GZMA', 'KLRD1'],
    'Neut':['CHIT1', 'S100A12', 'ORM1'],
    'Macrophage':['C1QA', 'C1QC', 'C1QB'],
    'MEMP':['CDK1', 'SPC25', 'CDCA3'],
    'MkP':['ITGB3', 'PLXDC2', 'PBX1'],
    'EP.1':['KLF1', 'TFR2', 'KCNH2'],
    'EP.2':['PKLR', 'TLCD4','SLC25A21'],
    'Erythroid':['SLC4A1', 'PHOSPHO1', 'TSPO2'],
    'Mast':['CMA1', 'MS4A2', 'HDC'],
    'HSC/CLP':['CCR9','FCMR', 'LCK'],
    'Pre Pro-B':['IL7R', 'RAG2', 'BLNK'],
    'Pro-B':['VPREB3', 'POU2AF1','PAX5'],
    'pDC':['LY6D', 'SFMBT2', 'IRF8'],
    'B':['TNFRSF13B', 'MS4A1', 'FCRLA'],
    'Endothelial':['CLDN5', 'COL4A1', 'STAB2'],
    'Fibroblast':['SOD3', 'PLAC9', 'ESM1'],
    'Prolif-unspecified':['H1-4', 'THBS1', 'SCD'],
}

desired_order = ['HSC/MPP','Prolif-MPP.1','Prolif-MPP.2','GMP', 'MoP.1','MoP.2','GrP','EoP','Prolif-GrP','NK','Neut','Macrophage','MEMP','MkP','EP.1','EP.2','Erythroid','Mast', 'HSC/CLP','Pre Pro-B','Pro-B','pDC', 'B','Endothelial','Fibroblast','Prolif-unspecified']
sc.pl.dotplot(adata, genes_to_show, groupby="cell_type", var_group_labels='cell_type', categories_order=desired_order, standard_scale='var', log=True, dendrogram=False,colorbar_title='Expression', size_title = '')

## barplot
combined_palette = sns.color_palette("Set1") + sns.color_palette("Set2") + sns.color_palette("Set3")

category_series = adata.obs['agetissue']
category_series = category_series.cat.reorder_categories(['Early2nd_FL','Late2nd_FL','Early3rd_FL','Late2nd_BM','Early3rd_BM','Adult_BM'])
freq_dist_by_leiden = pd.crosstab(category_series, adata.obs['cell_type'], normalize='index')

desired_order = ['HSC/MPP','Prolif-MPP.1','Prolif-MPP.2','GMP', 'MoP.1','MoP.2','GrP','EoP','Prolif-GrP','NK','Neut','Macrophage','MEMP','MkP','EP.1','EP.2','Erythroid','Mast', 'HSC/CLP','Pre Pro-B','Pro-B','pDC', 'B','Endothelial','Fibroblast','Prolif-unspecified']
freq_dist_by_leiden = freq_dist_by_leiden[desired_order]

# Create a color dictionary
color_dict={
     'HSC/MPP':combined_palette[0],
     'EP.1':combined_palette[1],
     'GMP':combined_palette[2],
     'EP.2':combined_palette[3],
     'GrP':combined_palette[4],
     'MoP.2':combined_palette[5],
     'Pre Pro-B':combined_palette[6],
     'pDC':combined_palette[7],
     'Pro-B':combined_palette[8],
     'MEMP':combined_palette[9],
     'MkP':combined_palette[10],
     'Mast':combined_palette[11],
     'Prolif-MPP.1':combined_palette[12],
     'Prolif-unspecified':combined_palette[13],
     'Prolif-GrP':combined_palette[14],
     'MoP.1':combined_palette[15],
     'HSC/CLP':combined_palette[16],
     'Prolif-MPP.2':combined_palette[17],
     'EoP':combined_palette[18],
     'Macrophage':combined_palette[19],
     'Endothelial':combined_palette[20],
     'Fibroblast':combined_palette[21],
     'B':combined_palette[22],
     'Neut':combined_palette[24],
     'NK':combined_palette[25],
     'Erythroid':combined_palette[26]}
desired_order.reverse()
ax = freq_dist_by_leiden.plot.bar(stacked=True, legend=None, color=color_dict)
legend_handles = [plt.Line2D([0], [0], color=color_dict[label], lw=4) for label in desired_order]
plt.legend(legend_handles, desired_order, loc='right', bbox_to_anchor=(1.5, 0.4))
plt.show()

## save
adata.write('./anndataobjs/Mf.HSPC.scVI.processed.h5ad')
