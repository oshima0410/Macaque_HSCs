#!/usr/bin/env python
# coding: utf-8

# Subsetting to purify HSC 

import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
from igraph import *
from fa2 import ForceAtlas2
from scipy.sparse import csr_matrix
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
adata = sc.read_h5ad('./anndataobjs/Mf.HSC.scVI.processed.h5ad')


## subsetting HSC cluster
adata_hsc = adata[adata.obs.cell_type.isin(['HSC']),:]

## wilcoxon rank sum test among each stage/tissue
sc.tl.rank_genes_groups(adata_hsc, groupby = 'agetissue', method = 'wilcoxon', tie_correct = True, pts = True, key_added = 'wilcoxon_hsc_agetissue')
agetissue_marker = pd.DataFrame(adata_hsc.uns['wilcoxon_hsc_agetissue']['names'])
agetissue_marker.head(100).to_csv(r'./results/240303.agetissue_hsc.wlcxn.csv', index=False)

agetissue_markers = sc.get.rank_genes_groups_df(adata_hsc, None, key='wilcoxon_hsc_agetissue')
agetissue_markers.to_csv("./results/hsc_agetissue_markers_allgenes_240223.csv")
agetissue_markers = agetissue_markers[(agetissue_markers.pvals_adj < 0.05) & (agetissue_markers.logfoldchanges > 0.5) & (agetissue_markers.pct_nz_group > 0.01)]
agetissue_markers = agetissue_markers[~agetissue_markers['names'].str.startswith('LOC')]

## dotplot
### select marker gene to show
agetissue_markers[agetissue_markers['group']=='Early2nd_FL'][:5].names.tolist()

genes_to_show = {
'Early2nd_FL':['TMSB10', 'RPL29', 'LDHA', 'FUOM', 'HMGA1','CD33', 'RPL22L1', 'ENO1', 'NPM1', 'CFL1'],
'Late2nd_FL':['JUND', 'DUSP1', 'SRSF5', 'ZFP36', 'FOS','RPL38', 'SMAD9', 'KEG98_p09', 'JUN', 'DDX5'],
'Early3rd_FL':['VCAN', 'HTR2A', 'CLEC10A', 'KCNH7', 'CRISP3','NKAIN2', 'AFF3', 'BRINP3', 'KEG98_p09', 'AHNAK'],
'Late2nd_BM':['RPS4X', 'APOD', 'RPL17', 'GSTM5', 'CACNA1B','RSRP1', 'LST1', 'LYN', 'ITM2C', 'TKT'],
'Early3rd_BM':['TSPO', 'RASD1', 'CLC', 'NKAIN2', 'NEGR1','AFF3', 'VCAN', 'PRSS57', 'PCDH9', 'HTR2A'],
'Adult_BM':['AR', 'JAKMIP1', 'ANXA1', 'RNF17', 'MAFA-F','MAMDC2', 'CACNB4', 'CLIC2', 'ADAM28', 'ANK3'],
}

sc.tl.dendrogram(adata_hsc, groupby="agetissue", use_rep="X_scVI")
desired_order = ['Early2nd_FL','Late2nd_FL','Early3rd_FL','Late2nd_BM','Early3rd_BM','Adult_BM']
sc.pl.dotplot(adata_hsc, genes_to_show, groupby="agetissue", var_group_labels='agetissue',  categories_order=desired_order, standard_scale='var', log=True,dendrogram=False,colorbar_title='Expression', size_title = '')

## violin plot
desired_order = ['Early2nd_FL','Late2nd_FL','Early3rd_FL','Late2nd_BM','Early3rd_BM','Adult_BM']
ax = sc.pl.violin(adata_hsc, ['TMSB10'], groupby='agetissue', order=desired_order, stripplot=False, rotation=90, show=False)
ax.yaxis.tick_right()
for line in ax.yaxis.get_gridlines():
    y_value = line.get_ydata()[0]
    if y_value == 2.5:
        line.set_visible(False)
plt.show()

## heatmap
all_hvg = {'Adult_BM':agetissue_marker['Adult_BM'].head(10),
           'Early2nd_FL':agetissue_marker['Early2nd_FL'].head(10),
           'Early3rd_BM':agetissue_marker['Early3rd_BM'].head(10),
           'Early3rd_FL':agetissue_marker['Early3rd_FL'].head(10),
           'Late2nd_BM':agetissue_marker['Late2nd_BM'].head(10),
           'Late2nd_FL':agetissue_marker['Late2nd_FL'].head(10),
}
sc.pl.heatmap(adata,all_hvg,groupby='agetissue',
              #layer = 'scaled',vmin=-2, vmax=2,
              standard_scale="var",
              swap_axes=True,
              dendrogram= True,
              #save ='heatmap231212.pdf',
              figsize=(8,12),
              cmap= 'coolwarm')

## subsampling for equalize the HSC cell count for each stage
adata_hsc.obs['agetissue'].value_counts()

a = adata[adata.obs['agetissue'] == 'Early2nd_FL']
b = adata[adata.obs['agetissue'] == 'Late2nd_FL']
c = adata[adata.obs['agetissue'] == 'Early3rd_BM']
d = adata[adata.obs['agetissue'] == 'Late2nd_BM']
e = adata[adata.obs['agetissue'] == 'Adult_BM']
f = adata[adata.obs['agetissue'] == 'Early3rd_FL']

sc.pp.subsample(a,n_obs = 210)
sc.pp.subsample(b,n_obs = 210)
sc.pp.subsample(c,n_obs = 210)
sc.pp.subsample(d,n_obs = 210)
sc.pp.subsample(e,n_obs = 210)
sc.pp.subsample(f,n_obs = 210)

adata_hsc = a.concatenate(b,c,d,e,f,
                      batch_categories= ['Early2nd_FL',
                                          'Late2nd_FL',
                                          'Early3rd_BM',
                                          'Late2nd_BM',
                                          'Adult_BM',
                                         'Early3rd_FL'])

sc.tl.rank_genes_groups(adata_hsc, groupby = 'agetissue', method = 'wilcoxon', tie_correct = True, pts = True, key_added = 'wilcoxon_hsc_agetissue')
agetissue_marker = pd.DataFrame(adata_hsc.uns['wilcoxon_hsc_agetissue']['names'])
agetissue_marker.head(100).to_csv(r'./results/240303.subset_agetissue_hsc.wlcxn.csv', index=False)
agetissue_markers = sc.get.rank_genes_groups_df(adata_hsc, None, key='wilcoxon_hsc_agetissue')
agetissue_markers = agetissue_markers[(agetissue_markers.pvals_adj < 0.05) & (agetissue_markers.logfoldchanges > 0.5) & (agetissue_markers.pct_nz_group > 0.01)]
agetissue_markers = agetissue_markers[~agetissue_markers['names'].str.startswith('LOC')]
agetissue_markers = agetissue_markers[~agetissue_markers['names'].str.startswith('KEG98')]

## heatmap
all_hvg = {'Early2nd_FL':agetissue_markers[agetissue_markers['group']=='Early2nd_FL'][:10].names.tolist(),
           'Late2nd_FL':agetissue_markers[agetissue_markers['group']=='Late2nd_FL'][:10].names.tolist(),
           'Early3rd_FL':agetissue_markers[agetissue_markers['group']=='Early3rd_FL'][:10].names.tolist(),
           'Late2nd_BM':agetissue_markers[agetissue_markers['group']=='Late2nd_BM'][:10].names.tolist(),
           'Early3rd_BL':agetissue_markers[agetissue_markers['group']=='Early3rd_BM'][:10].names.tolist(),
           'Adult_BM':agetissue_markers[agetissue_markers['group']=='Adult_BM'][:10].names.tolist(),
}

desired_order = ['Early2nd_FL','Late2nd_FL','Early3rd_FL','Late2nd_BM','Early3rd_BM','Adult_BM']
adata_hsc.obs['agetissue'] = adata_hsc.obs['agetissue'].astype('category')
adata_hsc.obs['agetissue'] = adata_hsc.obs['agetissue'].cat.reorder_categories(desired_order)
sc.pl.heatmap(adata_hsc, all_hvg, groupby='agetissue',
              standard_scale="var",
              swap_axes=True,
              dendrogram=False,
              figsize=(6,12),
              cmap='coolwarm',
             show_gene_labels=True)
