#!/usr/bin/env python
# coding: utf-8

# Subsetting to purify HSC and MPPs

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

## subset to purify HSC and MPPs
## select HSC, MPP.1-5
adata_hsc = adata[adata.obs.leiden_sub_R.isin(['0','1','2','3','5','11']),:]


## wilcoxon rank sum test among HSC/MPPs
sc.tl.rank_genes_groups(adata_hsc, groupby = 'cell_type', method = 'wilcoxon', tie_correct = True, pts = True, key_added = 'wilcoxon_hsc_celltype')
hspc_marker = pd.DataFrame(adata_hsc.uns['wilcoxon_hsc_celltype']['names'])
hspc_marker.head(100).to_csv(r'./results/240311.HSC-MPPs.wlcxn.csv', index=False)

hspc_markers = sc.get.rank_genes_groups_df(adata_hsc, None, key='wilcoxon_hsc_celltype')
hspc_markers.to_csv("./results/hspc_markers_allgenes_240130.csv")
hspc_markers = sc.get.rank_genes_groups_df(adata_hsc, None, key='wilcoxon_hsc_celltype')
hspc_markers = hspc_markers[(hspc_markers.pvals_adj < 0.05) & (hspc_markers.logfoldchanges > 1.0) & (hspc_markers.pct_nz_group > 0.05)]
hspc_markers = hspc_markers[~hspc_markers['names'].str.startswith('LOC')]

## dotplot
### select marker gene to show
hspc_markers[hspc_markers['group']=='MPP.5'][:5].names.tolist()

genes_to_show = {
'HSC':['MLLT3', 'MECOM', 'TENM4', 'HLF', 'NRIP1'],
'MPP.1':['MAFA-F', 'FOS', 'ANXA1', 'B2M','DUSP1'],
'MPP.2':['CDK6', 'HCST', 'LAT2', 'RFLNB', 'SELL'],
'MPP.3':['VIM', 'LMNA', 'TPPP3', 'TSPAN2', 'CRIP1'],
'MPP.4':['ZBTB16', 'PKIG', 'GATA2', 'PBX1', 'GATA1'],
'MPP.5':['TUBA1B', 'TYMS', 'MCM5', 'PCLAF', 'MCM7'],
}

sc.tl.dendrogram(adata_hsc, groupby="cell_type", use_rep="X_scVI")
desired_order = ['HSC','MPP.1','MPP.2','MPP.3','MPP.4','MPP.5']
sc.pl.dotplot(adata_hsc, genes_to_show, groupby="cell_type", var_group_labels='cell_type', categories_order=desired_order, standard_scale='var', log=True, dendrogram=False,colorbar_title='Expression', size_title = '')

## barplot
combined_palette = sns.color_palette("Set1") + sns.color_palette("Set2") + sns.color_palette("Set3") + sns.color_palette("Dark2")

category_series = adata_hsc.obs['agetissue']
category_series = category_series.cat.reorder_categories(['Early2nd_FL','Late2nd_FL','Early3rd_FL','Late2nd_BM','Early3rd_BM','Adult_BM'])
freq_dist_by_leiden = pd.crosstab(category_series, adata_hsc.obs['cell_type'], normalize='index')

desired_order = ['HSC','MPP.1','MPP.2','MPP.3','MPP.4','MPP.5']
freq_dist_by_leiden = freq_dist_by_leiden[desired_order]

color_dict={
     'HSC':combined_palette[30],
     'MPP.2':combined_palette[31],
     'MPP.3':combined_palette[32],
     'MPP.5':combined_palette[33],
     'MPP.4':combined_palette[34],
     'MPP.1':combined_palette[36],}

desired_order.reverse()
ax = freq_dist_by_leiden.plot.bar(stacked=True, legend=None, color=color_dict)
legend_handles = [plt.Line2D([0], [0], color=color_dict[label], lw=4) for label in desired_order]
plt.legend(legend_handles, desired_order, loc='right', bbox_to_anchor=(1.6, 0.4),fontsize='large')
plt.show()
