#!/usr/bin/env python
# coding: utf-8

# Data exploration

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
adata = sc.read_h5ad('./anndataobjs/Mf.HSPC.scVI.combined.processed.h5ad')

## PCA
sc.tl.pca(adata, svd_solver = 'arpack', n_comps=120, use_highly_variable=True)
sc.pl.pca(adata, color=['batch'], wspace=0.5)
sc.pl.pca_scatter(adata, color=['total_counts', 'n_genes_by_counts'], cmap = 'Reds')
sc.pl.pca_variance_ratio(adata, log=False, show=120, n_pcs=120)

## UMAP and leiden clustering
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors = 30, n_pcs = 60)
sc.tl.umap(adata, min_dist=0.5, spread = 1.0)
sc.tl.leiden(adata, key_added="leiden_scVI", resolution=1.0)
sc.pl.umap(adata, color=["leiden_scVI"], size = 10,frameon=True, add_outline= True)

## re-clustering
### #6 into three sub-clusters(later annotated 'Pre Pro-B','pDC','Pro-B')
sc.tl.leiden(adata, resolution=0.2, restrict_to=('leiden_scVI', ['6']), key_added='leiden_sub_R2')

### #16 into three sub-clusters(later annotated 'Endothelial','Macrophage','Fibroblast' )
sc.tl.leiden(adata, resolution=0.1, restrict_to=('leiden_sub_R2', ['16']), key_added='leiden_sub_R2')

### #17 into three sub-clusters(later annotated 'B','Neut','NK' )
sc.tl.leiden(adata, resolution=0.1, restrict_to=('leiden_sub_R2', ['17']), key_added='leiden_sub_R2')

### #11 into two sub-clusters(later annotated 'Prolif-unspecified','Prolif-GrP' )
sc.tl.leiden(adata, resolution=0.2, restrict_to=('leiden_sub_R2', ['11']), key_added='leiden_sub_R2')

categories = [str(i) for i in range(len(adata.obs['leiden_sub_R2'].cat.categories))]
adata.obs['leiden_sub_R2'] = adata.obs['leiden_sub_R2'].cat.rename_categories(categories)
sc.pl.umap(adata, color=["leiden_sub_R2"], size = 10,frameon=True, add_outline= True)

## PAGA & FDG using PAGA initialization
sc.tl.diffmap(adata)
sc.tl.paga(adata, groups = 'leiden_sub_R2')
sc.pl.paga(adata, color = ['leiden_sub_R2'], edge_width_scale=1, layout = 'fa',fontoutline=True, save=True)
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color=['leiden_sub_R2'], size=10, title ='')

## visualize the results
sc.pl.paga_compare(
    adata, threshold=0.1, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=False)

# Find marker genes

## Wilcoxon rank-sum
sc.tl.rank_genes_groups(adata, groupby = 'leiden_sub_R2', method = 'wilcoxon', tie_correct = True, pts = True, key_added = 'wilcoxon')
MF_CD34_cluster_genes = pd.DataFrame(adata.uns['wilcoxon']['names'])
MF_CD34_cluster_genes.head(100).to_csv(r'./results/240114.scVI.ALL_CD34.wlcxn_ver2.csv', index=False)

sc.tl.dendrogram(adata, groupby="leiden_sub_R2", use_rep="X_scVI")
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key="wilcoxon", groupby="leiden_sub_R2", var_group_labels= 'leiden_sub_R2', save= '.dotplot.scVI.pdf')

markers = sc.get.rank_genes_groups_df(adata, None, key='wilcoxon')
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
adata.uns['markers'] = markers

## scVI DE (1-vs-all)
model  = scvi.model.SCVI.load('my_model_240104/', adata)
markers_scvi = model.differential_expression(groupby = 'leiden_sub_R2')
markers_scvi = markers_scvi[(markers_scvi['is_de_fdr_0.05']) & (markers_scvi.lfc_mean > .5)]
adata.uns['scvi_markers'] = markers_scvi

# save
adata.write('./anndataobjs/Mf.HSPC.scVI.processed.h5ad')
