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
adata = sc.read_h5ad('./anndataobjs/Mf.HSC.scVI.combined.processed.h5ad')

## UMAP and leiden clustering
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors = 40, n_pcs = 60)
sc.tl.umap(adata, min_dist=1.0, spread = 1.5)
sc.tl.leiden(adata, key_added="leiden_scVI", resolution=1.2)
sc.pl.umap(adata, color=["leiden_scVI"], size = 10,frameon=True, add_outline= True)

## re-clustering
### #12 into 4 clusters (Later annotated as 'Mono/Macrophage','Endothelial','Fibroblast','B')
sc.tl.leiden(adata, resolution=0.1, restrict_to=('leiden_scVI', ['12']), key_added='leiden_sub_R')

### #1 into 2 clusters (Later annotated as 'MPP.3','MPP.4')
sc.tl.leiden(adata, resolution=0.3, restrict_to=('leiden_sub_R', ['1']), key_added='leiden_sub_R')

categories = [str(i) for i in range(len(adata.obs['leiden_sub_R'].cat.categories))]
adata.obs['leiden_sub_R'] = adata.obs['leiden_sub_R'].cat.rename_categories(categories)
sc.pl.umap(adata, color=["leiden_sub_R"], size = 10,frameon=True, add_outline= True)


## PAGA & FDG using PAGA initialization
sc.tl.diffmap(adata)
sc.tl.paga(adata, groups = 'leiden_sub_R')
sc.pl.paga(adata, color = ['leiden_sub_R'], node_size_scale=3.0, edge_width_scale=0.4, fontoutline=True)
sc.tl.draw_graph(adata, init_pos='paga')#, layout ='grid_fr')
sc.pl.draw_graph(adata, color=['leiden_sub_R'], size=10, title ='')

## visualize the results
sc.pl.paga_compare(
    adata, threshold=0.1, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=False)

# Find marker genes

## Wilcoxon rank-sum
sc.tl.rank_genes_groups(adata, groupby = 'leiden_sub_R', method = 'wilcoxon', tie_correct = True, pts = True, key_added = 'wilcoxon')
MF_HSC_cluster_genes.head(100).to_csv(r'./results/scVI.ALL_HSC.wlcxn.csv', index=False)
sc.get.rank_genes_groups_df(adata, group='0', key='wilcoxon').to_csv(r'./results/pHSC_cluser0_HSC.wlcxn.csv', index=False)

MF_HSC_cluster_genes = pd.DataFrame(adata.uns['wilcoxon']['names'])
MF_HSC_cluster_genes.head(60)

sc.tl.dendrogram(adata, groupby="leiden_sub_R", use_rep="X_scVI")
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key="wilcoxon", groupby="leiden_sub_R", var_group_labels= 'leiden_sub_R', save= '.dotplot.scVI.pdf')

markers = sc.get.rank_genes_groups_df(adata, None, key='wilcoxon')
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
adata.uns['markers'] = markers

## scVI DE (1-vs-all)
model  = scvi.model.SCVI.load('my_model_240109/', adata)
markers_scvi = model.differential_expression(groupby = 'leiden_sub_R')
markers_scvi = markers_scvi[(markers_scvi['is_de_fdr_0.05']) & (markers_scvi.lfc_mean > .5)]
adata.uns['scvi_markers'] = markers_scvi

# save 
adata.write('./anndataobjs/Mf.HSC.scVI.processed.h5ad')
