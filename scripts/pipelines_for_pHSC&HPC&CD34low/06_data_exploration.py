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
adata = sc.read_h5ad('./anndataobjs/Mf.E58CD34.scVI.combined.processed.h5ad')

## PCA
sc.tl.pca(adata, svd_solver = 'arpack', n_comps=120, use_highly_variable=True)
sc.pl.pca(adata, color=['batch'], wspace=0.5)
sc.pl.pca_scatter(adata, color=['total_counts', 'n_genes_by_counts'], cmap = 'Reds')
sc.pl.pca_variance_ratio(adata, log=False, show=120, n_pcs=120)

## UMAP and leiden clustering
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors = 30, n_pcs = 80)
sc.tl.umap(adata, min_dist=0.5, spread = 1.0)
sc.tl.leiden(adata, key_added="leiden_scVI", resolution=1.0)
sc.pl.umap(adata, color=["leiden_scVI"], size = 10,frameon=True, add_outline= True)

desired_order = ['HSC','HPC','CD34low']
adata.obs['population'] = pd.Categorical(adata.obs['population'], categories=desired_order, ordered=True)
sc.pl.umap(adata, color="population", add_outline=True, size=10, frameon=False,title='')

sc.pl.umap(adata,color = ['HLF','CALR','KLF1','LYZ','AZU1','EPX',
                        'ITGB3','HDC','VPREB3','IRF8'], 
           cmap = 'Reds', size = 10, ncols=5,frameon=False,add_outline=True,
           title = ['HSC/MPP\n(HLF)','GMP\n(CALR)','EP\n(KLF1)','MoP\n(LYZ)','GrP\n(AZU1)','EoP\n(EPX)',
                    'MkP\n(ITGB3)','Mast\n(HDC)','Pre Pro-B\n(VPREB3)','pDC\n(IRF8)'])

sc.pl.umap(adata,color = ['CD34','CD38','THY1'],
           cmap = 'Reds', size = 10, ncols=5,frameon=False,add_outline=True,
           title = ['CD34','CD38','THY1(CD90)'])

# save
adata.write('./anndataobjs/Mf.E58CD34.scVI.processed.h5ad')
