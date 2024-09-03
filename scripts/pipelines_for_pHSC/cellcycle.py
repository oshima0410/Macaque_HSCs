#!/usr/bin/env python
# coding: utf-8

#　Cell cycle

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

## import cell cycle related genes list
cell_cycle_genes = pd.read_csv('../regev_lab_cell_cycle_genes2.txt', sep='\t', header = None)

cell_cycle_genes[0]
cell_cycle_genes[0] = cell_cycle_genes[0].str.upper()
s_genes = pd.DataFrame(cell_cycle_genes[0][:43].to_list())
g2m_genes = pd.DataFrame(cell_cycle_genes[0][43:].to_list())
gene_list = adata.var_names
def find_genes(*genes):
    list = []
    for g in genes:
        list.extend([i for i in gene_list if i.split('.')[0].endswith(tuple('' + g))])
    return(list)

s_genes = find_genes(s_genes[0])
g2m_genes = find_genes(g2m_genes[0])
cell_cycle_genes = find_genes(cell_cycle_genes[0])

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
adata_ccgenes = adata[:, cell_cycle_genes]
sc.tl.umap(adata_ccgenes)

# visualization
sc.pl.umap(adata_ccgenes, color='phase',add_outline=True, size=30)

esired_order = ['G1','S','G2M']
adata_ccgenes.obs['phase'] = pd.Categorical(adata_ccgenes.obs['phase'], categories=desired_order, ordered=True)
combined_palette = sns.color_palette("RdBu")
color_dict={
    'G1': "#0173b2",  # Early phase
    'S': "#de8f05",  # Synthesis phase
    'G2M': "#029e73"  # Late phase
}
sc.pl.draw_graph(adata_ccgenes, color="phase", add_outline=True, size=10, frameon=False, palette=color_dict, title='')
sc.pl.draw_graph(adata_ccgenes, color="phase", groups='G1', add_outline=True, size=10, frameon=False,na_color='white', palette=color_dict, title='')
sc.pl.draw_graph(adata_ccgenes, color="phase", groups='S', add_outline=True, size=10, frameon=False,na_color='white', palette=color_dict, title='')
sc.pl.draw_graph(adata_ccgenes, color="phase", groups='G2M', add_outline=True, size=10, frameon=False,na_color='white', palette=color_dict, title='')
