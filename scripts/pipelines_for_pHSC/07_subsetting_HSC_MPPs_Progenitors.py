#!/usr/bin/env python
# coding: utf-8

# Subsetting to remove contaminated mature cells

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

## subset to remove mature cells
## remove Mast,Pre Pro-B/Pro-B, B,Endothelial, Macrophage
adata_hspc = adata[adata.obs.leiden_sub_R.isin(['0','1','2','3','4','5','6','7','8','9','11','12','17','20']),:]
adata_hspc

## paga
sc.tl.paga(adata_hspc, groups = 'cell_type')
sc.pl.paga(adata_hspc, color = ['cell_type'], edge_width_scale=1, layout = 'fa',fontoutline=True, save=True)


## plot paga
combined_palette = sns.color_palette("Set1") + sns.color_palette("Set2") + sns.color_palette("Set3")
sns.set_palette(combined_palette)

pos=pd.DataFrame(adata_hspc.obsm['X_draw_graph_fa'],index=adata_hspc.obs_names)
pos['group']=adata_hspc.obs[adata_hspc.uns['paga']['groups']]
pos=pos.groupby('group').mean()

ax=sc.pl.draw_graph(adata,add_outline = False, frameon=False, color = 'cell_type', show=False, title='',legend_loc='none')
sc.pl.paga(adata_hspc, color='cell_type',
           node_size_scale=1, edge_width_scale=0.3,
           threshold=0.5,
           pos=pos.values,
           random_state=0, ax=ax)

## paga connectivity
cell_cats = list(adata_hspc.obs["cell_type"].cat.categories)
population_size = adata_hspc.obs['cell_type'].value_counts()
connectivities = np.array(csr_matrix.todense(adata_hspc.uns["paga"]["connectivities"]), dtype = "float64")
connectivities[connectivities < .05] = 0.0
connectivities = pd.DataFrame(connectivities, columns = cell_cats, index = cell_cats)
order = ['HSC','MPP.1','MPP.2','MPP.3','MPP.4','MPP.5','Prolif-MPP.1','Prolif-MPP.2','MPP/GMP', 'GMP/GrP','MEMP','MkP','EP.1/EP.2/Erythroid']
connectivities_order = connectivities[order].T[order]
kws = dict(cbar_kws=dict(label='PAGA connectivity'))#,orientation='horizontal'))
sns.heatmap(connectivities_order,**kws,annot =  False)
plt.show()


## violinplot comparing each cell types
desired_order = ['HSC','MPP.1','MPP.2','MPP.3','MPP.4','MPP.5','MPP/GMP','Prolif-MPP.1','Prolif-MPP.2','MEMP','MkP','EP.1/EP.2/Erythroid','GMP/GrP']
ax = sc.pl.violin(adata_hspc, ['HLF'], groupby='cell_type', order=desired_order, stripplot=False, rotation=90, show=False)
ax.set_yticks([])
ax.set_ylim(bottom=0)  # This sets the bottom of the y-axis to 0
plt.show()

combined_palette = sns.color_palette("Set1",24)
desired_order = ['HSC','MPP.1','MPP.2','MPP.3','MPP.4','MPP.5','MPP/GMP','Prolif-MPP.1','Prolif-MPP.2','MEMP','MkP','EP.1/EP.2/Erythroid','GMP/GrP']
genes_to_show = ['HLF','MLLT3','MECOM','IGLL1','SELL','BTK','LAT2','LAG3','GATA2','MPIG6B','MPO','CLEC11A','MKI67','PCNA']
sc.pl.stacked_violin(adata_hspc, genes_to_show, groupby='cell_type', categories_order =desired_order,
                     figsize= [10,8], standard_scale ='var', swap_axes =True, cmap='Reds',
                     colorbar_title = 'expression',row_palette= combined_palette)
