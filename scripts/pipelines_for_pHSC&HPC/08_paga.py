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
adata = sc.read_h5ad('./anndataobjs/Mf.HSPC.scVI.processed.h5ad')


## paga
sc.tl.paga(adata, groups = 'cell_type')
sc.pl.paga(adata, color = ['cell_type'], edge_width_scale=1, layout = 'fa',fontoutline=True, save=True)

## plot paga
combined_palette = sns.color_palette("Set1") + sns.color_palette("Set2") + sns.color_palette("Set3")
sns.set_palette(combined_palette)

pos=pd.DataFrame(adata.obsm['X_draw_graph_fa'],index=adata.obs_names)
pos['group']=adata.obs[adata.uns['paga']['groups']]
pos=pos.groupby('group').mean()

ax=sc.pl.draw_graph(adata,add_outline = False, frameon=False, color = 'cell_type', show=False)
sc.pl.paga(adata, color='cell_type',
           node_size_scale=1, edge_width_scale=0.3,
           threshold=0.1,
           pos=pos.values,
           random_state=0, ax=ax)

## paga connectivity
cell_cats = list(adata.obs["cell_type"].cat.categories)
population_size = adata.obs['cell_type'].value_counts()
connectivities = np.array(csr_matrix.todense(adata.uns["paga"]["connectivities"]), dtype = "float64")
connectivities[connectivities < .05] = 0.0
connectivities = pd.DataFrame(connectivities, columns = cell_cats, index = cell_cats)
order = ['HSC/MPP','Prolif-MPP.1','Prolif-MPP.2','GMP', 'MoP.1','MoP.2','GrP','EoP','Prolif-GrP','NK','Neut','Macrophage','MEMP','MkP','EP.1','EP.2','Erythroid','Mast', 'HSC/CLP','Pre Pro-B','Pro-B','pDC', 'B','T','Endothelial','Fibroblast','Prolif-unspecified']
connectivities_order = connectivities[order].T[order]
kws = dict(cbar_kws=dict(label='PAGA connectivity'))
sns.heatmap(connectivities_order,**kws,annot =  False)
plt.show()

## removing mature clusters /clusters with few cells ('Macrophage','Endothelial','Fibroblast', 'B','Neut','NK','Erythroid',)
adata_sub = adata[adata.obs.leiden_sub_R2.isin(
    ['0','1','2','3','4','5','6','7','8','9','10','11','12','15','16','17','18']),:]

### paga
sc.tl.paga(adata_sub, groups = 'cell_type')
sc.pl.paga(adata_sub, color = ['cell_type'], edge_width_scale=1, layout = 'fa',fontoutline=True)

combined_palette = sns.color_palette("Set1") + sns.color_palette("Set2") + sns.color_palette("Set3")
sns.set_palette(combined_palette)
pos=pd.DataFrame(adata_sub.obsm['X_draw_graph_fa'],index=adata_sub.obs_names)
pos['group']=adata_sub.obs[adata_sub.uns['paga']['groups']]
pos=pos.groupby('group').mean()
ax=sc.pl.draw_graph(adata,add_outline = False, frameon=False, color = 'cell_type', palette = color_dict, show=False, title='', legend_loc='none')
sc.pl.paga(adata_sub, color='cell_type',
           node_size_scale=1, edge_width_scale=0.3,
           threshold=0.2,
           pos=pos.values,
           random_state=0, ax=ax)

### paga connectivity
cell_cats = list(adata_sub.obs["cell_type"].cat.categories)
population_size = adata_sub.obs['cell_type'].value_counts()
connectivities = np.array(csr_matrix.todense(adata_sub.uns["paga"]["connectivities"]), dtype = "float64")
connectivities[connectivities < .05] = 0.0
connectivities = pd.DataFrame(connectivities, columns = cell_cats, index = cell_cats)
order = ['HSC/MPP','Prolif-MPP.1','Prolif-MPP.2','HSC/CLP','pDC','Pre Pro-B','Pro-B','GMP','MoP.1','MoP.2','GrP','EoP','MEMP','MkP','EP.1','EP.2','Mast']
connectivities_order = connectivities[order].T[order]
kws = dict(cbar_kws=dict(label='PAGA connectivity'))#,orientation='horizontal'))
sns.heatmap(connectivities_order,**kws,annot =  False)
plt.show()
