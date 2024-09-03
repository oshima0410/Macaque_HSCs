#!/usr/bin/env python
# coding: utf-8

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

# Trajectory analysis
## Pseudotime analysis

# MEMP pathway
adata = sc.read_h5ad('./anndataobjs/Mf.HSPC.scVI.processed.h5ad')

## subset MEMP related clusters(HSC/MPP,EP.1,EP.2,MEMP,MkP,Mast,Prolif-MPP.1)
adata_sub = adata[adata.obs.leiden_sub_R2.isin(
    ['0','1','3','9','10','11','12',]),:]

sc.pp.neighbors(adata_sub, use_rep="X_scVI", n_neighbors = 30, n_pcs = 60)
sc.tl.umap(adata_sub, min_dist=0.5, spread = 1.0)
sc.tl.leiden(adata_sub, key_added="leiden_scVI", resolution=1.0)
sc.pl.umap(adata_sub, color=["cell_type"], size = 10,frameon=True, add_outline= True)

sc.tl.diffmap(adata_sub)
sc.tl.paga(adata_sub, groups = 'cell_type')
sc.pl.paga(adata_sub, color = ['cell_type'], edge_width_scale=1, layout = 'fa',fontoutline=True, save=True)
sc.tl.draw_graph(adata_sub, init_pos='paga')#, layout ='grid_fr')
sc.pl.draw_graph(adata_sub, color=['cell_type'], size=10)

## visualize paga on FDG
pos=pd.DataFrame(adata_sub.obsm['X_draw_graph_fa'],index=adata_sub.obs_names)
pos['group']=adata_sub.obs[adata_sub.uns['paga']['groups']]
pos=pos.groupby('group').mean()
ax=sc.pl.draw_graph(adata_sub,add_outline = False, frameon=False, color = 'cell_type', show=False)
sc.pl.paga(adata_sub, color='cell_type',
           node_size_scale=1, edge_width_scale=0.3,
           threshold=0.05,
           pos=pos.values,
           random_state=0, ax=ax)

## compute diffusion pseudotime
adata_sub.uns['iroot'] = np.flatnonzero(adata_sub.obs['cell_type'] == 'HSC/MPP')[0]
sc.tl.dpt(adata_sub)
sc.pl.draw_graph(adata_sub, color=['dpt_pseudotime'],cmap = 'Reds', size =30,legend_loc='on data')

## save the fdg_coordinates and h5ad file
pd.DataFrame(adata_sub.obsm["X_draw_graph_fa"]).to_csv("./results/Pseudotime_fdg_subset_240201.csv")
adata_sub.write('./anndataobjs/Pseudotime_fdg_sub_240201.h5ad')



## Do the same thing for EP pathway (HSC/MPP,Prolif-MPP.1, MEMP, EP.1,EP.2) and Mast pathway(HSC/MPP,Prolif-MPP.1, MEMP, Mast).
## Output files for EP and Mast pathway are below, respectively
### EP pathway
# pd.DataFrame(adata_sub.obsm["X_draw_graph_fa"]).to_csv("./results/Pseudotime_fdg_EP_subset_240229.csv")
# adata_sub.write('./anndataobjs/Pseudotime_fdg_EP_sub_240229.h5ad')
### Mast pathway
# pd.DataFrame(adata_sub.obsm["X_draw_graph_fa"]).to_csv("./results/Pseudotime_fdg_Mast_subset_240229.csv")
# adata_sub.write('./anndataobjs/Pseudotime_fdg_Mast_sub_240229.h5ad')



# MkP pathway
adata= sc.read_h5ad('./anndataobjs/Mf.HSC.scVI.processed.h5ad')
## subset MkP related clusters(HSC/MPP,Prolif-MPP.1, MEMP, MkP)
adata_sub = adata[adata.obs.cell_type.isin(
    ['HSC/MPP','Prolif-MPP1','MEMP','MkP']),:]

sc.pp.neighbors(adata_sub, use_rep="X_scVI", n_neighbors = 30, n_pcs = 60)
sc.tl.umap(adata_sub, min_dist=0.5, spread = 1.0)
sc.tl.leiden(adata_sub, key_added="leiden_scVI", resolution=1.0)
sc.pl.umap(adata_sub, color=["cell_type"], size = 10,frameon=True, add_outline= True)

sc.tl.diffmap(adata_sub)
sc.tl.paga(adata_sub, groups = 'cell_type')
sc.pl.paga(adata_sub, color = ['cell_type'], edge_width_scale=1, layout = 'fa',fontoutline=True, save=True)
sc.tl.draw_graph(adata_sub, init_pos='paga')
sc.pl.draw_graph(adata_sub, color=['cell_type'], size=10)

sc.tl.leiden(adata_sub, resolution=0.3, restrict_to=('leiden_scVI', ['3']), key_added='leiden_sub_R2')
sc.pl.draw_graph(adata_sub, color=["leiden_sub_R2"], size = 10,frameon=True, add_outline= True)

pos=pd.DataFrame(adata_sub.obsm['X_draw_graph_fa'],index=adata_sub.obs_names)
pos['group']=adata_sub.obs[adata_sub.uns['paga']['groups']]
pos=pos.groupby('group').mean()
ax=sc.pl.draw_graph(adata_sub,add_outline = False, frameon=False, color = 'cell_type', show=False)
sc.pl.paga(adata_sub, color='cell_type',
           node_size_scale=1, edge_width_scale=0.3,
           threshold=0.05,
           pos=pos.values,
           random_state=0, ax=ax)

## Further subset MkP indirect pathway
adata_sub = adata_sub[adata_sub.obs.leiden_sub_R2.isin(["1","2","3,0","3,1","4","5","6","7","8"]),:]

sc.pp.neighbors(adata_sub, use_rep="X_scVI", n_neighbors = 30, n_pcs = 60)
sc.tl.umap(adata_sub, min_dist=0.5, spread = 1.0)
sc.tl.leiden(adata_sub, key_added="leiden_scVI", resolution=1.0)
sc.pl.umap(adata_sub, color=["cell_type"], size = 10,frameon=True, add_outline= True)

sc.tl.diffmap(adata_sub)
sc.tl.paga(adata_sub, groups = 'leiden_scVI')
sc.pl.paga(adata_sub, color = ['leiden_scVI'], edge_width_scale=1, layout = 'fa',fontoutline=True, save=True)
sc.tl.draw_graph(adata_sub, init_pos='paga')
sc.pl.draw_graph(adata_sub, color=['leiden_scVI'], size=10)
sc.pl.draw_graph(adata_sub, color=["cell_type"], size = 10,frameon=True, add_outline= True)

### compute diffusion pseudotime
adata_sub.uns['iroot'] = np.flatnonzero(adata_sub.obs['leiden_scVI'] == '0')[0]
sc.tl.dpt(adata_sub)
sc.pl.draw_graph(adata_sub, color=['dpt_pseudotime'],cmap = 'Reds', size =30,legend_loc='on data')

### save the fdg_coordinates and h5ad file
pd.DataFrame(adata_sub.obsm["X_draw_graph_fa"]).to_csv("./results/Pseudotime_fdg_MkP_indirect_subset_240206.csv")
adata_sub.write('./anndataobjs/Pseudotime_fdg_MkP_indirect_sub_240206.h5ad')



## Further subset MkP direct pathway
adata_sub = adata_sub[adata_sub.obs.leiden_sub_R2.isin(["0","2","9"]),:]
sc.pp.neighbors(adata_sub, use_rep="X_scVI", n_neighbors = 30, n_pcs = 60)
sc.tl.umap(adata_sub, min_dist=0.5, spread = 1.0)
sc.tl.leiden(adata_sub, key_added="leiden_scVI", resolution=1.0)
sc.pl.umap(adata_sub, color=["cell_type"], size = 10,frameon=True, add_outline= True)

sc.tl.diffmap(adata_sub)
sc.tl.paga(adata_sub, groups = 'leiden_scVI')
sc.pl.paga(adata_sub, color = ['leiden_scVI'], edge_width_scale=1, layout = 'fa',fontoutline=True, save=True)
sc.tl.draw_graph(adata_sub, init_pos='paga')
sc.pl.draw_graph(adata_sub, color=['leiden_scVI'], size=10)

### further subset related clusters
adata_sub = adata_sub[adata_sub.obs.leiden_scVI.isin(["0","1","2","3","4","5",]),:]
adata_sub = adata_sub[adata_sub.obs.cell_type.isin(["HSC/MPP","MkP",]),:]

sc.pl.draw_graph(adata_sub, color=["cell_type"], size = 10,frameon=True, add_outline= True)
sc.pl.draw_graph(adata_sub, color=["leiden_scVI"], size = 10,frameon=True, add_outline= True)

### compute diffusion pseudotime
adata_sub.uns['iroot'] = np.flatnonzero(adata_sub.obs['leiden_scVI'] == '1')[0]
sc.tl.dpt(adata_sub)
sc.pl.draw_graph(adata_sub, color=['dpt_pseudotime'],cmap = 'Reds', size =30,legend_loc='on data')

### save the fdg_coordinates and h5ad file
pd.DataFrame(adata_sub.obsm["X_draw_graph_fa"]).to_csv("./results/Pseudotime_fdg_MkP_direct_subset_240206.csv")
adata_sub.write('./anndataobjs/Pseudotime_fdg_MkP_direct_sub_240206.h5ad')
