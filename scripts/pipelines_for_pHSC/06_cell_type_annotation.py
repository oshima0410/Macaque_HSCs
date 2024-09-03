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
adata = sc.read_h5ad('./anndataobjs/Mf.HSC.scVI.processed.h5ad')

## import cell_type label from pHSC+pHPC
metadata = pd.read_csv("./results/240119_barcodelist_pHSC_label_transfered.csv")
metadata.index=metadata["Unnamed: 0"]
adata.obs["cell_type_transfered"] = metadata["cell_type_transfer"]


## plot cell label on pHSC derived cells distribution
sc.pl.draw_graph(adata_sub,color = ['cell_type_transfered'],groups= 'HSC/MPP',title='HSC/MPP',palette=color_dict,size = 30,frameon=False,add_outline=True, legend_loc='none')
sc.pl.draw_graph(adata_sub,color = ['cell_type_transfered'],groups= 'MEMP',title='MEMP',palette=color_dict,size = 30,frameon=False,add_outline=True, legend_loc='none')
sc.pl.draw_graph(adata_sub,color = ['cell_type_transfered'],groups= 'MkP',title='MkP',palette=color_dict,size = 30,frameon=False,add_outline=True, legend_loc='none')
sc.pl.draw_graph(adata_sub,color = ['cell_type_transfered'],groups= 'HSC/CLP',title='HSC/CLP',palette=color_dict,size = 30,frameon=False,add_outline=True, legend_loc='none')
sc.pl.draw_graph(adata_sub,color = ['cell_type_transfered'],groups= 'GMP',title='GMP',palette=color_dict,size = 30,frameon=False,add_outline=True, legend_loc='none')
sc.pl.draw_graph(adata_sub,color = ['cell_type_transfered'],groups= 'Prolif-MPP.1',title='Prolif-MPP.1',palette=color_dict,size = 30,frameon=False,add_outline=True, legend_loc='none')
sc.pl.draw_graph(adata_sub,color = ['cell_type_transfered'],groups= 'Prolif-MPP.2',title='Prolif-MPP.2',palette=color_dict,size = 30,frameon=False,add_outline=True, legend_loc='none')
sc.pl.draw_graph(adata_sub,color = ['cell_type_transfered'],groups= 'Prolif-unspecified',title='Prolif-unspecified',palette=color_dict,size = 30,frameon=False,add_outline=True, legend_loc='none')


## Cluster annotation

### create a dictionary to map cluster to annotation label
cluster2annotation = {
     '0': 'HSC',
     '1': 'MPP.2',
     '2': 'MPP.3',
     '3': 'MPP.5',
     '4': 'MkP',
     '5': 'MPP.4',
     '6': 'MPP/GMP',
     '7': 'Prolif-MPP.1',
     '8': 'MEMP',
     '9': 'Prolif-MPP.2',
     '10': 'Prolif-unspecified',
     '11': 'MPP.1',
     '12': 'GMP/GrP',
     '13': 'Macrophage',
     '14': 'Endothelial',
     '15': 'Fibroblast',
     '16': 'B',
     '17': 'EP.1/EP.2/Erythroid',
     '18': 'Mast',
     '19': 'Pre Pro-B/Pro-B',
     '20': 'MkP',
}
adata.obs['cell_type'] = adata.obs['leiden_sub_R'].map(cluster2annotation).astype('category')

### visualiztion
combined_palette = sns.color_palette("Set1") + sns.color_palette("Set2") + sns.color_palette("Set3") + sns.color_palette("Dark2")
color_dict={
     'HSC':combined_palette[30],
     'MPP.2':combined_palette[31],
     'MPP.3':combined_palette[32],
     'MPP.5':combined_palette[33],
     'MkP':combined_palette[10],
     'MPP.4':combined_palette[34],
     'MPP/GMP':combined_palette[35],
     'Prolif-MPP.1':combined_palette[12],
     'MEMP':combined_palette[9],
     'Prolif-MPP.2':combined_palette[17],
     'Prolif-unspecified':combined_palette[13],
     'MPP.1':combined_palette[36],
     'GMP/GrP':combined_palette[4],
     'Macrophage':combined_palette[19],
     'Endothelial':combined_palette[20],
     'Fibroblast':combined_palette[21],
     'B':combined_palette[23],
     'EP.1/EP.2/Erythroid':combined_palette[1],
     'Mast':combined_palette[11],
     'Pre Pro-B/Pro-B':combined_palette[6],
     'MkP':combined_palette[10]
}
sns.set_palette(combined_palette)
sc.pl.draw_graph(adata,color = ['cell_type'],size = 10,legend_loc='on data',
                 palette = color_dict,frameon=False,add_outline=True, title='')

## Find marker genes
sc.tl.rank_genes_groups(adata, groupby = 'cell_type', method = 'wilcoxon', tie_correct = True, pts = True, key_added = 'wilcoxon_cell_type')
celltype_marker = pd.DataFrame(adata.uns['wilcoxon_cell_type']['names'])
celltype_marker.head(100).to_csv(r'./results/240229.cell_type.wlcxn.csv', index=False)

sc.tl.dendrogram(adata, groupby="cell_type", use_rep="X_scVI")
celltype_markers = celltype_markers[(celltype_markers.pvals_adj < 0.05) & (celltype_markers.logfoldchanges > 1.0) & (celltype_markers.pct_nz_group > 0.05)]
celltype_markers = celltype_markers[~celltype_markers['names'].str.startswith('LOC')]

### extract cell type specific genes
celltype_markers[celltype_markers['group']=='Pre Pro-B/Pro-B'][:5].names.tolist()

genes_to_show = {
'HSC':['TENM4', 'HLF', 'MLLT3'],
'MPP.1':['MAFA-F', 'FOS', 'CD74'],
'MPP.2':['GRID2', 'ITM2C', 'LAT2'],
'MPP.3':['VIM', 'TPPP3', 'CRIP2'],
'MPP.4':['SERPINA5', 'ZBTB16', 'HTR1F'],
'MPP.5':['UNG', 'MCM5', 'IGFBP2'],
'Prolif-MPP.1':['H1-1', 'SPC25', 'MXD3'],
'Prolif-MPP.2':['PIMREG', 'PLK1', 'CENPA'],
'MPP/GMP':['PTMS', 'IGLL1', 'MPO'],
'GMP/GrP':['AZU1', 'PRTN3', 'MS4A3'],
'MkP':['GATA1', 'TESC', 'ITGB3'],
'MEMP':['CA1', 'GATA1', 'RBPMS2'],
'EP.1/EP.2/Erythroid':['PKLR', 'SPTA1', 'ALAS2'],
'Mast':['HDC', 'PRG3', 'NTRK1'],
'Macrophage':['C1QA', 'C1QC', 'CRYBA2'],
'Pre Pro-B/Pro-B':['DNTT', 'ARPP21', 'IL7R'],
'B':['IL36A', 'TNFRSF13B', 'SLAMF7'],
'Endothelial':['COL4A1', 'NPR1', 'FCN3'],
'Fibroblast':['PLAC9', 'SOD3', 'RSPO1'],
'Prolif-unspecified':['H1-4', 'CLC', 'H1-5'],
}

desired_order = ['HSC','MPP.1','MPP.2','MPP.3','MPP.4','MPP.5','Prolif-MPP.1','Prolif-MPP.2','MPP/GMP','GMP/GrP', 'MkP','MEMP','EP.1/EP.2/Erythroid','Mast','Macrophage','Pre Pro-B/Pro-B','B','Endothelial','Fibroblast','Prolif-unspecified']
sc.pl.dotplot(adata, genes_to_show, groupby="cell_type", var_group_labels='cell_type', categories_order=desired_order, standard_scale='var', log=True, dendrogram=False,colorbar_title='Expression', size_title = '')
adata.write('./anndataobjs/Mf.HSC.scVI.processed.h5ad')
