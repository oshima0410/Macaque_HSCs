#!/usr/bin/env python
# coding: utf-8

#　Explore surface proteins detected in our scRNA-seq, referencing the human database


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

## import surface proteins list of human
hSurfaceome = pd.read_csv('./results/mart_export_human_surface_proteins.csv')
hSurfList = hSurfaceome['HGNC symbol'].unique().tolist()
hSurfList = [str(x) for x in hSurfList]
hSurfList = [x for x in hSurfList if x != 'nan']
len(hSurfList)


## extract Macaque surface proteins from adata
gene_list = adata.var_names
def find_genes_list(*genes):
    list = []
    for g in genes:
        list.extend([i for i in gene_list if i.split()[0].endswith(tuple('' + g))])
    return(list)
MfSurfList = find_genes_list(pd.DataFrame(hSurfList)[0])
macaque_surList = pd.DataFrame(find_genes_list(pd.DataFrame(hSurfList)[0]))
macaque_surList.to_csv('./MfSurfList.csv')

## Detecting HSC specific surface proteins, referencing diffxpy wald test
hsc_genes = pd.read_csv("./results/dedf_hsc_vs_mpp_240123.csv")
gene_list = hsc_genes['gene'].unique().tolist()
Mf_hsc_SurfList = find_genes_list(pd.DataFrame(hSurfList)[0])
macaque_surList = pd.DataFrame(find_genes_list(pd.DataFrame(hSurfList)[0]))
macaque_surList.to_csv('./Mf_HSC_amongAllgenes_SurfList.csv')

### matrixplot
cluster2annotation = {
     '0': 'HSC',
     '1': 'MPPs',
     '2': 'MPPs', 
     '3': 'MPPs',
     '4': '',
     '5': 'MPPs',
     '6': '',
     '7': '',
     '8': '',
     '9': '',
     '10': '',
     '11': 'MPPs',
     '12': '',
     '13': '',
     '14': '',
     '15': '',
     '16': '',
     '17': '',
     '18': '',
     '19': '',
     '20': '',
}
adata.obs['cell_type2'] = adata.obs['leiden_sub_R'].map(cluster2annotation).astype('category')

adata_sub = adata[adata.obs['cell_type2'].isin(['HSC', 'MPPs'])].copy()
tmp2 = pd.DataFrame(adata_sub.obs[["cell_type", "cell_type2","agetissue"]])
def create_column3(row):
    if row['cell_type2'] == 'HSC':
        return f'HSC_{row["agetissue"]}'
    else:
        return f'MPPs_{row["agetissue"]}'

tmp2['cell_type2_agetissue'] = tmp2.apply(create_column3, axis=1)
tmp2["cell_type_agetissue"] = tmp2["cell_type"].str.cat(tmp2["agetissue"], sep='_')

adata_sub.obs["cell_type_agetissue"] = tmp2["cell_type_agetissue"]
adata_sub.obs["cell_type2_agetissue"] = tmp2["cell_type2_agetissue"]

adata_sub.obs['cell_type_agetissue'] = pd.Categorical(adata_sub.obs['cell_type_agetissue'])
adata_sub.obs['cell_type2_agetissue'] = pd.Categorical(adata_sub.obs['cell_type2_agetissue'])

sc.tl.dendrogram(adata_sub, groupby="cell_type_agetissue", use_rep="X_scVI")

ene_names = {'HSC':['PTPRM','TEK','CNTNAP2','TENM4','SLC8A1','LRFN5','ABCB1','MDGA2','IL6ST','IL12RB2'],
              'MPP':['SELL']}

desired_order = ['HSC_Early2nd_FL', 'HSC_Late2nd_FL', 'HSC_Early3rd_FL', 'HSC_Late2nd_BM', 'HSC_Early3rd_BM', 'HSC_Adult_BM',
                 'MPP.1_Early2nd_FL', 'MPP.1_Late2nd_FL', 'MPP.1_Early3rd_FL', 'MPP.1_Late2nd_BM', 'MPP.1_Early3rd_BM', 'MPP.1_Adult_BM',
                'MPP.2_Early2nd_FL', 'MPP.2_Late2nd_FL', 'MPP.2_Early3rd_FL', 'MPP.2_Late2nd_BM', 'MPP.2_Early3rd_BM', 'MPP.2_Adult_BM',
                'MPP.3_Early2nd_FL', 'MPP.3_Late2nd_FL', 'MPP.3_Early3rd_FL', 'MPP.3_Late2nd_BM', 'MPP.3_Early3rd_BM', 'MPP.3_Adult_BM',
                'MPP.4_Early2nd_FL', 'MPP.4_Late2nd_FL', 'MPP.4_Early3rd_FL', 'MPP.4_Late2nd_BM', 'MPP.4_Early3rd_BM', 'MPP.4_Adult_BM',
                'MPP.5_Early2nd_FL', 'MPP.5_Late2nd_FL', 'MPP.5_Early3rd_FL', 'MPP.5_Late2nd_BM', 'MPP.5_Early3rd_BM', 'MPP.5_Adult_BM']

sc.pl.matrixplot(adata_sub, gene_names, 'cell_type_agetissue', dendrogram=False, cmap='Blues',
                 categories_order = desired_order,
                 standard_scale='var', colorbar_title='column scaled\nexpression',swap_axes=True )
