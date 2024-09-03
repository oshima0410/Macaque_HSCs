#!/usr/bin/env python
# coding: utf-8

# Differentially expressed genes
# Diffxpy

import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
import diffxpy.api as de

warnings.simplefilter(action = "ignore", category = FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="scanpy")
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")

sc.logging.print_header()


## data import
adata= sc.read_h5ad('./anndataobjs/Mf.HSC.scVI.processed.h5ad')

## subset HSC and MPPs
cluster2annotation = {
     '0': 'HSC',
     '1': 'MPPs',
     '2': 'MPPs',
     '3': 'MPPs',
     '4': 'MkP',
     '5': 'MPPs',
     '6': 'MPP/GMP',
     '7': 'Prolif-MPP.1',
     '8': 'MEMP',
     '9': 'Prolif-MPP.2',
     '10': 'Prolif-unspecified',
     '11': 'MPPs',
     '12': 'GMP/GrP',
     '13': 'Macrophage',
     '14': 'Endothelial',
     '15': 'Fibroblast',
     '16': 'T',
     '17': 'EP.1/EP.2/Erythroid',
     '18': 'Mast',
     '19': 'Pre Pro-B/Pro-B',
     '20': 'MkP',
}
adata.obs['cell_type2'] = adata.obs['leiden_sub_R'].map(cluster2annotation).astype('category')

subset = adata[adata.obs['cell_type2'].isin(['HSC', 'MPPs'])].copy()
subset.X = subset.X.toarray()
len(subset.var)

## wald test
res = de.test.wald(data=subset,
             formula_loc= '~ 1 + cell_type2',
             factor_loc_totest='cell_type2'
                  )

dedf_hsc_mpp = res.summary().sort_values('log2fc', ascending = False).reset_index(drop = True)
most_up = dedf_hsc_mpp.iloc[0].gene
i = np.where(subset.var_names == most_up)[0][0]

a = subset[subset.obs.cell_type2 == 'HSC'].X[:, i]
b = subset[subset.obs.cell_type2 == 'MPPs'].X[:, i]
print(f"{most_up} expression:")
print(f"HSC: {a.mean()}")
print(f"MPPs: {b.mean()}")

dedf_hsc_mpp['log2fc'] = dedf_hsc_mpp['log2fc']*-1
dedf_hsc_mpp = dedf_hsc_mpp.sort_values('log2fc', ascending = False).reset_index(drop = True)
dedf_hsc_mpp = dedf_hsc_mpp[(dedf_hsc_mpp.qval < 0.05) & (abs(dedf_hsc_mpp.log2fc) > .5)]
dedf_hsc_mpp = dedf_hsc_mpp[dedf_hsc_mpp['mean'] > 0.15]
dedf_hsc_mpp.to_csv("./results/dedf_hsc_vs_mpp_240123.csv")
