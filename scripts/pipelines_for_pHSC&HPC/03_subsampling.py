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

# Subsampling based on FACS gate frequency
adata = sc.read_h5ad('../anndataobjs/Mf.HSPC.combined.pre-preprocessed.post-doublet.h5ad')
adata = adata[adata.obs.doublet == 0.0, :]
adata.obs['batch'].value_counts()

a = adata[adata.obs['batch'] == 'HSC_HSC_1_E58F_FL_ST_CY153']
b = adata[adata.obs['batch'] == 'HPC_HPC_1_E58F_FL_ST_CY153']
c = adata[adata.obs['batch'] == 'HSC_HSC_3_E58F_FL_HT_CY158']
d = adata[adata.obs['batch'] == 'HPC_HPC_3_E58F_FL_HT_CY158']
e = adata[adata.obs['batch'] == 'HSC_HSC_4_E58M_FL_HT_CY165']
f = adata[adata.obs['batch'] == 'HPC_HPC_4_E58M_FL_HT_CY165']
g = adata[adata.obs['batch'] == 'HSC_HSC_4_E58M_FL_ST_CY165']
h = adata[adata.obs['batch'] == 'HPC_HPC_4_E58M_FL_ST_CY165']
i = adata[adata.obs['batch'] == 'HSC_HSC_5_E58M_FL_HT_CY167']
j = adata[adata.obs['batch'] == 'HPC_HPC_5_E58M_FL_HT_CY167']
k = adata[adata.obs['batch'] == 'HSC_HSC_5_E58M_FL_HT_CY168']
l = adata[adata.obs['batch'] == 'HPC_HPC_5_E58M_FL_HT_CY168']
m = adata[adata.obs['batch'] == 'HSC_HSC_6_9yoF_BM_HT_CE2032']
n = adata[adata.obs['batch'] == 'HPC_HPC_6_9yoF_BM_HT_CE2032']
o = adata[adata.obs['batch'] == 'HSC_HSC_7_9yoF_BM_HT_CE1959']
p = adata[adata.obs['batch'] == 'HPC_HPC_7_9yoF_BM_HT_CE1959']
q = adata[adata.obs['batch'] == 'HSC_HSC_8_E96F_FL_HT_CY145']
r = adata[adata.obs['batch'] == 'HPC_HPC_8_E96F_FL_HT_CY145']
s = adata[adata.obs['batch'] == 'HSC_HSC_8_E96F_BM_HT_CY145']
t = adata[adata.obs['batch'] == 'HPC_HPC_8_E96F_BM_HT_CY145']
u = adata[adata.obs['batch'] == 'HSC_HSC_9_10yF_BM_HT_CE2278+CE2064']
v = adata[adata.obs['batch'] == 'HPC_RAN_9_10yF_BM_HT_CE2278+CE2064']
w = adata[adata.obs['batch'] == 'HPC_RAP_9_10yF_BM_HT_CE2278+CE2064']
x = adata[adata.obs['batch'] == 'HSC_HSC_10_E58MF_FL_HT_CY167M170F']
y = adata[adata.obs['batch'] == 'HPC_HPC_10_E58MF_FL_HT_CY167M170F']
z = adata[adata.obs['batch'] == 'CD34hig_11_E93M_FL_HT_CY189']
aa = adata[adata.obs['batch'] == 'CD34hig_11_E93M_BM_HT_CY189']
ab = adata[adata.obs['batch'] == 'HSC_HSC_12_E121M_FL_HT_CY195']
ac = adata[adata.obs['batch'] == 'HPC_HPC_12_E121M_FL_HT_CY195']
ad = adata[adata.obs['batch'] == 'HSC_HSC_12_E121M_BM_HT_CY195']
ae = adata[adata.obs['batch'] == 'HPC_HPC_12_E121M_BM_HT_CY195']

## This is reflecting the actual ratio of pHSC vs pHPC in the sample
## Please refer to supplementary table sheet 2
sc.pp.subsample(a,n_obs = 1638)
sc.pp.subsample(b,n_obs = 8786)
sc.pp.subsample(c,n_obs = 1217)
sc.pp.subsample(d,n_obs = 8428)
sc.pp.subsample(e,n_obs = 1048)
sc.pp.subsample(f,n_obs = 5421)
sc.pp.subsample(g,n_obs = 658)
sc.pp.subsample(h,n_obs = 3405)
sc.pp.subsample(i,n_obs = 3819)
sc.pp.subsample(j,n_obs = 11583)
sc.pp.subsample(k,n_obs = 3216)
sc.pp.subsample(l,n_obs = 9754)
sc.pp.subsample(m,n_obs = 632)
sc.pp.subsample(n,n_obs = 3357)
sc.pp.subsample(o,n_obs = 108)
sc.pp.subsample(p,n_obs = 839)
sc.pp.subsample(q,n_obs = 2984)
sc.pp.subsample(r,n_obs = 13311)
sc.pp.subsample(s,n_obs = 2806)
sc.pp.subsample(t,n_obs = 15203)
sc.pp.subsample(u,n_obs = 2792)
sc.pp.subsample(v,n_obs = 2849)
sc.pp.subsample(w,n_obs = 16842)
sc.pp.subsample(x,n_obs = 2317)
sc.pp.subsample(y,n_obs = 15727)
sc.pp.subsample(aa,n_obs = 8490)
sc.pp.subsample(cc,n_obs = 9691)
sc.pp.subsample(dd,n_obs = 861)
sc.pp.subsample(ee,n_obs = 3651)
sc.pp.subsample(ff,n_obs = 3316)
sc.pp.subsample(gg,n_obs = 14619)

adata = a.concatenate(b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,aa,ab,ac,ad,ae,
                      batch_categories= ['HSC_HSC_1_E58F_FL_ST_CY153',
                                          'HPC_HPC_1_E58F_FL_ST_CY153',
                                          'HSC_HSC_3_E58F_FL_HT_CY158',
                                          'HPC_HPC_3_E58F_FL_HT_CY158',
                                          'HSC_HSC_4_E58M_FL_HT_CY165',
                                          'HPC_HPC_4_E58M_FL_HT_CY165',
                                          'HSC_HSC_4_E58M_FL_ST_CY165',
                                          'HPC_HPC_4_E58M_FL_ST_CY165',
                                          'HSC_HSC_5_E58M_FL_HT_CY167',
                                          'HPC_HPC_5_E58M_FL_HT_CY167',
                                          'HSC_HSC_5_E58M_FL_HT_CY168',
                                          'HPC_HPC_5_E58M_FL_HT_CY168',
                                          'HSC_HSC_6_9yoF_BM_HT_CE2032',
                                          'HPC_HPC_6_9yoF_BM_HT_CE2032',
                                          'HSC_HSC_7_9yoF_BM_HT_CE1959',
                                          'HPC_HPC_7_9yoF_BM_HT_CE1959',
                                          'HSC_HSC_8_E96F_FL_HT_CY145',
                                          'HPC_HPC_8_E96F_FL_HT_CY145',
                                          'HSC_HSC_8_E96F_BM_HT_CY145',
                                          'HPC_HPC_8_E96F_BM_HT_CY145',
                                          'HSC_HSC_9_10yF_BM_HT_CE2278+CE2064',
                                          'HPC_RAN_9_10yF_BM_HT_CE2278+CE2064',
                                          'HPC_RAP_9_10yF_BM_HT_CE2278+CE2064',
                                          'HSC_HSC_10_E58MF_FL_HT_CY167M170F',
                                          'HPC_HPC_10_E58MF_FL_HT_CY167M170F',
                                          'CD34hig_11_E93M_FL_HT_CY189',
                                          'CD34hig_11_E93M_BM_HT_CY189',
                                          'HSC_HSC_12_E121M_FL_HT_CY195',
                                          'HPC_HPC_12_E121M_FL_HT_CY195',
                                          'HSC_HSC_12_E121M_BM_HT_CY195',
                                          'HPC_HPC_12_E121M_BM_HT_CY195'])

adata.write('./anndataobjs/Mf.HSPC.combined.sub-preprocessed.h5ad')
