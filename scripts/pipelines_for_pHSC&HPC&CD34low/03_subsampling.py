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
adata = sc.read_h5ad('../anndataobjs/Mf.E58CD34.combined.pre-preprocessed.post-doublet.h5ad')
adata = adata[adata.obs.doublet == 0.0, :]
adata.obs['batch'].value_counts()

tmp1 = adata[adata.obs['batch'] == 'HSC_HSC_4_E58M_FL_HT_CY165']
tmp2 = adata[adata.obs['batch'] == 'HPC_HPC_4_E58M_FL_HT_CY165']
tmp3 = adata[adata.obs['batch'] == 'CD34low_4_E58M_FL_HT_CY165']
tmp4 = adata[adata.obs['batch'] == 'HSC_HSC_4_E58M_FL_ST_CY165']
tmp5 = adata[adata.obs['batch'] == 'HPC_HPC_4_E58M_FL_ST_CY165']
tmp6 = adata[adata.obs['batch'] == 'CD34low_4_E58M_FL_ST_CY165']
tmp7 = adata[adata.obs['batch'] == 'HSC_HSC_10_E58MF_FL_HT_CY167M170F']
tmp8 = adata[adata.obs['batch'] == 'HPC_HPC_10_E58MF_FL_HT_CY167M170F']
tmp9 = adata[adata.obs['batch'] == 'CD34low_10_E58MF_FL_HT_CY167M170F']

## This is reflecting the actual ratio of pHSC vs pHPC in the sample
## Please refer to supplementary table sheet 2
sc.pp.subsample(tmp1,n_obs = 1039)
sc.pp.subsample(tmp2,n_obs = 5379)
sc.pp.subsample(tmp3,n_obs = 2077)
sc.pp.subsample(tmp4,n_obs = 648)
sc.pp.subsample(tmp5,n_obs = 3355)
sc.pp.subsample(tmp6,n_obs = 1295)
sc.pp.subsample(tmp7,n_obs = 2315)
sc.pp.subsample(tmp8,n_obs = 15719)
sc.pp.subsample(tmp9,n_obs = 7951)

adata = tmp1.concatenate(tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,batch_categories= ['HSC_HSC_4_E58M_FL_HT_CY165',
                                                                     'HPC_HPC_4_E58M_FL_HT_CY165',
                                                                     'CD34low_4_E58M_FL_HT_CY165',
                                                                     'HSC_HSC_4_E58M_FL_ST_CY165',
                                                                     'HPC_HPC_4_E58M_FL_ST_CY165',
                                                                     'CD34low_4_E58M_FL_ST_CY165',
                                                                     'HSC_HSC_10_E58MF_FL_HT_CY167M170F',
                                                                     'HPC_HPC_10_E58MF_FL_HT_CY167M170F',
                                                                     'CD34low_10_E58MF_FL_HT_CY167M170F'])

adata.write('./anndataobjs/Mf.E58CD34.combined.sub-preprocessed.h5ad')
