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
import networkx as nx
import seaborn as sns

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=120, dpi_save=300)
sc.settings.n_jobs = 30
pd.options.display.max_columns = None

warnings.simplefilter(action = "ignore", category = FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="scanpy")
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")

# Data ingest
data_path = '/home/oshima/10xAnalysis/cHSC_10x/10XResults/'
subdirectories = [
    'HSC_1_FL153/outs/filtered_feature_bc_matrix/',
    'HSPC_1_FL153/outs/filtered_feature_bc_matrix/',
    'HSC_3_FL158/outs/filtered_feature_bc_matrix/',
    'HSPC_3_FL158/outs/filtered_feature_bc_matrix/',
    'HSC_4_FL165HT/outs/filtered_feature_bc_matrix/',
    'HSPC_4_FL165HT/outs/filtered_feature_bc_matrix/',
    'HSC_4_FL165ST/outs/filtered_feature_bc_matrix/',
    'HSPC_4_FL165ST/outs/filtered_feature_bc_matrix/',
    'HSC_5_FL167/outs/filtered_feature_bc_matrix/',
    'HSPC_5_FL167/outs/filtered_feature_bc_matrix/',
    'HSC_5_FL168/outs/filtered_feature_bc_matrix/',
    'HSPC_5_FL168/outs/filtered_feature_bc_matrix/',
    'HSC_6_CE2032F/outs/filtered_feature_bc_matrix/',
    'HSPC_6_CE2032F/outs/filtered_feature_bc_matrix/',
    'HSC_7_CE1959F/outs/filtered_feature_bc_matrix/',
    'HPC_7_CE1959F/outs/filtered_feature_bc_matrix/',
    'HSC_8_FL145/outs/filtered_feature_bc_matrix/',
    'HPC_8_FL145/outs/filtered_feature_bc_matrix/',
    'HSC_8_FBM145/outs/filtered_feature_bc_matrix/',
    'HPC_8_FBM145/outs/filtered_feature_bc_matrix/',
    'HSC_9_ABM/outs/filtered_feature_bc_matrix/',
    'HPC_45RAn_9_ABM/outs/filtered_feature_bc_matrix/',
    'HPC_45RAp_ABM/outs/filtered_feature_bc_matrix/',
    'HSC_10_FL167M170F/outs/filtered_feature_bc_matrix/',
    'HPC_10_FL167M170F/outs/filtered_feature_bc_matrix/',
    'CD34hi_11_FL189/outs/filtered_feature_bc_matrix/',
    'CD34hi_11_FBM189/outs/filtered_feature_bc_matrix/',
    'HSC_12_FL195/outs/filtered_feature_bc_matrix/',
    'HPC_12_FL195/outs/filtered_feature_bc_matrix/',
    'HSC_12_FBM195/outs/filtered_feature_bc_matrix/',
    'HPC_12_FBM195/outs/filtered_feature_bc_matrix/',
    ]
filenames = [data_path + subdir for subdir in subdirectories]
adatas = [sc.read_10x_mtx(filename, var_names = 'gene_symbols', cache=True) for filename in filenames]
adata = adatas[0].concatenate(adatas[1:], batch_categories = [
    'HSC_HSC_1_E58F_FL_ST_CY153',
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
    'HPC_HPC_12_E121M_BM_HT_CY195',
])

## set and rename adata.obs
adata.obs['population'] = [i.split()[0][:3] for i in adata.obs.batch]
adata.obs['population'] = adata.obs['population'].replace(to_replace='CD3', value='HSC+HPC', regex=False)

adata.obs['subpopulation'] = [i.split()[0][:7] for i in adata.obs.batch]
adata.obs['subpopulation'] = adata.obs['subpopulation'].replace(to_replace='RAP', value='CD45RA_positive', regex=True)
adata.obs['subpopulation'] = adata.obs['subpopulation'].replace(to_replace='RAN', value='CD45RA_negative', regex=True)
adata.obs['subpopulation'] = adata.obs['subpopulation'].replace(to_replace='CD34hig', value='HSC+HPC', regex=False)
adata.obs['subpopulation'] = adata.obs['subpopulation'].replace(to_replace='_HSC', value='', regex=True)
adata.obs['subpopulation'] = adata.obs['subpopulation'].replace(to_replace='_HPC', value='', regex=True)

adata.obs['exp'] = [i.split()[0][8:10] for i in adata.obs.batch]
adata.obs['exp'] = adata.obs['exp'].replace(to_replace='_', value='', regex=True)

adata.obs['sex'] = [i.split()[0][13:17] for i in adata.obs.batch]
adata.obs['sex'] = adata.obs['sex'].replace(to_replace='_FL', value='', regex=True)
adata.obs['sex'] = adata.obs['sex'].replace(to_replace='_BM', value='', regex=True)
adata.obs['sex'] = adata.obs['sex'].replace(to_replace='8MF_', value='M+F', regex=True)
adata.obs['sex'] = adata.obs['sex'].replace(to_replace='3M_F', value='M', regex=True)
adata.obs['sex'] = adata.obs['sex'].replace(to_replace='3M_B', value='M', regex=True)
adata.obs['sex'] = adata.obs['sex'].replace(to_replace='21M_', value='M', regex=True)

adata.obs['age'] = [i.split()[0][10:15] for i in adata.obs.batch]
adata.obs['age'] = adata.obs['age'].replace(to_replace='F_', value='', regex=True)
adata.obs['age'] = adata.obs['age'].replace(to_replace='M_', value='', regex=True)
adata.obs['age'] = adata.obs['age'].replace(to_replace='_E58M', value='E58', regex=True)
adata.obs['age'] = adata.obs['age'].replace(to_replace='_E93M', value='E93', regex=True)
adata.obs['age'] = adata.obs['age'].replace(to_replace='10y', value='10yo', regex=True)
adata.obs['age'] = adata.obs['age'].replace(to_replace='_E121', value='E121', regex=True)

adata.obs['tissue'] = [i.split()[0][10:19] for i in adata.obs.batch]
adata.obs['tissue'] = adata.obs['tissue'].replace(to_replace='_S', value='', regex=True)
adata.obs['tissue'] = adata.obs['tissue'].replace(to_replace='_H', value='', regex=True)
adata.obs['tissue'] = adata.obs['tissue'].replace(to_replace='F_', value='', regex=True)
adata.obs['tissue'] = adata.obs['tissue'].replace(to_replace='E58MF_FL', value='E58FL', regex=True)
adata.obs['tissue'] = adata.obs['tissue'].replace(to_replace='_E93M_FL_', value='E93FL', regex=True)
adata.obs['tissue'] = adata.obs['tissue'].replace(to_replace='_E93M_BM_', value='E93BM', regex=True)
adata.obs['tissue'] = adata.obs['tissue'].replace(to_replace='_', value='', regex=True)
adata.obs['tissue'] = adata.obs['tissue'].replace(to_replace='E58MFL', value='E58FL', regex=True)
adata.obs['tissue'] = adata.obs['tissue'].replace(to_replace='E121MFL', value='E121FL', regex=True)
adata.obs['tissue'] = adata.obs['tissue'].replace(to_replace='E121MBM', value='E121BM', regex=True)

adata.obs['kit'] = [i.split()[0][18:22] for i in adata.obs.batch]
adata.obs['kit'] = adata.obs['kit'].replace(to_replace='_C', value='', regex=True)
adata.obs['kit'] = adata.obs['kit'].replace(to_replace='L_', value='', regex=True)
adata.obs['kit'] = adata.obs['kit'].replace(to_replace='_HT_', value='HT', regex=True)
adata.obs['kit'] = adata.obs['kit'].replace(to_replace='M_', value='', regex=True)

adata.obs['donor'] = [i.split()[0][21:] for i in adata.obs.batch]
adata.obs['donor'] = adata.obs['donor'].replace(to_replace='_', value='', regex=True)
adata.obs['donor'] = adata.obs['donor'].replace(to_replace='T', value='', regex=True)
adata.obs['donor'] = adata.obs['donor'].replace(to_replace='CY167M170F', value='CY167+CY170', regex=True)

adata.obs['agetissue'] = [i.split()[0][10:18] for i in adata.obs.batch]
adata.obs['agetissue'] = adata.obs['agetissue'].replace(to_replace='E58F_FL_', value='Early2nd_FL', regex=False)
adata.obs['agetissue'] = adata.obs['agetissue'].replace(to_replace='E58M_FL_', value='Early2nd_FL', regex=False)
adata.obs['agetissue'] = adata.obs['agetissue'].replace(to_replace='9yoF_BM_', value='Adult_BM', regex=False)
adata.obs['agetissue'] = adata.obs['agetissue'].replace(to_replace='E96F_FL_', value='Late2nd_FL', regex=False)
adata.obs['agetissue'] = adata.obs['agetissue'].replace(to_replace='E96F_BM_', value='Late2nd_BM', regex=False)
adata.obs['agetissue'] = adata.obs['agetissue'].replace(to_replace='10yF_BM_', value='Adult_BM', regex=False)
adata.obs['agetissue'] = adata.obs['agetissue'].replace(to_replace='_E58MF_F', value='Early2nd_FL', regex=False)
adata.obs['agetissue'] = adata.obs['agetissue'].replace(to_replace='_E93M_FL', value='Late2nd_FL', regex=False)
adata.obs['agetissue'] = adata.obs['agetissue'].replace(to_replace='_E93M_BM', value='Late2nd_BM', regex=False)
adata.obs['agetissue'] = adata.obs['agetissue'].replace(to_replace='_E121M_F', value='Early3rd_FL', regex=False)
adata.obs['agetissue'] = adata.obs['agetissue'].replace(to_replace='_E121M_B', value='Early3rd_BM', regex=False)

# Preprocessing
sc.pl.highest_expr_genes(adata, n_top=20)
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=3)

## annotate the group of mitochondrial, ribosoma, hemoglobin genes as 'mt', 'ribo', 'hb'.
adata.var['mt'] = adata.var_names.str.startswith('KEG98_')
adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL"))
adata.var['hb'] = adata.var_names.str.contains(("^HB[^(P)]"))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'], percent_top=None, log1p=False, inplace=True)

## visualize the metrics and set threshold
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo', 'pct_counts_hb'],
             stripplot=True, jitter=0.4, size= 0.5,
             multi_panel=True, rotation=90, groupby= 'batch')
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

adata = adata[adata.obs.pct_counts_mt < 1.0, :]

## save
adata.write('./anndataobjs/Mf.HSPC.combined.pre-preprocessed.h5ad')
