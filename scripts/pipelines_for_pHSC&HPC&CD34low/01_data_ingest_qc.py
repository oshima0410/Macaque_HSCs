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
    'HSC_4_FL165HT/outs/filtered_feature_bc_matrix/',
    'HSPC_4_FL165HT/outs/filtered_feature_bc_matrix/',
    'CD34low_4_FL165HT/outs/filtered_feature_bc_matrix/',
    'HSC_4_FL165ST/outs/filtered_feature_bc_matrix/',
    'HSPC_4_FL165ST/outs/filtered_feature_bc_matrix/',
    'CD34low_4_FL165ST/outs/filtered_feature_bc_matrix/',
    'HSC_10_FL167M170F/outs/filtered_feature_bc_matrix/',
    'HPC_10_FL167M170F/outs/filtered_feature_bc_matrix/',
    'CD34low_10_FL167M170F/outs/filtered_feature_bc_matrix/',
     ]
filenames = [data_path + subdir for subdir in subdirectories]
adatas = [sc.read_10x_mtx(filename, var_names = 'gene_symbols', cache=True) for filename in filenames]
adata = adatas[0].concatenate(adatas[1:], batch_categories = [
adata = adatas[0].concatenate(adatas[1:], batch_categories = [
    'HSC_HSC_4_E58M_FL_HT_CY165',
    'HPC_HPC_4_E58M_FL_HT_CY165',
    'CD34low_4_E58M_FL_HT_CY165',
    'HSC_HSC_4_E58M_FL_ST_CY165',
    'HPC_HPC_4_E58M_FL_ST_CY165',
    'CD34low_4_E58M_FL_ST_CY165',
    'HSC_HSC_10_E58MF_FL_HT_CY167M170F',
    'HPC_HPC_10_E58MF_FL_HT_CY167M170F',
    'CD34low_10_E58MF_FL_HT_CY167M170F',
])

## set and rename adata.obs
adata.obs['population'] = [i.split()[0][:3] for i in adata.obs.batch]
adata.obs['population'] = adata.obs['population'].replace(to_replace='CD3', value='CD34low', regex=False)

adata.obs['exp'] = [i.split()[0][8:10] for i in adata.obs.batch]
adata.obs['exp'] = adata.obs['exp'].replace(to_replace='_', value='', regex=True)

adata.obs['donor'] = [i.split()[0][21:] for i in adata.obs.batch]
adata.obs.donor.unique()
adata.obs['donor'] = adata.obs['donor'].replace(to_replace='T_', value='', regex=True)

adata.obs['freshthaw'] = [i.split()[0][8:10] for i in adata.obs.batch]
adata.obs['freshthaw'] = adata.obs['freshthaw'].replace(to_replace='4_', value='fresh', regex=True)
adata.obs['freshthaw'] = adata.obs['freshthaw'].replace(to_replace='10', value='thaw', regex=True)

adata.obs['kit'] = [i.split()[0][18:20] for i in adata.obs.batch]
adata.obs['kit'] = adata.obs['kit'].replace(to_replace='L_', value='HT', regex=True)

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
adata.write('./anndataobjs/Mf.E58CD34.combined.pre-preprocessed.h5ad')
