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
## Visualization of monocle3 metadata on FDG along with pseudotime

＃MkP indirect pathway

## import data
adata_sub= sc.read_h5ad('./anndataobjs/Pseudotime_fdg_MkP_indirect_sub_240206.h5ad')
metadata = pd.read_csv("./results/240206_MkP_indirect_monocle3_pst_metadata.csv")
metadata.index=metadata["Unnamed: 0"]

## add metadata into obj
adata_sub.obs["ds_monocle3_pst_val"] = metadata["cell_color"]
adata_sub.obs["ds_monocle3_pst_val"] = adata_sub.obs["ds_monocle3_pst_val"].replace(np.inf, np.nan)

## check the plots
sc.pl.draw_graph(adata_sub, color=["leiden_sub_R2"],size = 10,frameon=True, add_outline= True)
sc.pl.draw_graph(adata_sub, color=["leiden_sub_R2"], legend_loc='on data',size = 10,frameon=True, add_outline= True)
sc.pl.draw_graph(adata_sub, color=["cell_type"], size = 10,frameon=False, title='', add_outline= True)
sc.pl.draw_graph(adata_sub, color=["ds_monocle3_pst_val"], size = 10,frameon=True, add_outline= True)


## select genes to show (top genes for monocle3 modues)
path_genes =[
'TENM4', #module2
'HLF',
'GUCY1A1',
'SPINK2', #module8
'ATP8B4', #module6
'HTR1F',
'GRID2',
'MSI2', #module3
'ERG',
'TOX', #module10
'AFF1',
'GLS',
'CDK6', #module1
'CHST11',
'CUX1',
'PCLAF', #module7
'RRM2',
'DIAPH3',
'BRCA2',
'RFC3',
'SPC25', #module13
'CDK1',
'CCNA2',
'UBE2C',
'H1-5', #module14
'H1-1',
'GATA1',#module5
'TESC',
'SLC38A8',
'TFR2',
'ANK1',
'HBQ1',
'EHD3',#module12
'THBS1',
'PLXDC2',
'VWF',
'ITGB3',
]

## matrixplot
sc.set_figure_params(scanpy=True, fontsize=16)
sc.pl.matrixplot(adata_sub, var_names=path_genes,
                      num_categories=40,
                      groupby='ds_monocle3_pst_val', swap_axes=True, use_raw=False, figsize=[10,14],log=False,
                 standard_scale = 'var',)


# MkP direct pathway

## import data
adata_sub= sc.read_h5ad('./anndataobjs/Pseudotime_fdg_MkP_direct_sub_240206.h5ad')
metadata = pd.read_csv("./results/240206_MkP_direct_monocle3_pst_metadata.csv")
metadata.index=metadata["Unnamed: 0"]

## add metadata into obj
adata_sub.obs["ds_monocle3_pst_val"] = metadata["cell_color"]
adata_sub.obs["ds_monocle3_pst_val"] = adata_sub.obs["ds_monocle3_pst_val"].replace(np.inf, np.nan)

## check the plots
sc.pl.draw_graph(adata_sub, color=["leiden_sub_R2"],size = 10,frameon=True, add_outline= True)
sc.pl.draw_graph(adata_sub, color=["leiden_sub_R2"], legend_loc='on data',size = 10,frameon=True, add_outline= True)
sc.pl.draw_graph(adata_sub, color=["cell_type"], size = 10,frameon=False, title='', add_outline= True)
sc.pl.draw_graph(adata_sub, color=["ds_monocle3_pst_val"], size = 10,frameon=True, add_outline= True)

## select genes to show (top genes for monocle3 modues)
path_genes =[
'HLF',  #module1
'TENM4',
'ADGRB3',
'CD48',#module4
'RAMP1',
'SPINK2',
'ROBO2', #module2
'GUCY1A1',
'MCTP1',
'RAPGEF2',
'INPP4B', #module3
'MEIS1',
'TRAPPC9',
'MSI2',
'ELMO1',
'IMMP2L',
'CDKAL1', #module5
'AGAP1',
'ZBTB16',
'TFR2',
'STON2',
'MED12L',
'SLC24A3',
'MTSS1',
'PBX1',
'GATA1', #module6
'PDLIM1',
'TESC',
'PLXDC2',
'FADS2',
'EHD3',
'PCLAF',
'VWF',
'LTBP1',
'THBS1',
'RBPMS2',
'ITGB3',
'ITGA2B',
]

## matrixplot
sc.set_figure_params(scanpy=True, fontsize=16)
sc.pl.matrixplot(adata_sub, var_names=path_genes,
                      num_categories=30,
                      groupby='ds_monocle3_pst_val', swap_axes=True, use_raw=False, figsize=[10,14],log=False,
                 standard_scale = 'var',)


# EP pathway

## import data
adata_sub= sc.read_h5ad('./anndataobjs/Pseudotime_fdg_EP_sub_240229.h5ad')
metadata = pd.read_csv("./results/240229_EP_monocle3_pst_metadata.csv")
metadata.index=metadata["Unnamed: 0"]

## add metadata into obj
adata_sub.obs["ds_monocle3_pst_val"] = metadata["cell_color"]
adata_sub.obs["ds_monocle3_pst_val"] = adata_sub.obs["ds_monocle3_pst_val"].replace(np.inf, np.nan)

## check the plots
sc.pl.draw_graph(adata_sub, color=["leiden_sub_R2"],size = 10,frameon=True, add_outline= True)
sc.pl.draw_graph(adata_sub, color=["leiden_sub_R2"], legend_loc='on data',size = 10,frameon=True, add_outline= True)
sc.pl.draw_graph(adata_sub, color=["cell_type"], size = 10,frameon=False, title='', add_outline= True)
sc.pl.draw_graph(adata_sub, color=["ds_monocle3_pst_val"], size = 10,frameon=True, add_outline= True)

## select genes to show (top genes for monocle3 modues)
path_genes =[
'TENM4',#3
'HLF',
'MECOM',
'GNG11', #5
'SPINK2',
'AIF1',
'RCSD1',
'CD200',
'TRIM22',
'ATP8B4', #6
'CSF3R',
'CELF2',
'SPC25', #12
'CDK1',
'CCNA2',
'UBE2C',
'CDCA3',
'NCAPG',
'MKI67',
'TOP2A',
'H1-1', #13
'H1-2',
'H1-3',
'RRM2', #4
'PCLAF',
'PKIG',
'GATA1',
'KLF1',
'TFR2',
'ANK1',
'PRDX2',
'PKLR', #11
'MYL4',
'SPTA1',
'LMNA',
'ALAS2',
]

## matrixplot
sc.set_figure_params(scanpy=True, fontsize=16)
sc.pl.matrixplot(adata_sub, var_names=path_genes,
                      num_categories=30,
                      groupby='ds_monocle3_pst_val', swap_axes=True, use_raw=False, figsize=[10,14],log=False,
                 standard_scale = 'var',)


# Mast pathway
## import data
adata_sub= sc.read_h5ad('./anndataobjs/Pseudotime_fdg_Mast_sub_240229.h5ad')
metadata = pd.read_csv("./results/240229_Mast_monocle3_pst_metadata.csv")
metadata.index=metadata["Unnamed: 0"]

## add metadata into obj
adata_sub.obs["ds_monocle3_pst_val"] = metadata["cell_color"]
adata_sub.obs["ds_monocle3_pst_val"] = adata_sub.obs["ds_monocle3_pst_val"].replace(np.inf, np.nan)

## check the plots
sc.pl.draw_graph(adata_sub, color=["leiden_sub_R2"],size = 10,frameon=True, add_outline= True)
sc.pl.draw_graph(adata_sub, color=["leiden_sub_R2"], legend_loc='on data',size = 10,frameon=True, add_outline= True)
sc.pl.draw_graph(adata_sub, color=["cell_type"], size = 10,frameon=False, title='', add_outline= True)
sc.pl.draw_graph(adata_sub, color=["ds_monocle3_pst_val"], size = 10,frameon=True, add_outline= True)

## select genes to show (top genes for monocle3 modues)
path_genes =[
'HOPX', #module1
'SPINK2',
'IFITM1',
'EGFL7',
'CDK6', #module3
'DIAPH3',
'BRCA2',
'RFC3',
'SPC25', #module14
'CDK1',
'CCNA2',
'H1-1', #module15
'H1-2',
'H1-3',
'ZEB2', #module7
'FNDC3B',
'AGAP1',
'PLD1',
'SLC24A3',
'ALOX15', #module4
'PRG3',
'PRG2',
'HDC',
'CMA1',
'MS4A2',
]

## matrixplot
sc.set_figure_params(scanpy=True, fontsize=16)
sc.pl.matrixplot(adata_sub, var_names=path_genes,
                      num_categories=22,
                      groupby='ds_monocle3_pst_val', swap_axes=True, use_raw=False, figsize=[10,14],log=False,
                 standard_scale = 'var',)
