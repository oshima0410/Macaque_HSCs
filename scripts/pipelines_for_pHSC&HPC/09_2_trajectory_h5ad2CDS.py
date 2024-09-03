#!/usr/bin/env python
# coding: utf-8

# Trajectory analysis
# Convert Scanpy object (h5ad) to CellDataSet object (CDS)

import anndata2ri
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
from igraph import *
import networkx as nx
import seaborn as sns

anndata2ri.activate()

get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('R', 'd = as.data.frame( R.Version() )')
get_ipython().run_line_magic('R', "d = d['version.string']")

get_ipython().run_cell_magic('R', '', '\nlibrary(monocle3)\nlibrary(ggplot2)\nlibrary(dplyr)\n')

adata = sc.read_h5ad('./anndataobjs/Pseudotime_fdg_sub_240201.h5ad')
data_mat_memp = adata.X.T
data_mat_memp
var_memp=adata.var.copy()
var_memp['gene_short_name'] = var_memp[var_memp.columns[0]]
obs_memp=adata.obs.copy()

get_ipython().run_cell_magic('R', '-i adata', 'adata\n')
get_ipython().run_cell_magic('R', '-i data_mat_memp -i obs_memp -i var_memp', '\ncolnames(data_mat_memp) <- rownames(obs_memp) \nrownames(data_mat_memp) <- rownames(var_memp) \n\nie_regions_cds <- new_cell_data_set(data_mat_memp, \n                                    cell_metadata=obs_memp, \n                                    gene_metadata=var_memp)\n\nprint("printing ie_regions_cds")\nprint(ie_regions_cds)\n')
get_ipython().run_cell_magic('R', '', '\n#Filter highly variable genes from our analysis \nprint("filter hvgs")\nhvg_mask = fData(ie_regions_cds)$highly_variable\nie_regions_cds <- ie_regions_cds[hvg_mask,]\n\n# print ie_regions_cds\nprint("printing cds again")\nprint(ie_regions_cds)\n')
get_ipython().run_cell_magic('R', '', '# save MEMP-pathway dataset in .RDS format for monocle3 \nsaveRDS(ie_regions_cds, "./RDS/Pseudotime_fdg_subset_240201.RDS")\n')

## Do the same thing for EP, Mast, MkP-direct, MkP-indirect pathways.
## Output files for each are below, respectively
### EP pathway
#"./RDS/Pseudotime_fdg_EP_direct_subset_240229.RDS"

### Mast mathway
#"./RDS/Pseudotime_fdg_Mast_direct_subset_240229.RDS"

### MkP indirect pathway
#"./RDS/Pseudotime_fdg_MkP_indirect_subset_240206.RDS"

### MkP direct pathway
#"./RDS/Pseudotime_fdg_MkP_direct_subset_240206.RDS"
