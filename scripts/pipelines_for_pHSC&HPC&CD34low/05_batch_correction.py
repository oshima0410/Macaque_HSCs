#!/usr/bin/env python
# coding: utf-8

# Batch correction

import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import scvi
import torch
import jax
warnings.simplefilter(action = "ignore", category = FutureWarning)
torch.cuda.is_available()

## data import
adata = sc.read_h5ad('./anndataobjs/Mf.E58CD34.combined.preprocessed.h5ad')
adata = adata[:, adata.var.highly_variable]
adata = adata.copy()

## run batch correction
## in this analysis, batch means independent experimental day
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["exp"],
    continuous_covariate_keys=["pct_counts_mt", "pct_counts_ribo"])

model = scvi.model.SCVI(adata, n_hidden=256, n_latent=80, n_layers=2, dropout_rate=0.1)
torch.set_float32_matmul_precision("highest")
model.train()

latent = model.get_latent_representation()
adata.obsm["X_scVI"] = model.get_latent_representation()
adata.layers["scvi_normalized"] = model.get_normalized_expression(library_size=10e4)

## save
model.save("my_model_240304/")
adata.write('./anndataobjs/Mf.E58CD34.scVI.combined.processed.h5ad')
