#!/usr/bin/env python
# coding: utf-8

# Doubletdetection

import numpy as np
import doubletdetection
import scanpy as sc
import matplotlib.pyplot as plt

sc.settings.n_jobs=8
sc.set_figure_params()
get_ipython().run_line_magic('matplotlib', 'inline')

## data import
adata = sc.read_h5ad('./anndataobjs/Mf.E58CD34.combined.pre-preprocessed.h5ad')
adata.var_names_make_unique()

## calculate doublet score
clf = doubletdetection.BoostClassifier(
    n_iters=10,
    clustering_algorithm="leiden",
    standard_scaling=True,
    pseudocount=0.1,
    n_jobs=-1,
)
doublets = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
doublet_score = clf.doublet_score()

adata.obs["doublet"] = doublets
adata.obs["doublet_score"] = doublet_score
adata.obs['doublet'].value_counts()

## visualize the results
f = doubletdetection.plot.convergence(clf, save='convergence_test.pdf', show=True, p_thresh=1e-16, voter_thresh=0.5)
sc.pl.violin(adata, "doublet_score")

## save
adata.write('./anndataobjs/Mf.E58CD34.combined.pre-preprocessed.post-doublet.h5ad')
