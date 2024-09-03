# Atlas of Cynomolgus Macaque Hematopoiesis
This repository contains analysis scripts used to explore scRNA-seq dataset for
"Oshima S. et al., Atlas of Cynomolgus Macaque Hematopoiesis"
doi: https://doi.org/10.1101/2024.04.22.590220

## Bundles
The bundle includes:
- tools for building scanpy object from multiple CellRanger count tables
- doublet removal
- subsampling
- batch correction
- computing data reduction coordinates like UMAP, paga, force directed graph and diffusion map
- trajectory analysis
- exploration of marker genes
- cell type or age-tissue comparison metrics
- cell cycle

## Sample species and reference genomes
- cynomolgus macaque ( Macaca fascicularis) [MFA1912RKSv2]


## Data accessibility
There are no restrictions on data availability for data presented in this study. 
The scRNA-seq raw data have been deposited at Sequence Read Archive (SRA) : PRJNA1090143, and the processed datasets have been uploaded at Gene Expression Omnibus (GEO) : GSE262140.
They are publicly available and can be downloaded as of the publication.

Abbreviation: 
E; embryonic day, FL; fetal liver, FBM; fetal bone marrow, pHSC; phenotipic hematopoietic stem cell, pHPC; phenotipic hematopoietic progenitor cell, pHSPC; phenotipic hematopoietic stem and progenitor cell, pHSC/MPP; phenotipic hematopoietic stem cells and multipotent progenitors, CD34low; CD34 low population cells

Sample names at GEO stands for;  FACS sorted population - cynomolgus donor ID - embryonic day and tissue - 10X genomics kit name (HT; high troughput(Chromium Next GEM Single Cell 5’ HT Reagent
Kits v.2) or ST; standard kit(Chromium Next GEM Single Cell 5’ Reagent Kits v.2))

## Prerequisites
### Python version 3.8 (or later)
- Scanpy (version 1.9.5)
- DoubletDetection (version 4.2)
- scVI (version 1.0.4)
- pytorch (version 2.1.3)
- jax (version 0.4.23)
- fa2 (version 0.3.5)
- diffxpy (version 0.7.4)
- numpy (version 1.24.4)
- pandas (version 2.1.4)
- igraph (version 0.11.3)
- forceatlas2-python (version 1.1)
- matplotlib (version 3.8.0)

### R version 4.2.2 
- Seurat (version 5.0.3)
- Monocle 3 (version 1.3.1)
- ClusterProfiler (version 4.10.0)
- dplyr (version 1.1.4)
- ggplot2 (version 3.5.1)
- Dose (version 3.28.2)
- org.Mmu.eg.db (version 3.18.0)
- tidyverse (version 2.0.0)

