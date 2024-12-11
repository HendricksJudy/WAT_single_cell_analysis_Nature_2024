#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import image
from matplotlib import transforms
import seaborn as sns
import numpy as np
import scipy
import bbknn as bk

import os
import anndata
import logging

import math
import tifffile
import cv2
import gc
import harmonypy as hm





sc.settings.set_figure_params(dpi=300,fontsize=10)





path_main="path_in/"


# # Load and prep data




# load global object
adata=sc.read_h5ad(path_main+"xenium_merged_global_v3.h5ad")





# load global scaled expression
# scaling within a cell type creates too much noise and fragmentation
adata.X=adata.layers["scaled"]





# subset to cell type of interest
adata=adata[adata.obs.cell_type.isin(["Adipocytes"])]
adata.raw=adata.copy()
gc.collect()


# # Clustering




sc.pp.pca(adata,svd_solver='arpack',n_comps=25,random_state=555,use_highly_variable=False)





# use Harmony to correct batch effects in PCA space
pca = adata.obsm['X_pca']
batch = adata.obs['SampleID']
meta_data = adata.obs

ho = hm.run_harmony(pca, meta_data, ["SampleID"], epsilon_harmony = -float('Inf'), max_iter_kmeans = 10, tau=2,max_iter_harmony = 10,random_state=555)

res = pd.DataFrame(ho.Z_corr)
res=res.T
adata.obsm['X_pca'] = res.values





# use BBKNN to get batch corrected neighbors
bk.bbknn(adata,n_pcs=10,local_connectivity=10,batch_key="SampleID",neighbors_within_batch=5)





sc.tl.umap(adata,spread=0.8,min_dist=1.5,random_state=555)





sc.tl.leiden(adata, key_added="cluster",resolution=0.65,random_state=555,n_iterations=-1)





# run DGE to get marker genes
sc.tl.rank_genes_groups(adata, "cluster", method="t-test_overestim_var",use_raw=True)
sc.tl.dendrogram(adata,groupby="cluster")
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, groupby="cluster",standard_scale="var")





# plot select genes to aid annotation
sc.pl.dotplot(adata,["CDKN1A","NAMPT","GRIA4","PGAP1","AGMO","FHOD3","PRR5L","FOSB","JUN","ADIPOQ","PECAM1","PDGFRA","PDGFRB","PTPRC","VGLL3","PNPLA3","MOGAT1","SLC14A2"],groupby="cluster",swap_axes=False, standard_scale="var")





# annotate based on marker genes to relect single-nuc annotation.
# merge clusters driven by a single-gene
# clusters that cannot be confidently matched to a single-nuc cluster are called Unassigned
adata.obs["cell_state"]=adata.obs["cluster"].replace({
"0":"Unassigned",
    "1":"AD1",
    "2":"AD1",
    "3":"AD5",
    "4":"AD5",
    "5":"Unassigned",
    "6":"AD3",
    "7":"AD6",
    "8":"Unassigned",
    "9":"Unassigned",
    })





# save cell type object
adata.write(path_main+"xenium_merged_adipocytes_v3.h5ad")





# export cell state annotation as a csv to add information to global object
cell_state=pd.DataFrame(adata.obs["cell_state"])
cell_state.to_csv("xenium_adip_cell_state_v3.csv")

