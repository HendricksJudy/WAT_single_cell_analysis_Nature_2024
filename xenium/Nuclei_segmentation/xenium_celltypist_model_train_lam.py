#!/usr/bin/env python
# coding: utf-8


import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import celltypist as ct

import matplotlib.pyplot as plt

import gc





path_in="/path_in/"


# # Load and prep data




# load single nuc data
INPUT = path+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_cell_states.h5ad"
adata = sc.read(INPUT)
adata.X=adata.layers["raw_counts"]
del adata.layers
gc.collect()





#select only myeloid cells
cell_type="Myeloid"
adata=adata[(adata.obs.cell_type_am.isin([cell_type]))&(~adata.obs.cell_state_am.isin(["Unassigned"])),].copy()
coord1=np.asarray(adata.obs["UMAP1_clean"])
coord1=np.split(coord1,len(adata.obs.index))
coord2=np.asarray(adata.obs["UMAP2_clean"])
coord2=np.split(coord2,len(adata.obs.index))
obsm=np.concatenate((coord1,coord2),axis=1).astype(float)
adata.obsm["X_umap"]=obsm
gc.collect()





# load cell state information with LAM subclustering if not in object
trem2=pd.read_csv("trem_subclustering.csv",index_col=0)





adata.obs["cell_state_2"]=trem2["cell_state_2"]





# Subset only to LAM subtypes
adata=adata[adata.obs.cell_state_2.isin(['LAM_ST1','LAM_ST2','LAM_ST3','LAM_ST4'])]





# load spatial data
spatial=sc.read_h5ad(path+"xenium_merged_global_v3.h5ad")





# select only genes that are both in the spatial and single-nuc data
common_genes=[x for x in spatial.var.index if x in adata.var.index]





#subset single nuc to common genes only. This prevents the model from including genes that cannot be used to predict in the spatial data
adata=adata[:,common_genes]





# normalize and log single-nuc data as we loaded the raw data above.
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)





adata.raw=adata.copy()


# # Train the model - 1st round




new_model = ct.train(adata, labels = 'cell_state_2', n_jobs = 8, feature_selection = True,balance_cell_type=True,use_SGD=True,alpha=0.001)
new_model.write('human_adip_lams_celltypist_fine_SGD_alpha001_spatial.pkl')





# use the model on the training dataset with probability match
predictions = ct.annotate(adata, model = "human_adip_lams_celltypist_fine_SGD_alpha001_spatial.pkl", majority_voting = True, mode = 'prob match', p_thres = 0.50)
adata = predictions.to_adata(insert_prob = False)
# create an object with only cells that have an exact match
clusters=list(set(adata.obs.cell_state_2))
adata_clean=adata[adata.obs.predicted_labels.isin(clusters)]
gc.collect()


# # Train the model - 2nd round




# train a 2nd model on the predicted labels from 1st round. This creates a model based only the cells that have the most robust identity, eliminating some noise.
new_model = ct.train(adata_clean, labels = 'predicted_labels', n_jobs = 8, feature_selection = True,balance_cell_type=True,use_SGD=True,alpha=0.01)
new_model.write('human_adip_lams_celltypist_fine_SGD_alpha001_r2_prob05_01_spatial.pkl')







