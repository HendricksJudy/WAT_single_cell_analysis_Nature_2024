import scanpy as sc
import pandas as pd
import numpy as np

import os
import anndata
import logging

import gc
import sys
import scenvi

#get condition from PBS submission
condition = "none"
if("--condition" in  sys.argv):
    condition = sys.argv[sys.argv.index("--condition") + 1]       
    
#get round number from PBS submission. Each condition was ran 10 times    
nround = "none"
if("--round" in  sys.argv):
    nround = sys.argv[sys.argv.index("--round") + 1]       

path="path_nuc/"
path_spatial="path_spatial/"

# Read and prep single-nuc data
INPUT = path+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_cell_states.h5ad"
adata = sc.read(INPUT)
adata.X=adata.layers["log1p_counts"]
adata.X=adata.X.todense()
adata=adata[adata.obs.condition2.isin([condition])] # subset 
del adata.layers
gc.collect()

# read and prep spatial data
spatial=sc.read_h5ad(path_spatial+"xenium_merged_global_v3.h5ad")
spatial=spatial[spatial.obs.condition.isin([condition])]

#Run ENVI
adata.var.highly_variable=True # set everything to Highly Variable so that ENVI does not subset any gene
envi_model = scenvi.ENVI(spatial_data = spatial, sc_data = adata, batch_key="SampleID_simple",num_cov_genes=-1,sc_genes=adata.var.index,num_HVG=len(adata.var.index))
envi_model.train()
envi_model.impute_genes()
envi_model.infer_niche_covet()
envi_model.infer_niche_celltype(cell_type_key='cell_type')

spatial.obsm['envi_latent'] = envi_model.spatial_data.obsm['envi_latent']
spatial.obsm['COVET'] = envi_model.spatial_data.obsm['COVET']
spatial.obsm['COVET_SQRT'] = envi_model.spatial_data.obsm['COVET_SQRT']
spatial.uns['COVET_genes'] =  envi_model.CovGenes
spatial.obsm['imputation'] = envi_model.spatial_data.obsm['imputation']
spatial.obsm['cell_type_niche'] = envi_model.spatial_data.obsm['cell_type_niche']

adata.obsm['envi_latent'] = envi_model.sc_data.obsm['envi_latent']
adata.obsm['COVET'] = envi_model.sc_data.obsm['COVET']
adata.obsm['COVET_SQRT'] = envi_model.sc_data.obsm['COVET_SQRT']
adata.obsm['cell_type_niche'] = envi_model.sc_data.obsm['cell_type_niche']
adata.uns['COVET_genes'] =  envi_model.CovGenes

envi_model.infer_niche_celltype(cell_type_key='cell_state_simple')
adata.obsm['cell_state_niche'] = envi_model.sc_data.obsm['cell_type_niche']
spatial.obsm['cell_state_niche'] = envi_model.spatial_data.obsm['cell_type_niche']


print(adata)
print(spatial)

adata.X=scipy.sparse.csr_matrix(adata.X)
spatial.X=scipy.sparse.csr_matrix(spatial.X)
spatial.obsm["imputation"].to_csv(path_spatial+"spatial_envi_imputation_"+condition+"_"+nround+".csv") # save the imputed matrix
del spatial.obsm["imputation"]

#spatial.write(path_spatial+"spatial_envi_trained_imput_"+condition+"_"+nround+".h5ad") # save imputed spatial object (optional)