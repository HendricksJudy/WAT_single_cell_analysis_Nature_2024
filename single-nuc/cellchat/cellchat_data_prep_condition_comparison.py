#!/usr/bin/env python
# coding: utf-8



import pandas as pd
import numpy as np
import gc
import scanpy as sc
import re





path_in="path_in/"
path_out="/path_out/cellchat/input_files/"
INPUT = path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_cell_states.h5ad"

adata = sc.read(INPUT)

adata.X=adata.layers["log1p_counts"]





#load a dataframe with the mean expression and percent expression of each gene for each cell state per sample 
df_expression_summary_type=pd.read_csv(path_in+"expression_matrix_state_am2.csv",index_col=0)





df_expression_summary_type





#create percentage counts column if not present yet, but count information exists
df_expression_summary_type["per_count"]=df_expression_summary_type["count_non_zeroes"]/df_expression_summary_type["count"]





df_expression_summary_type





# create a matrix with the mean percent exprepression per cell type
matrix=pd.crosstab(index=df_expression_summary_type.index, columns=df_expression_summary_type["cell_type_am_fine"], 
                  values=df_expression_summary_type["per_count"], aggfunc='mean')





matrix





# get the maximum percentage expression of each gene
genes_max_per=pd.DataFrame(matrix.T.max(),columns=["max_per_count"])





# subset to genes that have a maximum expression between 5 and 100%
genes_to_keep_05=np.asarray(genes_max_per[(genes_max_per.max_per_count>0.05) & (genes_max_per.max_per_count<1)].index)





adata = adata[:,genes_to_keep_05]





for i in set(adata.obs.condition2):
    gc.collect()
    sub=adata[adata.obs.condition2.isin([i])]
    sub=sub[~sub.obs.cell_type_am_fine.isin(["Unassigned","Lymphatic","B-cells","Lymphoid-Kit+"])] # remove very small clusters because they skew differences and can be driven by a small number of samples
    sub=sub[sub.obs.dataset.isin(["scott"])]
    sc.pp.subsample(sub,n_obs=20000,random_state=0) # subsample to ensure each condition has the same number of cells
    sub.raw=sub
    sub.write(path_out+"global_"+i+"_for_cellchat_no_lymphatic_bcells_kit_scott_only_stringent_subsampled.h5ad")
    gc.collect()







