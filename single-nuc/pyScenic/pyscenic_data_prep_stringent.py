#!/usr/bin/env python
# coding: utf-8



import pandas as pd
import numpy as np
import gc
import scanpy as sc
import loompy as lp
import re





path_in="input_path/"

path_out="output_path/"





INPUT = path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_cell_states.h5ad"

adata = sc.read(INPUT)

adata.X=adata.layers["raw_counts"]





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





adata.shape





#remove other genes that are not of interest
mito_genes = adata.var_names.str.startswith('MT-')
hb_genes = adata.var_names.str.contains('^HB[^(P)]')
ribo_genes = adata.var_names.str.startswith(("RPS","RPL"))
non_conding= adata.var_names.str.contains('LINC',"orf")
micro= adata.var_names.str.contains("MIR[0-9]")
contig_genes = adata.var_names.str.contains('^A[A-Z][0-9][0-9][0-9][0-9][0-9][0-9]')
as_genes = adata.var.index.str.endswith("-AS1","-AS2")
remove = np.add(mito_genes,hb_genes)
remove = np.add(remove,ribo_genes)
remove = np.add(remove,non_conding)
remove = np.add(remove,contig_genes)
remove = np.add(remove,as_genes)
remove = np.add(remove,micro)
keep = np.invert(remove)
adata = adata[:,keep]





adata.raw=adata.copy()


# # Subsample and export loom objects




#subsample each condition to have the same number of cells
for i in range(11):# do it 10 times due to stochastic nature of subsampling
    wl=adata[adata.obs.condition2.isin(["Weightloss"])]
    sc.pp.subsample(wl,n_obs=20000,random_state=i)
    obese=adata[adata.obs.condition2.isin(["Obese"])]
    sc.pp.subsample(obese,n_obs=20000,random_state=i)
    lean=adata[adata.obs.condition2.isin(["Lean"])]
    sc.pp.subsample(lean,n_obs=20000,random_state=i)
    adata_new=lean.concatenate(wl,obese)
    row_attrs = { 
        "Gene": np.array(adata_new.var.index),
        }
    col_attrs = { 
        "CellID":  np.array(adata_new.obs.index),
        "nGene": np.array( np.sum(adata_new.X.transpose()>0, axis=0)).flatten(),
        "nUMI": np.array( np.sum(adata_new.X.transpose(), axis=0)).flatten(),
        }

    lp.create(path_out+"global_all_stringent_sub"+str(i)+".loom", adata_new.X.transpose(), row_attrs, col_attrs)
    gc.collect()
#export full global as well
row_attrs = { 
    "Gene": np.array(adata.var.index),
    }
col_attrs = { 
    "CellID":  np.array(adata.obs.index),
    "nGene": np.array( np.sum(adata.X.transpose()>0, axis=0)).flatten(),
    "nUMI": np.array( np.sum(adata.X.transpose(), axis=0)).flatten(),
    }

lp.create("global_all_stringent.loom", adata.X.transpose(), row_attrs, col_attrs)







