#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import gc
import scanpy as sc
import loompy as lp
import re


# # Prep adipocyte data
# ### Macrophage metabolic networks were performed similarly




path_in="input_path/"

path_out="output_path/"





INPUT = path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_cell_states.h5ad"

adata = sc.read(INPUT)

adata.X=adata.layers["raw_counts"]





# load compass genes
df=pd.read_csv("compass_genes.csv",index_col=0)





df=df[(df.cell_type.isin(["Adipocytes"]))]





genes_to_keep=[x for x in df.associated_genes.unique() if x in adata.var.index]





tf=pd.read_csv("allTFs_hg38.txt",header=None)





tf=np.asarray(tf[0])





# only compass genes and TFs will be kept
genes_to_keep=np.append(genes_to_keep,tf)





genes_to_keep=np.unique(genes_to_keep)





genes_to_keep=[x for x in genes_to_keep if x in adata.var.index]





adata = adata[:,genes_to_keep]





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


# # Create loom objects




#subsample each condition to have the same number of cells
for i in range(11): # do it 10 times due to stochastic nature of subsampling
    wl=adata[(adata.obs.condition2.isin(["Weightloss"]))&(adata.obs.cell_type_am_fine.isin(["Adipocytes"]))] # macrophages were subsampled to 3000 instead
    sc.pp.subsample(wl,n_obs=9000,random_state=i)
    obese=adata[(adata.obs.condition2.isin(["Obese"]))&(adata.obs.cell_type_am_fine.isin(["Adipocytes"]))]# macrophages were subsampled to 3000 instead
    sc.pp.subsample(obese,n_obs=9000,random_state=i)
    lean=adata[(adata.obs.condition2.isin(["Lean"]))&(adata.obs.cell_type_am_fine.isin(["Adipocytes"]))]# macrophages were subsampled to 3000 instead
    sc.pp.subsample(lean,n_obs=9000,random_state=i)
    adata_new=lean.concatenate(wl,obese)
    row_attrs = { 
        "Gene": np.array(adata_new.var.index),
        }
    col_attrs = { 
        "CellID":  np.array(adata_new.obs.index),
        "nGene": np.array( np.sum(adata_new.X.transpose()>0, axis=0)).flatten(),
        "nUMI": np.array( np.sum(adata_new.X.transpose(), axis=0)).flatten(),
        }

    lp.create(path_out+"global_all_stringent_sub_compass_adipocytes"+str(i)+".loom", adata_new.X.transpose(), row_attrs, col_attrs)
    gc.collect()





#export full global as well
adata=adata[adata.obs.cell_type_am_fine.isin(["Adipocytes"])]
row_attrs = { 
    "Gene": np.array(adata.var.index),
    }
col_attrs = { 
    "CellID":  np.array(adata.obs.index),
    "nGene": np.array( np.sum(adata.X.transpose()>0, axis=0)).flatten(),
    "nUMI": np.array( np.sum(adata.X.transpose(), axis=0)).flatten(),
    }

lp.create(path_out+"global_all_stringent_compass_adipocytes.loom", adata.X.transpose(), row_attrs, col_attrs)







