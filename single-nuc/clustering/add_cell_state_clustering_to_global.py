#!/usr/bin/env python
# coding: utf-8



#import all tools  - if not installed please install them using "pip install toolname --user" or following the instructions on their GitHub page.
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import matplotlib.pyplot as plt
import os
import anndata 

get_ipython().run_line_magic('matplotlib', 'inline')





sc.settings.set_figure_params(dpi=300,fontsize=10,dpi_save=300)




# # Add all cell state information to the Global object




path_in="path_in/"





# load global object
adata=sc.read_h5ad(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated.h5ad")





#load individual cell type object
#alternatively export adata.obs for each object, if memory is limited
apc=sc.read_h5ad(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_APC.h5ad")
endo=sc.read_h5ad(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_Endothelial.h5ad")
mye=sc.read_h5ad(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_Myeloid.h5ad")
mur=sc.read_h5ad(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_Mural.h5ad")
lym=sc.read_h5ad(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_Lymphoid.h5ad")
adip=sc.read_h5ad(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_Adipocytes.h5ad")


# ### Add cell state




df1=pd.concat([apc.obs.cell_state_am.astype(str),endo.obs.cell_state_am.astype(str),mye.obs.cell_state_am.astype(str),mur.obs.cell_state_am.astype(str),lym.obs.cell_state_am.astype(str),adip.obs.cell_state_am.astype(str)])





df2=pd.concat([adata.obs.cell_type_am.astype(str),df1],axis=1)





df2[df2["cell_type_am"].isin(["Mast"])]="Mast"
df2[df2["cell_type_am"].isin(["Lymphatic"])]="Lymphatic"





df2





adata.obs["cell_state_am"]=df2["cell_state_am"]





sc.pl.umap(adata,color="cell_state_am")


# ### Add cell state long




df1=pd.concat([apc.obs.cell_state_am_long.astype(str),endo.obs.cell_state_am_long.astype(str),mye.obs.cell_state_am_long.astype(str),mur.obs.cell_state_am_long.astype(str),lym.obs.cell_state_am_long.astype(str),adip.obs.cell_state_am_long.astype(str)])





df1





df2=pd.concat([adata.obs.cell_type_am.astype(str),df1],axis=1)





df2





df2[df2["cell_type_am"].isin(["Mast"])]="Mast"
df2[df2["cell_type_am"].isin(["Lymphatic"])]="Lymphatic"





df2





adata.obs["cell_state_am_long"]=df2["cell_state_am_long"]





sc.pl.umap(adata,color="cell_state_am_long")





#save object
adata.write(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_cell_states.h5ad")





adata


# # Create annotation layer with more specific cell types, using cell state information

# ### This layer is used for all the cell type level analyses




adata=sc.read_h5ad(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_cell_states.h5ad")





df=pd.DataFrame(adata.obs["bmi"].astype(float))
df["condition2"]="Obese"
df["condition2"][df["bmi"]<(25.5)]="Lean"
adata.obs["condition2"]=adata.obs["condition"].astype(str)
adata.obs["condition2"][adata.obs.condition.isin(["NA"])]=df["condition2"]





cell_state_rename={
    "AD1":"Adipocytes",
    "AD2":"Adipocytes",
    "AD3":"Adipocytes",
    "AD4":"Adipocytes",
    "AD5":"Adipocytes",
    "AD6":"Adipocytes",
    "AD7":"Adipocytes",
    "AD8":"Adipocytes",
    "EC1":"Endothelial",
    "EC2":"Endothelial",
    "EC3":"Endothelial",
    "EC4":"Endothelial",
    "EC5":"Endothelial",
    "EC6":"Endothelial",
    "EC7":"Endothelial",
    "APC1":"ASC",
    "APC2":"APC",
    "APC3":"APC",
    "APC4":"APC",
    "APC5":"APC",
    "MYE1":"Macrophages",
    "MYE2":"Macrophages",
    "MYE3":"Macrophages",
    "MYE4":"Mono_DC",
    "MYE5":"Mono_DC",
    "MYE6":"Macrophages",
    "MYE7":"Mono_DC",
    "MYE8":"Macrophages",
    "MYE9":"Mono_DC",
    "MYE10":"Mono_DC",
    "Mu1":"Mural",
    "Mu2":"Mural",
    "Mu3":"Mural",
    "Mu4":"Mural",
    "Mu5":"Mural",
    "NK1":"NK_cells",
    "NK2":"NK_cells",
    "T1":"T-cells_CD4+",
    "T2":"T-cells_CD8+",
    "T3":"T-cells_CD4+",
    "T4":"T-cells_CD8+",
    "T5":"T-cells_CD4+",
    "T6":"T-cells_CD4+",
    "T7":"ILC-Kit+",
}





adata.obs["cell_type_am_fine"]=adata.obs.cell_state_am.replace(cell_state_rename)





sc.pl.umap(adata,color="cell_type_am_fine")


# # Add Subclustering UMAP information to global adata.obs

# ### umap with doublets included




df1=pd.concat([apc.obs.UMAP1_doublets.astype(str),
               endo.obs.UMAP1_doublets.astype(str),
               mye.obs.UMAP1_doublets.astype(str),
               mur.obs.UMAP1_doublets.astype(str),
               lym.obs.UMAP1_doublets.astype(str),
               adip.obs.UMAP1_doublets.astype(str)])





df1





df2=pd.concat([apc.obs.UMAP2_doublets.astype(str),
               endo.obs.UMAP2_doublets.astype(str),
               mye.obs.UMAP2_doublets.astype(str),
               mur.obs.UMAP2_doublets.astype(str),
               lym.obs.UMAP2_doublets.astype(str),
               adip.obs.UMAP2_doublets.astype(str)])





df2





df3=pd.concat([adata.obs.cell_type_am.astype(str),df1,df2],axis=1)





df3


# ### Add umap with doublets removed




df1=pd.concat([apc.obs.UMAP1_clean.astype(str),
               endo.obs.UMAP1_clean.astype(str),
               mye.obs.UMAP1_clean.astype(str),
               mur.obs.UMAP1_clean.astype(str),
               lym.obs.UMAP1_clean.astype(str),
               adip.obs.UMAP1_clean.astype(str)])





df1





df2=pd.concat([apc.obs.UMAP2_clean.astype(str),
               endo.obs.UMAP2_clean.astype(str),
               mye.obs.UMAP2_clean.astype(str),
               mur.obs.UMAP2_clean.astype(str),
               lym.obs.UMAP2_clean.astype(str),
               adip.obs.UMAP2_clean.astype(str)])





df2


# ### create merged dataframe




df3=pd.concat([df3,df1,df2],axis=1)





df3





# add to adata.obs
adata.obs["UMAP1_doublets"]=df3["UMAP1_doublets"] # use for umaps with doublet included
adata.obs["UMAP2_doublets"]=df3["UMAP2_doublets"]# use for umaps with doublet included
adata.obs["UMAP1_clean"]=df3["UMAP1_clean"]# use for umaps with doublet removed
adata.obs["UMAP2_clean"]=df3["UMAP2_clean"]# use for umaps with doublet removed





adata





#save object
adata.write(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_cell_states.h5ad")


# ### Plot specific cell type UMAPs from Global




#save the original coordinates
adata.obsm["X_pca_global"]=adata.obsm["X_pca"].copy()
adata.obsm["X_umap_global"]=adata.obsm["X_umap"].copy()





adata





#for the non doublet umap
cell_type="Adipocytes" # choose which cell type to plot
subset=adata[(adata.obs.cell_type_am.isin([cell_type]))
             &(~adata.obs.cell_state_am.isin(["Unassigned"]))# for the umap with doublets included, do not remove Unassigned
             ].copy() 
coord1=np.asarray(subset.obs["UMAP1_clean"])
coord1=np.split(coord1,len(subset.obs.index))
coord2=np.asarray(subset.obs["UMAP2_clean"])
coord2=np.split(coord2,len(subset.obs.index))
obsm=np.concatenate((coord1,coord2),axis=1).astype(float)
subset.obsm["X_umap"]=obsm





sc.pl.umap(subset,color=["cell_state_am_long"])







