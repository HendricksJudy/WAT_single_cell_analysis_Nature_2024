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
from scipy import spatial

import os
import anndata
import logging

import math
import tifffile
import cv2
import gc
import sys
import harmonypy as hm
import bbknn as bk





sc.settings.set_figure_params(dpi=150,fontsize=10)





path_main="path_in/"


# # Load and prep data




# load nuc segmentation object
adata=sc.read_h5ad(path_main+"xenium_merged_global_v3.h5ad")





# load and merge cell distance dataframes
df_final=None
for i in adata.obs.SampleID_simple.unique():
    df=pd.read_csv(path_main+i+"_all_distance.csv",index_col=0)
    if df_final is None:
        df_final=df.copy()
    else:
        df_final=pd.concat([df_final,df])





# create columns in distance dataframe that will be replaced with values from a dictionary
df_final["sample_target"]=df_final["sample"].astype(str)+"_"+df_final["target"].astype(str)
df_final["sample_query"]=df_final["sample"].astype(str)+"_"+df_final["query"].astype(str)





# create a column in adata.obs that will be used as dictionary index
adata.obs["sample_cells"]=adata.obs["SampleID_simple"].astype(str)+"_"+adata.obs["cell_ids"].astype(str)





# create dictionary with annotation of interest.
cstate=pd.Series(adata.obs["cell_state"].values,index=adata.obs["sample_cells"]).to_dict()





df_final=df_final[(df_final.distance<300)] # set maximum distance of interest





# create arrary to add annotation to distance dataframe
sample_target = [cstate[x] for x in df_final.sample_target] 
sample_query = [cstate[x] for x in df_final.sample_query]





# add annotation to distance daframe
df_final["target"]=sample_target
df_final["query"]=sample_query





# remove instance where query and target cell are the same (skip to account self for the niche rather than just the neighbors)
df_final["self"]=np.where(df_final["sample_target"]==df_final["sample_query"],"Yes","No")
df_final=df_final[df_final["self"].isin(["No"])]





# Remove Unassigned clusters from target cells
df_final=df_final[~df_final.target.isin(["Unassigned"])]





# create neighborhood matrix
test_df=pd.crosstab(df_final["sample_query"],df_final["target"])





#change matrix index to match adata.obs index
index_dict=pd.Series(adata.obs.index.values,index=adata.obs.sample_cells).to_dict()
test_df["ids"]=test_df.index
test_df["ids"]=test_df["ids"].replace(index_dict)
test_df.index=test_df["ids"]
del test_df["ids"]





# subset adata to only cells in neighbourhood matrix
# some cells might be lost due to not having neighbors within the set distance
adata=adata[adata.obs.index.isin(test_df.index)]





# normalize neighborhood matrix
# this is more important with large distance and adipose tissue due to adipocyte size expansion. 
#Also accounts for edges, which in adipose are just from tissue collection, not real boundaries
test_df=test_df.div(test_df.T.sum(),axis=0)





# reindex neighborhood matrix to match adata index
test_df=test_df.reindex(adata.obs.index)





# create anndata object with neighborhood matrix, and metadata from adata





niche_adata=anndata.AnnData(X=test_df,obs=adata.obs)





niche_adata.raw=niche_adata.copy()


# # Clustering




sc.pp.pca(niche_adata,svd_solver='arpack',n_comps=30,random_state=555,use_highly_variable=False) # PCA here is ran with maximum number of components





# Correct PCA for batch effects with Harmony
pca = niche_adata.obsm['X_pca']
batch = niche_adata.obs['SampleID']
meta_data = niche_adata.obs

ho = hm.run_harmony(pca, meta_data, ["SampleID"], epsilon_harmony = -float('Inf'), max_iter_kmeans = 20, tau=5,max_iter_harmony = 20,random_state=555)

res = pd.DataFrame(ho.Z_corr)
res=res.T
niche_adata.obsm['X_pca'] = res.values





# use BBKNN to get batch corrected neighbors
bk.bbknn(niche_adata,n_pcs=10,local_connectivity=10,batch_key="SampleID",neighbors_within_batch=5)





sc.tl.umap(niche_adata,spread=0.8,min_dist=1.5,random_state=555)





# leiden clustering
sc.tl.leiden(niche_adata, key_added="niche",resolution=0.3,random_state=555,n_iterations=-1)





gc.collect()





# get which cell states are enriched in each niche
sc.tl.rank_genes_groups(niche_adata, "niche", method="t-test_overestim_var",use_raw=True)
sc.tl.dendrogram(niche_adata,groupby="niche")
sc.pl.rank_genes_groups_dotplot(niche_adata, n_genes=5, groupby="niche",standard_scale="var")





#annotate niches. Niches with minor fluctions were merged
niche_adata.obs["niche_final"]=niche_adata.obs["niche"].replace({
    "0":"ad_niche",
    "1":"ad_niche",
    "2":"venous_niche",
    "3":"arterial_niche",
    "4":"stem_niche",
    "5":"stress_niche",
    "6":"ad_niche",
    "7":"ad_niche",
    "8":"venous_niche",
    "9":"arterial_niche"
})





# copy plotting information from adata.obs to be able to plot onto spatial data
niche_adata.obsm["spatial"]=adata.obsm["spatial"]
niche_adata.uns["spatial"]=adata.uns["spatial"]





# plotting example
sub_spatial=niche_adata[niche_adata.obs.SampleID_simple.isin(["A13"])]
sc.pl.spatial(
        sub_spatial,
        img_key="WGA", #whatever image you want as background
        library_id="0014271__A13__20240222__114657", #don't think you need to change this.
        color="niche_final",
        size=10,
        legend_loc=None,
        show=True,
        title=sub_spatial.obs.condition.unique(),
)
plt.show()





# save neighborhood object
niche_adata.write(path_main+"xenium_merged_niches_v3.h5ad")





# add niche clustering to nuc segmentation object and ave as a separate object
adata.obs["niche_raw"]=niche_adata.obs["niche"].copy()
adata.obs["niche_final"]=niche_adata.obs["niche_final"].copy()
adata.write(path_main+"xenium_global_with_niches_v4.h5ad")

