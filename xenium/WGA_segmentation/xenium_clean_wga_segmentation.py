#!/usr/bin/env python
# coding: utf-8


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns
import numpy as np
import scipy
import bbknn as bk

import os
import anndata
import logging

import math
import gc
import harmonypy as hm





sc.settings.set_figure_params(dpi=300,fontsize=10)





path_main="path_in/"


# # Load and prep data




names=["output-XETG00207__0014085__A1__20240208__122452",
"output-XETG00207__0014085__A5__20240208__122452",
"output-XETG00207__0014085__A11__20240208__122452",
"output-XETG00207__0014253__A4__20240208__122452",
"output-XETG00207__0014253__A8__20240208__122452",
"output-XETG00207__0014253__A16__20240208__122452",
"output-XETG00207__0014271__A2__20240222__114657",
"output-XETG00207__0014271__A6__20240222__114657",
"output-XETG00207__0014271__A13__20240222__114657",
"output-XETG00207__0014373__A3__20240222__114657",
"output-XETG00207__0014373__A7__20240222__114657",
"output-XETG00207__0014373__A14__20240222__114657",
]





# load nuc segmentation object
adata_orig=sc.read_h5ad(path_main+"xenium_merged_global.h5ad")





orig_ids=adata_orig.obs.cell_ids.unique()





# load and merge wga segmentation objects
adata=None
for i in names:
    x,y=i.split("XETG00207__")
    adata_tmp=sc.read_h5ad(path_main+i+"/adip_seg.h5ad")
    adata_tmp.var_names_make_unique()
    adata_tmp.obs["SampleID"]=y
    if adata is None:
        adata=adata_tmp.copy()
    else:
        adata=sc.concat([adata,adata_tmp],uns_merge="unique",keys=names,index_unique="-",label="library_id")
    del adata_tmp
    gc.collect()





# simplify sample ID
sample_dict={"0014085__A1__20240208__122452":"A1",
"0014085__A5__20240208__122452":"A5",
"0014085__A11__20240208__122452":"A11",
"0014253__A4__20240208__122452":"A4",
"0014253__A8__20240208__122452":"A8",
"0014253__A16__20240208__122452":"A16",
"0014271__A2__20240222__114657":"A2",
"0014271__A6__20240222__114657":"A6",
"0014271__A13__20240222__114657":"A13",
"0014373__A3__20240222__114657":"A3",
"0014373__A7__20240222__114657":"A7",
"0014373__A14__20240222__114657":"A14",
            }
adata.obs["SampleID_simple"]=adata.obs.SampleID.replace(sample_dict)





adata.layers["raw"]=adata.X.copy()





df=pd.DataFrame(adata.X,columns=adata.var.index,index=adata.obs.index) 





# normalize count values based on total number of counts in the cell. This has to be done because large adipocytes have a higher likelyhood of being assigned noise
adata.X=df.div(adata.obs["n_counts"],axis=0)*5000





adata.layers["norm"]=adata.X.copy()





sc.pp.log1p(adata)





adata.layers["norm_log1p"]=adata.X.copy()
adata.raw=adata.copy()





adata.obs["Area"]=adata.obs["Area"].str.replace("s_","")





# Keep only objects obtained from WGA segmentation
adata.obs["Area"]=np.where(adata.obs["Area"].isin(orig_ids),np.nan,adata.obs["Area"])
adata.obs["Area"]=adata.obs["Area"][pd.to_numeric(adata.obs['Area'], errors='coerce').notnull()]
adata=adata[~adata.obs["Area"].isin([np.nan])]





adata.obs["log10_Area"]=np.log10(adata.obs["Area"])





# add condition information
sample_dict={
"A1":"Obese",
"A5":"Weightloss",
"A11":"Lean",
"A4":"Obese",
"A8":"Weightloss",
"A16":"Lean",
"A2":"Obese",
"A3":"Obese",
"A6":"Weightloss",
"A7":"Weightloss",
"A13":"Lean",
"A14":"Lean"
}
adata.obs["condition"]=adata.obs["SampleID_simple"].replace(sample_dict)





sc.pl.violin(adata,["log10_Area"],groupby="condition",rotation=90)


# # Clustering to identify cells with high noise content




adata_all=adata.copy()





sc.pp.scale(adata_all)





adata_all.layers["scaled"]=adata_all.X.copy()





adata_all.obs[["SampleID_simple","cell_ids_new"]]=adata_all.obs["cell_ids"].str.split("_",expand=True)





adata=adata_all[~adata_all.obs.cell_ids_new.isin(orig_ids)]





adata.X=adata.layers["norm_log1p"]





adata=adata[(adata[:,["ADIPOQ"]].X>0),:] # remove any cells that has no ADIPOQ expression (pan-Adipocyte marker)





adata.X=adata.layers["scaled"]





sc.pp.pca(adata,svd_solver='arpack',n_comps=50,random_state=555,use_highly_variable=False)





# batch correct PCA with Harmony
pca = adata.obsm['X_pca']
batch = adata.obs['SampleID']
meta_data = adata.obs

ho = hm.run_harmony(pca, meta_data, ["SampleID"], epsilon_harmony = -float('Inf'), max_iter_kmeans = 10, tau=2,max_iter_harmony = 10,random_state=555)

res = pd.DataFrame(ho.Z_corr)
res=res.T
adata.obsm['X_pca'] = res.values





# get batch corrected neighbors
bk.bbknn(adata,n_pcs=40,local_connectivity=10,batch_key="SampleID",neighbors_within_batch=5)





sc.tl.umap(adata,spread=0.8,min_dist=1.5,random_state=555)





sc.tl.leiden(adata, key_added="cluster",resolution=0.15,random_state=555,n_iterations=-1)





# run DGE to identify marker genes
sc.tl.rank_genes_groups(adata, "cluster", method="t-test_overestim_var",use_raw=True)
sc.tl.dendrogram(adata,groupby="cluster")
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, groupby="cluster",standard_scale="var")





# hand selected genes to aid identification
sc.pl.dotplot(adata,["ADIPOQ","PDGFRA","VWF","PTPRC","MYH11","POSTN","JUN","CD55","PNPLA3","SLC14A2"],groupby="cluster",swap_axes=False, standard_scale="var")





gc.collect()





# remove clusters that have a high content of genes specific to non-adipocyte linegages
adata=adata[~adata.obs.cluster.isin(["1"])]





gc.collect()





# save clean object
adata.write(path_main+"xenium_wga_seg_adip.h5ad")







