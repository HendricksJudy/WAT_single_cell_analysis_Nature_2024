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





path_main="/path_in/"


# # Load and merge objects




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





adata=None
for i in names:
    x,y=i.split("XETG00207__")
    adata_tmp=sc.read_h5ad(path_main+i+"/"+y+"cell_ids_nuc_only.h5ad")
    adata_tmp.var_names_make_unique()
    adata_tmp.obs["SampleID"]=y
    excluded=pd.read_csv(path_main+i+"/excluded_cells.csv") # csv with coordinates to exclude due to auto-fluorescence (optional)
    adata_tmp.obs["Exclude"]=np.where((adata_tmp.obs.cell_ids.isin(np.asarray(excluded["Cell ID"]))),"Yes","No") # optional from above
    if adata is None:
        adata=adata_tmp.copy()
    else:
        adata=sc.concat([adata,adata_tmp],uns_merge="unique",keys=names,index_unique="-",label="library_id")
    del adata_tmp
    del excluded
    gc.collect()





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





adata=adata[adata.obs.Exclude.isin(["No"])] # keep only cells from regions with no auto-fluorescence (optional - see above)
gc.collect()





# remove cells with low number of transcripts. Value will depend on panel used, how many genes it has, and how well each cell type is represented.
# determine empyrically. 
adata=adata[adata.obs.n_counts>40] 





adata.layers["raw"]=adata.X.copy()





sc.pp.log1p(adata)





adata.layers["norm_log1p"]=adata.X.copy()
adata.raw=adata.copy()


# # Clustering




sc.pp.scale(adata)





adata.layers["scaled"]=adata.X.copy()





sc.pp.pca(adata,svd_solver='arpack',n_comps=50,random_state=555,use_highly_variable=False)





# run harmony to correct batch effects on PCA space
pca = adata.obsm['X_pca']
batch = adata.obs['SampleID']
meta_data = adata.obs

ho = hm.run_harmony(pca, meta_data, ["SampleID"], epsilon_harmony = -float('Inf'), max_iter_kmeans = 10, tau=2,max_iter_harmony = 10,random_state=555)

res = pd.DataFrame(ho.Z_corr)
res=res.T
adata.obsm['X_pca'] = res.values





# run bbknn to get batch-corrected neighbors
bk.bbknn(adata,n_pcs=40,local_connectivity=10,batch_key="SampleID",neighbors_within_batch=5)





sc.tl.umap(adata,spread=0.8,min_dist=1.5,random_state=555)





sc.tl.leiden(adata, key_added="cluster",resolution=0.2,random_state=555,n_iterations=-1) # low resolution Leiden to get general cell types





# run DGE to get marker genes
sc.tl.rank_genes_groups(adata, "cluster", method="t-test_overestim_var",use_raw=True)
sc.tl.dendrogram(adata,groupby="cluster")
sc.pl.rank_genes_groups_dotplot(adata, n_genes=3, groupby="cluster",standard_scale="var")





# annotate based on marker genes
adata.obs["cell_type"]=adata.obs["cluster"].replace({
"0":"Endothelial",
    "1":"Adipocytes",
    "2":"APC-ASC",
    "3":"Mural",
    "4":"Myeloid",
    "5":"Endothelial",
    "6":"Lymphoid",
    "7":"Mast"
    })





# revert X to non-scaled matrix prior to saving
adata.X=adata.layers["norm_log1p"].copy()





#save object
adata.write(path_main+"xenium_merged_global_V3.h5ad")





gc.collect()

